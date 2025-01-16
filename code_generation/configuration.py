from __future__ import annotations  # needed for type annotations in > python 3.7

import logging
import copy
from typing import Any, Dict, List, Set, Union

from code_generation.exceptions import (
    ScopeConfigurationError,
    ConfigurationError,
    EraConfigurationError,
    InvalidOutputError,
    SampleConfigurationError,
    InvalidShiftError,
)
from code_generation.modifiers import EraModifier, SampleModifier
from code_generation.optimizer import ProducerOrdering
from code_generation.producer import (
    ProducerGroup,
    CollectProducersOutput,
    TProducerInput,
    TProducerStore,
)
from code_generation.quantity import (
    NanoAODQuantity,
    QuantitiesInput,
    QuantitiesStore,
    QuantityGroup,
)
from code_generation.rules import ProducerRule, RemoveProducer
from code_generation.systematics import SystematicShift, SystematicShiftByQuantity
from code_generation.helpers import is_empty

log = logging.getLogger(__name__)
# type aliases
TConfiguration = Dict[
    str,
    Union[
        List[Union[str, int, float, bool, Dict[Any, Any]]],
        str,
        int,
        float,
        bool,
        EraModifier,
        SampleModifier,
        Dict[Any, Any],
    ],
]


class Configuration(object):
    """
    Configuration class for for the CROWN configuration. This class
    holds all parts of the configuration, from the sample, era, scope,
    and systematics, to the output. All modifications to the configuration should be done through this class.
    """

    def __init__(
        self,
        era: str,
        sample: str,
        scopes: Union[str, List[str]],
        shifts: Union[str, List[str], Set[str]],
        available_sample_types: Union[str, List[str]],
        available_eras: Union[str, List[str]],
        available_scopes: Union[str, List[str]],
        global_scope: bool = True,
    ):
        """

        Initiation of a new configuration.

        Args:
            era: The era of the sample.
            sample: The sample type of the sample.
            scopes: The scopes to be used in the configuration.
                By default, the global scope will always be added and run first, all other scopes will be run in
                parallel as children of the global scope.
            shifts: The systematics to be used in the configuration.
            available_sample_types: The available sample types.
            available_eras: The available eras.
            available_scopes: The available scopes.

        """
        self.era = era
        self.sample = sample
        self.selected_scopes = set(scopes)
        self.selected_shifts = shifts
        self.available_sample_types = set(available_sample_types)
        self.available_eras = set(available_eras)
        self.available_scopes = set(available_scopes)
        self.available_outputs: Dict[str, QuantitiesStore] = {}
        self.available_shifts: Dict[str, Set[str]] = {}
        if global_scope:
            self.global_scope = "global"
        else:
            self.global_scope = None

        self.producers: TProducerStore = {}
        self.unpacked_producers: TProducerStore = {}
        self.scopes: List[str] = []
        self.outputs: QuantitiesStore = {}
        self.shifts: Dict[str, Dict[str, TConfiguration]] = {}
        self.rules: Set[ProducerRule] = set()
        self.config_parameters: Dict[str, TConfiguration] = {}

        self.setup_defaults()

    def _validate_scopes(self) -> None:
        """
        Function to validate the selected scopes. If the scope is not available, an error is raised.

        Args:
            None

        Returns:
            None
        """
        missing_scopes = self.selected_scopes - self.available_scopes
        if len(missing_scopes) > 0:
            raise ScopeConfigurationError(missing_scopes, self.available_scopes)

    def _validate_sample_type(self) -> None:
        """
        Function to validate the selected sample type. If the sample type is not available, an error is raised.

        Args:
            None

        Returns:
            None
        """
        if self.sample not in self.available_sample_types:
            raise SampleConfigurationError(self.sample, self.available_sample_types)

    def _validate_era(self) -> None:
        """
        Function to validate the selected era. If the era is not available, an error is raised.

        Args:
            None

        Returns:
            None
        """
        if self.era not in self.available_eras:
            raise EraConfigurationError(self.era, self.available_eras)

    def _set_sample_parameters(self) -> None:
        """
        Helper function to add sample type variables to the configuration.
        The variables look like ``is_${sampletype}`` and can be used in all producer
        calls to check the type of the sample.

        Args:
            None

        Returns:
            None
        """
        sample_parameters: Dict[str, bool] = {}
        for sampletype in self.available_sample_types:
            if self.sample == sampletype:
                sample_parameters["is_{}".format(sampletype)] = True
            else:
                sample_parameters["is_{}".format(sampletype)] = False
        for scope in self.scopes:
            self.config_parameters[scope].update(sample_parameters)

    def setup_defaults(self) -> None:
        """
        Function used to add some defaults to the configuration. This function is called by the ``__init__`` function.
        For all configured scopes, the nessessay variables are added to the configuration.
        The validation of the initial settings is also done here.

        Args:
            None

        Returns:
            None
        """
        self._validate_scopes()
        self._validate_sample_type()
        self._validate_era()
        self.scopes = [self.global_scope]
        for scope in self.available_scopes:
            self.scopes.append(scope)
        for scope in self.scopes:
            self.producers[scope] = []
            self.unpacked_producers[scope] = {}
            self.outputs[scope] = set()
            self.shifts[scope] = {}
            self.available_outputs[scope] = {}
            for sampletype in self.available_sample_types:
                self.available_outputs[scope][sampletype] = set()
            self.config_parameters[scope] = {}
            self.available_shifts[scope] = set()
        self._set_sample_parameters()

    def add_config_parameters(
        self, scopes: Union[str, List[str]], parameters: TConfiguration
    ) -> None:
        """
        Function to add new config parameters to the configuration. Modifiers are used to
        modify the configuration, are directly resolved using the resolve_modifiers function.

        Args:
            scopes: The scopes to which the parameters should be added. This can be a list of scopes or a single scope.
            parameters: The parameters to be added. This must be a dictionary of parameters.
                If multiple scopes are given, the parameters are added to all scopes.

        Returns:
            None
        """
        if not isinstance(scopes, list):
            scopes = [scopes]
        for scope in scopes:
            self.config_parameters[scope].update(self.resolve_modifiers(parameters))

    def add_producers(
        self, scopes: Union[str, List[str]], producers: TProducerInput
    ) -> None:
        """
        Function used to add producers to the configuration. Internally, a set of all
        available outputs is updated, which is later used to check if all required ouputs are available.

        Args:
            scopes: The scopes to which the producers should be added. This can be a list of scopes or a single scope.
            producers: The producers to be added. If multiple scopes are given, the producers are added to all scopes.

        Returns:
            None
        """
        if not isinstance(scopes, list):
            scopes = [scopes]
        if not isinstance(producers, list):
            producers = [producers]
        for scope in scopes:
            self.producers[scope].extend(producers)
            for sampletype in self.available_sample_types:
                self.available_outputs[scope][sampletype].update(
                    CollectProducersOutput(producers, scope)
                )
            self.unpack_producergroups(scope, producers)

    def unpack_producergroups(
        self,
        scope: str,
        producers: TProducerInput,
        parent: Union[TProducerInput, None] = None,
        depth: int = 0,
    ) -> None:
        """
        Function used to add producers to the configuration. Internally, a set of all
        available outputs is updated, which is later used to check if all required ouputs are available.

        Args:
            scope: The scope to which the producers should be added.
            producers: The producers to be added.

        Returns:
            None
        """

        if isinstance(producers, list):
            # we always want to know the toplevel producergroup, so if the parent evaluates to false, we set it to the first producer.
            # If a prent is given, we set it to the parent, since this means we are in a producergroup. This is important if we
            # have nested producergroups, this way every producer is assigned to the outermost producergroup, which is important for the
            # potential removal of a single producer.
            for producer in producers:
                if is_empty(parent):
                    parent_producer = producer
                else:
                    parent_producer = parent
                self.unpack_producergroups(
                    scope=scope,
                    producers=producer,
                    parent=parent_producer,
                    depth=depth + 1,
                )
        else:
            if isinstance(producers, ProducerGroup):
                log.debug("{} Unpacking ".format("    " * depth))
                for sub_producer in producers.producers[scope]:
                    if is_empty(parent):
                        parent_producer = producers
                    else:
                        parent_producer = parent
                        self.unpack_producergroups(
                            scope=scope,
                            producers=sub_producer,
                            parent=parent_producer,
                            depth=depth + 1,
                        )
            else:
                if is_empty(parent):
                    log.debug("{} {}".format("    " * depth, producers))
                    self.unpacked_producers[scope][producers] = producers
                else:
                    log.debug("{} {}".format("    " * depth, parent))
                    self.unpacked_producers[scope][producers] = parent

    def add_outputs(
        self, scopes: Union[str, List[str]], output: QuantitiesInput
    ) -> None:
        """
        Function used to add outputs to the configuration.

        Args:
            scopes: The scopes to which the outputs should be added.
                This can be a list of scopes or a single scope.
            output: The outputs to be added. If multiple scopes are given, the outputs are added to all scopes.

        Returns:
            None
        """
        if not isinstance(scopes, list):
            scopes = [scopes]
        if not isinstance(output, list):
            output = [output]
        for scope in scopes:
            self.outputs[scope].update(output)

    def add_shift(
        self,
        shift: Union[SystematicShift, SystematicShiftByQuantity],
        samples: Union[str, List[str], None] = None,
        exclude_samples: Union[str, List[str], None] = None,
    ) -> None:
        """
        Function used to add a systematics shift to the configuration. During this step, the shift is validated and applied.

        Args:
            shift: The shift to be added. This must be a SystematicShift object.
            samples: The samples to which the shift should be applied. This can be a list of samples or a single sample.
                If ths option is not set, the shift is applied, regardless of the sample type.
            exclude_samples: The samples to which the shift should not be applied. This can be a list of samples or a single sample. if exclude_samples is set, samples cannot be used.

        Returns:
            None
        """
        if not is_empty(exclude_samples) and not is_empty(samples):
            raise ConfigurationError(
                f"You cannot use samples and exclude_samples at the same time -> Shift {shift}, samples {samples}, exclude_samples {exclude_samples}"
            )
        if not is_empty(samples):
            if isinstance(samples, str):
                samples = [samples]
            for sample in samples:
                if sample not in self.available_sample_types:
                    raise ConfigurationError(
                        f"Sampletype {sample} is not available -> Shift {shift}, available_sample_types {self.available_sample_types}, sample_types {samples}"
                    )
        if not is_empty(exclude_samples):
            if isinstance(exclude_samples, str):
                exclude_samples = [exclude_samples]
            for excluded_sample in exclude_samples:
                if excluded_sample not in self.available_sample_types:
                    raise ConfigurationError(
                        f"Sampletype {excluded_sample} is not available for exclusion -> Shift {shift}, available_sample_types {self.available_sample_types}, excluded_sample_types {exclude_samples}"
                    )
            samples = list(self.available_sample_types - set(exclude_samples))
        if self._is_valid_shift(shift):
            log.debug("Shift {} is valid".format(shift.shiftname))
            if not isinstance(shift, SystematicShift):
                raise TypeError("shift must be of type SystematicShift")
            if isinstance(samples, str):
                samples = [samples]
            if is_empty(samples) or self.sample in samples:
                scopes_to_shift = [
                    scope for scope in shift.get_scopes() if scope in self.scopes
                ]
                # if a modifier is used within the shift config, we have to resolve it
                #  here using the resolve_modifiers function
                if self.global_scope in scopes_to_shift:
                    for scope in self.scopes:
                        if scope in shift.get_scopes():
                            self._add_available_shift(shift, scope)
                            shift.apply(scope)
                            self.shifts[scope][
                                shift.shiftname
                            ] = self.resolve_modifiers(shift.get_shift_config(scope))
                        else:
                            self._add_available_shift(shift, scope)
                            shift.apply(self.global_scope)
                            self.shifts[scope][
                                shift.shiftname
                            ] = self.resolve_modifiers(
                                shift.get_shift_config(self.global_scope)
                            )
                else:
                    for scope in scopes_to_shift:
                        self._add_available_shift(shift, scope)
                        shift.apply(scope)
                        self.shifts[scope][shift.shiftname] = self.resolve_modifiers(
                            shift.get_shift_config(scope)
                        )

    def _is_valid_shift(
        self, shift: Union[SystematicShift, SystematicShiftByQuantity]
    ) -> bool:
        """
        Function to check if a shift is valid. A shift is condisered valid, if its name
        matches the name of the shifts prodivded in self.selected_shifts.
        If none is selected, all shifts are invalid, if all is selected, all shifts are valid.

        Args:
            shift: The shift to be checked.

        Returns:
            bool: True if the shift is valid, False otherwise.
        """
        # first check if the scopes are correct
        for scope in shift.get_scopes():
            if scope not in list(self.available_scopes) + [self.global_scope]:
                raise ConfigurationError(
                    "Shift {} has scope {} which is not in the list of avialble scopes {}".format(
                        shift.shiftname, scope, self.available_scopes
                    )
                )
        if len(self.selected_shifts) == 1 and "all" in self.selected_shifts:
            return True
        elif len(self.selected_shifts) == 1 and "none" in self.selected_shifts:
            return False
        else:
            return any(
                [
                    shiftname.lower() in shift.shiftname.lower()
                    for shiftname in self.selected_shifts
                ]
            )

    def _add_available_shift(
        self, shift: Union[SystematicShift, SystematicShiftByQuantity], scope
    ) -> None:
        """Add a shift to the set of available shifts

        Args:
            shift: The shift to be added.
            scope: The scope to which the shift should be added.
        """
        log.debug("Adding shift {} to scope {}".format(shift.shiftname, scope))
        self.available_shifts[scope].add(shift.shiftname.lower())

    def add_modification_rule(
        self, scopes: Union[str, List[str]], rule: ProducerRule
    ) -> None:
        """
        Function used to add a rule to the configuration.

        Args:
            scopes: The scopes to which the rule should be added. This can be a list of scopes or a single scope.
            rule: The rule to be added. This must be a ProducerRule object.

        Returns:
            None

        """
        if not isinstance(rule, ProducerRule):
            raise TypeError("Rule must be of type ProducerRule")
        if not isinstance(scopes, list):
            scopes = [scopes]
        rule.set_available_sampletypes(self.available_sample_types)
        rule.set_scopes(scopes)
        rule.set_global_scope(self.global_scope)
        self.rules.add(rule)

    def resolve_modifiers(self, configuration_dict: Dict[Any, Any]) -> TConfiguration:
        """
        Function used to resolve mofifiers used in the configuration. This function is called by the add_config_parameters function.

        Args:
            configuration_dict: The configuration dictionary to be resolved.

        Returns:
            dict: The resolved configuration dictionary.
        """
        resolved_dict: TConfiguration = {}
        for key in configuration_dict:
            resolved_value = None
            log.debug("Resolved: {}".format(resolved_dict))
            log.debug("Parameter: {}".format(key))
            log.debug("Paramters: {}".format(configuration_dict))
            if isinstance(configuration_dict[key], dict):
                resolved_dict.update(self.resolve_modifiers(configuration_dict[key]))
            elif isinstance(configuration_dict[key], SampleModifier):
                log.debug("Applying Samplemodifier")
                resolved_value = configuration_dict[key].apply(self.sample)
            elif isinstance(configuration_dict[key], EraModifier):
                log.debug("Applying Eramodifier")
                resolved_value = configuration_dict[key].apply(self.era)
            else:
                resolved_value = configuration_dict[key]
            # after resolving, step down iterately, if nested Modifiers are used
            if isinstance(resolved_value, list):
                for i, value in enumerate(resolved_value):
                    if isinstance(value, dict):
                        resolved_value[i] = self.resolve_modifiers(value)
                    elif isinstance(value, SampleModifier):
                        resolved_value[i] = value.apply(self.sample)
                    elif isinstance(value, EraModifier):
                        resolved_value[i] = value.apply(self.era)
                    else:
                        resolved_value[i] = value
            resolved_dict.update({key: resolved_value})
        return resolved_dict

    def _remove_empty_scopes(self) -> None:
        """
        Internal function used to remove empty scopes from the configuration. This function is called by the optimize function,
        which should be called after all configuration changes have been made.

        Args:
            None

        Returns:
            None
        """
        # we have to use a seperate list, because we cannot modify the list while iterating over it without breaking stuff
        scopes_to_test = [scope for scope in self.scopes]
        for scope in scopes_to_test:
            if (
                len(self.producers[scope]) == 0 or scope not in self.selected_scopes
            ) and scope is not self.global_scope:
                log.warning("Removing unrequested / empty scope {}".format(scope))
                self.scopes.remove(scope)
                del self.producers[scope]
                del self.outputs[scope]
                del self.shifts[scope]
                del self.config_parameters[scope]
                for sampletype in self.available_sample_types:
                    del self.available_outputs[scope][sampletype]

    def _apply_rules(self) -> None:
        """
        Internal function used to apply all rules to the configuration.
        This function is called by the optimize function.

        """
        for rule in self.rules:
            rule.apply(
                self.sample, self.producers, self.unpacked_producers, self.outputs
            )
            # also update the set of available outputs if the affected sample is the current sample
            if self.sample in rule.samples:
                for scope in rule.affected_scopes():
                    if isinstance(rule, RemoveProducer):
                        self.available_outputs[scope][
                            self.sample
                        ] -= CollectProducersOutput(rule.affected_producers(), scope)
                    else:
                        self.available_outputs[scope][self.sample].update(
                            CollectProducersOutput(rule.affected_producers(), scope)
                        )

    def optimize(self) -> None:
        """
        Function used to optimize the configuration. Optimizaion steps are:

            1. Remove empty scopes
            2. Apply rules
            3. Optimizing producer ordering (this does not change the configuration, but only the order of producers)

        Args:
            None

        Returns:
            None
        """
        self._apply_rules()
        self._remove_empty_scopes()
        for scope in self.scopes:
            log.debug("Optimizing Producer Ordering in scope {}".format(scope))
            ordering = ProducerOrdering(
                global_producers=list(self.producers[self.global_scope]),
                scope=scope,
                producer_ordering=list(self.producers[scope]),
            )
            ordering.Optimize()
            self.producers[scope] = ordering.optimized_ordering

    def _validate_outputs(self) -> None:
        """
        Function used to validate the defined outputs. If an output is requested in the configuration,
        but is not available, since no producer will be able to produce it, an error is raised.

        Args:
            None

        Returns:
            None
        """
        for scope in [scope for scope in self.scopes if scope != self.global_scope]:
            required_outputs = set(
                output
                for output in self.outputs[scope] | self.outputs[self.global_scope]
                if not isinstance(output, NanoAODQuantity)
            )
            # merge the two sets of outputs
            provided_outputs = (
                self.available_outputs[scope][self.sample]
                | self.available_outputs[self.global_scope][self.sample]
            )
            missing_outputs = required_outputs - provided_outputs
            if len(missing_outputs) > 0:
                raise InvalidOutputError(scope, missing_outputs)

    def _remove_empty_configkeys(self, config) -> None:
        """
        Function used to remove empty configuration parameters from the configuration. Empty parameters are
        * empty list: []
        * empty dict: {}
        * empty string: ""
        * Nonetype: None

            Args:
                None

            Returns:
                None
        """
        for key in list(config):
            if isinstance(config[key], dict):
                self._remove_empty_configkeys(config[key])
            # special case for extended vector producers, here we can have a list, that contains empty dicts
            elif isinstance(config[key], list):
                subdict = copy.deepcopy(config[key])
                for i, value in enumerate(subdict):
                    if value == {}:
                        log.info(
                            "Removing {}, (from {}) since it is an empty configuration parameter".format(
                                value, key
                            )
                        )
                        config[key].remove(value)
                    # does this work ?
                    if isinstance(value, dict):
                        self._remove_empty_configkeys(value)

            elif is_empty(config[key]):
                log.info(
                    "Removing {} since it is an empty configuration parameter".format(
                        key
                    )
                )
                del config[key]

    def _validate_all_shifts(self) -> None:
        """
        Function to validate the set of selected shifts against the set of available shifts.
        If a shift is required, that is not set, an error is raised.

        Args:
            None

        Returns:
            None
        """
        if len(self.selected_shifts) == 0:
            raise ConfigurationError("No shifts selected")
        elif len(self.selected_shifts) == 1 and "all" in self.selected_shifts:
            log.info("All shifts are selected")
        elif len(self.selected_shifts) == 1 and "none" in self.selected_shifts:
            log.info("Nominal is selected, no shifts will be applied")
            for scope in self.scopes:
                self.shifts[scope] = {}
        else:
            # check if all selected shifts are available
            for scope in self.available_shifts:
                if len(self.available_shifts[scope]) == 0:
                    continue
                # we do not need to check the global scope, since shifts from
                # the global scope are always propagated down to all scopes
                if scope in self.selected_scopes:
                    for shift in self.shifts[scope].keys():
                        log.debug(
                            "Validating shift {} in scope {}".format(shift, scope)
                        )
                        if not any(
                            [
                                shift.lower() in available_shift
                                for available_shift in self.available_shifts[scope]
                            ]
                        ):
                            raise InvalidShiftError(shift, self.sample, scope)
        log.info("Shift configuration is valid")

    def _validate_parameters(self) -> None:
        """
        Function used to check, if parameters are set in the configuration, that are not used by any producer.
        """
        # first collect all parameters that are used by any producer
        log.info("------------------------------------")
        log.info("Checking for unused parameters")
        log.info(
            "Unused parameters are not an error, but can be a sign of a misconfiguration"
        )

        required_parameters = {}
        for scope in self.scopes:
            log.debug("  Collecting required parameters for scope {}".format(scope))
            required_parameters[scope] = set()
            for producer in self.producers[scope]:
                log.debug("    Collecting parameters for producer {}".format(producer))
                try:
                    required_parameters[scope] |= producer.parameters[scope]
                except KeyError:
                    log.error(
                        f"Producer {producer} is not correctly setup for scope {scope}"
                    )
                    raise ConfigurationError(
                        "Producer {} is not correctly setup for scope {}".format(
                            producer, scope
                        )
                    )
        # now check, which parameters are set in the configuration, but not used by any producer
        for scope in self.scopes:
            log.info("------------------------------------")
            log.info("  Validating parameters in scope {}".format(scope))
            log.info("------------------------------------")
            for parameter in self.config_parameters[scope]:
                # the sample parameters like is_data are skipped here, since they are added by default to every scope
                sampletype_parameters = [
                    f"is_{sampletype}" for sampletype in self.available_sample_types
                ]
                if (
                    parameter not in required_parameters[scope]
                    and parameter not in sampletype_parameters
                ):
                    log.info(
                        "    [{} Scope] Parameter {} is set in the configuration, but not used by any requested producer".format(
                            scope, parameter
                        )
                    )
            log.info("------------------------------------")

    def validate(self) -> None:
        """
        Function used to validate the configuration. During the validation, the following steps are performed:

            - Validate the outputs
            - Validate the shifts

        Args:
            None

        Returns:
            None
        """
        self._validate_outputs()
        self._validate_parameters()
        self._remove_empty_configkeys(self.config_parameters)
        self._validate_all_shifts()

    def report(self) -> None:
        """
        Function used to log a summary of the configuration in form of a report.

        Args:
            None

        Returns:
            None
        """
        running_scopes = self.scopes
        total_producers = sum([len(self.producers[scope]) for scope in running_scopes])
        # if a ExtendedVectorProducer is used, count the correct number of output quantities to be written out
        total_quantities = [
            sum(
                [
                    (
                        len(self.config_parameters[scope][output.vec_config])
                        if isinstance(output, QuantityGroup)
                        else 1
                    )
                    for output in self.outputs[scope]
                ]
            )
            for scope in running_scopes
        ]
        total_shifts = sum([len(self.shifts[scope]) for scope in running_scopes])
        log.info("------------------------------------")
        log.info("Configuration Report")
        log.info("------------------------------------")
        log.info("  Sample: {}".format(self.sample))
        log.info("  Era: {}".format(self.era))
        log.info("  Scopes: {}".format(self.scopes))
        log.info("  Total number of producers: {}".format(total_producers))
        for scope in running_scopes:
            log.info("       {}: {}".format(scope, len(self.producers[scope])))
        log.info("  Total number of quantities: {}".format(sum(total_quantities)))
        for i, scope in enumerate(running_scopes):
            log.info("       {}: {}".format(scope, total_quantities[i]))
        log.info("  Total number of shifts: {}".format(total_shifts))
        for scope in running_scopes:
            log.info("       {}: {}".format(scope, len(self.shifts[scope])))
        log.info("------------------------------------")

    def __str__(self) -> str:
        """
        String representation of the configuration.
        """
        returnstr = "Configuration: \n"
        returnstr += "  Era: {}\n".format(self.era)
        returnstr += "  Sample: {}\n".format(self.sample)
        returnstr += "  Scopes: {}\n".format(self.scopes)
        returnstr += "  Shifts: {}\n".format(self.shifts)
        returnstr += "  Rules:  {}\n".format(self.rules)
        returnstr += "  Outputs:\n"
        for scope in self.scopes:
            returnstr += "    {}: {}\n".format(scope, self.outputs[scope])
        returnstr += "  Producers:\n"
        for scope in self.scopes:
            returnstr += "    {}: {}\n".format(scope, self.producers[scope])
        returnstr += "  Config Parameters:\n"
        for scope in self.scopes:
            returnstr += "    {}: {}\n".format(scope, self.config_parameters[scope])
        return returnstr

    def expanded_configuration(self) -> Configuration:
        """Function used to generate an expanded version of the configuration, where all shifts are applied.
        This expanded configuration is used by the code generator to generate the C++ code.

        Returns:
            Configuration: Expanded configuration
        """
        expanded_configuration = {}
        for scope in self.scopes:
            expanded_configuration[scope] = {}
            expanded_configuration[scope]["nominal"] = self.config_parameters[scope]
            if len(self.shifts[scope]) > 0:
                for shift in self.shifts[scope]:
                    log.debug("Adding shift {} in scope {}".format(shift, scope))
                    log.debug("  {}".format(self.shifts[scope][shift]))
                    try:
                        expanded_configuration[scope][shift] = (
                            self.config_parameters[scope] | self.shifts[scope][shift]
                        )
                    except KeyError:
                        expanded_configuration[scope][shift] = {}
                        expanded_configuration[scope][shift] = (
                            self.config_parameters[scope] | self.shifts[scope][shift]
                        )
        self.config_parameters = expanded_configuration
        return self
