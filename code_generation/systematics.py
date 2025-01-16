from __future__ import annotations  # needed for type annotations in > python 3.7

import logging
from typing import Any, Dict, List, Set, Tuple, Union

from code_generation.modifiers import EraModifier, SampleModifier
from code_generation.producer import (
    Producer,
    ProducerGroup,
    TProducerInput,
    TProducerStore,
)
from code_generation.quantity import NanoAODQuantity
from code_generation.helpers import is_empty

log = logging.getLogger(__name__)

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


class SystematicShift(object):
    """Class containing systematic shifts. A systematic shift is a variation of a set of configuration parameters.

    A dummy shift looks like this::

        SystematicShift(
                name="shiftname",
                shift_config={
                    ("scope1", "scope2"): {
                        "paramter1": True,
                        "paramter2": 42.0,
                    },
                },
                producers={
                    "scope1": modified_producer1,
                    "scope2": modified_producer2,
                },
                ignore_producers={
                    "scope1": ignored_producer1,
                    "scope2": ignored_producer2,
                },
            )

    Args:
        name (str): Name of the systematic shift.

        shift_config: Dictionary containing the configuration parameters. The dictionary keys have to be strings, or tuples, in case the same configuration change is used for multiple scopes. The dictionary values have to be dictionaries containing the changed configuration parameters.

        producers: Dictionary containing the producers that are affected by the systematic shift. The dictionary keys have to be strings, or tuples, in case the same producer is affected by the systematic shift in multiple scopes. The dictionary values have to be the modified producers.

        ignore_producers: Dictionary containing the producers that are not affected by the systematic shift. The dictionary keys have to be strings, or tuples, in case the same producer is not affected by the systematic shift in multiple scopes. The dictionary values have to be the ignored producers.

        scopes (List[str], optional): List of scopes that are affected by the systematic shift. If not given, all scopes are affected.

    """

    def __init__(
        self,
        name: str,
        shift_config: Dict[Union[str, Tuple[str, ...]], TConfiguration],
        producers: Dict[Union[str, Tuple[str, ...]], TProducerInput],
        scopes: Union[List[str], None] = None,
        ignore_producers: Dict[Union[str, Tuple[str, ...]], TProducerInput] = {},
    ):
        self.shiftname: str = "__" + name
        self.input_producers: Dict[
            str, Union[Producer, ProducerGroup, List[Producer | ProducerGroup]]
        ] = self.expand_producer_dict_keys(producers)
        self.input_ignore_producers: Dict[
            str, Union[Producer, ProducerGroup, List[Producer | ProducerGroup]]
        ] = self.expand_producer_dict_keys(ignore_producers)
        self.producers: TProducerStore = {}
        self.ignore_producers: TProducerStore = {}
        self.shift_config: Dict[
            str, TConfiguration
        ] = self.expand_configuration_dict_keys(shift_config)
        self.scopes: Set[str] = self.determine_scopes(scopes)
        self.validate()

    def expand_producer_dict_keys(
        self,
        dict_to_expand: Dict[Union[str, Tuple[str, ...]], TProducerInput],
    ) -> Dict[str, Union[Producer, ProducerGroup, List[Producer | ProducerGroup]]]:
        """
        Function used to expand dictionaries. If the key is a string,
        it is returned as is. If the key is a tuple, the tuple is expanded
        and the resulting dict with one key per tuple entry is returned.

        Args:
            dict_to_expand (dict): Dictionary to expand.

        Returns:
            dict: Expanded dictionary.
        TProducerListStore = Dict[str, Dict[str, List[Union[Producer, ProducerGroup]]]]
        """
        exanded_dict: Dict[
            str, Union[Producer, ProducerGroup, List[Producer | ProducerGroup]]
        ] = {}
        for key in dict_to_expand.keys():
            # if not (isinstance(key, str) or isinstance(key, tuple)):
            #     errormsg = "Producer dict key {} must be a string or tuple for shift {}".format(
            #         key, self.shiftname
            #     )
            #     raise ValueError(errormsg)
            if isinstance(key, str):
                scopes = (key,)
            else:
                scopes = key
            for scope in scopes:
                exanded_dict[scope] = dict_to_expand[key]
        return exanded_dict

    def expand_configuration_dict_keys(
        self,
        dict_to_expand: Dict[Union[str, Tuple[str, ...]], TConfiguration],
    ) -> Dict[str, TConfiguration]:
        """
        Function used to expand the configuration dictionary.
        If the key is a string, it is returned as is. If the key is a tuple,
        the tuple is expanded and the resulting dict with one key per tuple entry is returned.

        Args:
            dict_to_expand (dict): Dictionary to expand.

        Returns:
            dict: Expanded dictionary.
        """
        exanded_dict: Dict[str, TConfiguration] = {}
        for key in dict_to_expand:
            # if not (isinstance(key, str) or isinstance(key, tuple)):
            #     errormsg = "Configuration dict key {} must be a string or tuple for shift {}".format(
            #         key, self.shiftname
            #     )
            #     raise ValueError(errormsg)
            if isinstance(key, str):
                scopes = (key,)
            else:
                scopes = key
            for scope in scopes:
                exanded_dict[scope] = dict_to_expand[key]
        return exanded_dict

    def __str__(self) -> str:
        returnstr = "SystematicShift: {}\n".format(self.shiftname)
        returnstr += "  Scopes: {}\n".format(self.scopes)
        returnstr += "  Producers: {}\n".format(self.producers)
        returnstr += "  Ignore producers: {}\n".format(self.ignore_producers)
        returnstr += "  Shift config: {}\n".format(self.shift_config)
        returnstr += "  ------------ \n"
        return returnstr

    def __repr__(self) -> str:
        returnstr = "SystematicShift: {}\n".format(self.shiftname)
        returnstr += "  Scopes: {}\n".format(self.scopes)
        returnstr += "  Producers: {}\n".format(self.producers)
        returnstr += "  Ignore producers: {}\n".format(self.ignore_producers)
        returnstr += "  Shift config: {}\n".format(self.shift_config)
        returnstr += "  ------------ \n"
        return returnstr

    def determine_scopes(self, scopes: Union[List[str], str, None]) -> Set[str]:
        """
        Function used to determine the scopes that are affected by the systematic shift.
        If no scope is specified by the user, the scopes are determined from the shift_config,
        producers and ignore_producers.

        Args:
            scopes (set): Set of scopes that are affected by the systematic shift.

        Returns:
            set: Set of scopes that are affected by the systematic shift.
        """
        if is_empty(scopes):
            scope_set: Set[str] = (
                set(self.shift_config.keys())
                | set(self.input_producers.keys())
                | set(self.input_ignore_producers.keys())
            )
        elif isinstance(scopes, str):
            scope_set = set(scopes)
        else:
            scope_set = set(scopes)
        return scope_set

    def normalize_inputs(
        self,
        input_producers: Dict[
            str, Union[Producer, ProducerGroup, List[Producer | ProducerGroup]]
        ],
    ) -> TProducerStore:
        """

        Helper function used in the validation of the systematic shift. If makes sure,
        that the provided dict of producers is in a defined format, meaning a dict with one key per
        scope and a list of producers as value.

        Args:
            producers: Producers can be a single Producer object, a list of Producer objects or a dict
                with one key per scope and a list of producers as value.

        Returns:
            dict: Dictionary containing the producers with one key per scope and a list of producers as value.

        """
        temp_dict: Dict[str, List[Producer | ProducerGroup]] = {}
        resolved_dict: TProducerStore = {}
        scopes = input_producers.keys()
        for scope in scopes:
            temp_input = input_producers[scope]
            if isinstance(temp_input, Producer) or isinstance(
                temp_input, ProducerGroup
            ):
                temp_dict[scope] = [temp_input]
            else:
                temp_dict[scope] = temp_input
            resolved_dict[scope] = []
            try:
                temp = temp_dict[scope]
                if isinstance(temp, Producer):
                    resolved_dict[scope] = [temp]
                else:
                    resolved_dict[scope] = temp
            except KeyError:
                resolved_dict[scope] = []
                log.debug(
                    "Setting empty Producer set for shift {} in scope {}".format(
                        self.shiftname, scope
                    )
                )
        return resolved_dict

    def validate(self) -> None:
        """
        Function used to validate the systematic shift. The provided producers
        and ignore_producers are converted into a defined format. If no
        configuration is provided for a scope, an empty dict is used.

        Args:
            None

        Returns:
            None
        """
        log.debug("Validating systematic shift {}".format(self.shiftname))
        log.debug(" input Producers: {}".format(self.input_producers))
        self.producers = self.normalize_inputs(self.input_producers)
        self.ignore_producers = self.normalize_inputs(self.input_ignore_producers)
        unconfigured_scopes = set(self.scopes) - set(self.shift_config.keys())
        # for uncongured scopes, set an empty config
        if len(unconfigured_scopes) > 0:
            for scope in unconfigured_scopes:
                self.shift_config[scope] = {}

    def add_producer(
        self,
        producer: Producer,
        input_scope: Union[List[str], None] = None,
    ) -> None:
        """
        Function used to add a producer to the list of producers affected by the
        systematic shift. After adding the producer, the shift is validated.

        Args:
            producer: Producer to add.
            scope: Scope to which the producer should be added. If no sopce is
                provided, the producer is added to all scopes of the systematic.

        Returns:
            None
        """
        scopes: Set[str] = self.scopes
        if isinstance(input_scope, str):
            scopes = set([input_scope])
        for scope in scopes:
            self.producers[scope].append(producer)

    def add_ignore_producer(
        self,
        producer: Producer,
        scopes: Union[str, None, Set[str], List[str]] = None,
    ) -> None:
        """
        Function used to add an ignored producer to the list of ignored producers,
        which are untouched by the systematic shift. After adding the ignored producer, the shift is validated.

        Args:
            producer: Producer to ignore.
            scope: Scope to which the ignored producer should be added. If no sopce is
                provided, the ignored producer is added to all scopes of the systematic.

        Returns:
            None
        """
        if is_empty(scopes):
            scopes = self.scopes
        if isinstance(scopes, str):
            scopes = set(scopes)
        for scope in scopes:
            self.ignore_producers[scope].append(producer)

    def add_scope(self, scope: str) -> None:
        """
        Function used to add a scope to the list of scopes affected
        by the systematic shift. After adding the scope, the shift is validated.

        Args:
            scope: Scope to add.

        Returns:
            None
        """
        self.scopes.add(scope)
        self.input_producers[scope] = []
        self.input_ignore_producers[scope] = []
        self.validate()

    def add_config(self, config: TConfiguration, scope: str) -> None:
        """
        Function used to add a configuration to the list of configurations for a scope. After adding the configuration, the shift is validated.

        Args:
            config: Configuration to add.
            scope: Scope to which the configuration should be added.

        Returns:
            None
        """
        self.shift_config[scope] = config
        self.validate()

    def get_scopes(self) -> Set[str]:
        """
        This function returns a set of all scopes that are affected by the systematic shift.

        Args:
            None

        Returns:
            set: Set of scopes that are affected by the systematic shift.
        """
        return self.scopes

    def get_shift_config(self, scope: str) -> TConfiguration:
        """
        This function returns the configuration for a given scope.

        Args:
            scope: Scope for which the configuration should be returned.

        Returns:
            dict: Configuration for the given scope.
        """
        return self.shift_config[scope]

    def apply(self, scope: str) -> None:
        """
        Function used to apply the systematic shift to the given producers. For the given scope, all producers aer shifted using producer.shift, while, for all ignored producers, the producer.ignore_shift function is called. If the scope is not defined in the shift, no shift is applied.

        Args:
            scope: Scope for which the shift should be applied.

        Returns:
            None
        """
        log.debug("Applying systematic shift \n{}".format(self))
        if scope in self.scopes:
            if scope in self.ignore_producers.keys():
                for producer in self.ignore_producers[scope]:
                    producer.ignore_shift(self.shiftname, scope)
            if scope in self.producers.keys():
                for producer in self.producers[scope]:
                    producer.shift(self.shiftname, scope)


class SystematicShiftByQuantity(SystematicShift):
    """
    Class used to define a systematic shift that is defined by a quantity.
    """

    def __init__(
        self,
        name: str,
        quantity_change: Dict[NanoAODQuantity, Union[str, NanoAODQuantity]],
        scopes: Union[List[str], None] = None,
    ):
        """
        Constructor for the SystematicShiftByQuantity class.

        Args:
            name: Name of the systematic shift.
            quantity_change: Dictionary of quantities that should be changed.
            scopes: List of scopes that are affected by the systematic shift.
        """
        super().__init__(name, {}, {}, scopes, {})
        self.quantity_change: Dict[
            NanoAODQuantity, Union[str, NanoAODQuantity]
        ] = quantity_change
        self.quantities: Set[NanoAODQuantity] = set(quantity_change.keys())

    def apply(self, scope: str) -> None:
        """
        Function used to apply the systematic shift to the given producers. For the given scope, all producers aer shifted using producer.shift, while, for all ignored producers, the producer.ignore_shift function is called. If the scope is not defined in the shift, no shift is applied.

        Args:
            scope: Scope for which the shift should be applied.

        Returns:
            None
        """
        for quantity in self.quantities:
            quantity.register_external_shift(
                shift_name=self.shiftname, external_name=self.quantity_change[quantity]
            )
