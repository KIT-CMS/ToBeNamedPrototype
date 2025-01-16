from __future__ import annotations  # needed for type annotations in > python 3.7

import logging
import re
from typing import Any, Dict, List, Set, Union

from code_generation.exceptions import (
    InvalidProducerConfigurationError,
    ConfigurationError,
)
from code_generation.helpers import is_empty

import code_generation.quantity as q

log = logging.getLogger(__name__)


class SafeDict(Dict[Any, Any]):
    def __missing__(self, key: Any) -> Any:
        return "{" + key + "}"


class Producer:
    def __init__(
        self,
        name: str,
        call: str,
        input: Union[List[q.Quantity], Dict[str, List[q.Quantity]]],
        output: Union[List[q.Quantity], None],
        scopes: List[str],
    ):
        """
        A Producer is a class that holds all information about a producer. Input quantities are

        Args:
            name: Name of the producer
            call: The call of the producer. This is the C++ function call of the producer
            input: A list of input quantities or a dict with scope specific input quantities
            output: A list of output quantities
            scopes: A list of scopes in which the producer is used

        """
        log.debug("Setting up a new producer {}".format(name))

        # sanity checks
        if not isinstance(input, list) and not isinstance(input, dict):
            log.error(
                "Exception (%s): Argument 'input' must be a list or a dict!" % name
            )
            raise Exception
        if not isinstance(output, list) and output is not None:
            log.error(
                "Exception (%s): Argument 'output' must be a list or None!" % name
            )
            raise Exception
        self.name: str = name
        self.call: str = call
        self.output: Union[List[q.Quantity], None] = output
        self.scopes = scopes
        self.parameters: Dict[str, Set[str]] = self.extract_parameters()
        # if input not given as dict and therfore not scope specific transform into dict with all scopes
        if not isinstance(input, dict):
            inputdict = {}
            for scope in self.scopes:
                inputdict[scope] = input.copy() if isinstance(input, list) else input
        else:
            inputdict = input
        self.input: Dict[str, List[q.Quantity]] = inputdict
        # keep track of variable dependencies
        if not is_empty(self.output):
            for scope in self.scopes:
                for input_quantity in self.input[scope]:
                    for output_quantity in self.output:
                        input_quantity.adopt(output_quantity, scope)
        log.debug("-----------------------------------------")
        log.debug("| Producer: {}".format(self.name))
        log.debug("| Call: {}".format(self.call))
        for scope in self.scopes:
            if is_empty(self.input[scope]):
                log.debug("| Inputs ({}): None".format(scope))
            else:
                log.debug(
                    "| Inputs ({}): {}".format(
                        scope, [input.name for input in self.input[scope]]
                    )
                )
        if is_empty(self.output):
            log.debug("| Output: None")
        else:
            log.debug("| Outputs: {}".format([output.name for output in self.output]))
        log.debug("| scopes: {}".format(self.scopes))
        log.debug("-----------------------------------------")

    def __str__(self) -> str:
        return "Producer: {}".format(self.name)

    def __repr__(self) -> str:
        return "Producer: {}".format(self.name)

    def extract_parameters(self) -> Dict[str, Set[str]]:
        """
        Function used to extract all parameters from a producer call. Parameters are enclosed in curly brackets.
        Reserved parameters are:
        - {output} : will be replaced by the output quantities
        - {input} : will be replaced by the input quantities
        - {output_vec} : will be replaced by the output quantities in vector form
        - {input_vec} : will be replaced by the input quantities in vector form
        - {df} : will be replaced by the dataframe name
        - {vec_open} : will be replaced by the opening bracket of the vector
        - {vec_close} : will be replaced by the closing bracket of the vector
        """
        regex = "\{([A-Za-z0-9_]+)\}"
        reserved_parameters = [
            "output",
            "input",
            "output_vec",
            "input_vec",
            "df",
            "vec_open",
            "vec_close",
        ]
        parameters = {}
        for scope in self.scopes:
            parameters[scope] = set(
                [
                    x
                    for x in re.findall(regex, self.call)
                    if x not in reserved_parameters
                ]
            )
        return parameters

    # Check if a output_quantity is already used as an output by
    # another producer within the same scope.
    # If this occurs, a Exception is thrown, since this is not possible with dataframes
    def reserve_output(self, scope: str) -> None:
        """
        Check if a output_quantity is already used as an output by another producer within the same scope.
        This is an internal function and should not be called by the user.

        Args:
            scope: The scope in which the output is reserved

        """

        if not is_empty(self.output):
            for output_quantity in self.output:
                output_quantity.reserve_scope(scope)

    def shift(self, name: str, scope: str = "global") -> None:
        """Add a shift to the producer. This is an internal function and should not be called by the user.

        Args:
            name (str): Name of the shift
            scope (str, optional): Name of the scope where the shift is to be applied. Defaults to "global".

        Raises:
            Exception: If the producer does not have any output, or if the producer does not exist in the given scope, an exception is thrown
        """
        if scope not in self.scopes:
            log.error(
                "Trying to add shift %s to producer %s in scope %s, but producer does not exist in given scope!"
                % (name, self.name, scope)
            )
            raise Exception
        if is_empty(self.output):
            log.error(
                "Exception (%s): output %s cannot be shifted ! How did you end up here ?"
                % (name, self.output)
            )
            raise Exception
        for entry in self.output:
            entry.shift(name, scope)

    def shift_inputs(self, name: str, scope: str, inputs_to_shift: List[str]) -> None:
        """
        Shift all inputs of a producer. This is an internal function and should not be called by the user.

        Args:
            name: Name of the shift
            scope: Name of the scope where the shift is to be applied
        """
        if scope not in self.scopes:
            log.error(
                "Trying to add shift %s to producer %s in scope %s, but producer does not exist in given scope!"
                % (name, self.name, scope)
            )
            raise Exception
        for entry in self.input[scope]:
            if entry in inputs_to_shift:
                log.debug(f"Shifting {entry.name}")
                entry.shift(name, scope)

    def ignore_shift(self, name: str, scope: str = "global") -> None:
        """Ingore a given shift for a producer. This in an internal function and should not be called by the user.

        Args:
            name (str): Name of the shift
            scope (str, optional): Name of the scope where the shift is to be ignored. Defaults to "global".

        Raises:
            Exception: If the producer does not have any output, or if the producer does not exist in the given scope, an exception is thrown
        """
        if scope not in self.scopes:
            log.error(
                "Trying to add shift %s to producer %s in scope %s, but producer does not exist in given scope!"
                % (name, self.name, scope)
            )
            raise Exception
        if is_empty(self.output):
            log.error(
                "Exception (%s): output %s cannot be shifted ! How did you end up here ?"
                % (name, self.output)
            )
            raise Exception
        for entry in self.output:
            entry.ignore_shift(name, scope)

    def writecall(
        self, config: Dict[str, Dict[str, Any]], scope: str, shift: str = "nominal"
    ) -> str:
        """Function to generate the nessessary C++ calls for the code generation

        Args:
            config (Dict[str, Dict[str, Any]]): Configuration dict containing the parameter and input / output information for the producer
            scope (str): The scope in which the producer is to be called
            shift (str, optional): The shift, for which the function call should be generated. Defaults to "nominal".

        Raises:
            Exception: Raises an expection, the the requested shift is not available in the configuration

        Returns:
            str: The generated C++ call
        """
        if is_empty(self.output):
            config[shift]["output"] = ""
            config[shift]["output_vec"] = ""
        else:
            # log.warning(f"Available shifts: {config.keys()}")
            # log.warning(f"Configuration: {config[shift]}")
            # log.warning("Writing call for {}".format(self.name))
            # log.warning("output: {}".format(self.output))
            # log.warning("name: {}".format(self.name))
            # log.warning("shift: {}".format(shift))
            # log.warning("Scopes: {}".format(config[shift].keys()))
            config[shift]["output"] = (
                '"' + '","'.join([x.get_leaf(shift, scope) for x in self.output]) + '"'
            )
            config[shift]["output_vec"] = (
                '{"'
                + '","'.join([x.get_leaf(shift, scope) for x in self.output])
                + '"}'
            )
        config[shift]["input"] = (
            '"'
            + '", "'.join([x.get_leaf(shift, scope) for x in self.input[scope]])
            + '"'
        )
        config[shift]["input_vec"] = (
            '{"'
            + '","'.join([x.get_leaf(shift, scope) for x in self.input[scope]])
            + '"}'
        )
        config[shift]["df"] = "{df}"
        config[shift]["vec_open"] = "{vec_open}"
        config[shift]["vec_close"] = "{vec_close}"
        # if a bool is used in the python configuration, convert it to a c++ bool value
        # True -> true, False -> false
        for para in config[shift]:
            if isinstance(config[shift][para], bool):
                if config[shift][para]:
                    log.debug("Found a boolean True ! - converting to C++ syntax")
                    config[shift][para] = "true"
                else:
                    log.debug("Found a boolean False ! - converting to C++ syntax")
                    config[shift][para] = "false"
        try:
            return self.call.format(
                **config[shift]
            )  # use format (not format_map here) such that missing config entries cause an error
        except KeyError as e:
            log.error(
                "Error in {} Producer, key {} is not found in configuration".format(
                    self.name, e
                )
            )
            log.error(config[shift])
            log.error("Call: {}".format(self.call))
            raise Exception

    def writecalls(
        self, config: Dict[str, Dict[str, Dict[str, str]]], scope: str
    ) -> List[str]:
        """Function to generate calls for all shifts of a producer, wraping the writecall function

        Args:
            Configuration dict containing the parameter and input / output information for the producer
            scope (str): The scope in which the producer is to be called

        Raises:
            Exception: Raises an Exception, of the requested scope is not valid for the producer

        Returns:
            List[str]: Returns a list of C++ calls for all shifts of the producer
        """
        if scope not in self.scopes:
            log.error(
                "Exception ({}): Tried to use producer in scope {}, which the producer is not forseen for!".format(
                    self.name, scope
                )
            )
            raise Exception
        calls = [self.writecall(config, scope)]
        if not is_empty(self.output):
            list_of_shifts = self.output[0].get_shifts(
                scope
            )  # all entries must have same shifts
            for shift in list_of_shifts:
                calls.append(self.writecall(config, scope, shift))
        return calls

    def get_inputs(self, scope: str) -> List[q.Quantity]:
        """Get a list of all inputs of the producer in a given scope

        Args:
            scope (str): The scope in which the inputs are requested

        Raises:
            Exception: Raises an Exception, of the requested scope is not valid for the producer

        Returns:
            List[str]: Returns a list of Quantity objects, which are the inputs of the producer
        """
        if scope not in self.scopes:
            log.error(
                "Exception ({}): Tried to get producer inputs in scope {}, which the producer is not forseen for!".format(
                    self.name, scope
                )
            )
            raise Exception
        return self.input[scope]

    def get_outputs(self, scope: str) -> List[Union[q.QuantityGroup, q.Quantity]]:
        """Get a list of all outputs of the producer in a given scope

        Args:
            scope (str): The scope in which the outputs are requested

        Raises:
            Exception: Raises an Exception, of the requested scope is not valid for the producer

        Returns:
            List[Union[q.QuantityGroup, q.Quantity]]: Returns a list of Quantity objects, which are the outputs of the producer
        """
        if scope not in self.scopes:
            log.error(
                "Exception ({}): Tried to get producer outputs in scope {}, which the producer is not forseen for!".format(
                    self.name, scope
                )
            )
            raise Exception
        if is_empty(self.output):
            return []
        else:
            return self.output


class VectorProducer(Producer):
    def __init__(
        self,
        name: str,
        call: str,
        input: Union[List[q.Quantity], Dict[str, List[q.Quantity]]],
        output: Union[List[q.Quantity], None],
        scopes: List[str],
        vec_configs: List[str],
    ):
        """A Vector Producer is a Producer which can be configured to produce multiple calls and outputs at once, deprecated in favor of the ExtendedVectorProducer

        Args:
            name (str): Name of the producer
            call (str): The call to be made in C++, with placeholders for the parameters
            input (Union[List[q.Quantity], Dict[str, List[q.Quantity]]]): The inputs of the producer, either a list of Quantity objects, or a dict with the scope as key and a list of Quantity objects as value
            output (Union[List[q.Quantity], None]): The outputs of the producer, either a list of Quantity objects, or None if the producer does not produce any output
            scopes (List[str]): The scopes in which the producer is to be called
            vec_configs (List[str]): A list of strings, which are the names of the parameters to be varied in the vectorized call
        """
        self.name = name
        super().__init__(name, call, input, output, scopes)
        self.vec_configs = vec_configs

    def __str__(self) -> str:
        return "VectorProducer: {}".format(self.name)

    def __repr__(self) -> str:
        return "VectorProducer: {}".format(self.name)

    def writecalls(
        self, config: Dict[str, Dict[str, Dict[str, str]]], scope: str
    ) -> List[str]:
        """Function to generate calls for all shifts of a producer, wraping the writecall function

        Args:
            Configuration dict containing the parameter and input / output information for the producer
            scope (str): The scope in which the producer is to be called

        Raises:
            Exception: Raises an Exception, of the requested scope is not valid for the producer

        Returns:
            List[str]: Returns a list of C++ calls for all shifts of the producer
        """
        basecall = self.call
        calls: List[str] = []
        shifts = ["nominal"]
        if self.output:
            shifts.extend(self.output[0].get_shifts(scope))
        for shift in shifts:
            # check that all config lists (and output if applicable) have same length
            log.debug("self.vec_configs[0]: {}".format(self.vec_configs[0]))
            log.debug("len(self.vec_configs): {}".format(len(self.vec_configs)))
            log.debug("available shifts: {}".format(config.keys()))
            n_versions = len(config[shift][self.vec_configs[0]])
            for key in self.vec_configs:
                if n_versions != len(config[shift][key]):
                    log.error(
                        "Following lists in config must have same length: %s, %s"
                        % (self.vec_configs[0], key)
                    )
                    raise Exception
            if not is_empty(self.output) and len(self.output) != n_versions:
                log.error(
                    "{} expects either no output or same amount as entries in config lists !".format(
                        self
                    )
                )
                log.error("Number of expected outputs: {}".format(n_versions))
                log.error("List of outputs: {}".format(self.output))
                raise Exception
            for i in range(n_versions):
                helper_dict: Dict[Any, Any] = {}
                for key in self.vec_configs:
                    helper_dict[key] = config[shift][key][i]
                if not is_empty(self.output):
                    helper_dict["output"] = (
                        '"' + self.output[i].get_leaf(shift, scope) + '"'
                    )
                    helper_dict["output_vec"] = (
                        '{"' + self.output[i].get_leaf(shift, scope) + '"}'
                    )
                self.call = basecall.format_map(SafeDict(helper_dict))
                calls.append(self.writecall(config, scope, shift))
        self.call = basecall
        return calls


class ExtendedVectorProducer(Producer):
    def __init__(
        self,
        name: str,
        call: str,
        input: Union[List[q.Quantity], Dict[str, List[q.Quantity]]],
        output: str,
        scope: Union[List[str], str],
        vec_config: str,
    ):
        """A ExtendedVectorProducer is a Producer which can be configured to produce multiple calls and outputs at once

        Args:
            name (str): Name of the producer
            call (str): The call to be made in C++, with placeholders for the parameters
            input (Union[List[q.Quantity], Dict[str, List[q.Quantity]]]): The inputs of the producer, either a list of Quantity objects, or a dict with the scope as key and a list of Quantity objects as value
            output (Union[List[q.Quantity], None]): The outputs of the producer, either a list of Quantity objects, or None if the producer does not produce any output
            scopes (List[str]): The scopes in which the producer is to be called
            vec_configs (List[str]): The key of the vec config in the config dict
        """
        # we create a Quantity Group, which is updated during the writecalls() step
        self.outputname = output
        self.vec_config = vec_config
        if not isinstance(scope, list):
            scope = [scope]
        quantity_group = q.QuantityGroup(name)
        # set the vec config key of the quantity group
        quantity_group.set_vec_config(vec_config)
        super().__init__(name, call, input, [quantity_group], scope)
        if is_empty(self.output):
            raise InvalidProducerConfigurationError(self.name)
        # add the vec config to the parameters of the producer
        for scope in self.scopes:
            self.parameters[scope].add(self.vec_config)

    def __str__(self) -> str:
        return "ExtendedVectorProducer: {}".format(self.name)

    def __repr__(self) -> str:
        return "ExtendedVectorProducer: {}".format(self.name)

    @property
    def output_group(self) -> q.QuantityGroup:
        if is_empty(self.output):
            raise Exception("ExtendedVectorProducer has no output!")
        if not isinstance(self.output[0], q.QuantityGroup):
            log.error("ExtendedVectorProducer expects a QuantityGroup as output!")
            raise Exception
        return self.output[0]

    def writecalls(
        self, config: Dict[str, Dict[str, Dict[str, Any]]], scope: str
    ) -> List[str]:
        """Function to generate all calls for all shifts of the ExtendedVectorProducer, wraping the writecall function

        Args:
            Configuration dict containing the parameter and input / output information for the ExtendedVectorProducer
            scope (str): The scope in which the ExtendedVectorProducer is to be called

        Raises:
            Exception: Raises an Exception, of the requested scope is not valid for the ExtendedVectorProducer

        Returns:
            List[str]: Returns a list of C++ calls for all shifts of the ExtendedVectorProducer
        """
        n_versions = len(config["nominal"][self.vec_config])
        log.debug("Number of extended producers to be created {}".format(n_versions))
        if is_empty(self.output):
            raise InvalidProducerConfigurationError(self.name)
        if not isinstance(self.output[0], q.QuantityGroup):
            log.error("ExtendedVectorProducer expects a QuantityGroup as output!")
            raise Exception
        for i in range(n_versions):
            self.output[0].add(config["nominal"][self.vec_config][i][self.outputname])
        basecall = self.call
        calls: List[str] = []
        shifts = ["nominal"]
        shifts.extend(self.output[0].get_shifts(scope))
        for shift in shifts:
            for i in range(n_versions):
                # the information for the producer is directly read from the configuration
                helper_dict = config[shift][self.vec_config][i]
                helper_dict["output"] = (
                    '"' + self.output[0].quantities[i].get_leaf(shift, scope) + '"'
                )
                helper_dict["output_vec"] = (
                    '{"' + self.output[0].quantities[i].get_leaf(shift, scope) + '"}'
                )
                self.call = basecall.format_map(SafeDict(helper_dict))
                calls.append(self.writecall(config, scope, shift))
        self.call = basecall
        return calls


class BaseFilter(Producer):
    def __init__(
        self,
        name: str,
        call: str,
        input: Union[List[q.Quantity], Dict[str, List[q.Quantity]]],
        scopes: List[str],
    ):
        """A BaseFilter is a Producer which does not produce any output, but is used to filter events.
        This class should not be used by the user, use the Filter class instead.

        Args:
            name (str): name of the filter
            call (str): The call to be made in C++, with placeholders for the parameters
            input (Union[List[q.Quantity], Dict[str, List[q.Quantity]]]): The inputs of the filter,
            either a list of Quantity objects, or a dict with the scope as key and a list of Quantity objects as value
            scopes (List[str]): The scopes in which the filter is to be called
        """
        super().__init__(name, call, input, None, scopes)

    def __str__(self) -> str:
        return "BaseFilter: {}".format(self.name)

    def __repr__(self) -> str:
        return "BaseFilter: {}".format(self.name)

    def writecall(
        self, config: Dict[str, Dict[str, Dict[str, str]]], scope: str, shift: str = ""
    ) -> str:
        """
        Do not use, filters do not support method writecall!
        """
        log.critical("{}: Filters do not support method writecall!".format(self.name))
        raise Exception

    def writecalls(
        self, config: Dict[str, Dict[str, Dict[str, str]]], scope: str
    ) -> List[str]:
        """The writecalls function of a BaseFilter is used to generate the C++ calls for the filter

        Args:
            config (Dict[str, Dict[str, Dict[str, str]]]): The configuration dict containing the
            parameters and input / output information for the filter
            scope (str): The scope in which the filter is to be called

        Raises:
            Exception: Raises an Exception, if the requested scope is not valid for the filter

        Returns:
            List[str]: Returns a list of C++ calls for the filter
        """
        inputs: List[str] = []
        for quantity in self.input[scope]:
            inputs.extend(quantity.get_leaves_of_scope(scope))
        config["nominal"]["input"] = '"' + '", "'.join(inputs) + '"'
        config["nominal"]["input_vec"] = '{"' + '","'.join(inputs) + '"}'
        config["nominal"]["df"] = "{df}"
        try:
            return [
                self.call.format(**config["nominal"])
            ]  # use format (not format_map here) such that missing config entries cause an error
        except KeyError as e:
            log.error(
                "Error in {} Basefilter, key {} is not found in configuration".format(
                    self.name, e
                )
            )
            log.error("Call: {}".format(self.call))
            raise Exception


class ProducerGroup:
    PG_count = 1  # counter for internal quantities used by ProducerGroups

    def __init__(
        self,
        name: str,
        call: Union[str, None],
        input: Union[List[q.Quantity], Dict[str, List[q.Quantity]], None],
        output: Union[List[q.Quantity], None],
        scopes: List[str],
        subproducers: Union[
            List[Producer | ProducerGroup],
            Dict[str, List[Producer | ProducerGroup]],
        ],
    ):
        """A ProducerGroup can be used to group multiple producers. This is useful to keep the configuration simpler and to ensure that the producers are called in the correct order. ProducerGroups can be nested.

        Args:
            name (str): Name of the ProducerGroup
            call (Union[str, None]): Typically, this is None
            input (Union[List[q.Quantity], Dict[str, List[q.Quantity]], None]): The inputs of the ProducerGroup
            output (Union[List[q.Quantity], None]): Output quantities of the Producer Group
            scopes (List[str]): The scopes in which the ProducerGroup is to be called
            subproducers (Union[ List[Producer  |  ProducerGroup], Dict[str, List[Producer  |  ProducerGroup]], ]): The subproducers of the ProducerGroup, either a list of Producer or ProducerGroup objects, or a dict with the scope as key and a list of Producer or ProducerGroup objects as value
        """
        self.name = name
        self.call = call
        self.output = output
        self.scopes = scopes
        self.producers: Dict[str, List[Producer | ProducerGroup]] = {}
        # if subproducers are given as dict and therefore scope specific transform into dict with all scopes
        if not isinstance(subproducers, dict):
            log.debug("Converting subproducer list to dictionary")
            for scope in self.scopes:
                self.producers[scope] = [producer for producer in subproducers]
        else:
            self.producers = subproducers
        # do a consistency check for the scopes
        self.check_producer_scopes()
        # if input not given as dict and therefore not scope specific transform into dict with all scopes
        if not isinstance(input, dict):
            inputdict = {}
            for scope in self.scopes:
                inputdict[scope] = input.copy() if isinstance(input, list) else input
            self.input = inputdict
        else:
            self.input = dict(input)
        # If call is provided, this is supposed to consume output of subproducers. Creating these internal products below:
        if not is_empty(self.call):
            log.debug("Constructing {}".format(self.name))
            log.debug(" --> Scopes: {}".format(self.scopes))
            for scope in self.scopes:
                log.debug("Adding for scope: {}".format(scope))
                log.debug(" --> Producers {}".format(self.producers[scope]))
                for subproducer in self.producers[scope]:
                    # skip producers without output
                    log.debug(
                        "Adding {} / {} : output: {}".format(
                            subproducer, scope, subproducer.output
                        )
                    )
                    if subproducer.output is None:
                        continue
                    # if the subproducer does not have an output quantity, we assign an internally tracked quantity
                    if subproducer.output == []:
                        # create quantities that are produced by subproducers and then collected by the final call of the producer group
                        subproducer.output = [
                            q.Quantity(
                                "PG_internal_quantity_%i" % self.__class__.PG_count
                            )
                        ]  # quantities of vector producers will be duplicated later on when config is known
                        log.debug(
                            "Added new internal quantitiy: {}".format(
                                subproducer.output
                            )
                        )
                        self.__class__.PG_count += 1
                        if isinstance(subproducer, ProducerGroup):
                            subproducer.producers[scope][-1].output = subproducer.output
                        for subproducer_scope in subproducer.scopes:
                            for quantity in subproducer.input[subproducer_scope]:
                                quantity.adopt(subproducer.output[0], subproducer_scope)
                    for output_quantity in subproducer.output:
                        if output_quantity not in self.input[scope]:
                            self.input[scope].append(output_quantity)
            # treat own collection function as subproducer
            self.setup_own_producer()
        log.debug("-----------------------------------------")
        log.debug("| ProducerGroup: {}".format(self.name))
        log.debug("| Call: {}".format(self.call))
        for scope in self.scopes:
            if is_empty(self.input[scope]):
                log.debug("| Inputs ({}): None".format(scope))
            else:
                log.debug(
                    "| Inputs ({}): {}".format(
                        scope, [input.name for input in self.input[scope]]
                    )
                )
        if is_empty(self.output):
            log.debug("| Output: None")
        else:
            log.debug("| Outputs: {}".format([output.name for output in self.output]))
        log.debug("| scopes: {}".format(self.scopes))
        for scope in self.scopes:
            log.debug(
                "| Producers ({}): {}".format(
                    scope, [producer.name for producer in self.producers[scope]]
                )
            )

        log.debug("-----------------------------------------")
        # in the end, determine all parameters used by the producer group
        self.parameters = self.extract_parameters()

    def __str__(self) -> str:
        return "ProducerGroup: {}".format(self.name)

    def __repr__(self) -> str:
        return "ProducerGroup: {}".format(self.name)

    def extract_parameters(self) -> Dict[str, Set[str]]:
        parameters = {}
        for scope in self.scopes:
            parameters[scope] = set()
            for producer in self.producers[scope]:
                if scope in producer.parameters.keys():
                    parameters[scope] = parameters[scope].union(
                        producer.parameters[scope]
                    )
                else:
                    log.warn(
                        f"ProducerGroup {self} is setup for scope {scope}, but producer {producer} is not configured for this scope."
                    )
        return parameters

    def check_producer_scopes(self) -> None:
        """Function to validate the scopes of the subproducers. Internal function.

        Raises:
            Exception: If a scope is not found in the subproducer configuration, an exception is raised
        """
        for scope in self.scopes:
            if scope not in self.producers.keys():
                raise Exception(
                    "ProducerGroup {}: scope {} not found in subproducer configuration: {}".format(
                        self.name, scope, self.producers
                    )
                )

    def setup_own_producer(self) -> None:
        """
        Function to setup a new producer within the ProducerGroup. Internal function.
        """
        producer = Producer(self.name, self.call, self.input, self.output, self.scopes)
        for scope in self.scopes:
            self.producers[scope].append(producer)

    # for a producer group, step iteratively
    # through the subproducers and reserve the output there
    def reserve_output(self, scope: str) -> None:
        """Function used to reserve the output for every producer in the group. Internal function.

        Args:
            scope (str): Scope for which the output is reserved
        """
        for subproducer in self.producers[scope]:
            subproducer.reserve_output(scope)

    def shift(self, name: str, scope: str = "global") -> None:
        """Function used to add a shift for every producer in the group. Wraps the shift function of a producer. Internal function.

        Args:
            name (str): name of the shift
            scope (str, optional): name of the scope. Defaults to "global".
        """
        for producer in self.producers[scope]:
            producer.shift(name, scope)

    def ignore_shift(self, name: str, scope: str = "global") -> None:
        """Function used to ignore a shift for every producer in the group. Wraps the ignore_shift function of a producer. Internal function.

        Args:
            name (str): name of the shift to be ingored
            scope (str, optional): name of the scope. Defaults to "global".
        """
        for producer in self.producers[scope]:
            producer.ignore_shift(name, scope)

    def writecall(
        self, config: Dict[str, Dict[str, Dict[str, str]]], scope: str, shift: str = ""
    ) -> str:
        raise NotImplementedError("This function is not supported to a ProducerGroup!")

    def writecalls(
        self, config: Dict[str, Dict[str, Dict[str, str]]], scope: str
    ) -> List[str]:
        """Function to generate all calls for all shifts and all producer in the group, wraping the writecall function

        Args:
            Configuration dict containing the parameter and input / output information for the ProducerGroup
            scope (str): The scope in which the ProducerGroup is to be called

        Raises:
            Exception: Raises an Exception, of the requested scope is not valid for the ProducerGroup

        Returns:
            List[str]: Returns a list of C++ calls for all shifts of the ProducerGroup
        """
        calls: List[str] = []
        for producer in self.producers[scope]:
            # duplicate outputs of vector subproducers if they were generated automatically
            if (
                not is_empty(self.call)
                and isinstance(producer, VectorProducer)
                and not is_empty(producer.output)
            ):
                for i in range(len(config["nominal"][producer.vec_configs[0]]) - 1):
                    producer.output.append(
                        producer.output[0].copy(
                            "PG_internal_quantity_%i" % self.__class__.PG_count
                        )
                    )
                    self.__class__.PG_count += 1
                    if producer.output[-1] not in self.input[scope]:
                        self.input[scope].append(producer.output[-1])
            # retrieve calls of subproducers
            calls.extend(producer.writecalls(config, scope))
        return calls

    def get_inputs(self, scope: str) -> List[q.Quantity]:
        """Get a list of all inputs of the ProducerGroup in a given scope

        Args:
            scope (str): The scope in which the inputs are requested

        Raises:
            Exception: Raises an Exception, of the requested scope is not valid for the ProducerGroup

        Returns:
            List[str]: Returns a list of Quantity objects, which are the inputs of the ProducerGroup
        """

        inputs: List[q.Quantity] = []
        log.debug("Getting inputs for {}".format(self.name))
        for subproducer in self.producers[scope]:
            log.debug("  --> {} {}".format(subproducer, subproducer.get_inputs(scope)))
            inputs.extend(subproducer.get_inputs(scope))
        return inputs

    def get_outputs(self, scope: str) -> List[q.Quantity]:
        """Get a list of all outputs of the ProducerGroup in a given scope

        Args:
            scope (str): The scope in which the outputs are requested

        Raises:
            Exception: Raises an Exception, of the requested scope is not valid for the ProducerGroup

        Returns:
            List[Union[q.QuantityGroup, q.Quantity]]: Returns a list of Quantity objects, which are the outputs of the ProducerGroup
        """
        outputs: List[q.Quantity] = []
        log.debug("Getting outputs for {}".format(self.name))
        for subproducer in self.producers[scope]:
            log.debug("  --> {} {}".format(subproducer, subproducer.get_outputs(scope)))
            outputs.extend(subproducer.get_outputs(scope))
        return outputs


class Filter(ProducerGroup):
    def __init__(
        self,
        name: str,
        call: str,
        input: Union[List[q.Quantity], Dict[str, List[q.Quantity]]],
        scopes: List[str],
        subproducers: Union[
            List[Producer | ProducerGroup],
            Dict[str, List[Producer | ProducerGroup]],
        ],
    ):
        """A Filter is used to filter events. Wraps the BaseFilter class, and is a ProducerGroup.

        Args:
            name (str): name of the filter
            call (str): the C++ function call to be used for the filter
            input (Union[List[q.Quantity], Dict[str, List[q.Quantity]]]): The input quantities for the filter
            scopes (List[str]): The scopes in which the filter is to be called
            subproducers (Union[ List[Producer  |  ProducerGroup], Dict[str, List[Producer  |  ProducerGroup]], ]): The subproducers of the filter
        """
        self.__class__.PG_count = ProducerGroup.PG_count
        super().__init__(name, call, input, None, scopes, subproducers)
        ProducerGroup.PG_count = self.__class__.PG_count
        self.call: str = call

    def __str__(self) -> str:
        return "Filter: {}".format(self.name)

    def __repr__(self) -> str:
        return "Filter: {}".format(self.name)

    def setup_own_producer(self) -> None:
        for scope in self.scopes:
            self.producers[scope].append(
                BaseFilter(self.name, self.call, self.input, self.scopes)
            )


def CollectProducersOutput(
    producers: List[Union[ProducerGroup, Producer]], scope: str
) -> Set[q.Quantity]:
    output: Set[q.Quantity] = set()
    for producer in producers:
        if not is_empty(producer.output):
            output |= set(producer.output)
        if isinstance(producer, ProducerGroup):
            try:
                for prod in producer.producers[scope]:
                    output |= CollectProducerOutput(prod, scope)
            except KeyError:
                raise ConfigurationError(
                    "Scope {} not found in subproducer configuration: {}".format(
                        scope, producer.producers
                    )
                )
    return set(output)


def CollectProducerOutput(
    producer: Union[ProducerGroup, Producer], scope: str
) -> Set[q.Quantity]:
    output: Set[q.Quantity] = set()
    if producer.output:
        output |= set(producer.output)
    if isinstance(producer, ProducerGroup):
        try:
            for prod in producer.producers[scope]:
                output |= CollectProducerOutput(prod, scope)
        except KeyError:
            raise ConfigurationError(
                "Scope {} not found in subproducer configuration: {}".format(
                    scope, producer.producers
                )
            )
    return set(output)


TProducerInput = Union[Producer, ProducerGroup, List[Union[Producer, ProducerGroup]]]
TProducerSet = Union[Producer, ProducerGroup, Set[Union[Producer, ProducerGroup]]]
TProducerStore = Dict[str, List[Union[Producer, ProducerGroup]]]
TProducerListStore = Dict[str, TProducerStore]
