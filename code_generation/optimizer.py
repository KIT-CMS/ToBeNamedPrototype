from __future__ import annotations  # needed for type annotations in > python 3.7
from code_generation.quantity import NanoAODQuantity, Quantity
from code_generation.producer import Filter, BaseFilter, Producer, ProducerGroup
from typing import Set, Tuple, Union, List
import logging

log = logging.getLogger(__name__)


class ProducerOrdering:
    """
    Class used to check if the producer ordering is correct,
    If it is not, the Optimize function will auto-correct it.
    Additionally, the optimize attempts to put filters at the top of the list,
    or as far up as possible. A wrong configuration due to missing inputs will also be caught here.

    If the scope is not global, the outputs generated by producers in the global scope are also considered.
    """

    def __init__(
        self,
        global_producers: List[Producer | ProducerGroup],
        scope: str,
        producer_ordering: List[Producer | ProducerGroup],
    ):
        """
        Init function

        Args:
            config: The configuration dictionary
            global_producers: The list of producers in the global scope
            scope: The scope of the producer ordering
            producer_ordering: The producer ordering to be optimized
        """
        self.global_producers: List[Producer | ProducerGroup] = global_producers
        self.ordering: List[Producer | ProducerGroup] = producer_ordering
        self.size = len(self.ordering)
        self.scope = scope
        self.optimized: bool = False
        self.optimized_ordering: List[Producer | ProducerGroup] = []
        self.global_outputs = self.get_global_outputs()

    def get_position(self, producer: Producer | ProducerGroup) -> int:
        """
        Helper Function to get the position of a producer in the ordering list

        Args:
            producer: The producer to get the position of

        Returns:
            The position of the producer in the current ordering
        """
        for i, prod in enumerate(self.ordering):
            if prod == producer:
                return i
        raise Exception("Producer not in ordering")
        return -1

    def get_producer(self, position: int) -> Producer | ProducerGroup:
        """
        Helper function to get the producer at a given position

        Args:
            position: The position of the producer

        Returns:
            The producer at the given position
        """
        return self.ordering[position]

    def get_global_outputs(self) -> List[Quantity]:
        """
        Function used to generate a list of all outputs generated by the global scope.

        Args:
            None

        Returns:
            A list of all outputs generated by the global scope
        """
        outputs: List[Quantity] = []
        for producer in self.global_producers:
            if producer.get_outputs("global") is not None:
                outputs.extend(
                    [
                        quantity
                        for quantity in producer.get_outputs("global")
                        if not isinstance(quantity, NanoAODQuantity)
                    ]
                )
        return outputs

    def MoveFiltersUp(self) -> None:
        """
        Function used to relocate all filters to the top of the ordering

        Args:
            None

        Returns:
            None
        """
        new_ordering: List[Producer | ProducerGroup] = []
        for producer in self.ordering:
            if isinstance(producer, Filter) or isinstance(producer, BaseFilter):
                new_ordering.insert(0, producer)
            else:
                new_ordering.append(producer)
        for i, prod in enumerate(self.ordering):
            log.debug(" --> {}. : {}".format(i, prod))
        for i, prod in enumerate(new_ordering):
            log.debug(" --> {}. : {}".format(i, prod))
        self.ordering = new_ordering

    def Optimize(self) -> None:
        """
        The main function of this class. During the optimization,
        finding a correct ordering is attempted. This is done as follows:

        1. Bring all filters to the beginning of the ordering.

        2. Check if the ordering is already correct. The ordering is correct,
           if, for all producers in the ordering, all inputs can be found in
           the outputs of preceding producers. If the scope is not global,
           all outputs from producers in the global scope are also considered.

        3. If the ordering is correct, return.

        4. If the ordering is not correct,

            1. find all inputs, that have to be produced before the wrong producer
            2. put one producer, which is responsible for creating the input, in front of the wrong producer
            3. repeat from step 2

        The sorting algorithm should take at most ``2*(number of producers)`` steps.
        If this limit is reached, the optimization is
        considered to be failed and an Exception is raised.
        If a missing input cant be found in all outputs,
        the Optimize function will raise an Exception.

        Args:
            None

        Returns:
            None
        """
        # first bring filters to the top
        self.MoveFiltersUp()
        # run optimization in a while loop
        counter = 0
        while not self.optimized:
            if counter > 2 * self.size + 1:
                log.error("Could not optimize ordering")
                log.error("Please check, if all needed producers are activated")
                raise Exception
            wrongProducer, wrong_inputs = self.check_ordering()
            if wrongProducer is not None:
                producers_to_relocate = self.find_inputs(wrongProducer, wrong_inputs)
                # if len(producers_to_relocate) == 0:
                #     self.optimized = True
                #     break
                # else:
                for producer_to_relocate in producers_to_relocate:
                    counter += 1
                    self.relocate_producer(
                        producer_to_relocate,
                        self.get_position(producer_to_relocate),
                        self.get_position(wrongProducer),
                    )
        self.optimized_ordering = self.ordering
        log.info(
            "Optimization for scope {} done after {} steps: {}".format(
                self.scope, counter, self.optimized_ordering
            )
        )

    def check_ordering(
        self,
    ) -> Tuple[Union[Producer, ProducerGroup, None], List[Quantity]]:
        """
        Function used to check the ordering.
        If at least of one the inputs of a producer cannot be found in
        the list of all preceding outputs, the ordering is not correct.
        If the whole odering is correct, the optimized flag is set to
        true and the ordering is considered to be correct.

        Args:
            None

        Returns:
            A tuple of the wrong producer and a list of the inputs
            that are not found in the outputs
        """
        outputs = []
        if self.scope != "global":
            outputs = self.global_outputs
        for producer_to_check in self.ordering:
            temp_outputs = producer_to_check.get_outputs(self.scope)
            if temp_outputs is not None:
                outputs.extend(
                    [
                        quantity
                        for quantity in temp_outputs
                        if not isinstance(quantity, NanoAODQuantity)
                    ]
                )
            inputs = [
                quantity
                for quantity in producer_to_check.get_inputs(self.scope)
                if not isinstance(quantity, NanoAODQuantity)
            ]
            invalid_inputs = self.invalid_inputs(inputs, outputs)
            if len(invalid_inputs) > 0:
                return producer_to_check, invalid_inputs
        self.optimized = True
        return None, []

    def invalid_inputs(
        self, inputs: List[Quantity], outputs: List[Quantity]
    ) -> List[Quantity]:
        """
        Helper function used to check, if for a list of inputs, a match in a list of outputs can be found.

        Args:
            inputs: The list of quantity inputs to check
            outputs: The list of quantity outputs to check

        Returns:
            wrong_inputs: A list of inputs, that are not in the outputs list
        """
        wrong_inputs: List[Quantity] = []
        for input in inputs:
            if input not in outputs:
                wrong_inputs.append(input)
                log.debug("Input {} not in outputs".format(input))
        return wrong_inputs

    def find_inputs(
        self, producer: Producer | ProducerGroup, inputs: List[Quantity]
    ) -> List[Producer | ProducerGroup]:
        """
        Function used to locate the producers responsible for creating the given inputs.
        The function return a list of producers, that have to be run before the tested producer.
        If a needed input is not found, the function raise an Exception.

        Args:
            producer: The producer to check
            inputs: The list of inputs to check

        Returns:
            producers_to_relocate: A list of producers, that have to be run before the given producer
        """
        producers_to_relocate: Set[Producer | ProducerGroup] = set()
        log.debug("Trying to find inputs {}".format(inputs))
        for input in inputs:
            found = False
            for producer in self.ordering:
                if input in producer.get_outputs(self.scope):
                    found = True
                    log.debug(
                        "found {} in outputs from producer {} ({}) --> Rank {}".format(
                            input,
                            producer,
                            producer.get_outputs(self.scope),
                            self.ordering.index(producer),
                        )
                    )
                    producers_to_relocate.add(producer)
            if self.scope != "global":
                for producer in self.global_producers:
                    if input in producer.get_outputs("global"):
                        found = True
                        log.debug("found {} in global scope ..".format(input))
            if not found:
                log.error(
                    "{} (required by {}) was not found in any producer output!".format(
                        input, producer
                    )
                )
                log.error("Please check if all needed producers are activated !")
                raise Exception
        return list(producers_to_relocate)

    def relocate_producer(
        self,
        producer: Producer | ProducerGroup,
        old_position: int,
        new_position: int,
    ) -> None:
        """
        Function used to relocate a producer to a given position.

        Args:
            producer: The producer to relocate
            old_position: The old position of the producer in the ordering
            new_position: The new position of the producer in the ordering
        """
        log.debug(
            "Relocating Producer {} from rank {} to rank {}".format(
                producer, old_position, new_position
            )
        )
        updated_ordering = list(range(self.size))
        if old_position > new_position:
            for position in updated_ordering:
                if position <= old_position and position > new_position:
                    updated_ordering[position] = position - 1
                if position == new_position:
                    updated_ordering[position] = old_position
        if old_position < new_position:
            for position in updated_ordering:
                if position >= old_position and position < new_position:
                    updated_ordering[position] = position + 1
                if position == new_position:
                    updated_ordering[position] = old_position
        if old_position == new_position:
            log.debug("How did we get here ??")
        new_ordering = [self.ordering[i] for i in updated_ordering]
        log.debug(
            "New ordering - ",
        )
        for i, prod in enumerate(new_ordering):
            log.debug(" --> {}. : {}".format(i, prod))
        self.ordering = new_ordering
