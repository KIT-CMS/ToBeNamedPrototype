import logging
from typing import List, Union

from code_generation.producer import (
    CollectProducersOutput,
    Producer,
    ProducerGroup,
    TProducerInput,
    TProducerStore,
)
from code_generation.quantity import QuantitiesStore
from code_generation.exceptions import SampleRuleConfigurationError

log = logging.getLogger(__name__)


class ProducerRule:
    def __init__(
        self,
        producers: TProducerInput,
        samples: Union[str, List[str]] = [],
        exclude_samples: Union[str, List[str]] = [],
        scopes: Union[str, List[str]] = "global",
        invert: bool = False,
        update_output: bool = True,
    ):
        """ProducerRule is a base class for all rules that modify producers.

        Args:
            producers: A list of producers or producer groups to be modified.
            samples: A list of samples, for which the rule should be applied. Only one of samples and exclude_samples can be defined.
            exclude_samples: A list of samples, for which the rule should not be applied. Only one of samples and exclude_samples can be defined.
            scopes: The scopes, in which the rule should be applied. Defaults to "global".
            invert: If set, the invert of the rule is applied. Defaults to False.
            update_output: If set, the output quantities are updated. Defaults to True.
        """
        if isinstance(producers, ProducerGroup) or isinstance(producers, Producer):
            producers = [producers]
        self.producers = producers
        if isinstance(scopes, str):
            scopes = [scopes]
        self.scopes = scopes
        self.invert = invert
        self.update_output = update_output
        self.global_scope = "global"
        if isinstance(exclude_samples, str):
            self.exclude_samples = [exclude_samples]
        else:
            self.exclude_samples = exclude_samples
        if isinstance(samples, str):
            self.samples = [samples]
        else:
            self.samples = samples

    def set_available_sampletypes(self, available_samples) -> None:
        # sanitize input
        if isinstance(available_samples, str):
            self.available_samples = [available_samples]
        else:
            self.available_samples = available_samples
        # make sure that either samples or exclude_samples are defined
        if self.exclude_samples == [] and self.samples == []:
            raise ValueError(
                f"ProducerRule: Either samples or exclude_samples have to be defined!: (Rule: {self}, Samples: {self.samples}, Excluded Samples: {self.exclude_samples})"
            )
        if self.exclude_samples != [] and self.samples != []:
            raise ValueError(
                f"ProducerRule: Both samples and are exclude_samples are defined, pick one!: (Rule: {self}, Samples: {self.samples}, Excluded Samples: {self.exclude_samples})"
            )
        # make sure that the sampletypes are valid
        self.validate_sampletypes(self.samples)
        self.validate_sampletypes(self.exclude_samples)
        # if exclude_samples are defined, we have to contstruct the list of samples them from the list of available samples
        if self.exclude_samples != []:
            self.samples = list(set(self.available_samples) - set(self.exclude_samples))

    def set_scopes(self, scopes: List[str]) -> None:
        if isinstance(scopes, str):
            scopes = [scopes]
        self.scopes = scopes

    def affected_scopes(self) -> List[str]:
        return self.scopes

    def affected_producers(self) -> List[Union[Producer, ProducerGroup]]:
        return self.producers

    def set_global_scope(self, global_scope: str) -> None:
        self.global_scope = global_scope

    def validate_sampletypes(self, sampletypes: List[str]) -> None:
        """Function to check, if a rule is applicable, or if one of the defined samples is not available.

        Args:
            available_samples (List[str]): List of available samples.

        Returns:
            None
        """
        for sample in sampletypes:
            if sample not in self.available_samples:
                raise SampleRuleConfigurationError(sample, self, self.available_samples)

    # Evaluate whether modification should be applied depending on sample and inversion flag
    def is_applicable(self, sample: str) -> bool:
        applicable = sample in self.samples
        if self.invert:
            applicable = not applicable
        return applicable

    # Placeholder for the actual operation on a list. To be overwritten by inheriting classes
    def update_producers(
        self,
        producers_to_be_updated: TProducerStore,
        unpacked_producers: TProducerStore,
    ) -> None:
        log.error("Operation not implemented for ProducerRule base class!")

    def update_outputs(self, outputs_to_be_updated: QuantitiesStore) -> None:
        log.error("Operation not implemented for ProducerRule base class!")

    def apply(
        self,
        sample: str,
        producers_to_be_updated: TProducerStore,
        unpacked_producers: TProducerStore,
        outputs_to_be_updated: QuantitiesStore,
    ) -> None:
        if self.is_applicable(sample):
            log.critical(f"Applying rule {self} for sample {sample}")
            log.debug("For sample {}, applying >> {} ".format(sample, self))
            self.update_producers(producers_to_be_updated, unpacked_producers)
            self.update_outputs(outputs_to_be_updated)


# Modifier class that can remove producers from lists
class RemoveProducer(ProducerRule):
    def __str__(self) -> str:
        return "ProducerRule - remove {} for {} in scopes {}".format(
            self.producers, self.samples, self.scopes
        )

    def __repr__(self) -> str:
        return "ProducerRule - remove {} for {} in scopes {}".format(
            self.producers, self.samples, self.scopes
        )

    def update_producers(
        self,
        producers_to_be_updated: TProducerStore,
        unpacked_producers: TProducerStore,
    ) -> None:
        log.debug("Producers to be updated: {}".format(producers_to_be_updated))
        log.debug("scopes: {}".format(self.scopes))
        log.debug("Producers to be removed: {}".format(self.producers))
        for scope in self.scopes:
            for producer in self.producers:
                if producer in producers_to_be_updated[scope]:
                    log.debug(
                        "RemoveProducer: Removing {} from producer in scope {}".format(
                            producer, scope
                        )
                    )
                    producers_to_be_updated[scope].remove(producer)
                else:
                    # in this case, the producer does not exist, possibly because it is part of a producer group,
                    # so we have to check this further
                    if producer in unpacked_producers[scope].keys():
                        # if the producer is part of a producer group, we have to remove the whole group
                        # and add all remaining producers from the group to the list of producers to be updated
                        log.debug("Found {} within a producer group".format(producer))
                        corresponding_producer_group = unpacked_producers[scope][
                            producer
                        ]
                        log.debug(
                            "Removing {} from producer group {}".format(
                                producer, corresponding_producer_group
                            )
                        )
                        log.debug(
                            "Replacing {} with its unpacked producers".format(
                                corresponding_producer_group
                            )
                        )
                        producers_to_be_updated[scope].remove(
                            corresponding_producer_group
                        )
                        for unpacked_producer in unpacked_producers[scope]:
                            if (
                                unpacked_producers[scope][unpacked_producer]
                                == corresponding_producer_group
                            ):
                                producers_to_be_updated[scope].append(unpacked_producer)
                        # now remove the producer we initially wanted to remove
                        log.debug(producers_to_be_updated[scope])
                        producers_to_be_updated[scope].remove(producer)

                    else:
                        raise ConnectionError(
                            "Producer {} not found in scope {}, cannot apply \n {}".format(
                                producer, scope, self
                            )
                        )

    def update_outputs(self, outputs_to_be_updated: QuantitiesStore) -> None:
        if self.update_output:
            outputs: QuantitiesStore = {}
            # if the producer is in the global scope, we add the output to all running scopes
            if self.scopes == [self.global_scope]:
                log.debug(
                    "Updating outputs for producer in global scope --> adding output to all scopes"
                )
                scopes = outputs_to_be_updated.keys()
                for scope in scopes:
                    outputs[scope] = CollectProducersOutput(
                        self.producers, self.global_scope
                    )
            else:
                scopes = self.scopes
                for scope in scopes:
                    outputs[scope] = CollectProducersOutput(self.producers, scope)
            for scope in scopes:
                for output in outputs[scope]:
                    if output in outputs_to_be_updated[scope]:
                        log.debug(
                            "RemoveProducer: Removing {} from outputs in scope {}".format(
                                output, scope
                            )
                        )
                        outputs_to_be_updated[scope].remove(output)


# Modifier class that can append producers to lists
class AppendProducer(ProducerRule):
    def __str__(self) -> str:
        return "ProducerRule - add {} for {} in scopes {}".format(
            self.producers, self.samples, self.scopes
        )

    def __repr__(self) -> str:
        return "ProducerRule - add {} for {} in scopes {}".format(
            self.producers, self.samples, self.scopes
        )

    def update_producers(
        self,
        producers_to_be_updated: TProducerStore,
        unpacked_producers: TProducerStore,
    ) -> None:
        for scope in self.scopes:
            for producer in self.producers:
                log.debug(
                    "AppendProducer: Adding {} to producers in scope {}".format(
                        producer, scope
                    )
                )
                producers_to_be_updated[scope].append(producer)

    def update_outputs(self, outputs_to_be_updated: QuantitiesStore) -> None:
        if self.update_output:
            outputs: QuantitiesStore = {}
            # if the producer is in the global scope, we add the output to all running scopes
            if self.scopes == [self.global_scope]:
                log.debug(
                    "Updating outputs for producer in global scope --> adding output to all scopes"
                )
                scopes = outputs_to_be_updated.keys()
                for scope in scopes:
                    outputs[scope] = CollectProducersOutput(
                        self.producers, self.global_scope
                    )
            else:
                scopes = self.scopes
                for scope in scopes:
                    outputs[scope] = CollectProducersOutput(self.producers, scope)
            for scope in scopes:
                for output in outputs[scope]:
                    log.debug(
                        "AppendProducer: Adding {} to outputs in scope {}".format(
                            output, scope
                        )
                    )
                    outputs_to_be_updated[scope].add(output)


class ReplaceProducer(ProducerRule):
    def __str__(self) -> str:
        return "ProducerRule - replace {} with {} for {} in scopes {}".format(
            self.producers[0], self.producers[1], self.samples, self.scopes
        )

    def __repr__(self) -> str:
        return "ProducerRule - replace {} for {} in scopes {}".format(
            self.producers, self.samples, self.scopes
        )

    def update_producers(
        self,
        producers_to_be_updated: TProducerStore,
        unpacked_producers: TProducerStore,
    ) -> None:
        log.debug("Producers to be updated: {}".format(producers_to_be_updated))
        log.debug("scopes: {}".format(self.scopes))
        log.debug("Producers to be replaced: {}".format(self.producers))
        producer = self.producers[0]
        new_producer = self.producers[1]
        for scope in self.scopes:
            if producer in producers_to_be_updated[scope]:
                log.debug(
                    "ReplaceProducer: Replace {} from producer in scope {} with {}".format(
                        producer, scope, new_producer
                    )
                )
                producers_to_be_updated[scope].remove(producer)
                producers_to_be_updated[scope].append(new_producer)
            else:
                # in this case, the producer does not exist, possibly because it is part of a producer group,
                # so we have to check this further
                if producer in unpacked_producers[scope].keys():
                    # if the producer is part of a producer group, we have to remove the whole group
                    # and add all remaining producers from the group to the list of producers to be updated
                    log.debug("Found {} within a producer group".format(producer))
                    corresponding_producer_group = unpacked_producers[scope][producer]
                    log.debug(
                        "Replacing {} from producer group {} with {}".format(
                            producer, corresponding_producer_group, new_producer
                        )
                    )
                    log.debug(
                        "Replacing {} with its unpacked producers".format(
                            corresponding_producer_group
                        )
                    )
                    producers_to_be_updated[scope].remove(corresponding_producer_group)
                    for unpacked_producer in unpacked_producers[scope]:
                        if (
                            unpacked_producers[scope][unpacked_producer]
                            == corresponding_producer_group
                        ):
                            producers_to_be_updated[scope].append(unpacked_producer)
                    # now remove the producer we initially wanted to remove
                    log.debug(producers_to_be_updated[scope])
                    producers_to_be_updated[scope].remove(producer)
                    producers_to_be_updated[scope].append(new_producer)

                else:
                    raise ConnectionError(
                        "Producer {} not found in scope {}, cannot apply \n {}".format(
                            producer, scope, self
                        )
                    )

    def update_outputs(self, outputs_to_be_updated: QuantitiesStore) -> None:
        if self.update_output:
            removed_outputs: QuantitiesStore = {}
            added_outputs: QuantitiesStore = {}
            # if the producer is in the global scope, we add the output to all running scopes
            if self.scopes == [self.global_scope]:
                log.debug(
                    "Updating outputs for producer in global scope --> adding output to all scopes"
                )
                scopes = outputs_to_be_updated.keys()
                for scope in scopes:
                    removed_outputs[scope] = CollectProducersOutput(
                        [self.producers[0]], self.global_scope
                    )
                    added_outputs[scope] = CollectProducersOutput(
                        [self.producers[1]], self.global_scope
                    )
            else:
                scopes = self.scopes
                for scope in scopes:
                    removed_outputs[scope] = CollectProducersOutput(
                        [self.producers[0]], scope
                    )
                    added_outputs[scope] = CollectProducersOutput(
                        [self.producers[1]], scope
                    )
            if added_outputs == removed_outputs:
                log.debug(
                    f"Outputs {added_outputs} are identical, no need to update outputs"
                )
            else:
                for scope in scopes:
                    for removed_output in removed_outputs[scope]:
                        if removed_output in outputs_to_be_updated[scope]:
                            log.debug(
                                "ReplaceProducer: Removing {} from outputs in scope {}".format(
                                    removed_output, scope
                                )
                            )
                            outputs_to_be_updated[scope].remove(removed_output)
                    for added_output in added_outputs[scope]:
                        if added_output not in outputs_to_be_updated[scope]:
                            log.debug(
                                "ReplaceProducer: Adding {} from outputs in scope {}".format(
                                    added_output, scope
                                )
                            )
                            outputs_to_be_updated[scope].add(added_output)
