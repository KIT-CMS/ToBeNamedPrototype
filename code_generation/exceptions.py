from __future__ import annotations  # needed for type annotations in > python 3.7
from typing import List, Set, Union
from code_generation.quantity import Quantity
from code_generation.helpers import is_empty


class ConfigurationError(Exception):
    """
    Exception raised when the configuration provided by the user is not valid.
    """

    def __init__(self, message: str):
        self.message = message
        super().__init__(self.message)


class InvalidOutputError(ConfigurationError):
    """
    Exception raised when the list of output provided by the user is not valid.
    """

    def __init__(self, scope: str, outputs: Union[Set[Quantity], List[Quantity]]):
        self.message = "The required outputs {} for the scope '{}' are not provided by any producer \n Please modify the list of producers or the list of required outputs".format(
            outputs, scope
        )
        super().__init__(self.message)


class InvalidInputError(ConfigurationError):
    """
    Exception raised when the list of avialable inputs does not cover all quantities required.
    """

    def __init__(self, scope: str, outputs: Union[Set[str], List[str]]):
        self.message = "The required inputs {} for the scope '{}' are not provided by any inputfile or producer \n Please check the error message above to find all misconfigured producers".format(
            outputs, scope
        )
        super().__init__(self.message)


class ScopeConfigurationError(ConfigurationError):
    """
    Exception raised when the scope configuration provided by the user is not valid.
    """

    def __init__(self, scopes: Set[str], available_scopes: Union[Set[str], List[str]]):
        self.message = "Scopes {} cannot be used in the configuration since it is not setup properly. Available scopes are {}".format(
            scopes, available_scopes
        )
        super().__init__(self.message)


class SampleConfigurationError(ConfigurationError):
    """
    Exception raised when the sample type configuration provided by the user is not valid.
    """

    def __init__(self, sample: str, available_samples: Union[Set[str], List[str]]):
        self.message = "Sampletype {} cannot be used in the configuration since it is not setup properly. Available samples types are {}".format(
            sample, available_samples
        )
        super().__init__(self.message)


class SampleRuleConfigurationError(ConfigurationError):
    """
    Exception raised when the sample type used for a Rule provided by the user is not valid.
    """

    def __init__(
        self,
        sample: str,
        rule,
        available_samples: Union[Set[str], List[str]],
    ):
        self.message = "Sampletype {} cannot be used in Rule {} since the type is not defined. Available samples types are {}".format(
            sample, rule, available_samples
        )
        super().__init__(self.message)


class EraConfigurationError(ConfigurationError):
    """
    Exception raised when the era configuration provided by the user is not valid.
    """

    def __init__(self, era: str, available_eras: Union[Set[str], List[str]]):
        self.message = "Era {} cannot be used in the configuration since it is not setup properly. Available eras are {}".format(
            era, available_eras
        )
        super().__init__(self.message)


class InvalidProducerConfigurationError(ConfigurationError):
    """
    Exception raised when the producer configuration provided by the user is not valid.
    """

    def __init__(self, producer: str):
        self.message = "Producer {} is not setup properly".format(
            producer,
        )
        super().__init__(self.message)


class InvalidShiftError(ConfigurationError):
    """
    Exception raised when the shift configuration provided by the user is not valid.
    """

    def __init__(self, shift: str, sample: str, scope: Union[str, None] = None):
        if is_empty(scope):
            self.message = "Shift {} is not setup properly or not available for sampletype {}".format(
                shift, sample
            )
        else:
            self.message = "Shift {} is not setup for scope {} or not available for sampletype {}".format(
                shift, scope, sample
            )
        super().__init__(self.message)


class InsufficientShiftInformationError(ConfigurationError):
    """
    Exception raised if "all" is used for the shift settings for a FriendTree
    """

    def __init__(self, shift: Union[str, List[str]], available_shifts: List[str]):
        self.message = "Shift(s) {} cannot be used for the FriendTreeConfiguration, it is not found in the provided ntuples \n".format(
            shift
        )
        self.message += "Available shifts are: {}".format(available_shifts)
        super().__init__(self.message)
