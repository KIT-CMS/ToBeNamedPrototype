from __future__ import annotations  # needed for type annotations in > python 3.7
from typing import List, Union, Dict, Any
from code_generation.exceptions import (
    SampleConfigurationError,
    EraConfigurationError,
)
from code_generation.helpers import is_empty

ConfigurationParameters = Union[str, int, float, bool]

ModifierDict = Dict[
    str,
    Union[
        List[ConfigurationParameters],
        List[Dict[str, ConfigurationParameters]],
        Dict[str, ConfigurationParameters],
        ConfigurationParameters,
        List[Any],
    ],
]
ModifierResolved = Union[
    List[ConfigurationParameters],
    Dict[str, ConfigurationParameters],
    List[Dict[str, ConfigurationParameters]],
    ConfigurationParameters,
    None,
]


class Modifier(object):
    def __init__(self, modifier_dict: ModifierDict):
        """Genric Modifier class.

        Args:
            modifier_dict : Dict of parameters to be modified.
        """
        self.modifier_dict = modifier_dict

    def __str__(self) -> str:
        return "Modifier: {}".format(self.modifier_dict)

    def __repr__(self) -> str:
        return "Modifier: {}".format(self.modifier_dict)


class SampleModifier(Modifier):
    def __init__(
        self,
        modifier_dict: ModifierDict,
        default: Union[str, int, float, bool, Dict, None] = None,
    ):
        """A Sample Modifier is a Modifier, that modifies the configuration based on the given sample

        Args:
            modifier_dict : A dict containing the information, how a parameter should be modified based on the sample.
            default: If set, the default is used for all sample not specified in the modifier dict. Defaults to None.
        """
        super(SampleModifier, self).__init__(modifier_dict)
        self.modifier_dict = modifier_dict
        self.default = default
        self.samples: List[str] = list(self.modifier_dict.keys())

    def apply(self, sample: str) -> ModifierResolved:
        """
        Applies the modifier to the given sample.

        Args:
            sample: The sample to apply the modifier to

        Returns:
            A modified configuration
        """
        if sample in self.samples:
            return self.modifier_dict[sample]
        elif not is_empty(self.default):
            return self.default
        else:
            raise SampleConfigurationError(sample, self.samples)


class EraModifier(Modifier):
    def __init__(
        self,
        modifier_dict: ModifierDict,
        default: Union[str, int, float, bool, Dict, None] = None,
    ):
        """A Era Modifier is a Modifier, that modifies the configuration based on the given era

        Args:
            modifier_dict : A dict containing the information, how a parameter should be modified based on the sample.
            default: If set, the default is used for all sample not specified in the modifier dict. Defaults to None.
        """
        super(EraModifier, self).__init__(modifier_dict)
        self.modifier_dict = modifier_dict
        self.default = default
        self.eras: List[str] = list(self.modifier_dict.keys())

    def apply(self, era: str) -> ModifierResolved:
        """
        Applies the modifier to the given era.

        Args:
            era: The era to apply the modifier to

        Returns:
            A modified configuration
        """
        if era in self.eras:
            return self.modifier_dict[era]
        elif not is_empty(self.default):
            return self.default
        else:
            raise EraConfigurationError(era, self.eras)
