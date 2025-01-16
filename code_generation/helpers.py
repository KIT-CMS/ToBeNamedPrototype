from __future__ import annotations  # needed for type annotations in > python 3.7

# File with helper functions for the CROWN code generation


def is_empty(value):
    """
    Check if a value is empty.

    Args:
        value: The value that should be checked.

    Returns:
        bool: Whether the input value is considered 'empty'
    """
    # List of all values that should be considered empty despite not having a length.
    empty_values = [None]

    try:
        length = len(value)
    except TypeError:
        length = -1
    bool_val = value in empty_values or length == 0
    return bool_val
