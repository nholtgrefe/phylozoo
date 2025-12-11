"""
Tools module.

This module provides various utility functions.
"""

import random
import string



def id_generator(size: int = 6, chars: str | None = None) -> str:
    """
    Generate a random identifier string.

    Parameters
    ----------
    size : int, optional
        Length of the identifier, by default 6
    chars : str | None, optional
        Characters to use for generation. If None, uses ASCII letters and digits, by default None

    Returns
    -------
    str
        Random identifier string
    """
    if chars is None:
        chars = string.ascii_letters + string.digits
    return "".join(random.choice(chars) for _ in range(size))


def normalized_l_p_norm(vector: list[float], p: float = 2.0) -> float:
    """
    Calculate the normalized L-p norm of a vector.

    Parameters
    ----------
    vector : list[float]
        Input vector
    p : float, optional
        Norm parameter (p=2 is Euclidean norm), by default 2.0

    Returns
    -------
    float
        Normalized L-p norm value

    Notes
    -----
    This is a placeholder implementation. Implement actual norm calculation here.
    """
    if not vector:
        return 0.0
    # Placeholder implementation
    return sum(abs(x) ** p for x in vector) ** (1.0 / p) / len(vector) if p > 0 else 0.0
