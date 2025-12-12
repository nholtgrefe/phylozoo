"""
Distance matrix I/O module.

This module provides functions for reading and writing distance matrices to/from files.

TODO: Add reading from Nexus files, and string?
    Add other nexus options (upper, etc.)
"""

from __future__ import annotations

import os
from typing import TypeVar

from .base import DistanceMatrix

T = TypeVar('T')


def to_nexus(distance_matrix: DistanceMatrix) -> str:
    """
    Convert a distance matrix to a NEXUS format string.
    
    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to convert.
    
    Returns
    -------
    str
        The NEXUS format string representation of the distance matrix.
    
    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.io import to_nexus
    >>> 
    >>> matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
    >>> dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
    >>> nexus_str = to_nexus(dm)
    >>> print(nexus_str[:50])
    #NEXUS
    
    BEGIN Taxa;
        DIMENSIONS ntax=3;
    
    Notes
    -----
    The NEXUS format includes:
    - Taxa block with label names
    - Distances block with lower triangular matrix
    - Format: triangle=LOWER diagonal LABELS
    """
    n = len(distance_matrix)
    nexus_string = "#NEXUS\n\nBEGIN Taxa;\n"
    nexus_string += f"    DIMENSIONS ntax={n};\n"
    nexus_string += "    TAXLABELS\n"
    nexus_string += "\n".join(f"        {taxon}" for taxon in distance_matrix.labels)
    nexus_string += ";\nEND;\n\n"
    
    nexus_string += "BEGIN Distances;\n"
    nexus_string += f"    DIMENSIONS ntax={n};\n"
    nexus_string += "    FORMAT triangle=LOWER diagonal LABELS;\n"
    nexus_string += "    MATRIX\n"
    
    matrix = distance_matrix._matrix
    for i, taxon in enumerate(distance_matrix.labels):
        row = " ".join(f"{matrix[i, j]:.6f}" for j in range(i + 1))
        nexus_string += f"    {taxon} {row}\n"
    
    nexus_string += ";\nEND;\n"
    
    return nexus_string


def save_nexus(
    distance_matrix: DistanceMatrix,
    file_path: str,
    overwrite: bool = False
) -> None:
    """
    Save a distance matrix to a NEXUS format file.
    
    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to save.
    file_path : str
        Path to the output file. Must end with '.nexus' or '.nex'.
    overwrite : bool, optional
        If True, overwrite existing file. If False, raise error if file exists.
        By default False.
    
    Raises
    ------
    ValueError
        If file_path does not end with '.nexus' or '.nex'.
    FileExistsError
        If file already exists and overwrite is False.
    OSError
        If the file cannot be written (e.g., permission error).
    
    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.io import save_nexus
    >>> 
    >>> matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
    >>> dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
    >>> save_nexus(dm, 'distances.nexus', overwrite=True)
    >>> 
    >>> # Read back (future functionality)
    >>> # dm2 = load_nexus('distances.nexus')
    
    Notes
    -----
    This function uses `to_nexus()` to generate the NEXUS string and then
    writes it to a file. See `to_nexus()` for details on the NEXUS format.
    """
    # Validate file extension
    if not (file_path.endswith('.nexus') or file_path.endswith('.nex')):
        raise ValueError("File path must end with '.nexus' or '.nex'")
    
    # Check if file exists
    if not overwrite and os.path.exists(file_path):
        raise FileExistsError(
            f"File {file_path} already exists. Use overwrite=True to overwrite it."
        )
    
    # Ensure directory exists
    directory = os.path.dirname(file_path)
    if directory:
        os.makedirs(directory, exist_ok=True)
    
    # Generate NEXUS string and write to file
    nexus_string = to_nexus(distance_matrix)
    with open(file_path, "w") as file:
        file.write(nexus_string)

