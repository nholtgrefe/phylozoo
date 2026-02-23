"""
Distance matrix I/O module.

This module provides format handlers for reading and writing distance matrices
to/from files. Format handlers are registered with FormatRegistry for use with
the IOMixin system.

The following format handlers are defined and registered:
- **nexus**: NEXUS format for distance matrices (extensions: .nexus, .nex, .nxs)
  - Writer: `to_nexus()` - Converts DistanceMatrix to NEXUS string
  - Reader: `from_nexus()` - Parses NEXUS string to DistanceMatrix
- **phylip**: PHYLIP format for distance matrices (extensions: .phy, .phylip)
  - Writer: `to_phylip()` - Converts DistanceMatrix to PHYLIP string
  - Reader: `from_phylip()` - Parses PHYLIP string to DistanceMatrix
- **csv**: CSV format for distance matrices (extensions: .csv)
  - Writer: `to_csv()` - Converts DistanceMatrix to CSV string
  - Reader: `from_csv()` - Parses CSV string to DistanceMatrix

These handlers are automatically registered when this module is imported.
DistanceMatrix inherits from IOMixin, so you can use:
- `dm.save('file.nexus')` - Save to file (auto-detects format)
- `dm.load('file.nexus')` - Load from file (auto-detects format)
- `dm.to_string(format='phylip')` - Convert to string
- `dm.from_string(string, format='csv')` - Parse from string
- `DistanceMatrix.convert('in.nexus', 'out.phy')` - Convert between formats
- `DistanceMatrix.convert_string(str1, 'nexus', 'phylip')` - Convert strings
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np

from phylozoo.utils.io import FormatRegistry
from phylozoo.utils.io.formats import nexus as nexus_fmt
from phylozoo.utils.io.formats import phylip as phylip_fmt
from phylozoo.utils.exceptions import PhyloZooParseError, PhyloZooValueError

from .base import DistanceMatrix


def to_nexus(distance_matrix: DistanceMatrix, **kwargs: Any) -> str:
    """
    Convert a distance matrix to a NEXUS format string.
    
    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to convert.
    **kwargs
        Additional arguments:
        - triangle (str): Triangle format - 'LOWER', 'UPPER', or 'BOTH' (default: 'LOWER')
    
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
    >>> '#NEXUS' in nexus_str
    True
    >>> 'triangle=LOWER' in nexus_str
    True
    >>> 
    >>> # Upper triangular format
    >>> nexus_str_upper = to_nexus(dm, triangle='UPPER')
    >>> 'triangle=UPPER' in nexus_str_upper
    True
    
    Notes
    -----
    The NEXUS format includes:
    - Taxa block with label names
    - Distances block with matrix in specified triangle format
    - Format options: triangle=LOWER, triangle=UPPER, or triangle=BOTH
    """
    triangle = kwargs.get('triangle', 'LOWER').upper()
    if triangle not in ('LOWER', 'UPPER', 'BOTH'):
        raise PhyloZooValueError(f"triangle must be 'LOWER', 'UPPER', or 'BOTH', got '{triangle}'")

    n = len(distance_matrix)
    matrix = distance_matrix._matrix
    body = f"    DIMENSIONS ntax={n};\n"
    body += f"    FORMAT triangle={triangle} diagonal LABELS;\n"
    body += "    MATRIX\n"
    if triangle == 'LOWER':
        for i, taxon in enumerate(distance_matrix.labels):
            row = " ".join(f"{matrix[i, j]:.6f}" for j in range(i + 1))
            body += f"    {taxon} {row}\n"
    elif triangle == 'UPPER':
        for i, taxon in enumerate(distance_matrix.labels):
            row = " ".join(f"{matrix[i, j]:.6f}" for j in range(i, n))
            body += f"    {taxon} {row}\n"
    else:  # BOTH
        for i, taxon in enumerate(distance_matrix.labels):
            row = " ".join(f"{matrix[i, j]:.6f}" for j in range(n))
            body += f"    {taxon} {row}\n"
    return (
        nexus_fmt.nexus_header()
        + nexus_fmt.write_taxa_block(distance_matrix.labels)
        + nexus_fmt.write_block("Distances", body)
    )


def from_nexus(nexus_string: str, **kwargs: Any) -> DistanceMatrix:
    """
    Parse a NEXUS format string and create a DistanceMatrix.
    
    Parameters
    ----------
    nexus_string : str
        NEXUS format string containing distance matrix data.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    DistanceMatrix
        Parsed distance matrix.
    
    Raises
    ------
    PhyloZooParseError
        If the NEXUS string is malformed or cannot be parsed (e.g., missing Taxa or Distances blocks,
        mismatched number of taxa and matrix rows, invalid matrix format, invalid distance values).
    
    Examples
    --------
    >>> from phylozoo.core.distance.io import from_nexus
    >>> 
    >>> # Lower triangular format
    >>> nexus_str = '''#NEXUS
    ... 
    ... BEGIN Taxa;
    ...     DIMENSIONS ntax=3;
    ...     TAXLABELS
    ...         A
    ...         B
    ...         C
    ...     ;
    ... END;
    ... 
    ... BEGIN Distances;
    ...     DIMENSIONS ntax=3;
    ...     FORMAT triangle=LOWER diagonal LABELS;
    ...     MATRIX
    ...     A 0.000000
    ...     B 1.000000 0.000000
    ...     C 2.000000 1.000000 0.000000
    ...     ;
    ... END;'''
    >>> 
    >>> dm = from_nexus(nexus_str)
    >>> len(dm)
    3
    >>> dm.get_distance('A', 'B')
    1.0
    
    Notes
    -----
    This parser supports:
    - A Taxa block with TAXLABELS
    - A Distances block with FORMAT triangle=LOWER/UPPER/BOTH diagonal LABELS
    - Lower triangular, upper triangular, or full matrix formats
    """
    labels, blocks = nexus_fmt.parse_nexus(nexus_string)
    content = blocks.get("Distances")
    if content is None:
        raise PhyloZooParseError(
            f"NEXUS file contains no Distances block (found: {list(blocks.keys())})"
        )

    n = len(labels)
    if not labels:
        raise PhyloZooParseError("No taxa labels found in NEXUS string")

    # Extract FORMAT line to determine triangle type from block content
    format_match = re.search(r'FORMAT\s+(.*?);', content, re.IGNORECASE)
    triangle = 'LOWER'
    if format_match:
        format_section = format_match.group(1).upper()
        if 'TRIANGLE=UPPER' in format_section:
            triangle = 'UPPER'
        elif 'TRIANGLE=BOTH' in format_section:
            triangle = 'BOTH'
        elif 'TRIANGLE=LOWER' in format_section:
            triangle = 'LOWER'

    # Extract MATRIX section from block content
    matrix_match = re.search(r'MATRIX\s+(.*?);', content, re.DOTALL | re.IGNORECASE)
    if not matrix_match:
        raise PhyloZooParseError("Could not find MATRIX in Distances block")
    matrix_section = matrix_match.group(1)
    matrix_lines = [line.strip() for line in matrix_section.strip().split('\n') if line.strip()]
    
    if len(matrix_lines) != n:
        raise PhyloZooParseError(
            f"Number of matrix rows ({len(matrix_lines)}) does not match "
            f"number of taxa ({n})"
        )
    
    # Parse matrix based on triangle format
    matrix = np.zeros((n, n), dtype=np.float64)
    
    for i, line in enumerate(matrix_lines):
        parts = line.split()
        
        # First part is the label (should match labels[i])
        label = parts[0]
        if label != labels[i]:
            raise PhyloZooParseError(
                f"Matrix row {i+1} label '{label}' does not match taxa label '{labels[i]}'"
            )
        
        if triangle == 'LOWER':
            # Lower triangular: row i has values for columns 0..i
            if len(parts) < i + 2:  # label + (i+1) values
                raise PhyloZooParseError(
                    f"Matrix row {i+1} has insufficient values. "
                    f"Expected {i+2} values (label + {i+1} distances), got {len(parts)}"
                )
            
            for j in range(i + 1):
                try:
                    value = float(parts[j + 1])
                    matrix[i, j] = value
                    matrix[j, i] = value  # Symmetric
                except (ValueError, IndexError) as e:
                    raise PhyloZooParseError(
                        f"Could not parse distance value at row {i+1}, column {j+1}: {e}"
                    ) from e
        
        elif triangle == 'UPPER':
            # Upper triangular: row i has values for columns i..n-1
            expected_values = n - i
            if len(parts) < expected_values + 1:  # label + (n-i) values
                raise PhyloZooParseError(
                    f"Matrix row {i+1} has insufficient values. "
                    f"Expected {expected_values + 1} values (label + {expected_values} distances), got {len(parts)}"
                )
            
            for j in range(i, n):
                try:
                    value = float(parts[j - i + 1])  # Offset by 1 for label
                    matrix[i, j] = value
                    matrix[j, i] = value  # Symmetric
                except (ValueError, IndexError) as e:
                    raise PhyloZooParseError(
                        f"Could not parse distance value at row {i+1}, column {j+1}: {e}"
                    ) from e
        
        else:  # BOTH
            # Full matrix: row i has all n values
            if len(parts) < n + 1:  # label + n values
                raise PhyloZooParseError(
                    f"Matrix row {i+1} has insufficient values. "
                    f"Expected {n + 1} values (label + {n} distances), got {len(parts)}"
                )
            
            for j in range(n):
                try:
                    value = float(parts[j + 1])
                    matrix[i, j] = value
                except (ValueError, IndexError) as e:
                    raise PhyloZooParseError(
                        f"Could not parse distance value at row {i+1}, column {j+1}: {e}"
                    ) from e
    
    # Create DistanceMatrix
    return DistanceMatrix(matrix, labels=labels)


def to_phylip(distance_matrix: DistanceMatrix, **kwargs: Any) -> str:
    """
    Convert a distance matrix to PHYLIP format string.
    
    PHYLIP format consists of:
    - First line: number of taxa
    - Subsequent lines: taxon name (padded to 10 chars) followed by all distances
    
    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to convert.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    str
        The PHYLIP format string representation of the distance matrix.
    
    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.io import to_phylip
    >>> 
    >>> matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
    >>> dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
    >>> phylip_str = to_phylip(dm)
    >>> print(phylip_str[:30])
    3
    A          0.00000 1.00000
    
    Notes
    -----
    Taxon names are padded to 10 characters (standard PHYLIP format).
    Distances are formatted with 5 decimal places.
    """
    n = len(distance_matrix)
    matrix = distance_matrix._matrix
    rows = [
        (str(taxon), ' '.join(f"{matrix[i, j]:.5f}" for j in range(n)))
        for i, taxon in enumerate(distance_matrix.labels)
    ]
    return phylip_fmt.write_phylip_matrix(n, rows)


def from_phylip(phylip_string: str, **kwargs: Any) -> DistanceMatrix:
    """
    Parse a PHYLIP format string and create a DistanceMatrix.
    
    Parameters
    ----------
    phylip_string : str
        PHYLIP format string containing distance matrix data.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    DistanceMatrix
        Parsed distance matrix.
    
    Raises
    ------
    PhyloZooParseError
        If the PHYLIP string is malformed or cannot be parsed (e.g., empty string, invalid number of taxa,
        insufficient values, invalid distance values, non-symmetric matrix).
    PhyloZooValueError
        If the number of taxa is not positive.
    
    Examples
    --------
    >>> from phylozoo.core.distance.io import from_phylip
    >>> 
    >>> phylip_str = '''3
    ... A          0.00000 1.00000 2.00000
    ... B          1.00000 0.00000 1.00000
    ... C          2.00000 1.00000 0.00000
    ... '''
    >>> 
    >>> dm = from_phylip(phylip_str)
    >>> len(dm)
    3
    >>> dm.get_distance('A', 'B')
    1.0
    
    Notes
    -----
    This parser expects:
    - First line: number of taxa
    - Subsequent lines: taxon name (first 10 chars or until whitespace) followed by distances
    - Full matrix format (not just lower triangle)
    """
    n, rows = phylip_fmt.parse_phylip_matrix(phylip_string)
    labels = [r[0] for r in rows]
    matrix = np.zeros((n, n), dtype=np.float64)
    for i, (_, rest) in enumerate(rows):
        parts = rest.split()
        if len(parts) < n:
            raise PhyloZooParseError(
                f"Line {i + 2} has insufficient values. "
                f"Expected {n} distances, got {len(parts)}"
            )
        for j in range(n):
            try:
                value = float(parts[j])
                matrix[i, j] = value
            except (ValueError, IndexError) as e:
                raise PhyloZooParseError(
                    f"Could not parse distance value at row {i + 1}, column {j + 1}: {e}"
                ) from e
    if not np.allclose(matrix, matrix.T, rtol=1e-10, atol=1e-10):
        raise PhyloZooParseError("PHYLIP matrix is not symmetric")
    return DistanceMatrix(matrix, labels=labels)


def to_csv(distance_matrix: DistanceMatrix, **kwargs: Any) -> str:
    """
    Convert a distance matrix to CSV format string.
    
    CSV format consists of:
    - First row: header with empty first cell, then taxon labels
    - Subsequent rows: taxon label in first column, then distances
    
    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to convert.
    **kwargs
        Additional arguments:
        - delimiter (str): Field delimiter (default: ',')
        - include_header (bool): Include header row (default: True)
    
    Returns
    -------
    str
        The CSV format string representation of the distance matrix.
    
    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.io import to_csv
    >>> 
    >>> matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
    >>> dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
    >>> csv_str = to_csv(dm)
    >>> print(csv_str[:20])
    ,A,B,C
    A,0.0,1.0
    
    Notes
    -----
    Default delimiter is comma. Use delimiter='\t' for tab-separated values.
    """
    delimiter = kwargs.get('delimiter', ',')
    include_header = kwargs.get('include_header', True)
    
    n = len(distance_matrix)
    csv_lines = []
    
    matrix = distance_matrix._matrix
    
    # Header row
    if include_header:
        header = delimiter.join([''] + [str(label) for label in distance_matrix.labels])
        csv_lines.append(header)
    
    # Data rows
    for i, taxon in enumerate(distance_matrix.labels):
        row_values = [str(taxon)] + [f"{matrix[i, j]:.6f}" for j in range(n)]
        csv_lines.append(delimiter.join(row_values))
    
    return '\n'.join(csv_lines) + '\n'


def from_csv(csv_string: str, **kwargs: Any) -> DistanceMatrix:
    """
    Parse a CSV format string and create a DistanceMatrix.
    
    Parameters
    ----------
    csv_string : str
        CSV format string containing distance matrix data.
    **kwargs
        Additional arguments:
        - delimiter (str): Field delimiter (default: ','). Can be ',' or '\t' or whitespace
        - has_header (bool): Whether first row is a header (default: True)
    
    Returns
    -------
    DistanceMatrix
        Parsed distance matrix.
    
    Raises
    ------
    PhyloZooParseError
        If the CSV string is malformed or cannot be parsed (e.g., empty string, no data rows,
        invalid distance values, mismatched dimensions, non-symmetric matrix).
    
    Examples
    --------
    >>> from phylozoo.core.distance.io import from_csv
    >>> 
    >>> csv_str = ''',A,B,C
    ... A,0.0,1.0,2.0
    ... B,1.0,0.0,1.0
    ... C,2.0,1.0,0.0
    ... '''
    >>> 
    >>> dm = from_csv(csv_str)
    >>> len(dm)
    3
    >>> dm.get_distance('A', 'B')
    1.0
    
    Notes
    -----
    This parser expects:
    - First row (if has_header=True): empty first cell, then taxon labels
    - Subsequent rows: taxon label in first column, then distances
    - Delimiter can be comma, tab, or whitespace
    """
    delimiter = kwargs.get('delimiter', ',')
    has_header = kwargs.get('has_header', True)
    
    lines = [line.strip() for line in csv_string.strip().split('\n') if line.strip()]
    
    if not lines:
        raise PhyloZooParseError("CSV string is empty")
    
    # Parse header if present
    start_idx = 0
    labels_from_header = None
    
    if has_header:
        if not lines:
            raise PhyloZooParseError("CSV string has header flag but is empty")
        
        header_parts = [p.strip() for p in lines[0].split(delimiter)]
        if header_parts[0] == '':  # First cell is empty (standard CSV format)
            labels_from_header = [label for label in header_parts[1:] if label]
            start_idx = 1
        else:
            # First cell is not empty - might be no header or wrong delimiter
            # Try to auto-detect: if first value is numeric, it's probably data
            try:
                float(header_parts[1])
                has_header = False  # First line is actually data
                start_idx = 0
            except (ValueError, IndexError):
                # Assume it's a header without empty first cell
                labels_from_header = [label for label in header_parts if label]
                start_idx = 1
    else:
        start_idx = 0
    
    if start_idx >= len(lines):
        raise PhyloZooParseError("No data rows found in CSV string")
    
    # Parse data rows
    labels: list[str] = []
    matrix_rows: list[list[float]] = []
    
    for line in lines[start_idx:]:
        parts = [p.strip() for p in line.split(delimiter) if p.strip()]
        if not parts:
            continue
        
        taxon = parts[0]
        labels.append(taxon)
        
        try:
            distances = [float(part) for part in parts[1:]]
            matrix_rows.append(distances)
        except ValueError as e:
            raise PhyloZooParseError(f"Could not parse distances for taxon '{taxon}': {e}") from e
    
    if not labels:
        raise PhyloZooParseError("No taxa found in CSV string")
    
    n = len(labels)
    
    # Verify dimensions
    if labels_from_header:
        if len(labels_from_header) != n:
            raise PhyloZooParseError(
                f"Header has {len(labels_from_header)} labels but data has {n} rows"
            )
        # Use labels from header (they should match, but header is authoritative)
        labels = labels_from_header
    
    # Check that all rows have the same number of distances
    expected_distances = n
    for i, row in enumerate(matrix_rows):
        if len(row) != expected_distances:
            raise PhyloZooParseError(
                f"Row {i + 1} (taxon '{labels[i]}') has {len(row)} distances, "
                f"expected {expected_distances}"
            )
    
    # Build matrix
    matrix = np.zeros((n, n), dtype=np.float64)
    for i, row in enumerate(matrix_rows):
        for j, dist in enumerate(row):
            matrix[i, j] = dist
    
    # Verify symmetry
    if not np.allclose(matrix, matrix.T, rtol=1e-10, atol=1e-10):
        raise PhyloZooParseError("CSV matrix is not symmetric")
    
    # Create DistanceMatrix
    return DistanceMatrix(matrix, labels=labels)


# Register format handlers with FormatRegistry
# Register format handlers with FormatRegistry
FormatRegistry.register(
    DistanceMatrix, 'nexus',
    reader=from_nexus,
    writer=to_nexus,
    extensions=['.nexus', '.nex', '.nxs'],
    default=True
)

FormatRegistry.register(
    DistanceMatrix, 'phylip',
    reader=from_phylip,
    writer=to_phylip,
    extensions=['.phy', '.phylip']
)

FormatRegistry.register(
    DistanceMatrix, 'csv',
    reader=from_csv,
    writer=to_csv,
    extensions=['.csv']
)

