"""
Shared PHYLIP matrix format structure.

Generic layout: first line = n, then n lines of (label, rest).
Class-specific parsing (e.g. float matrix, symmetry) stays in core/distance/io.
"""

from __future__ import annotations

from phylozoo.utils.exceptions import PhyloZooParseError

# PHYLIP: label is first 10 characters (or until whitespace)
_PHYLIP_LABEL_LEN = 10


def parse_phylip_matrix(phylip_string: str) -> tuple[int, list[tuple[str, str]]]:
    """
    Parse PHYLIP matrix layout into n and rows (label, rest of line).

    Parameters
    ----------
    phylip_string : str
        Full PHYLIP file content.

    Returns
    -------
    tuple[int, list[tuple[str, str]]]
        (n, rows). n = number of taxa; rows = list of (label, rest) per line.
        Label is first 10 chars or until whitespace; rest is the remainder.

    Raises
    ------
    PhyloZooParseError
        If string is empty or first line is not an integer.
    """
    lines = [line.rstrip('\n\r') for line in phylip_string.strip().split('\n') if line.strip()]

    if not lines:
        raise PhyloZooParseError("PHYLIP string is empty")

    try:
        n = int(lines[0])
    except ValueError as e:
        raise PhyloZooParseError(
            f"Could not parse number of taxa from first line: {e}"
        ) from e

    if n <= 0:
        raise PhyloZooParseError(
            f"Number of taxa must be positive, got {n}"
        )

    if len(lines) < n + 1:
        raise PhyloZooParseError(
            f"PHYLIP string has {len(lines) - 1} data lines, expected {n}"
        )

    rows: list[tuple[str, str]] = []
    for i in range(1, n + 1):
        line = lines[i]
        # Label: first 10 characters (PHYLIP standard) or until first whitespace
        if len(line) <= _PHYLIP_LABEL_LEN:
            label = line.strip()
            rest = ''
        else:
            label = line[:_PHYLIP_LABEL_LEN].strip()
            rest = line[_PHYLIP_LABEL_LEN:].strip()
        rows.append((label, rest))

    return n, rows


def write_phylip_matrix(n: int, rows: list[tuple[str, str]]) -> str:
    """
    Build PHYLIP matrix string (first line n, then n lines of label + rest).

    Parameters
    ----------
    n : int
        Number of taxa.
    rows : list[tuple[str, str]]
        Each element is (label, rest) for one row. Label is padded to 10 chars.

    Returns
    -------
    str
        Full PHYLIP matrix content (no trailing newline required by callers).
    """
    lines = [str(n)]
    for label, rest in rows:
        label_padded = str(label).ljust(_PHYLIP_LABEL_LEN)
        lines.append(f"{label_padded}{rest}")
    return '\n'.join(lines) + '\n'
