"""
Split system I/O module.

This module provides format handlers for reading and writing split systems
to/from files. Format handlers are registered with FormatRegistry for use with
the IOMixin system.

The following format handlers are defined and registered:
- **nexus**: NEXUS format for split systems (extensions: .nexus, .nex, .nxs)
  - Writer: `to_nexus_split_system()` - Converts SplitSystem to NEXUS string
  - Reader: `from_nexus_split_system()` - Parses NEXUS string to SplitSystem
  - Writer: `to_nexus_weighted_split_system()` - Converts WeightedSplitSystem to NEXUS string
  - Reader: `from_nexus_weighted_split_system()` - Parses NEXUS string to WeightedSplitSystem

These handlers are automatically registered when this module is imported.
SplitSystem and WeightedSplitSystem inherit from IOMixin, so you can use:
- `system.save('file.nexus')` - Save to file (auto-detects format)
- `system.load('file.nexus')` - Load from file (auto-detects format)
- `system.to_string(format='nexus')` - Convert to string
- `system.from_string(string, format='nexus')` - Parse from string
- `SplitSystem.convert('in.nexus', 'out.nexus')` - Convert between formats

Notes
-----
The NEXUS format for splits supports weights via the FORMAT WEIGHTS=YES option
in the SPLITS block. WeightedSplitSystem uses this feature, while SplitSystem
writes splits without weights.
"""

from __future__ import annotations

import re
from typing import Any

from ...utils.io import FormatRegistry
from .base import Split
from .splitsystem import SplitSystem
from .weighted_splitsystem import WeightedSplitSystem


def _get_taxon_labels(elements: frozenset) -> list[str]:
    """
    Get sorted list of taxon labels from elements.
    
    Parameters
    ----------
    elements : frozenset
        Set of elements (taxa).
    
    Returns
    -------
    list[str]
        Sorted list of taxon labels as strings.
    """
    return sorted(str(elem) for elem in elements)


def _split_to_nexus_format(split: Split, labels: list[str]) -> tuple[str, str]:
    """
    Convert a split to NEXUS format representation.
    
    Parameters
    ----------
    split : Split
        The split to convert.
    labels : list[str]
        Sorted list of all taxon labels.
    
    Returns
    -------
    tuple[str, str]
        Tuple of (set1_str, set2_str) where each is a space-separated list of labels.
        The smaller set comes first, or lexicographically first if equal size.
    """
    set1_labels = sorted(str(elem) for elem in split.set1)
    set2_labels = sorted(str(elem) for elem in split.set2)
    
    # Ensure canonical ordering (smaller set first, or lexicographically first)
    if len(set1_labels) < len(set2_labels):
        return (' '.join(set1_labels), ' '.join(set2_labels))
    elif len(set2_labels) < len(set1_labels):
        return (' '.join(set2_labels), ' '.join(set1_labels))
    else:
        # Same size, use lexicographic ordering
        if set1_labels < set2_labels:
            return (' '.join(set1_labels), ' '.join(set2_labels))
        else:
            return (' '.join(set2_labels), ' '.join(set1_labels))


def to_nexus_split_system(split_system: SplitSystem, **kwargs: Any) -> str:
    """
    Convert a SplitSystem to a NEXUS format string.
    
    Parameters
    ----------
    split_system : SplitSystem
        The split system to convert.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    str
        The NEXUS format string representation of the split system.
    
    Examples
    --------
    >>> from phylozoo.core.split import Split, SplitSystem
    >>> from phylozoo.core.split.io import to_nexus_split_system
    >>> 
    >>> split1 = Split({1, 2}, {3, 4})
    >>> split2 = Split({1, 3}, {2, 4})
    >>> system = SplitSystem([split1, split2])
    >>> nexus_str = to_nexus_split_system(system)
    >>> print(nexus_str)
    #NEXUS
    <BLANKLINE>
    BEGIN TAXA;
        DIMENSIONS NTAX=4;
        TAXLABELS
            1
            2
            3
            4
        ;
    END;
    <BLANKLINE>
    BEGIN SPLITS;
        DIMENSIONS NSPLITS=2;
        FORMAT LABELS=YES;
        MATRIX
            [1] (1 2) (3 4)
            [2] (1 3) (2 4)
        ;
    END;
    
    Notes
    -----
    The NEXUS format includes:
    - TAXA block with taxon labels
    - SPLITS block with split definitions (no weights for unweighted systems)
    """
    if len(split_system.elements) == 0:
        # Empty split system
        return (
            "#NEXUS\n\nBEGIN TAXA;\n"
            "    DIMENSIONS NTAX=0;\n"
            "    TAXLABELS\n"
            "    ;\n"
            "END;\n\n"
            "BEGIN SPLITS;\n"
            "    DIMENSIONS NSPLITS=0;\n"
            "    FORMAT LABELS=YES;\n"
            "    MATRIX\n"
            "    ;\n"
            "END;\n"
        )
    
    labels = _get_taxon_labels(split_system.elements)
    n = len(labels)
    
    nexus_string = "#NEXUS\n\nBEGIN TAXA;\n"
    nexus_string += f"    DIMENSIONS NTAX={n};\n"
    nexus_string += "    TAXLABELS\n"
    nexus_string += "\n".join(f"        {label}" for label in labels)
    nexus_string += "\n    ;\nEND;\n\n"
    
    nexus_string += "BEGIN SPLITS;\n"
    nexus_string += f"    DIMENSIONS NSPLITS={len(split_system)};\n"
    nexus_string += "    FORMAT LABELS=YES;\n"
    nexus_string += "    MATRIX\n"
    
    # Write each split
    # Sort splits for consistent output
    sorted_splits = sorted(
        split_system.splits,
        key=lambda s: (len(s.set1), sorted(str(e) for e in s.set1))
    )
    
    for i, split in enumerate(sorted_splits, start=1):
        set1_str, set2_str = _split_to_nexus_format(split, labels)
        nexus_string += f"        [{i}] ({set1_str}) ({set2_str})\n"
    
    nexus_string += "    ;\nEND;\n"
    
    return nexus_string


def to_nexus_weighted_split_system(weighted_system: WeightedSplitSystem, **kwargs: Any) -> str:
    """
    Convert a WeightedSplitSystem to a NEXUS format string.
    
    Parameters
    ----------
    weighted_system : WeightedSplitSystem
        The weighted split system to convert.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    str
        The NEXUS format string representation of the weighted split system.
    
    Examples
    --------
    >>> from phylozoo.core.split import Split, WeightedSplitSystem
    >>> from phylozoo.core.split.io import to_nexus_weighted_split_system
    >>> 
    >>> split1 = Split({1, 2}, {3, 4})
    >>> split2 = Split({1, 3}, {2, 4})
    >>> system = WeightedSplitSystem({split1: 0.8, split2: 0.6})
    >>> nexus_str = to_nexus_weighted_split_system(system)
    >>> print(nexus_str)
    #NEXUS
    <BLANKLINE>
    BEGIN TAXA;
        DIMENSIONS NTAX=4;
        TAXLABELS
            1
            2
            3
            4
        ;
    END;
    <BLANKLINE>
    BEGIN SPLITS;
        DIMENSIONS NSPLITS=2;
        FORMAT LABELS=YES WEIGHTS=YES;
        MATRIX
            [1] (1 2) (3 4) 0.800000
            [2] (1 3) (2 4) 0.600000
        ;
    END;
    
    Notes
    -----
    The NEXUS format includes:
    - TAXA block with taxon labels
    - SPLITS block with FORMAT WEIGHTS=YES and split definitions with weights
    """
    if len(weighted_system.elements) == 0:
        # Empty split system
        return (
            "#NEXUS\n\nBEGIN TAXA;\n"
            "    DIMENSIONS NTAX=0;\n"
            "    TAXLABELS\n"
            "    ;\n"
            "END;\n\n"
            "BEGIN SPLITS;\n"
            "    DIMENSIONS NSPLITS=0;\n"
            "    FORMAT LABELS=YES WEIGHTS=YES;\n"
            "    MATRIX\n"
            "    ;\n"
            "END;\n"
        )
    
    labels = _get_taxon_labels(weighted_system.elements)
    n = len(labels)
    
    nexus_string = "#NEXUS\n\nBEGIN TAXA;\n"
    nexus_string += f"    DIMENSIONS NTAX={n};\n"
    nexus_string += "    TAXLABELS\n"
    nexus_string += "\n".join(f"        {label}" for label in labels)
    nexus_string += "\n    ;\nEND;\n\n"
    
    nexus_string += "BEGIN SPLITS;\n"
    nexus_string += f"    DIMENSIONS NSPLITS={len(weighted_system)};\n"
    nexus_string += "    FORMAT LABELS=YES WEIGHTS=YES;\n"
    nexus_string += "    MATRIX\n"
    
    # Write each split with its weight
    # Sort splits for consistent output
    sorted_splits = sorted(
        weighted_system.splits,
        key=lambda s: (len(s.set1), sorted(str(e) for e in s.set1))
    )
    
    for i, split in enumerate(sorted_splits, start=1):
        set1_str, set2_str = _split_to_nexus_format(split, labels)
        weight = weighted_system.get_weight(split)
        nexus_string += f"        [{i}] ({set1_str}) ({set2_str}) {weight:.6f}\n"
    
    nexus_string += "    ;\nEND;\n"
    
    return nexus_string


def from_nexus_split_system(nexus_string: str, **kwargs: Any) -> SplitSystem:
    """
    Parse a NEXUS format string and create a SplitSystem.
    
    Parameters
    ----------
    nexus_string : str
        NEXUS format string containing split system data.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    SplitSystem
        Parsed split system.
    
    Raises
    ------
    ValueError
        If the NEXUS string is malformed or cannot be parsed.
    
    Examples
    --------
    >>> from phylozoo.core.split.io import from_nexus_split_system
    >>> 
    >>> nexus_str = '''#NEXUS
    ... 
    ... BEGIN TAXA;
    ...     DIMENSIONS NTAX=4;
    ...     TAXLABELS
    ...         1
    ...         2
    ...         3
    ...         4
    ...     ;
    ... END;
    ... 
    ... BEGIN SPLITS;
    ...     DIMENSIONS NSPLITS=2;
    ...     FORMAT LABELS=YES;
    ...     MATRIX
    ...         [1] (1 2) (3 4)
    ...         [2] (1 3) (2 4)
    ...     ;
    ... END;'''
    >>> 
    >>> system = from_nexus_split_system(nexus_str)
    >>> len(system)
    2
    
    Notes
    -----
    This parser expects:
    - A TAXA block with TAXLABELS
    - A SPLITS block with FORMAT LABELS=YES (weights optional, ignored if present)
    - Split definitions in format: [n] (taxa1 taxa2 ...) (taxa3 taxa4 ...) [weight]
    """
    # Extract taxa labels
    taxa_match = re.search(
        r'BEGIN\s+TAXA;.*?TAXLABELS\s+(.*?);\s*END;',
        nexus_string,
        re.DOTALL | re.IGNORECASE
    )
    if not taxa_match:
        raise ValueError("Could not find TAXA block with TAXLABELS in NEXUS string")
    
    taxa_section = taxa_match.group(1)
    labels = [line.strip() for line in taxa_section.strip().split('\n') if line.strip()]
    
    if not labels:
        # Empty split system
        return SplitSystem()
    
    # Create mapping from label string to element (try to preserve original type)
    # Try to convert labels back to their original types (int, float, etc.)
    def try_convert(label: str) -> Any:
        """Try to convert label to original type."""
        # Try int first
        try:
            return int(label)
        except ValueError:
            pass
        # Try float
        try:
            return float(label)
        except ValueError:
            pass
        # Keep as string
        return label
    
    label_to_elem = {label: try_convert(label) for label in labels}
    elements_set = set(label_to_elem.values())
    
    # Extract splits block
    splits_match = re.search(
        r'BEGIN\s+SPLITS;.*?MATRIX\s+(.*?);\s*END;',
        nexus_string,
        re.DOTALL | re.IGNORECASE
    )
    if not splits_match:
        raise ValueError("Could not find SPLITS block with MATRIX in NEXUS string")
    
    matrix_section = splits_match.group(1)
    split_lines = [line.strip() for line in matrix_section.strip().split('\n') if line.strip()]
    
    splits: list[Split] = []
    
    # Parse each split line: [n] (taxa1 taxa2 ...) (taxa3 taxa4 ...) [weight]
    # Weight is optional for unweighted systems
    split_pattern = r'\[\d+\]\s+\(([^)]+)\)\s+\(([^)]+)\)(?:\s+([\d.]+))?'
    
    for line in split_lines:
        match = re.match(split_pattern, line)
        if not match:
            raise ValueError(f"Could not parse split line: {line}")
        
        set1_str = match.group(1).strip()
        set2_str = match.group(2).strip()
        # Weight is optional - ignore it for unweighted systems
        
        # Parse taxa in each set
        set1_elems = {label_to_elem[label.strip()] for label in set1_str.split() if label.strip()}
        set2_elems = {label_to_elem[label.strip()] for label in set2_str.split() if label.strip()}
        
        if not set1_elems or not set2_elems:
            raise ValueError(f"Empty set in split: {line}")
        
        # Verify all elements are in the taxa set
        if not set1_elems.issubset(elements_set) or not set2_elems.issubset(elements_set):
            raise ValueError(f"Split contains elements not in taxa: {line}")
        
        # Verify sets are disjoint and cover all elements
        if set1_elems & set2_elems:
            raise ValueError(f"Split sets overlap: {line}")
        
        if set1_elems | set2_elems != elements_set:
            raise ValueError(f"Split does not cover all taxa: {line}")
        
        splits.append(Split(set1_elems, set2_elems))
    
    return SplitSystem(splits)


def from_nexus_weighted_split_system(nexus_string: str, **kwargs: Any) -> WeightedSplitSystem:
    """
    Parse a NEXUS format string and create a WeightedSplitSystem.
    
    Parameters
    ----------
    nexus_string : str
        NEXUS format string containing weighted split system data.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    WeightedSplitSystem
        Parsed weighted split system.
    
    Raises
    ------
    ValueError
        If the NEXUS string is malformed, cannot be parsed, or weights are missing
        when FORMAT WEIGHTS=YES is specified.
    
    Examples
    --------
    >>> from phylozoo.core.split.io import from_nexus_weighted_split_system
    >>> 
    >>> nexus_str = '''#NEXUS
    ... 
    ... BEGIN TAXA;
    ...     DIMENSIONS NTAX=4;
    ...     TAXLABELS
    ...         1
    ...         2
    ...         3
    ...         4
    ...     ;
    ... END;
    ... 
    ... BEGIN SPLITS;
    ...     DIMENSIONS NSPLITS=2;
    ...     FORMAT LABELS=YES WEIGHTS=YES;
    ...     MATRIX
    ...         [1] (1 2) (3 4) 0.8
    ...         [2] (1 3) (2 4) 0.6
    ...     ;
    ... END;'''
    >>> 
    >>> system = from_nexus_weighted_split_system(nexus_str)
    >>> len(system)
    2
    
    Notes
    -----
    This parser expects:
    - A TAXA block with TAXLABELS
    - A SPLITS block with FORMAT LABELS=YES WEIGHTS=YES (or just LABELS=YES if weights optional)
    - Split definitions in format: [n] (taxa1 taxa2 ...) (taxa3 taxa4 ...) weight
    - If FORMAT WEIGHTS=YES is specified, all splits must have weights
    """
    # Extract taxa labels (same as unweighted)
    taxa_match = re.search(
        r'BEGIN\s+TAXA;.*?TAXLABELS\s+(.*?);\s*END;',
        nexus_string,
        re.DOTALL | re.IGNORECASE
    )
    if not taxa_match:
        raise ValueError("Could not find TAXA block with TAXLABELS in NEXUS string")
    
    taxa_section = taxa_match.group(1)
    labels = [line.strip() for line in taxa_section.strip().split('\n') if line.strip()]
    
    if not labels:
        # Empty split system
        return WeightedSplitSystem()
    
    # Create mapping from label string to element (try to preserve original type)
    # Try to convert labels back to their original types (int, float, etc.)
    def try_convert(label: str) -> Any:
        """Try to convert label to original type."""
        # Try int first
        try:
            return int(label)
        except ValueError:
            pass
        # Try float
        try:
            return float(label)
        except ValueError:
            pass
        # Keep as string
        return label
    
    label_to_elem = {label: try_convert(label) for label in labels}
    elements_set = set(label_to_elem.values())
    
    # Extract splits block and check for WEIGHTS=YES
    splits_match = re.search(
        r'BEGIN\s+SPLITS;.*?FORMAT\s+(.*?);.*?MATRIX\s+(.*?);\s*END;',
        nexus_string,
        re.DOTALL | re.IGNORECASE
    )
    if not splits_match:
        raise ValueError("Could not find SPLITS block with MATRIX in NEXUS string")
    
    format_section = splits_match.group(1)
    matrix_section = splits_match.group(2)
    
    # Check if weights are required
    has_weights = 'WEIGHTS=YES' in format_section.upper()
    
    split_lines = [line.strip() for line in matrix_section.strip().split('\n') if line.strip()]
    
    splits: list[Split] = []
    weights: dict[Split, float] = {}
    
    # Parse each split line: [n] (taxa1 taxa2 ...) (taxa3 taxa4 ...) [weight]
    # Use a pattern that captures any text after the second parentheses as the weight
    split_pattern = r'\[\d+\]\s+\(([^)]+)\)\s+\(([^)]+)\)(?:\s+(.+))?$'
    
    for line in split_lines:
        match = re.match(split_pattern, line)
        if not match:
            raise ValueError(f"Could not parse weighted split line: {line}")
        
        set1_str = match.group(1).strip()
        set2_str = match.group(2).strip()
        weight_str = match.group(3)
        
        # If weights are required, they must be present and valid
        if has_weights:
            if weight_str is None:
                raise ValueError(f"Weight required but missing in line: {line}")
            try:
                weight = float(weight_str)
                if weight <= 0:
                    raise ValueError(f"Weight must be positive, got {weight} in line: {line}")
            except ValueError as e:
                # Check if it's our custom error about positive weights
                if "must be positive" in str(e):
                    raise
                # Otherwise it's a parsing error
                raise ValueError(f"Could not parse weight '{weight_str}' in line: {line}") from e
        else:
            # Weights optional - use default weight of 1.0 if not present
            if weight_str is not None:
                try:
                    weight = float(weight_str)
                    if weight <= 0:
                        raise ValueError(f"Weight must be positive, got {weight} in line: {line}")
                except ValueError as e:
                    # Check if it's our custom error about positive weights
                    if "must be positive" in str(e):
                        raise
                    # Otherwise it's a parsing error
                    raise ValueError(f"Could not parse weight '{weight_str}' in line: {line}") from e
            else:
                weight = 1.0
        
        # Parse taxa in each set
        set1_elems = {label_to_elem[label.strip()] for label in set1_str.split() if label.strip()}
        set2_elems = {label_to_elem[label.strip()] for label in set2_str.split() if label.strip()}
        
        if not set1_elems or not set2_elems:
            raise ValueError(f"Empty set in split: {line}")
        
        # Verify all elements are in the taxa set
        if not set1_elems.issubset(elements_set) or not set2_elems.issubset(elements_set):
            raise ValueError(f"Split contains elements not in taxa: {line}")
        
        # Verify sets are disjoint and cover all elements
        if set1_elems & set2_elems:
            raise ValueError(f"Split sets overlap: {line}")
        
        if set1_elems | set2_elems != elements_set:
            raise ValueError(f"Split does not cover all taxa: {line}")
        
        split = Split(set1_elems, set2_elems)
        splits.append(split)
        weights[split] = weight
    
    return WeightedSplitSystem(weights)


# Register format handlers with FormatRegistry
FormatRegistry.register(
    SplitSystem, 'nexus',
    reader=from_nexus_split_system,
    writer=to_nexus_split_system,
    extensions=['.nexus', '.nex', '.nxs'],
    default=True
)

FormatRegistry.register(
    WeightedSplitSystem, 'nexus',
    reader=from_nexus_weighted_split_system,
    writer=to_nexus_weighted_split_system,
    extensions=['.nexus', '.nex', '.nxs'],
    default=True
)

