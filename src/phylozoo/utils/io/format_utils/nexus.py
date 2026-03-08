"""
Shared NEXUS format structure.

Provides parse_nexus (labels + all blocks), write helpers (header, TAXA block, block).
Class-specific parsing of block content stays in core/distance/io, core/sequence/io, core/split/io.
"""

from __future__ import annotations

import re

from phylozoo.utils.exceptions import PhyloZooParseError

# Regex for TAXA block: BEGIN TAXA; ... TAXLABELS ... ; END;
_RE_TAXA = re.compile(
    r'BEGIN\s+Taxa;.*?TAXLABELS\s+(.*?);\s*END;',
    re.DOTALL | re.IGNORECASE,
)
# Regex for any block: BEGIN Name; content ; END;
_RE_BLOCK = re.compile(
    r'BEGIN\s+(\w+);\s*([\s\S]*?);\s*END;',
    re.IGNORECASE,
)


def parse_nexus(nexus_string: str) -> tuple[list[str], dict[str, str]]:
    """
    Parse a NEXUS string into labels (from TAXA block) and data blocks.

    Parameters
    ----------
    nexus_string : str
        Full NEXUS file content.

    Returns
    -------
    tuple[list[str], dict[str, str]]
        (labels, blocks). labels from TAXA block; blocks maps canonical block name
        (Distances, CHARACTERS, SPLITS) to full block content (including trailing ";").

    Raises
    ------
    PhyloZooParseError
        If no TAXA block with TAXLABELS is found.

    Notes
    -----
    A file may contain multiple data blocks (e.g. Distances and SPLITS).
    Callers use only the block they need (e.g. DistanceMatrix uses "Distances").
    """
    # TAXA block
    taxa_match = _RE_TAXA.search(nexus_string)
    if not taxa_match:
        raise PhyloZooParseError(
            "Could not find Taxa block with TAXLABELS in NEXUS string"
        )
    taxa_section = taxa_match.group(1)
    labels = [line.strip() for line in taxa_section.strip().split('\n') if line.strip()]

    # Data blocks (exclude TAXA)
    blocks: dict[str, str] = {}
    for match in _RE_BLOCK.finditer(nexus_string):
        name = match.group(1)
        if name.upper() == 'TAXA':
            continue
        # Canonical names: Distances, CHARACTERS, SPLITS
        canonical = (
            name.upper()
            if name.upper() in ('CHARACTERS', 'SPLITS')
            else name.capitalize()
        )
        # Include trailing ";" so block content can be parsed with MATRIX\s+(.*?);
        body = match.group(2).strip()
        if not body.endswith(';'):
            body += ';'
        blocks[canonical] = body

    return labels, blocks


def nexus_header() -> str:
    """
    Return the NEXUS file header (#NEXUS and blank line).

    Returns
    -------
    str
        The NEXUS file header.
    """

    return "#NEXUS\n\n"


def write_taxa_block(labels: list[str]) -> str:
    """
    Build the TAXA block string (BEGIN TAXA; DIMENSIONS ntax=N; TAXLABELS ... ; END;).

    Parameters
    ----------
    labels : list[str]
        Taxon labels.

    Returns
    -------
    str
        Full TAXA block including BEGIN/END.
    """
    n = len(labels)
    out = "BEGIN TAXA;\n"
    out += f"    DIMENSIONS ntax={n};\n"
    out += "    TAXLABELS\n"
    out += "\n".join(f"        {lab}" for lab in labels)
    out += "\n    ;\nEND;\n\n"
    return out


def write_block(block_name: str, body: str) -> str:
    """
    Build a NEXUS data block (BEGIN name; body END;).

    Parameters
    ----------
    block_name : str
        Block name (e.g. Distances, CHARACTERS, SPLITS).
    body : str
        Block body (commands and MATRIX etc.); need not include trailing ";".

    Returns
    -------
    str
        Full block including BEGIN/END.
    """
    if body and not body.rstrip().endswith(';'):
        body = body.rstrip() + "\n    ;\n"
    else:
        body = body.rstrip() + "\n"
    return f"BEGIN {block_name};\n{body}END;\n"
