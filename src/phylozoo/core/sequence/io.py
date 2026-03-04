"""
MSA I/O module.

This module provides format handlers for reading and writing MSAs
to/from files. Format handlers are registered with FormatRegistry for use with
the IOMixin system.

The following format handlers are defined and registered:

- **fasta**: FASTA format for sequence alignments (extensions: .fasta, .fa, .fas)

  - Writer: `to_fasta()` - Converts MSA to FASTA string
  - Reader: `from_fasta()` - Parses FASTA string to MSA

- **nexus**: NEXUS format for sequence alignments (extensions: .nexus, .nex, .nxs)

  - Writer: `to_nexus()` - Converts MSA to NEXUS string
  - Reader: `from_nexus()` - Parses NEXUS string to MSA

These handlers are automatically registered when this module is imported.
MSA inherits from IOMixin, so you can use:

- `msa.save('file.fasta')` - Save to file (auto-detects format)
- `msa.load('file.fasta')` - Load from file (auto-detects format)
- `msa.to_string(format='fasta')` - Convert to string
- `MSA.from_string(string, format='fasta')` - Parse from string
- `MSA.convert('in.fasta', 'out.nexus')` - Convert between formats
"""

from __future__ import annotations

import re
from typing import Any

from phylozoo.utils.io import FormatRegistry
from phylozoo.utils.io.formats import nexus as nexus_fmt
from phylozoo.utils.exceptions import PhyloZooParseError

from .base import MSA


def to_fasta(msa: MSA, **kwargs: Any) -> str:
    """
    Convert an MSA to FASTA format string.
    
    FASTA format consists of:

    - Each sequence starts with a '>' followed by the taxon identifier
    - The sequence follows on subsequent lines (can be wrapped)
    - Sequences are separated by newlines
    
    Parameters
    ----------
    msa : MSA
        The MSA to convert.
    **kwargs
        Additional arguments:

        - line_length (int): Maximum line length for sequences (default: 80).
          Set to 0 or None for no wrapping.
    
    Returns
    -------
    str
        The FASTA format string representation of the MSA.
    
    Examples
    --------
    >>> from phylozoo.core.sequence import MSA
    >>> from phylozoo.core.sequence.io import to_fasta
    >>> 
    >>> sequences = {"taxon1": "ACGTACGT", "taxon2": "TGCAACGT"}
    >>> msa = MSA(sequences)
    >>> fasta_str = to_fasta(msa)
    >>> '>taxon1' in fasta_str
    True
    >>> 'ACGTACGT' in fasta_str
    True
    
    Notes
    -----
    Taxon identifiers are converted to strings for output.
    Sequences are wrapped to 80 characters per line by default.
    """
    line_length = kwargs.get('line_length', 80)
    if line_length is None or line_length <= 0:
        line_length = None  # No wrapping
    
    fasta_lines = []
    
    for taxon in msa.taxa_order:
        # Header line: >taxon_identifier
        fasta_lines.append(f">{taxon}")
        
        # Get sequence
        seq = msa.get_sequence(taxon)
        if seq is None:
            seq = ''
        
        # Wrap sequence if line_length is specified
        if line_length and len(seq) > line_length:
            for i in range(0, len(seq), line_length):
                fasta_lines.append(seq[i:i + line_length])
        else:
            fasta_lines.append(seq)
    
    return '\n'.join(fasta_lines) + '\n'


def from_fasta(fasta_string: str, **kwargs: Any) -> MSA:
    """
    Parse a FASTA format string and create an MSA.
    
    Parameters
    ----------
    fasta_string : str
        FASTA format string containing sequence alignment data.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    MSA
        Parsed MSA.
    
    Raises
    ------
    PhyloZooParseError
        If the FASTA string is malformed or cannot be parsed (e.g., empty taxon identifier,
        sequence data before header, no sequences found).
    
    Examples
    --------
    >>> from phylozoo.core.sequence.io import from_fasta
    >>> fasta_str = ">taxon1\\nACGTACGT\\n>taxon2\\nTGCAACGT\\n"
    >>> msa = from_fasta(fasta_str)
    >>> len(msa)
    2
    >>> msa.get_sequence("taxon1")
    'ACGTACGT'
    
    Notes
    -----
    This parser:

    - Handles multi-line sequences (concatenates lines between headers)
    - Strips whitespace from sequence lines
    - Converts taxon identifiers to strings (FASTA format uses strings)
    - Raises error if sequences have different lengths
    """
    sequences: dict[str, str] = {}
    current_taxon: str | None = None
    current_seq: list[str] = []
    
    lines = fasta_string.split('\n')
    
    for line in lines:
        line = line.strip()
        
        if not line:
            continue
        
        # Check if this is a header line (starts with '>')
        if line.startswith('>'):
            # Save previous sequence if any
            if current_taxon is not None:
                seq_str = ''.join(current_seq)
                if seq_str:  # Only add if sequence is not empty
                    sequences[current_taxon] = seq_str
                current_seq = []
            
            # Extract taxon identifier (everything after '>', stripped)
            current_taxon = line[1:].strip()
            if not current_taxon:
                raise PhyloZooParseError("FASTA header line has empty taxon identifier")
        
        elif current_taxon is not None:
            # This is a sequence line
            # Remove whitespace and add to current sequence
            seq_line = line.replace(' ', '').replace('\t', '')
            if seq_line:
                current_seq.append(seq_line)
        else:
            # Sequence line before any header
            raise PhyloZooParseError("FASTA string has sequence data before first header line")
    
    # Don't forget the last sequence
    if current_taxon is not None:
        seq_str = ''.join(current_seq)
        if seq_str:
            sequences[current_taxon] = seq_str
    
    if not sequences:
        raise PhyloZooParseError("No sequences found in FASTA string")
    
    # Create MSA (will validate that all sequences have same length)
    return MSA(sequences)


def to_nexus(msa: MSA, **kwargs: Any) -> str:
    """
    Convert an MSA to NEXUS format string.
    
    NEXUS format for sequences consists of:

    - #NEXUS header
    - TAXA block with taxon labels
    - CHARACTERS block with sequence data
    
    Parameters
    ----------
    msa : MSA
        The MSA to convert.
    **kwargs
        Additional arguments:

        - datatype (str): Data type (default: 'DNA'). Can be 'DNA', 'RNA', 'PROTEIN', etc.
        - missing (str): Missing data character (default: 'N')
        - gap (str): Gap character (default: '-')
    
    Returns
    -------
    str
        The NEXUS format string representation of the MSA.
    
    Examples
    --------
    >>> from phylozoo.core.sequence import MSA
    >>> from phylozoo.core.sequence.io import to_nexus
    >>> 
    >>> sequences = {"taxon1": "ACGTACGT", "taxon2": "TGCAACGT"}
    >>> msa = MSA(sequences)
    >>> nexus_str = to_nexus(msa)
    >>> '#NEXUS' in nexus_str
    True
    >>> 'BEGIN TAXA' in nexus_str
    True
    >>> 'BEGIN CHARACTERS' in nexus_str
    True
    
    Notes
    -----
    The NEXUS format includes:

    - Taxa block with DIMENSIONS and TAXLABELS
    - Characters block with DIMENSIONS, FORMAT, and MATRIX
    - Taxon identifiers are converted to strings
    """
    datatype = kwargs.get('datatype', 'DNA').upper()
    missing = kwargs.get('missing', 'N')
    gap = kwargs.get('gap', '-')

    n = msa.num_taxa
    seq_length = msa.sequence_length
    body = f"    DIMENSIONS nchar={seq_length};\n"
    body += f"    FORMAT datatype={datatype} missing={missing} gap={gap};\n"
    body += "    MATRIX\n"
    for taxon in msa.taxa_order:
        seq = msa.get_sequence(taxon) or ''
        body += f"        {taxon}    {seq}\n"
    return (
        nexus_fmt.nexus_header()
        + nexus_fmt.write_taxa_block(msa.taxa_order)
        + nexus_fmt.write_block("CHARACTERS", body)
    )


def from_nexus(nexus_string: str, **kwargs: Any) -> MSA:
    """
    Parse a NEXUS format string and create an MSA.
    
    Parameters
    ----------
    nexus_string : str
        NEXUS format string containing sequence alignment data.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    MSA
        Parsed MSA.
    
    Raises
    ------
    PhyloZooParseError
        If the NEXUS string is malformed or cannot be parsed (e.g., missing TAXA or CHARACTERS blocks,
        mismatched number of taxa and matrix rows, invalid matrix format).
    
    Examples
    --------
    >>> from phylozoo.core.sequence.io import from_nexus
    >>> 
    >>> nexus_str = '''#NEXUS
    ... 
    ... BEGIN TAXA;
    ...     DIMENSIONS ntax=2;
    ...     TAXLABELS
    ...         taxon1
    ...         taxon2
    ...     ;
    ... END;
    ... 
    ... BEGIN CHARACTERS;
    ...     DIMENSIONS nchar=8;
    ...     FORMAT datatype=DNA missing=N gap=-;
    ...     MATRIX
    ...         taxon1    ACGTACGT
    ...         taxon2    TGCAACGT
    ...     ;
    ... END;'''
    >>> 
    >>> msa = from_nexus(nexus_str)
    >>> len(msa)
    2
    >>> msa.get_sequence("taxon1")
    'ACGTACGT'
    
    Notes
    -----
    This parser supports:

    - TAXA block with TAXLABELS
    - CHARACTERS block with MATRIX
    - Handles missing and gap characters as specified in FORMAT
    - Converts taxon identifiers to strings
    """
    labels, blocks = nexus_fmt.parse_nexus(nexus_string)
    content = blocks.get("CHARACTERS")
    if content is None:
        raise PhyloZooParseError(
            f"NEXUS file contains no CHARACTERS block (found: {list(blocks.keys())})"
        )

    n = len(labels)
    if not labels:
        raise PhyloZooParseError("No taxa labels found in NEXUS string")

    matrix_match = re.search(r'MATRIX\s+(.*?);', content, re.DOTALL | re.IGNORECASE)
    if not matrix_match:
        raise PhyloZooParseError("Could not find MATRIX in CHARACTERS block")
    matrix_section = matrix_match.group(1)
    matrix_lines = [line.strip() for line in matrix_section.strip().split('\n') if line.strip()]
    
    if len(matrix_lines) != n:
        raise PhyloZooParseError(
            f"Number of matrix rows ({len(matrix_lines)}) does not match "
            f"number of taxa ({n})"
        )
    
    # Parse sequences
    sequences: dict[str, str] = {}
    
    for i, line in enumerate(matrix_lines):
        # Split line: taxon identifier and sequence
        # Taxon identifier is typically separated by whitespace from sequence
        parts = line.split(None, 1)  # Split on whitespace, max 1 split
        
        if len(parts) < 2:
            raise PhyloZooParseError(
                f"Matrix row {i+1} does not have both taxon identifier and sequence"
            )
        
        taxon = parts[0].strip()
        seq = parts[1].strip().replace(' ', '').replace('\t', '')  # Remove whitespace
        
        # Verify taxon matches expected label
        if taxon != labels[i]:
            raise PhyloZooParseError(
                f"Matrix row {i+1} taxon '{taxon}' does not match taxa label '{labels[i]}'"
            )
        
        sequences[taxon] = seq
    
    # Create MSA (will validate that all sequences have same length)
    return MSA(sequences)


# Register format handlers with FormatRegistry
FormatRegistry.register(
    MSA, 'fasta',
    reader=from_fasta,
    writer=to_fasta,
    extensions=['.fasta', '.fa', '.fas'],
    default=True
)

FormatRegistry.register(
    MSA, 'nexus',
    reader=from_nexus,
    writer=to_nexus,
    extensions=['.nexus', '.nex', '.nxs']
)

