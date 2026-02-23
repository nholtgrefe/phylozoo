"""
Tests for phylozoo.utils.io.formats.nexus.

Covers parse_nexus, nexus_header, write_taxa_block, write_block.
"""

import pytest

from phylozoo.utils.io.formats import nexus as nexus_fmt
from phylozoo.utils.exceptions import PhyloZooParseError


class TestParseNexus:
    """Tests for parse_nexus."""

    def test_parse_nexus_single_block(self) -> None:
        """Parse NEXUS with TAXA and one Distances block."""
        s = """#NEXUS

BEGIN TAXA;
    DIMENSIONS ntax=3;
    TAXLABELS
        A
        B
        C
    ;
END;

BEGIN Distances;
    DIMENSIONS ntax=3;
    FORMAT triangle=LOWER;
    MATRIX
        A 0
        B 1 0
        C 2 1 0
    ;
END;
"""
        labels, blocks = nexus_fmt.parse_nexus(s)
        assert labels == ['A', 'B', 'C']
        assert 'Distances' in blocks
        assert 'MATRIX' in blocks['Distances']

    def test_parse_nexus_multi_block(self) -> None:
        """Parse NEXUS with TAXA, Distances, and SPLITS blocks."""
        s = """#NEXUS

BEGIN TAXA;
    DIMENSIONS ntax=4;
    TAXLABELS
        1
        2
        3
        4
    ;
END;

BEGIN Distances;
    DIMENSIONS ntax=4;
    FORMAT triangle=LOWER;
    MATRIX
        1 0
        2 1 0
        3 2 1 0
        4 3 2 1 0
    ;
END;

BEGIN SPLITS;
    DIMENSIONS NSPLITS=2;
    FORMAT LABELS=YES;
    MATRIX
        [1] (1 2) (3 4)
        [2] (1 3) (2 4)
    ;
END;
"""
        labels, blocks = nexus_fmt.parse_nexus(s)
        assert labels == ['1', '2', '3', '4']
        assert 'Distances' in blocks
        assert 'SPLITS' in blocks
        assert 'MATRIX' in blocks['Distances']
        assert 'MATRIX' in blocks['SPLITS']

    def test_parse_nexus_no_taxa_raises(self) -> None:
        """parse_nexus raises when no TAXA block is present."""
        s = "#NEXUS\n\nBEGIN Distances;\nMATRIX\n;\nEND;\n"
        with pytest.raises(PhyloZooParseError, match="Taxa block"):
            nexus_fmt.parse_nexus(s)


class TestNexusHeader:
    """Tests for nexus_header."""

    def test_nexus_header(self) -> None:
        """nexus_header returns #NEXUS and blank line."""
        assert nexus_fmt.nexus_header() == "#NEXUS\n\n"


class TestWriteTaxaBlock:
    """Tests for write_taxa_block."""

    def test_write_taxa_block_empty(self) -> None:
        """write_taxa_block with empty labels."""
        out = nexus_fmt.write_taxa_block([])
        assert "BEGIN TAXA" in out
        assert "ntax=0" in out
        assert "END;" in out

    def test_write_taxa_block_labels(self) -> None:
        """write_taxa_block with labels."""
        out = nexus_fmt.write_taxa_block(['A', 'B', 'C'])
        assert "BEGIN TAXA" in out
        assert "ntax=3" in out
        assert "A" in out and "B" in out and "C" in out
        assert "END;" in out


class TestWriteBlock:
    """Tests for write_block."""

    def test_write_block(self) -> None:
        """write_block produces BEGIN name; body END;."""
        body = "    DIMENSIONS ntax=2;\n    MATRIX\n    "
        out = nexus_fmt.write_block("Distances", body)
        assert out.startswith("BEGIN Distances;\n")
        assert "DIMENSIONS ntax=2" in out
        assert "MATRIX" in out
        assert out.rstrip().endswith("END;")
