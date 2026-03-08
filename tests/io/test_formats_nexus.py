"""
Tests for phylozoo.utils.io.format_utils.nexus.

Covers parse_nexus, nexus_header, write_taxa_block, write_block.
"""

import pytest

from phylozoo.utils.io.format_utils import nexus as nexus_fmt
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

    @pytest.mark.parametrize(
        "labels,expected_ntax,expected_substrings",
        [
            ([], "0", ["BEGIN TAXA", "ntax=0", "END;"]),
            (
                ['A', 'B', 'C'],
                "3",
                ["BEGIN TAXA", "ntax=3", "A", "B", "C", "END;"],
            ),
        ],
    )
    def test_write_taxa_block(
        self,
        labels: list[str],
        expected_ntax: str,
        expected_substrings: list[str],
    ) -> None:
        """write_taxa_block with empty or non-empty labels."""
        out = nexus_fmt.write_taxa_block(labels)
        assert f"ntax={expected_ntax}" in out
        for sub in expected_substrings:
            assert sub in out


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
