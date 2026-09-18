"""
Tests for the text (ASCII-art) drawing of networks.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooTypeError, PhyloZooValueError
from phylozoo.viz import to_ascii
from tests.fixtures import directed_networks, sd_networks
from tests.fixtures.directed_networks import (
    DTREE_LARGE_BINARY,
    LEVEL_1_DNETWORK_SINGLE_HYBRID,
    LEVEL_2_DNETWORK_MULTIPLE_BLOBS,
    LEVEL_3_DNETWORK_LARGE_MANY_HYBRIDS,
)
from tests.fixtures.sd_networks import (
    LEVEL_1_SDNETWORK_TWO_BLOBS,
    SDTREE_EMPTY,
    SDTREE_LARGE_BINARY,
    SDTREE_SINGLE_NODE,
)

DOTTED = set("┄┊")
HEADS = set("><v^")
GLYPHS = set("─│┌┐└┘├┤┬┴┼┄┊><v^○●◆ ")


class TestToAscii:
    """Structural properties of the drawing."""

    def test_snapshot_single_hybrid(self) -> None:
        """Exact drawing of the small single-hybrid example (regression)."""
        net = DirectedPhyNetwork(
            edges=[(5, 3), (5, 4), (3, 2), (4, 2), (2, 1), (3, 6), (4, 7)],
            nodes=[(1, {"label": "A"}), (6, {"label": "B"}), (7, {"label": "C"})],
        )
        expected = "\n".join(
            [
                "┌───────────────●───────● C",
                "│               v",
                "○       ┌┄┄┄┄┄┄>◆───────● A",
                "└───────●",
                "        └───────────────● B",
            ]
        )
        assert to_ascii(net) == expected

    @pytest.mark.parametrize(
        "net",
        [LEVEL_1_DNETWORK_SINGLE_HYBRID, LEVEL_2_DNETWORK_MULTIPLE_BLOBS, DTREE_LARGE_BINARY],
    )
    def test_counts(self, net: DirectedPhyNetwork) -> None:
        """One root, one glyph per node, every leaf label once, one arrowhead per hybrid edge."""
        text = to_ascii(net)
        drawing = text.replace("\n", "")
        assert drawing.count("○") == 1
        assert drawing.count("◆") == len(net.hybrid_nodes)
        assert drawing.count("●") == net.number_of_nodes() - 1 - len(net.hybrid_nodes)
        for leaf in net.leaves:
            assert text.count(" " + net.get_label(leaf)) >= 1
        assert sum(drawing.count(h) for h in HEADS) == len(net.hybrid_edges)

    def test_parallel_edges_stay_visible(self) -> None:
        """Every copy of a parallel edge gets its own arrowhead (drawn as a detour loop)."""
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_PARALLEL_EDGES

        net = LEVEL_1_DNETWORK_PARALLEL_EDGES
        text = to_ascii(net)
        assert sum(text.count(h) for h in HEADS) == len(net.hybrid_edges)
        assert set("┌┐└┘") & set(text)  # the detour has corners

    def test_tree_has_no_dotted_lines(self) -> None:
        """A tree is drawn with solid lines only."""
        text = to_ascii(DTREE_LARGE_BINARY)
        assert not (set(text) & (DOTTED | HEADS))
        assert "┼" not in text

    def test_only_known_glyphs(self) -> None:
        """Apart from labels, the drawing uses the documented character set."""
        net = LEVEL_2_DNETWORK_MULTIPLE_BLOBS
        text = to_ascii(net)
        labels = {net.get_label(leaf) for leaf in net.leaves}
        for line in text.splitlines():
            body = line
            for label in labels:
                body = body.replace(" " + label, "")
            assert set(body) <= GLYPHS, line

    def test_one_row_per_leaf_when_compact(self) -> None:
        """rows_per_leaf=1 gives exactly one line per leaf; 2 gives twice as many minus one."""
        net = DTREE_LARGE_BINARY
        assert len(to_ascii(net, rows_per_leaf=1).splitlines()) == len(net.leaves)
        assert len(to_ascii(net, rows_per_leaf=2).splitlines()) == 2 * len(net.leaves) - 1

    def test_auto_compaction(self) -> None:
        """Large networks get a narrower column width and no blank rows, and fit max_width."""
        net = LEVEL_3_DNETWORK_LARGE_MANY_HYBRIDS  # > 25 leaves
        text = to_ascii(net, max_width=60)
        assert len(text.splitlines()) == len(net.leaves)
        assert max(len(line) for line in text.splitlines()) <= 60 + 3  # col_width floor of 3
        wide = to_ascii(net, max_width=400)
        assert max(len(line) for line in wide.splitlines()) > max(
            len(line) for line in text.splitlines()
        )

    def test_too_many_leaves(self) -> None:
        """Networks above max_leaves are refused with a pointer to plot()."""
        with pytest.raises(PhyloZooValueError, match="max_leaves"):
            to_ascii(DTREE_LARGE_BINARY, max_leaves=5)
        assert to_ascii(DTREE_LARGE_BINARY, max_leaves=len(DTREE_LARGE_BINARY.leaves))

    def test_deterministic(self) -> None:
        """Repeated calls give the same text."""
        assert to_ascii(LEVEL_2_DNETWORK_MULTIPLE_BLOBS) == to_ascii(
            LEVEL_2_DNETWORK_MULTIPLE_BLOBS
        )

    def test_semidirected(self) -> None:
        """Semi-directed networks are rooted and drawn; root_location is honoured."""
        text = to_ascii(LEVEL_1_SDNETWORK_TWO_BLOBS)
        assert text.count("◆") == len(LEVEL_1_SDNETWORK_TWO_BLOBS.hybrid_nodes)
        node = next(iter(LEVEL_1_SDNETWORK_TWO_BLOBS.internal_nodes))
        rooted = to_ascii(LEVEL_1_SDNETWORK_TWO_BLOBS, root_location=node)
        assert rooted.count("○") == 1

    def test_tiny_networks(self) -> None:
        """Empty and single-node networks do not crash."""
        assert to_ascii(SDTREE_EMPTY) == ""
        assert to_ascii(SDTREE_SINGLE_NODE) == "● A"

    def test_methods(self, capsys: pytest.CaptureFixture[str]) -> None:
        """to_ascii() / print_ascii() methods exist on both classes."""
        net = LEVEL_1_DNETWORK_SINGLE_HYBRID
        assert net.to_ascii() == to_ascii(net)
        net.print_ascii()
        assert capsys.readouterr().out.strip() == to_ascii(net)
        assert SDTREE_LARGE_BINARY.to_ascii() == to_ascii(SDTREE_LARGE_BINARY)

    def test_wrong_type(self) -> None:
        """Other objects are rejected."""
        with pytest.raises(PhyloZooTypeError):
            to_ascii("not a network")  # type: ignore[arg-type]

    def test_all_fixtures(self) -> None:
        """Every fixture draws without error and with the right hybrid count."""
        for module in (directed_networks, sd_networks):
            for name in module.NETWORK_METADATA:
                net = getattr(module, name)
                text = to_ascii(net, trials=1)
                if net.number_of_nodes() > 1:
                    assert text.count("◆") == len(net.hybrid_nodes), name

    def test_works_without_matplotlib(self) -> None:
        """The text drawing does not import matplotlib."""
        code = (
            "import sys; sys.modules['matplotlib'] = None\n"
            "from phylozoo import DirectedPhyNetwork\n"
            "from phylozoo.viz import to_ascii\n"
            "net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], "
            "nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])\n"
            "print(to_ascii(net))\n"
        )
        src = str(Path(__file__).resolve().parents[2] / "src")
        env = dict(os.environ, PYTHONPATH=src + os.pathsep + os.environ.get("PYTHONPATH", ""))
        result = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, check=True, env=env
        )
        assert "● A" in result.stdout and "● B" in result.stdout
