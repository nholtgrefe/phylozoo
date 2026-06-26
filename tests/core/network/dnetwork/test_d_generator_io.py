"""I/O tests for :class:`DirectedGenerator` (IOMixin via its DirectedMultiGraph)."""

import os
import tempfile

import pytest

from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
from phylozoo.core.network.dnetwork.generator.base import DirectedGenerator
from phylozoo.utils.io import IOMixin


def _gen() -> DirectedGenerator:
    """A level-1 generator (root with parallel edges to a hybrid node)."""
    return DirectedGenerator(DirectedMultiGraph(edges=[(8, 4), (8, 4)]))


def _sig(g: DirectedMultiGraph) -> list:
    return sorted(g.edges_iter(keys=True))


class TestDirectedGeneratorIO:
    def test_is_iomixin_with_directed_multigraph_formats(self) -> None:
        assert issubclass(DirectedGenerator, IOMixin)
        assert DirectedGenerator._default_format == "dot"
        assert set(DirectedGenerator._supported_formats) == {"dot", "edgelist"}

    @pytest.mark.parametrize("fmt", ["dot", "edgelist"])
    def test_round_trip_string(self, fmt: str) -> None:
        g = _gen()
        g2 = DirectedGenerator.from_string(g.to_string(fmt), fmt)
        assert isinstance(g2, DirectedGenerator)
        assert _sig(g.graph) == _sig(g2.graph)

    def test_default_format_is_dot_digraph(self) -> None:
        g = _gen()
        assert g.to_string() == g.to_string("dot")
        assert g.to_string("dot").lstrip().startswith("digraph")

    @pytest.mark.parametrize("ext", [".dot", ".el"])
    def test_save_load_file_autodetect(self, ext: str) -> None:
        g = _gen()
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "gen" + ext)
            g.save(path)
            g2 = DirectedGenerator.load(path)
            assert isinstance(g2, DirectedGenerator)
            assert _sig(g.graph) == _sig(g2.graph)

    def test_convert_string_between_formats(self) -> None:
        g = _gen()
        as_el = DirectedGenerator.convert_string(g.to_string("dot"), "dot", "edgelist")
        g2 = DirectedGenerator.from_string(as_el, "edgelist")
        assert _sig(g.graph) == _sig(g2.graph)
