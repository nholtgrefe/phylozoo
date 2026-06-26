"""I/O tests for :class:`SemiDirectedGenerator` (IOMixin via its MixedMultiGraph)."""

import os
import tempfile

import pytest

from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
from phylozoo.core.network.sdnetwork.generator.base import SemiDirectedGenerator
from phylozoo.core.network.sdnetwork.generator.attachment import attach_leaves_to_generator
from phylozoo.core.network.sdnetwork.generator.construction import all_level_k_generators
from phylozoo.core.network.dnetwork.generator.side import HybridSide
from phylozoo.utils.io import IOMixin

# A real level-3 generator that exercises parallel directed edges plus directed
# and undirected edges (built once at import; generators are immutable).
_GEN = next(
    g
    for g in all_level_k_generators(3)
    if g.parallel_directed_edge_sides and any(True for _ in g.graph.undirected_edges_iter())
)


def _gen() -> SemiDirectedGenerator:
    return _GEN


def _dsig(g: MixedMultiGraph) -> list:
    return sorted(g.directed_edges_iter(keys=True))


def _usig(g: MixedMultiGraph) -> list:
    return sorted((tuple(sorted((a, b))), k) for a, b, k in g.undirected_edges_iter(keys=True))


class TestSemiDirectedGeneratorIO:
    def test_is_iomixin_with_mixed_multigraph_formats(self) -> None:
        assert issubclass(SemiDirectedGenerator, IOMixin)
        assert SemiDirectedGenerator._default_format == "phylozoo-dot"
        assert set(SemiDirectedGenerator._supported_formats) == {"phylozoo-dot", "dot"}

    @pytest.mark.parametrize("fmt", ["phylozoo-dot", "dot"])
    def test_round_trip_string(self, fmt: str) -> None:
        g = _gen()
        g2 = SemiDirectedGenerator.from_string(g.to_string(fmt), fmt)
        assert isinstance(g2, SemiDirectedGenerator)
        assert _dsig(g.graph) == _dsig(g2.graph)
        assert _usig(g.graph) == _usig(g2.graph)

    def test_default_format_and_headers(self) -> None:
        g = _gen()
        assert g.to_string() == g.to_string("phylozoo-dot")
        assert g.to_string("phylozoo-dot").lstrip().startswith("graph {")
        assert g.to_string("dot").lstrip().startswith("digraph")

    @pytest.mark.parametrize("ext", [".pzdot", ".dot"])
    def test_save_load_file_autodetect(self, ext: str) -> None:
        g = _gen()
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "gen" + ext)
            g.save(path)
            g2 = SemiDirectedGenerator.load(path)
            assert isinstance(g2, SemiDirectedGenerator)
            assert _dsig(g.graph) == _dsig(g2.graph)
            assert _usig(g.graph) == _usig(g2.graph)

    def test_round_tripped_generator_is_usable(self) -> None:
        """A generator loaded directly attaches leaves without a strip-and-rebuild step."""
        g2 = SemiDirectedGenerator.from_string(_gen().to_string("dot"), "dot")
        side_taxa = {h: [f"x{i}"] for i, h in enumerate(g2.hybrid_sides)}
        edge_side = next(s for s in g2.sides if not isinstance(s, HybridSide))
        side_taxa[edge_side] = ["y0"]
        net = attach_leaves_to_generator(g2, side_taxa)
        assert net.number_of_nodes() > 0

    def test_convert_string_between_formats(self) -> None:
        g = _gen()
        as_dot = SemiDirectedGenerator.convert_string(
            g.to_string("phylozoo-dot"), "phylozoo-dot", "dot"
        )
        assert as_dot.lstrip().startswith("digraph")
        g2 = SemiDirectedGenerator.from_string(as_dot, "dot")
        assert _dsig(g.graph) == _dsig(g2.graph)
        assert _usig(g.graph) == _usig(g2.graph)
