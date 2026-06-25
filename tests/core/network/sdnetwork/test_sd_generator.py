"""
Tests for semi-directed level-k generators (SemiDirectedGenerator, all_level_k_generators).
"""

from __future__ import annotations

import pytest

from phylozoo.core.network.sdnetwork.generator import (
    SemiDirectedGenerator,
    all_level_k_generators,
    dgenerator_to_sdgenerator,
)
from phylozoo.core.network.sdnetwork.generator.attachment import attach_leaves_to_generator
from phylozoo.core.network.sdnetwork.generator.side import UndirEdgeSide
from phylozoo.core.network.dnetwork.generator.construction import (
    all_level_k_generators as all_level_k_dgenerators,
)
from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph


class TestAllLevelKGenerators:
    """Tests for all_level_k_generators (semi-directed)."""

    def test_level_0_returns_one(self) -> None:
        """Level-0 SD generators: exactly one (single node)."""
        gens = all_level_k_generators(0)
        assert len(gens) == 1
        gen = next(iter(gens))
        assert gen.level == 0
        assert gen.graph.number_of_nodes() == 1
        assert gen.graph.number_of_edges() == 0

    def test_level_1_returns_one(self) -> None:
        """Level-1 SD generators: exactly one (one node, one undirected self-loop)."""
        gens = all_level_k_generators(1)
        assert len(gens) == 1
        gen = next(iter(gens))
        assert gen.level == 1
        assert gen.graph.number_of_nodes() == 1
        assert gen.graph.number_of_edges() == 1

    def test_level_2_returns_two(self) -> None:
        """Level-2 SD generators: exactly 2 (known count)."""
        gens = all_level_k_generators(2)
        assert len(gens) == 2
        for gen in gens:
            assert gen.level == 2
            assert isinstance(gen, SemiDirectedGenerator)

    def test_level_3_returns_seventeen(self) -> None:
        """Level-3 SD generators: exactly 17 (known count)."""
        gens = all_level_k_generators(3)
        assert len(gens) == 17
        for gen in gens:
            assert gen.level == 3
            assert isinstance(gen, SemiDirectedGenerator)

    def test_negative_level_raises(self) -> None:
        """Negative level raises."""
        with pytest.raises(Exception):  # PhyloZooValueError or ValueError
            all_level_k_generators(-1)


class TestSemiDirectedGeneratorFromMixedGraph:
    """Tests that converted directed generators yield valid SD generators (validation on)."""

    def test_sd_gen_from_converted_d_gen_level_3(self) -> None:
        """Every level-3 directed generator converts to a valid SD generator (validation on)."""
        d_gens = list(all_level_k_dgenerators(3))
        for d_gen in d_gens:
            sd_gen = dgenerator_to_sdgenerator(d_gen)
            assert isinstance(sd_gen, SemiDirectedGenerator)
            assert sd_gen.level == 3

    def test_sd_gen_from_copy_of_converted_graph(self) -> None:
        """Building SemiDirectedGenerator from a copy of a converted graph still validates."""
        d_gens = list(all_level_k_dgenerators(3))
        d_gen = next(iter(d_gens))
        sd_gen = dgenerator_to_sdgenerator(d_gen)
        graph_copy = sd_gen.graph.copy()
        gen_from_copy = SemiDirectedGenerator(graph_copy)
        assert gen_from_copy.level == sd_gen.level
        assert graph_copy.number_of_nodes() == gen_from_copy.graph.number_of_nodes()
        assert graph_copy.number_of_edges() == gen_from_copy.graph.number_of_edges()


class TestHybridSidesOccupiedReticulation:
    """
    A hybrid node whose exit is already an undirected (backbone) edge is a hybrid
    *node* but not a hybrid *side*: attaching a pendant leaf there would make it an
    invalid hybrid (total degree = in-degree + 2). Regression test for that case.
    """

    @staticmethod
    def _occupied_generator() -> SemiDirectedGenerator:
        # Level-4 generator in which reticulation node 3 (in-edges 2->3, 5->3) also
        # has an undirected backbone edge 3--6, so its child slot is occupied.
        graph = MixedMultiGraph(
            directed_edges=[(2, 3), (5, 3), (5, 4), (6, 1), (6, 4), (8, 1), (8, 7), (9, 7)],
            undirected_edges=[(2, 5), (2, 9), (3, 6), (8, 9)],
        )
        return SemiDirectedGenerator(graph)

    def test_occupied_reticulation_is_hybrid_node_but_not_hybrid_side(self) -> None:
        gen = self._occupied_generator()
        assert gen.level == 4
        assert 3 in gen.hybrid_nodes  # still a reticulation
        side_nodes = {s.node for s in gen.hybrid_sides}
        assert 3 not in side_nodes  # but NOT an attachable side
        # the three free reticulations remain hybrid sides
        assert side_nodes == {1, 4, 7}

    def test_attach_leaves_skips_occupied_reticulation(self) -> None:
        """Attaching one leaf per (free) hybrid side + leaves on an edge side builds
        a valid network, instead of raising a degree error on the occupied node."""
        gen = self._occupied_generator()
        side_taxa: dict = {s: [f"r{i}"] for i, s in enumerate(gen.hybrid_sides)}
        undir_side = next(s for s in gen.sides if isinstance(s, UndirEdgeSide))
        side_taxa[undir_side] = ["a", "b"]
        net = attach_leaves_to_generator(gen, side_taxa)
        # 3 free reticulation leaves + 2 edge leaves
        assert {"r0", "r1", "r2", "a", "b"} == set(net.taxa)

    def test_level_1_node_is_still_a_hybrid_side(self) -> None:
        """The fix must not drop the level-1 bidirected self-loop hybrid side."""
        gen = next(iter(all_level_k_generators(1)))
        assert gen.level == 1
        assert len(gen.hybrid_sides) == 1
