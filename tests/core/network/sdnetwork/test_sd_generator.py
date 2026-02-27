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
from phylozoo.core.network.dnetwork.generator.construction import (
    all_level_k_generators as all_level_k_dgenerators,
)


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
