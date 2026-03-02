"""
Tests for attachment of leaves to semi-directed generators.

This module tests attaching taxa to ``SemiDirectedGenerator`` sides using the
``attach_leaves_to_generator`` utility.
"""

from __future__ import annotations

import pytest

from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
from phylozoo.core.network.sdnetwork.generator import (
    SemiDirectedGenerator,
    all_level_k_generators,
    attach_leaves_to_generator,
)
from phylozoo.core.network.sdnetwork.generator.side import (
    BidirectedEdgeSide,
    UndirEdgeSide,
)
from phylozoo.core.network.dnetwork.generator.side import (
    DirEdgeSide,
    HybridSide,
    IsolatedNodeSide,
    Side,
)
from phylozoo.utils.exceptions import PhyloZooValueError


class TestAttachLeavesNodeSides:
    """Tests for attaching leaves to node sides of SD generators."""

    def test_attach_three_taxa_to_level0_node_side(self) -> None:
        """
        Attach exactly three taxa to a level-0 SD generator's node side (binary).

        The single node receives three leaves via directed edges.
        """
        mm = MixedMultiGraph()
        mm.add_node(0)
        gen = SemiDirectedGenerator(mm)

        side: Side = gen.sides[0]
        assert isinstance(side, IsolatedNodeSide)
        taxa = ["A", "B", "C"]

        network = attach_leaves_to_generator(gen, {side: taxa})

        assert network.taxa == {"A", "B", "C"}
        assert len(network.leaves) == 3

    def test_level0_required_side_must_appear_in_mapping(self) -> None:
        """IsolatedNodeSide must appear in side_taxa; omitting it raises."""
        mm = MixedMultiGraph()
        mm.add_node(0)
        gen = SemiDirectedGenerator(mm)

        with pytest.raises(PhyloZooValueError, match="must appear in side_taxa"):
            attach_leaves_to_generator(gen, {})

    def test_level0_node_side_requires_exactly_three_taxa(self) -> None:
        """IsolatedNodeSide must receive exactly three taxa."""
        mm = MixedMultiGraph()
        mm.add_node(0)
        gen = SemiDirectedGenerator(mm)
        side = gen.sides[0]

        with pytest.raises(PhyloZooValueError, match="exactly three taxa"):
            attach_leaves_to_generator(gen, {side: ["A", "B"]})
        with pytest.raises(PhyloZooValueError, match="exactly three taxa"):
            attach_leaves_to_generator(gen, {side: ["A", "B", "C", "D"]})

    def test_attach_single_taxon_to_hybrid_side_level1(self) -> None:
        """Attach one taxon to the hybrid side and one along the bidirected edge (removes self-loop)."""
        mm = MixedMultiGraph(undirected_edges=[(1, 1)])
        gen = SemiDirectedGenerator(mm)

        hybrid_side = gen.hybrid_sides[0]
        bidir_side = gen.edge_sides[0]
        assert isinstance(bidir_side, BidirectedEdgeSide)

        network = attach_leaves_to_generator(
            gen, {hybrid_side: ["X"], bidir_side: ["Y"]}
        )

        assert "X" in network.taxa
        assert "Y" in network.taxa
        assert len(network.leaves) >= 2

    def test_hybrid_side_requires_exactly_one_taxon(self) -> None:
        """Hybrid side must receive exactly one taxon."""
        mm = MixedMultiGraph(undirected_edges=[(1, 1)])
        gen = SemiDirectedGenerator(mm)
        hybrid_side = gen.hybrid_sides[0]

        with pytest.raises(PhyloZooValueError, match="exactly one taxon"):
            attach_leaves_to_generator(gen, {hybrid_side: []})
        with pytest.raises(PhyloZooValueError, match="exactly one taxon"):
            attach_leaves_to_generator(gen, {hybrid_side: ["X", "Y"]})

    def test_level1_requires_bidirected_side_with_taxon(self) -> None:
        """
        Level-1 generators must provide at least one taxon on the bidirected edge side.

        Omitting the bidirected side or giving it an empty list should raise.
        """
        mm = MixedMultiGraph(undirected_edges=[(1, 1)])
        gen = SemiDirectedGenerator(mm)

        hybrid_side = gen.hybrid_sides[0]
        bidir_side = gen.edge_sides[0]
        assert isinstance(bidir_side, BidirectedEdgeSide)

        # Omit bidirected side entirely
        with pytest.raises(PhyloZooValueError, match="bidirected edge side"):
            attach_leaves_to_generator(gen, {hybrid_side: ["X"]})

        # Provide empty taxa list for bidirected side
        with pytest.raises(PhyloZooValueError, match="bidirected edge side"):
            attach_leaves_to_generator(gen, {hybrid_side: ["X"], bidir_side: []})


class TestAttachLeavesEdgeSides:
    """Tests for attaching leaves to edge sides (Dir, Undir, Bidirected)."""

    def test_edge_sides_may_be_omitted_or_empty(self) -> None:
        """Edge sides may be omitted or given an empty list."""
        gens = list(all_level_k_generators(2))
        gen = gens[0]
        required: dict[Side, list[str]] = {}
        for i, s in enumerate(gen.sides):
            if isinstance(s, HybridSide):
                required[s] = [f"H{i}"]
            elif isinstance(s, IsolatedNodeSide):
                required[s] = ["A", "B", "C"]
        # Ensure at least two taxa in total by attaching one extra leaf on a
        # directed edge side (if present). Other edge sides may still be
        # omitted or have empty taxa lists.
        dir_sides = [s for s in gen.sides if isinstance(s, DirEdgeSide)]
        if dir_sides:
            required[dir_sides[0]] = ["E1"]

        network = attach_leaves_to_generator(gen, required)
        assert network.number_of_nodes() >= 1

    def test_attach_along_bidirected_edge_side_level1(self) -> None:
        """Attach taxa along the bidirected edge side of a level-1 generator."""
        mm = MixedMultiGraph(undirected_edges=[(1, 1)])
        gen = SemiDirectedGenerator(mm)
        hybrid_side = gen.hybrid_sides[0]
        bidir_side = gen.edge_sides[0]
        assert isinstance(bidir_side, BidirectedEdgeSide)

        network = attach_leaves_to_generator(
            gen, {hybrid_side: ["H"], bidir_side: ["P", "Q"]}
        )

        assert network.taxa >= {"H", "P", "Q"}
        assert len(network.leaves) >= 3


class TestAttachLeavesLevel2:
    """Tests for level-2 SD generator attachment."""

    def test_attach_to_level2_generator(self) -> None:
        """Level-2 generator: required hybrid sides get one taxon each."""
        gens = list(all_level_k_generators(2))
        assert len(gens) >= 1
        gen = gens[0]

        side_taxa: dict[Side, list[str]] = {}
        for i, s in enumerate(gen.sides):
            if isinstance(s, HybridSide):
                side_taxa[s] = [f"Hybrid{i}"]
        assert len(side_taxa) >= 1, "Level-2 generator should have at least one hybrid side"

        # Add a single extra leaf on a directed edge side (if present) so that
        # at least two taxa are attached in total.
        dir_sides = [s for s in gen.sides if isinstance(s, DirEdgeSide)]
        if dir_sides:
            side_taxa[dir_sides[0]] = ["Extra"]

        network = attach_leaves_to_generator(gen, side_taxa)
        assert network.number_of_nodes() >= 1
        assert any(t.startswith("Hybrid") for t in network.taxa)

    def test_requires_at_least_two_taxa_total_level2(self) -> None:
        """
        Attaching fewer than two taxa in total on a level-2 generator should raise.
        """
        gens = list(all_level_k_generators(2))
        # Choose a generator with exactly one hybrid side so that the global
        # \"at least two taxa\" check is the first failing condition.
        gen = next(g for g in gens if len(g.hybrid_sides) == 1)

        hybrid_side = gen.hybrid_sides[0]

        with pytest.raises(PhyloZooValueError, match="At least two taxa"):
            attach_leaves_to_generator(gen, {hybrid_side: ["H"]})
