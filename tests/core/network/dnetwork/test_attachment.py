"""
Tests for attachment of leaves to generators.

This module tests attaching taxa to ``DirectedGenerator`` sides using the
``attach_leaves_to_generator`` utility.
"""

from __future__ import annotations

import pytest

from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
from phylozoo.core.network.dnetwork.generator import (
    DirectedGenerator,
    DirEdgeSide,
    HybridSide,
    IsolatedNodeSide,
    NodeSide,
    Side,
    attach_leaves_to_generator,
)
from phylozoo.utils.exceptions import PhyloZooValueError


class TestAttachLeavesNodeSides:
    """Tests for attaching leaves to node sides."""

    def test_attach_three_taxa_to_level0_node_side(self) -> None:
        """
        Attach exactly three taxa to a level-0 generator's node side (binary).

        The single node of the level-0 generator becomes the root of the
        resulting network and has one outgoing edge to each of the three leaves.
        """
        graph: DirectedMultiGraph[int] = DirectedMultiGraph()
        graph.add_node(0)
        gen = DirectedGenerator(graph)

        side: Side = gen.sides[0]
        assert isinstance(side, IsolatedNodeSide)
        taxa = ["A", "B", "C"]

        network = attach_leaves_to_generator(gen, {side: taxa})

        assert network.root_node == 0
        assert network.taxa == {"A", "B", "C"}
        assert len(network.leaves) == 3

    def test_level0_required_side_must_appear_in_mapping(self) -> None:
        """
        IsolatedNodeSide must appear in side_taxa; omitting it raises.
        """
        graph: DirectedMultiGraph[int] = DirectedMultiGraph()
        graph.add_node(0)
        gen = DirectedGenerator(graph)

        with pytest.raises(PhyloZooValueError, match="must appear in side_taxa"):
            attach_leaves_to_generator(gen, {})

    def test_level0_node_side_requires_exactly_three_taxa(self) -> None:
        """
        Non-hybrid node side (e.g. IsolatedNodeSide) must receive exactly three taxa.

        Providing fewer or more than three should raise PhyloZooValueError.
        """
        graph: DirectedMultiGraph[int] = DirectedMultiGraph()
        graph.add_node(0)
        gen = DirectedGenerator(graph)
        side = gen.sides[0]

        with pytest.raises(PhyloZooValueError, match="exactly three taxa"):
            attach_leaves_to_generator(gen, {side: ["A", "B"]})
        with pytest.raises(PhyloZooValueError, match="exactly three taxa"):
            attach_leaves_to_generator(gen, {side: ["A"]})
        with pytest.raises(PhyloZooValueError, match="exactly three taxa"):
            attach_leaves_to_generator(gen, {side: []})
        with pytest.raises(PhyloZooValueError, match="exactly three taxa"):
            attach_leaves_to_generator(gen, {side: ["A", "B", "C", "D"]})

    def test_attach_single_taxon_to_hybrid_side(self) -> None:
        """
        Attach a single taxon to a hybrid side of a level-1 generator, plus one
        additional leaf on an edge side so that at least two taxa are attached
        in total.

        The hybrid node obtains out-degree 1 (to the new leaf) and the taxon
        appears in the network's taxa.
        """
        graph: DirectedMultiGraph[int] = DirectedMultiGraph(edges=[(0, 1), (0, 1)])
        gen = DirectedGenerator(graph)

        hybrid_sides = gen.hybrid_sides
        assert len(hybrid_sides) == 1
        hybrid_side: HybridSide[int] = hybrid_sides[0]

        edge_side: DirEdgeSide[int] = gen.edge_sides[0]

        taxa_hybrid = ["X"]
        taxa_edge = ["Y"]
        network = attach_leaves_to_generator(
            gen, {hybrid_side: taxa_hybrid, edge_side: taxa_edge}
        )

        assert "X" in network.taxa
        outdeg = network._graph.outdegree(hybrid_side.node)  # type: ignore[attr-defined]
        assert outdeg == 1

    def test_hybrid_side_required_must_appear_in_mapping(self) -> None:
        """
        Every hybrid side must appear in side_taxa; omitting it raises.
        """
        graph: DirectedMultiGraph[int] = DirectedMultiGraph(edges=[(0, 1), (0, 1)])
        gen = DirectedGenerator(graph)

        with pytest.raises(PhyloZooValueError, match="must appear in side_taxa"):
            attach_leaves_to_generator(gen, {})

    def test_hybrid_side_requires_exactly_one_taxon(self) -> None:
        """
        Hybrid sides must receive exactly one taxon.

        Providing zero or more than one taxon should raise PhyloZooValueError.
        """
        graph: DirectedMultiGraph[int] = DirectedMultiGraph(edges=[(0, 1), (0, 1)])
        gen = DirectedGenerator(graph)

        hybrid_side = gen.hybrid_sides[0]

        with pytest.raises(PhyloZooValueError, match="must have exactly one taxon"):
            attach_leaves_to_generator(gen, {hybrid_side: []})

        with pytest.raises(PhyloZooValueError, match="must have exactly one taxon"):
            attach_leaves_to_generator(gen, {hybrid_side: ["X", "Y"]})


class TestAttachLeavesEdgeSides:
    """Tests for attaching leaves to directed edge sides."""

    def test_attach_taxa_along_edge_side(self) -> None:
        """
        Attach taxa along a DirEdgeSide by subdividing the edge.

        The original parallel edge is subdivided into a path and each
        subdivision node receives exactly one leaf with the correct taxon.
        """
        graph: DirectedMultiGraph[int] = DirectedMultiGraph(edges=[(0, 1), (0, 1)])
        gen = DirectedGenerator(graph)

        edge_sides = gen.edge_sides
        assert len(edge_sides) == 2
        edge_side: DirEdgeSide[int] = edge_sides[0]

        # Also attach exactly one taxon to the hybrid side to obtain
        # a valid network (hybrid node with out-degree 1).
        hybrid_side: HybridSide[int] = gen.hybrid_sides[0]

        taxa_edge = ["A", "B"]
        taxa_hybrid = ["H"]
        network = attach_leaves_to_generator(
            gen, {edge_side: taxa_edge, hybrid_side: taxa_hybrid}
        )

        # Taxa should match exactly the attached labels
        assert network.taxa == {"A", "B", "H"}

    def test_dir_edge_sides_may_be_omitted_or_empty(self) -> None:
        """
        DirEdgeSides may be omitted from side_taxa or given an empty list.
        """
        graph: DirectedMultiGraph[int] = DirectedMultiGraph(edges=[(0, 1), (0, 1)])
        gen = DirectedGenerator(graph)
        hybrid_side = gen.hybrid_sides[0]
        node_side: NodeSide[int] = NodeSide(node=gen.root_node)

        # Only required sides (hybrid and a level-0-like node side) present; no edge sides
        network = attach_leaves_to_generator(
            gen, {hybrid_side: ["H"], node_side: ["R1", "R2", "R3"]}
        )
        assert network.taxa == {"H", "R1", "R2", "R3"}

    def test_requires_at_least_two_taxa_total(self) -> None:
        """
        Attaching fewer than two taxa in total should raise PhyloZooValueError.
        """
        graph: DirectedMultiGraph[int] = DirectedMultiGraph(edges=[(0, 1), (0, 1)])
        gen = DirectedGenerator(graph)
        hybrid_side = gen.hybrid_sides[0]

        with pytest.raises(PhyloZooValueError, match="At least two taxa"):
            attach_leaves_to_generator(gen, {hybrid_side: ["H"]})


class TestAttachLeavesMixed:
    """Mixed tests for attaching leaves to multiple sides."""

    def test_attach_to_multiple_sides(self) -> None:
        """
        Attach leaves to root node side (three taxa), edge side, and hybrid side.

        Sides without entries in the mapping remain unchanged.
        """
        graph: DirectedMultiGraph[int] = DirectedMultiGraph(edges=[(0, 1), (0, 1)])
        gen = DirectedGenerator(graph)

        # Root node is not in generator.sides for level-1; use NodeSide for attachment
        node_side: NodeSide[int] = NodeSide(node=gen.root_node)
        edge_side: DirEdgeSide[int] = gen.edge_sides[0]
        hybrid_side: HybridSide[int] = gen.hybrid_sides[0]

        side_taxa: dict[Side, list[str]] = {
            node_side: ["R1", "R2", "R3"],
            edge_side: ["L1", "L2"],
            hybrid_side: ["H1"],
        }

        network = attach_leaves_to_generator(gen, side_taxa)

        assert network.taxa == {"R1", "R2", "R3", "L1", "L2", "H1"}

