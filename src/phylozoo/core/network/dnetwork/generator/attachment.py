"""
Attachment utilities for directed level-k generators.

This module provides functions for attaching leaves (taxa) to the sides of a
``DirectedGenerator`` to obtain a full ``DirectedPhyNetwork``.

Leaves are generated as new nodes in the underlying ``DirectedMultiGraph`` and
are attached either to node sides (``NodeSide`` / ``HybridSide``) or along
edge sides (``DirEdgeSide``), depending on the side type.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, TypeVar

from phylozoo.core.network.dnetwork.conversions import dnetwork_from_graph
from phylozoo.core.network.dnetwork.generator.base import DirectedGenerator
from phylozoo.core.network.dnetwork.generator.side import (
    DirEdgeSide,
    HybridSide,
    Level0NodeSide,
    NodeSide,
    Side,
)
from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
from phylozoo.core.network.dnetwork.base import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooValueError

T = TypeVar("T")


def _attach_leaves_to_node_side(
    graph: DirectedMultiGraph[T],
    side: NodeSide[T],
    taxa: Sequence[str],
) -> None:
    """
    Attach leaves to a node side.

    For a ``HybridSide``, exactly one taxon must be provided. For a non-hybrid
    node side (e.g. ``Level0NodeSide``), exactly three taxa must be provided
    (binary networks). Otherwise a :class:`PhyloZooValueError` is raised.

    Parameters
    ----------
    graph : DirectedMultiGraph[T]
        Graph to which new leaf nodes and edges are added.
    side : NodeSide[T]
        Node side (including ``HybridSide`` and ``Level0NodeSide``) that serves
        as attachment point.
    taxa : Sequence[str]
        Ordered sequence of taxon labels to attach as leaves.

    Raises
    ------
    PhyloZooValueError
        If ``side`` is a ``HybridSide`` and the number of taxa is not exactly
        one, or if ``side`` is a non-hybrid node side and the number of taxa
        is not exactly three.
    """
    if isinstance(side, HybridSide):
        if len(taxa) != 1:
            raise PhyloZooValueError(
                f"Hybrid side {side} must have exactly one taxon, got {len(taxa)}."
            )
    else:
        # Non-hybrid node side (e.g. Level0NodeSide): binary networks require 3 leaves
        if len(taxa) != 3:
            raise PhyloZooValueError(
                f"Non-hybrid node side {side} must have exactly three taxa (binary), got {len(taxa)}."
            )

    if side.node not in graph.nodes():
        raise PhyloZooValueError(
            f"Node {side.node!r} from side {side} does not exist in generator graph."
        )

    for taxon in taxa:
        leaf_id = next(graph.generate_node_ids(1))
        # Attach leaf node with label attribute
        graph.add_node(leaf_id, label=taxon)
        graph.add_edge(side.node, leaf_id)


def _attach_leaves_to_dir_edge_side(
    graph: DirectedMultiGraph[T],
    side: DirEdgeSide[T],
    taxa: Sequence[str],
) -> None:
    """
    Attach leaves along a directed edge side.

    The edge ``(u, v, key)`` corresponding to ``side`` is subdivided into a
    path and each taxon is attached at a subdivision node. If there are
    ``k`` taxa, the edge is transformed into::

        u -> w1 -> w2 -> ... -> wk -> v

    and taxon ``taxa[i]`` is attached as a leaf to node ``w_{i+1}``.

    Parameters
    ----------
    graph : DirectedMultiGraph[T]
        Graph to which new internal and leaf nodes are added.
    side : DirEdgeSide[T]
        Edge side identifying the edge to subdivide.
    taxa : Sequence[str]
        Ordered sequence of taxon labels to attach along the edge. The order
        determines the order of attachment from ``u`` towards ``v``.

    Raises
    ------
    PhyloZooValueError
        If the edge corresponding to ``side`` does not exist in the graph.
    """
    if not taxa:
        return

    u, v, key = side.u, side.v, side.key

    # Ensure the edge exists before attempting to remove it
    if not graph.has_edge(u, v, key=key):
        raise PhyloZooValueError(
            f"Edge ({u!r}, {v!r}, key={key!r}) from side {side} does not exist in generator graph."
        )

    # Remove the original edge and build a path u -> w1 -> ... -> wk -> v
    graph.remove_edge(u, v, key=key)

    previous: T = u
    attachment_nodes: list[T] = []

    for _ in taxa:
        w = next(graph.generate_node_ids(1))
        graph.add_node(w)
        graph.add_edge(previous, w)
        previous = w
        attachment_nodes.append(w)

    # Connect the last subdivision node to v
    graph.add_edge(previous, v)

    # Attach one leaf per subdivision node, preserving taxon order
    for attach_node, taxon in zip(attachment_nodes, taxa):
        leaf_id = next(graph.generate_node_ids(1))
        graph.add_node(leaf_id, label=taxon)
        graph.add_edge(attach_node, leaf_id)


def attach_leaves_to_generator(
    generator: DirectedGenerator[T],
    side_taxa: Mapping[Side, Sequence[str]],
) -> DirectedPhyNetwork[T]:
    """
    Attach leaves to a generator to build a binary directed phylogenetic network.

    This function takes a :class:`DirectedGenerator` and a mapping from sides to
    ordered lists of taxa. It returns a :class:`DirectedPhyNetwork` obtained by
    attaching new leaf nodes to the generator's underlying graph:

    - Every ``HybridSide`` of the generator must appear in ``side_taxa`` with
      exactly one taxon (cannot be omitted or empty).
    - Every ``Level0NodeSide`` (the single side when the generator has level 0)
      must appear in ``side_taxa`` with exactly three taxa (cannot be omitted or empty).
    - ``DirEdgeSide`` sides may be omitted or given an empty list; they receive no leaves.

    Parameters
    ----------
    generator : DirectedGenerator[T]
        The generator whose sides will receive attached leaves.
    side_taxa : Mapping[Side, Sequence[str]]
        Mapping from sides of ``generator`` to ordered sequences of taxon
        labels. The order of taxa for each side determines the attachment
        order along that side. The order of sides in the mapping is not
        important.

    Returns
    -------
    DirectedPhyNetwork[T]
        A directed phylogenetic network obtained by attaching leaves to the
        generator according to ``side_taxa``.

    Raises
    ------
    PhyloZooValueError
        If a required side (hybrid or level-0 node side) is missing from
        ``side_taxa`` or has the wrong number of taxa; or if a side refers to
        a node or edge not present in the generator.
    """
    sides = generator.sides
    required_node_sides: list[NodeSide[T]] = []
    for s in sides:
        if isinstance(s, HybridSide):
            required_node_sides.append(s)
        elif isinstance(s, Level0NodeSide):
            required_node_sides.append(s)

    for side in required_node_sides:
        if side not in side_taxa:
            raise PhyloZooValueError(
                f"Required side {side} must appear in side_taxa with "
                f"{1 if isinstance(side, HybridSide) else 3} taxon(s)."
            )
        taxa = side_taxa[side]
        if isinstance(side, HybridSide):
            if len(taxa) != 1:
                raise PhyloZooValueError(
                    f"Hybrid side {side} must have exactly one taxon, got {len(taxa)}."
                )
        else:
            if len(taxa) != 3:
                raise PhyloZooValueError(
                    f"Level0NodeSide {side} must have exactly three taxa (binary), got {len(taxa)}."
                )

    # Work on a copy of the generator's graph to avoid mutating the original
    graph: DirectedMultiGraph[T] = generator.graph.copy()

    for side, taxa in side_taxa.items():
        # Mapping does not need to contain all sides. For edge sides, an empty
        # taxa list is a no-op. For hybrid sides we still call the helper to
        # enforce the "exactly one taxon" requirement.
        if isinstance(side, DirEdgeSide):
            if not taxa:
                continue
            _attach_leaves_to_dir_edge_side(graph, side, taxa)
        elif isinstance(side, NodeSide):
            _attach_leaves_to_node_side(graph, side, taxa)
        else:
            raise PhyloZooValueError(
                f"Unsupported side type {type(side)} in side_taxa mapping."
            )

    # Convert the augmented graph into a full DirectedPhyNetwork
    return dnetwork_from_graph(graph)

