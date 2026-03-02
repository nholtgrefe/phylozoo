"""
Attachment utilities for semi-directed level-k generators.

This module provides functions for attaching leaves (taxa) to the sides of a
``SemiDirectedGenerator`` to obtain a full semi-directed phylogenetic network
(as a ``SemiDirectedPhyNetwork``).

Leaves are generated as new nodes in the underlying ``MixedMultiGraph`` and
are attached either to node sides (``NodeSide`` / ``HybridSide`` / ``IsolatedNodeSide``)
or along edge sides (``DirEdgeSide``, ``UndirEdgeSide``, ``BidirectedEdgeSide``),
depending on the side type.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TypeVar, TYPE_CHECKING

from phylozoo.core.network.dnetwork.generator.side import (
    DirEdgeSide,
    HybridSide,
    IsolatedNodeSide,
    NodeSide,
    Side,
)
from phylozoo.core.network.sdnetwork.conversions import sdnetwork_from_graph
from phylozoo.core.network.sdnetwork.generator.base import SemiDirectedGenerator
from phylozoo.core.network.sdnetwork.generator.side import (
    BidirectedEdgeSide,
    UndirEdgeSide,
)
from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph as MMG
from phylozoo.utils.exceptions import PhyloZooValueError

T = TypeVar("T")

if TYPE_CHECKING:
    from phylozoo.core.network.sdnetwork.base import SemiDirectedPhyNetwork


def _attach_leaves_to_node_side(
    graph: MixedMultiGraph[T],
    side: NodeSide[T],
    taxa: Sequence[str],
) -> None:
    """
    Attach leaves to a node side.

    For a ``HybridSide``, exactly one taxon must be provided. For a non-hybrid
    node side (e.g. ``IsolatedNodeSide``), exactly three taxa must be provided
    (binary networks). Otherwise a :class:`PhyloZooValueError` is raised.

    Parameters
    ----------
    graph : MixedMultiGraph[T]
        Graph to which new leaf nodes and edges are added.
    side : NodeSide[T]
        Node side (including ``HybridSide`` and ``IsolatedNodeSide``) that serves
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
        if len(taxa) != 3:
            raise PhyloZooValueError(
                f"Non-hybrid node side {side} must have exactly three taxa (binary), "
                f"got {len(taxa)}."
            )

    if side.node not in graph.nodes():
        raise PhyloZooValueError(
            f"Node {side.node!r} from side {side} does not exist in generator graph."
        )

    for taxon in taxa:
        leaf_id = next(graph.generate_node_ids(1))
        graph.add_node(leaf_id, label=taxon)
        graph.add_undirected_edge(side.node, leaf_id)


def _attach_leaves_to_edge_side(
    graph: MixedMultiGraph[T],
    side: DirEdgeSide[T] | UndirEdgeSide[T] | BidirectedEdgeSide[T],
    taxa: Sequence[str],
) -> None:
    """
    Attach leaves along an edge side (directed, undirected, or bidirected).

    For ``DirEdgeSide``: the directed edge is subdivided into undirected segments
    u - w1 - ... - wk and a final directed edge wk -> v; one leaf per wi.
    For ``UndirEdgeSide``: the undirected edge is subdivided into an undirected
    path u - w1 - ... - wk - v; one leaf per wi.
    For ``BidirectedEdgeSide``: the self-loop is replaced by directed w1 -> node,
    wk -> node and undirected path w1 - ... - wk; one leaf per wi.

    Parameters
    ----------
    graph : MixedMultiGraph[T]
        Graph to which new internal and leaf nodes are added.
    side : DirEdgeSide[T] | UndirEdgeSide[T] | BidirectedEdgeSide[T]
        Edge side identifying the edge to subdivide.
    taxa : Sequence[str]
        Ordered sequence of taxon labels to attach along the edge.

    Raises
    ------
    PhyloZooValueError
        If the edge corresponding to ``side`` does not exist in the graph.
    """
    if not taxa:
        return

    if isinstance(side, DirEdgeSide):
        u, v, key = side.u, side.v, side.key
        if not graph._directed.has_edge(u, v, key):
            raise PhyloZooValueError(
                f"Edge ({u!r}, {v!r}, key={key!r}) from side {side} does not exist "
                "in generator graph."
            )
        graph.remove_edge(u, v, key=key)
        previous = u
        attachment_nodes = []
        for _ in taxa:
            w = next(graph.generate_node_ids(1))
            graph.add_node(w)
            graph.add_undirected_edge(previous, w)
            previous = w
            attachment_nodes.append(w)
        graph.add_directed_edge(previous, v)
        for attach_node, taxon in zip(attachment_nodes, taxa):
            leaf_id = next(graph.generate_node_ids(1))
            graph.add_node(leaf_id, label=taxon)
            graph.add_undirected_edge(attach_node, leaf_id)

    elif isinstance(side, UndirEdgeSide):
        u, v, key = side.u, side.v, side.key
        norm_u, norm_v, _ = MMG.normalize_undirected_edge(u, v, key=key)
        if not graph._undirected.has_edge(norm_u, norm_v, key):
            raise PhyloZooValueError(
                f"Undirected edge ({u!r}, {v!r}, key={key!r}) from side {side} does "
                "not exist in generator graph."
            )
        graph.remove_edge(norm_u, norm_v, key=key)
        previous = norm_u
        attachment_nodes = []
        for _ in taxa:
            w = next(graph.generate_node_ids(1))
            graph.add_node(w)
            graph.add_undirected_edge(previous, w)
            previous = w
            attachment_nodes.append(w)
        graph.add_undirected_edge(previous, norm_v)
        for attach_node, taxon in zip(attachment_nodes, taxa):
            leaf_id = next(graph.generate_node_ids(1))
            graph.add_node(leaf_id, label=taxon)
            graph.add_undirected_edge(attach_node, leaf_id)

    else:
        # BidirectedEdgeSide: self-loop on side.node
        node, key = side.node, side.key
        if not graph._undirected.has_edge(node, node, key):
            raise PhyloZooValueError(
                f"Bidirected edge (node={node!r}, key={key!r}) from side {side} does "
                "not exist in generator graph."
            )
        graph.remove_edge(node, node, key=key)
        attachment_nodes = []
        for _ in taxa:
            w = next(graph.generate_node_ids(1))
            graph.add_node(w)
            attachment_nodes.append(w)
        # Directed: w1 -> node, wk -> node; undirected path w1 - w2 - ... - wk
        graph.add_directed_edge(attachment_nodes[0], node)
        graph.add_directed_edge(attachment_nodes[-1], node)
        for i in range(len(attachment_nodes) - 1):
            graph.add_undirected_edge(attachment_nodes[i], attachment_nodes[i + 1])
        for attach_node, taxon in zip(attachment_nodes, taxa):
            leaf_id = next(graph.generate_node_ids(1))
            graph.add_node(leaf_id, label=taxon)
            graph.add_undirected_edge(attach_node, leaf_id)


def attach_leaves_to_generator(
    generator: SemiDirectedGenerator[T],
    side_taxa: Mapping[Side, Sequence[str]],
) -> "SemiDirectedPhyNetwork[T]":
    """
    Attach leaves to a semi-directed generator to build a binary semi-directed network.

    Takes a :class:`SemiDirectedGenerator` and a mapping from sides to ordered
    lists of taxa. Returns a :class:`SemiDirectedPhyNetwork` by attaching new
    leaf nodes to the generator's underlying graph:

    - Every ``HybridSide`` must appear in ``side_taxa`` with exactly one taxon.
    - Every ``IsolatedNodeSide`` (level-0) must appear with exactly three taxa.
    - Edge sides (``DirEdgeSide``, ``UndirEdgeSide``, ``BidirectedEdgeSide``)
      may be omitted or given an empty list.

    Parameters
    ----------
    generator : SemiDirectedGenerator[T]
        The generator whose sides will receive attached leaves.
    side_taxa : Mapping[Side, Sequence[str]]
        Mapping from sides to ordered sequences of taxon labels.

    Returns
    -------
    SemiDirectedPhyNetwork[T]
        A semi-directed phylogenetic network with leaves attached per ``side_taxa``.

    Raises
    ------
    PhyloZooValueError
        If a required side is missing or has the wrong number of taxa; if fewer
        than two taxa are attached in total; or if a side refers to a node or
        edge not present in the generator. For a level-1 generator, the single
        ``BidirectedEdgeSide`` must also receive at least one taxon.

    Examples
    --------
    Build a level-1 semi-directed generator from a mixed multigraph (one node
    with an undirected self-loop), then attach leaves to the hybrid and
    bidirected edge sides:

    >>> from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
    >>> from phylozoo.core.network.sdnetwork.generator.base import SemiDirectedGenerator
    >>> mm = MixedMultiGraph(undirected_edges=[(1, 1)])
    >>> gen = SemiDirectedGenerator(mm)
    >>> hybrid_side = gen.hybrid_sides[0]
    >>> bidir_side = gen.edge_sides[0]
    >>> network = attach_leaves_to_generator(
    ...     gen, {hybrid_side: ["H"], bidir_side: ["P", "Q"]}
    ... )
    >>> sorted(network.taxa)
    ['H', 'P', 'Q']
    """
    sides = generator.sides
    required_node_sides: list[NodeSide[T]] = []
    level1_bidirected_sides: list[BidirectedEdgeSide[T]] = []

    for s in sides:
        if isinstance(s, (HybridSide, IsolatedNodeSide)):
            required_node_sides.append(s)
        if generator.level == 1 and isinstance(s, BidirectedEdgeSide):
            level1_bidirected_sides.append(s)

    for side in required_node_sides:
        if side not in side_taxa:
            n = 1 if isinstance(side, HybridSide) else 3
            raise PhyloZooValueError(
                f"Required side {side} must appear in side_taxa with {n} taxon(s)."
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
                    f"IsolatedNodeSide {side} must have exactly three taxa (binary), "
                    f"got {len(taxa)}."
                )

    # For level-1 generators, require that the bidirected edge side receives
    # at least one taxon, so that the self-loop is always subdivided away.
    if generator.level == 1:
        for side in level1_bidirected_sides:
            taxa = side_taxa.get(side)
            if not taxa:
                raise PhyloZooValueError(
                    f"Level-1 generator bidirected edge side {side} must have at least one taxon."
                )

    # Require at least two taxa in total across all sides
    total_taxa = sum(len(taxa) for taxa in side_taxa.values())
    if total_taxa < 2:
        raise PhyloZooValueError(
            f"At least two taxa must be attached in total, got {total_taxa}."
        )

    graph = generator.graph.copy()

    for side, taxa in side_taxa.items():
        if isinstance(side, (DirEdgeSide, UndirEdgeSide, BidirectedEdgeSide)):
            if not taxa:
                continue
            _attach_leaves_to_edge_side(graph, side, taxa)
        elif isinstance(side, NodeSide):
            _attach_leaves_to_node_side(graph, side, taxa)
        else:
            raise PhyloZooValueError(
                f"Unsupported side type {type(side)} in side_taxa mapping."
            )

    return sdnetwork_from_graph(graph, network_type="semi-directed")
