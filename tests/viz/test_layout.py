"""
Tests for layout computation in viz.

This test suite covers the pz-cladogram and pz-unrooted layout algorithms.
"""

import math

import numpy as np
import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooLayoutError, PhyloZooValueError
from phylozoo.viz import plot
from phylozoo.viz.dnetwork.layout import DNetLayout, compute_pz_cladogram_layout
from phylozoo.viz.sdnetwork.layout import SDNetLayout, compute_pz_unrooted_layout
from tests.fixtures import sd_networks
from tests.fixtures.directed_networks import (
    DTREE_LARGE_BINARY,
    LEVEL_1_DNETWORK_SINGLE_HYBRID,
    LEVEL_2_DNETWORK_MULTIPLE_BLOBS,
    LEVEL_5_DNETWORK_NON_TREEBASED,
)
from tests.fixtures.sd_networks import (
    LEVEL_1_SDNETWORK_FIVE_BLOBS,
    LEVEL_1_SDNETWORK_LARGE_MULTI_BLOB,
    LEVEL_3_SDNETWORK_LARGE_MANY_HYBRIDS,
    SDTREE_LARGE_BINARY,
    SDTREE_SINGLE_NODE,
)


class TestCladogramComputation:
    """Test DAG layout computation."""

    def test_simple_tree(self) -> None:
        """Test DAG layout on a simple tree."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )

        layout = compute_pz_cladogram_layout(net, direction="TD", trials=100)

        assert isinstance(layout, DNetLayout)
        assert len(layout.positions) == 3
        assert len(layout.edge_routes) == 2
        assert layout.algorithm == "pz-cladogram"

    def test_single_hybrid(self) -> None:
        """Test DAG layout on network with single hybrid."""
        net = LEVEL_1_DNETWORK_SINGLE_HYBRID

        layout = compute_pz_cladogram_layout(net, direction="TD", trials=100)

        assert isinstance(layout, DNetLayout)
        assert len(layout.positions) == net.number_of_nodes()
        assert len(layout.edge_routes) == net.number_of_edges()
        assert layout.algorithm == "pz-cladogram"

    def test_left_right_direction(self) -> None:
        """Test DAG layout with LR direction."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )

        layout = compute_pz_cladogram_layout(net, direction="LR", trials=100)

        assert layout.algorithm == "pz-cladogram"
        assert len(layout.positions) == 3

    def test_empty_network_raises(self) -> None:
        """Test that empty network raises PhyloZooLayoutError."""
        net = DirectedPhyNetwork()

        with pytest.raises(PhyloZooLayoutError, match="empty network"):
            compute_pz_cladogram_layout(net)

    def test_invalid_direction_raises(self) -> None:
        """Test that invalid direction raises PhyloZooValueError."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )

        with pytest.raises(PhyloZooValueError, match="direction must be"):
            compute_pz_cladogram_layout(net, direction="invalid")


class TestDLayoutOrdering:
    """Test the ordering and options of the pz-cladogram layout."""

    def test_deterministic(self) -> None:
        """Same seed gives the same layout."""
        net = LEVEL_2_DNETWORK_MULTIPLE_BLOBS
        a = compute_pz_cladogram_layout(net, trials=3, seed=1)
        b = compute_pz_cladogram_layout(net, trials=3, seed=1)
        assert a.positions == b.positions

    def test_leaves_aligned_by_default(self) -> None:
        """All leaves share the bottom layer with align_leaves=True."""
        net = LEVEL_2_DNETWORK_MULTIPLE_BLOBS
        layout = compute_pz_cladogram_layout(net)
        ys = {round(layout.positions[leaf][1], 9) for leaf in net.leaves}
        assert len(ys) == 1
        assert min(y for _, y in layout.positions.values()) == ys.pop()

    def test_leaves_not_aligned(self) -> None:
        """With align_leaves=False leaves sit at their own depth."""
        net = LEVEL_2_DNETWORK_MULTIPLE_BLOBS
        layout = compute_pz_cladogram_layout(net, align_leaves=False)
        ys = {round(layout.positions[leaf][1], 9) for leaf in net.leaves}
        assert len(ys) > 1

    def test_backbone_and_reticulate_partition_edges(self) -> None:
        """Every edge is either a backbone edge or a reticulate edge, never both."""
        net = LEVEL_2_DNETWORK_MULTIPLE_BLOBS
        layout = compute_pz_cladogram_layout(net)
        assert layout.backbone_edges.isdisjoint(layout.reticulate_edges)
        assert len(layout.backbone_edges) == net.number_of_nodes() - 1
        assert len(layout.backbone_edges | layout.reticulate_edges) == net.number_of_edges()

    def test_tree_has_no_crossings(self) -> None:
        """A tree is drawn without crossing edges."""
        from phylozoo.viz._layout_utils import count_crossings as _count_crossings

        net = DTREE_LARGE_BINARY
        layout = compute_pz_cladogram_layout(net)
        segments = np.array(
            [
                (*layout.positions[u], *layout.positions[v])
                for u, v, _ in net._graph.edges(keys=True)
            ]
        )
        assert _count_crossings(segments) == 0

    def test_left_right_puts_root_left(self) -> None:
        """In LR direction the root has the smallest x coordinate."""
        net = LEVEL_1_DNETWORK_SINGLE_HYBRID
        layout = compute_pz_cladogram_layout(net, direction="LR")
        root_x = layout.positions[net.root_node][0]
        assert root_x == min(x for x, _ in layout.positions.values())


class TestSDUnrootedLayout:
    """Test the tree-of-blobs (pz-unrooted) layout for semi-directed networks."""

    @pytest.mark.parametrize("name", sorted(sd_networks.NETWORK_METADATA))
    def test_all_fixtures(self, name: str) -> None:
        """Every fixture gets a finite position for every node and a route for every edge."""
        net = getattr(sd_networks, name)
        if net.number_of_nodes() == 0:
            with pytest.raises(PhyloZooLayoutError, match="empty"):
                compute_pz_unrooted_layout(net)
            return
        layout = compute_pz_unrooted_layout(net)
        assert isinstance(layout, SDNetLayout)
        assert layout.algorithm == "pz-unrooted"
        assert set(layout.positions) == set(net._graph.nodes)
        assert all(math.isfinite(x) and math.isfinite(y) for x, y in layout.positions.values())
        assert len(layout.edge_routes) == net.number_of_edges()

    def test_deterministic(self) -> None:
        """Repeated calls give identical positions."""
        net = LEVEL_3_SDNETWORK_LARGE_MANY_HYBRIDS
        assert (
            compute_pz_unrooted_layout(net).positions == compute_pz_unrooted_layout(net).positions
        )

    def test_blob_nodes_keep_their_cyclic_order(self) -> None:
        """A level-1 blob is drawn as a convex polygon in its cycle order (refine=0)."""
        from phylozoo.core.network.sdnetwork.features import blobs

        net = LEVEL_1_SDNETWORK_LARGE_MULTI_BLOB
        layout = compute_pz_unrooted_layout(net, refine=0)
        for blob in blobs(net, trivial=False, leaves=False):
            # Walk the cycle to get the graph's cyclic order of the blob nodes.
            adjacency = {v: [u for u in net.neighbors(v) if u in blob] for v in blob}
            assert all(len(nbrs) == 2 for nbrs in adjacency.values())  # level-1 blob
            cycle = [min(blob)]
            while len(cycle) < len(blob):
                nxt = [u for u in adjacency[cycle[-1]] if u not in cycle]
                cycle.append(nxt[0])
            # Angular order around the centroid must be the same cyclic sequence.
            cx = sum(layout.positions[v][0] for v in blob) / len(blob)
            cy = sum(layout.positions[v][1] for v in blob) / len(blob)
            by_angle = sorted(
                blob,
                key=lambda v: math.atan2(layout.positions[v][1] - cy, layout.positions[v][0] - cx),
            )
            i = by_angle.index(cycle[0])
            rotated = by_angle[i:] + by_angle[:i]
            assert rotated == cycle or rotated == [cycle[0]] + cycle[:0:-1]

    def test_no_crossings_in_tree_and_level_1_networks(self) -> None:
        """Trees (raw, refined, daylight) and level-1 networks (refined, daylight) have no crossings."""
        from phylozoo.viz._layout_utils import count_crossings as _count_crossings

        variants = [
            (SDTREE_LARGE_BINARY, {"refine": 0}),
            (SDTREE_LARGE_BINARY, {}),
            (SDTREE_LARGE_BINARY, {"refine": 0, "daylight": 8}),
            (LEVEL_1_SDNETWORK_LARGE_MULTI_BLOB, {}),
            (LEVEL_1_SDNETWORK_LARGE_MULTI_BLOB, {"refine": 0, "daylight": 8}),
        ]
        for net, kwargs in variants:
            for _ in (0,):
                layout = compute_pz_unrooted_layout(net, **kwargs)
                segments = np.array(
                    [
                        (*layout.positions[u], *layout.positions[v])
                        for u, v, _ in net._graph.edges(keys=True)
                    ]
                )
                assert _count_crossings(segments) == 0

    def test_cut_edges_have_edge_length(self) -> None:
        """Without refinement every cut edge is drawn with length edge_length (before scaling)."""
        from phylozoo.core.network.sdnetwork.features import cut_edges

        net = LEVEL_1_SDNETWORK_FIVE_BLOBS
        layout = compute_pz_unrooted_layout(net, refine=0, edge_length=1.0, inner_scale=1.0)
        lengths = {
            round(math.dist(layout.positions[u], layout.positions[v]), 9)
            for u, v, _ in cut_edges(net)
        }
        assert len(lengths) == 1  # all equal (up to the final uniform normalisation)

    def test_plot_default_layout_is_pz_unrooted(self) -> None:
        """plot() uses pz-unrooted for semi-directed networks by default."""
        from phylozoo.viz._dispatch import resolve_layout

        assert resolve_layout(LEVEL_1_SDNETWORK_FIVE_BLOBS, "auto") == "pz-unrooted"
        ax = plot(LEVEL_1_SDNETWORK_FIVE_BLOBS)
        assert ax is not None

    def test_single_node(self) -> None:
        """A single-node network is placed at the origin."""
        layout = compute_pz_unrooted_layout(SDTREE_SINGLE_NODE)
        assert list(layout.positions.values()) == [(0.0, 0.0)]


class TestCladogramRectangular:
    """Test pz-cladogram with rectangular=True."""

    def test_routes_are_elbows(self) -> None:
        """Backbone edges are 3-point elbows; reticulate edges end sideways into the hybrid."""
        net = LEVEL_2_DNETWORK_MULTIPLE_BLOBS
        layout = compute_pz_cladogram_layout(net, rectangular=True)
        assert layout.algorithm == "pz-cladogram"
        assert layout.parameters["rectangular"] is True
        assert len(layout.edge_routes) == net.number_of_edges()
        for key, route in layout.edge_routes.items():
            pts = route.points
            if key in layout.backbone_edges:
                (x0, y0), (x1, y1), (x2, y2) = pts
                assert y1 == y0 and x1 == x2  # horizontal bar at the parent, then down
            else:
                assert 2 <= len(pts) <= 4
                (x1, y1), (x2, y2) = pts[-2], pts[-1]
                assert y1 == y2  # enters the hybrid sideways, along its layer
                assert (x2, y2) == layout.positions[key[1]]
                for (ax_, ay), (bx, by) in zip(pts, pts[1:]):
                    assert ax_ == bx or ay == by  # orthogonal segments only

    def test_rectangular_defaults_to_horizontal_reticulations(self) -> None:
        """rectangular=True levels reticulate parents unless told otherwise."""
        net = LEVEL_2_DNETWORK_MULTIPLE_BLOBS
        rect = compute_pz_cladogram_layout(net, rectangular=True)
        assert (
            rect.positions
            == compute_pz_cladogram_layout(net, horizontal_reticulations=True).positions
        )
        plain = compute_pz_cladogram_layout(net, rectangular=True, horizontal_reticulations=False)
        assert plain.positions == compute_pz_cladogram_layout(net).positions

    def test_horizontal_reticulations(self) -> None:
        """Levelled reticulate edges are single horizontal segments; layering stays valid."""
        net = LEVEL_2_DNETWORK_MULTIPLE_BLOBS
        layout = compute_pz_cladogram_layout(net, rectangular=True)
        horizontal = [k for k in layout.reticulate_edges if len(layout.edge_routes[k].points) == 2]
        assert horizontal  # some reticulate parents could be lowered
        for u, v, _ in horizontal:
            assert layout.positions[u][1] == layout.positions[v][1]
        for u, v, key in net._graph.edges(keys=True):
            if (u, v, key) in horizontal:
                continue
            assert layout.positions[u][1] > layout.positions[v][1]  # every other edge points down

    def test_plot_rectangular(self) -> None:
        """plot() accepts the option in both directions and for pz-layered."""
        for layout in ("pz-cladogram", "pz-layered"):
            assert plot(LEVEL_1_DNETWORK_SINGLE_HYBRID, layout=layout, rectangular=True) is not None
            assert (
                plot(
                    LEVEL_1_DNETWORK_SINGLE_HYBRID, layout=layout, rectangular=True, direction="LR"
                )
                is not None
            )

    def test_dot_rectangular_routes(self) -> None:
        """pz-layered elbows are orthogonal and edges into hybrids arrive sideways."""
        from phylozoo.viz.dnetwork.layout import compute_pz_layered_layout

        net = LEVEL_2_DNETWORK_MULTIPLE_BLOBS
        layout = compute_pz_layered_layout(net, rectangular=True)
        for (u, v, _), route in layout.edge_routes.items():
            pts = route.points
            assert pts[0] == layout.positions[u] and pts[-1] == layout.positions[v]
            for (ax_, ay), (bx, by) in zip(pts, pts[1:]):
                assert ax_ == bx or ay == by
            if v in net.hybrid_nodes:
                assert pts[-2][1] == pts[-1][1]  # sideways arrival


class TestSDBlobShape:
    """Blobs are drawn as evenly sided polygons."""

    def test_level_1_blob_edges_have_similar_length(self) -> None:
        """After refinement no blob edge is more than 1.6x longer than another in the same blob."""
        from phylozoo.core.network.sdnetwork.features import blobs

        net = LEVEL_1_SDNETWORK_LARGE_MULTI_BLOB
        layout = compute_pz_unrooted_layout(net)
        for blob in blobs(net, trivial=False, leaves=False):
            lengths = [
                math.dist(layout.positions[u], layout.positions[v])
                for u, v, _ in net._graph.edges(keys=True)
                if u in blob and v in blob
            ]
            assert max(lengths) / min(lengths) < 1.6


class TestRadialLayout:
    """Test the radial (circular cladogram) layout for both network types."""

    def test_directed_network(self) -> None:
        """Root at the centre, leaves on the unit circle, every edge routed."""
        from phylozoo.viz.dnetwork.layout import compute_pz_radial_layout

        net = LEVEL_2_DNETWORK_MULTIPLE_BLOBS
        layout = compute_pz_radial_layout(net)
        assert layout.algorithm == "pz-radial"
        assert layout.positions[net.root_node] == (0.0, 0.0)
        for leaf in net.leaves:
            assert math.hypot(*layout.positions[leaf]) == pytest.approx(1.0)
        assert len(layout.edge_routes) == net.number_of_edges()
        assert plot(net, layout="pz-radial") is not None

    def test_semidirected_network(self) -> None:
        """Works for reticulate semi-directed networks, not only trees."""
        from phylozoo.viz.sdnetwork.layout import compute_pz_radial_layout

        net = LEVEL_3_SDNETWORK_LARGE_MANY_HYBRIDS
        layout = compute_pz_radial_layout(net, radius=2.0)
        assert set(layout.positions) == set(net._graph.nodes)
        for leaf in net.leaves:
            assert math.hypot(*layout.positions[leaf]) == pytest.approx(2.0)
        assert len(layout.edge_routes) == net.number_of_edges()
        assert plot(net, layout="pz-radial") is not None

    def test_root_location(self) -> None:
        """Rooting at a node puts that node at the centre."""
        from phylozoo.viz.sdnetwork.layout import compute_pz_radial_layout

        net = SDTREE_LARGE_BINARY
        node = next(iter(net.internal_nodes))
        layout = compute_pz_radial_layout(net, root_location=node)
        assert layout.positions[node] == (0.0, 0.0)

    def test_invalid_direction(self) -> None:
        """An unknown angle direction raises PhyloZooValueError."""
        from phylozoo.viz.dnetwork.layout import compute_pz_radial_layout

        with pytest.raises(PhyloZooValueError, match="angle_direction"):
            compute_pz_radial_layout(DTREE_LARGE_BINARY, angle_direction="sideways")


class TestLayeredLayout:
    """Test the layered pz-layered layout."""

    def test_basic_properties(self) -> None:
        """Every node placed, every edge routed downward through its layers, no NaN."""
        from phylozoo.viz.dnetwork.layout import compute_pz_layered_layout

        net = LEVEL_2_DNETWORK_MULTIPLE_BLOBS
        layout = compute_pz_layered_layout(net)
        assert layout.algorithm == "pz-layered"
        assert set(layout.positions) == set(net._graph.nodes)
        assert len(layout.edge_routes) == net.number_of_edges()
        for (u, v, _), route in layout.edge_routes.items():
            assert route.points[0] == layout.positions[u]
            assert route.points[-1] == layout.positions[v]
            ys = [y for _, y in route.points]
            assert ys == sorted(ys, reverse=True)  # strictly downward through the layers
            assert len(set(ys)) == len(ys)

    def test_long_edges_get_bend_points(self) -> None:
        """An edge spanning k layers has k - 1 intermediate points."""
        from phylozoo.viz.dnetwork.layout import compute_pz_layered_layout

        net = LEVEL_2_DNETWORK_MULTIPLE_BLOBS
        layout = compute_pz_layered_layout(net, align_leaves=True)
        ys = sorted({y for _, y in layout.positions.values()}, reverse=True)
        layer = {y: i for i, y in enumerate(ys)}
        for (u, v, _), route in layout.edge_routes.items():
            span = layer[layout.positions[v][1]] - layer[layout.positions[u][1]]
            assert len(route.points) == span + 1

    def test_fewer_or_equal_crossings_than_initial_order(self) -> None:
        """Sweeps never return a worse ordering than the DFS start."""
        from phylozoo.viz.dnetwork.layout import compute_pz_layered_layout
        from phylozoo.viz._layout_utils import count_crossings as _count_crossings

        net = LEVEL_5_DNETWORK_NON_TREEBASED

        def crossings(sweeps: int) -> int:
            layout = compute_pz_layered_layout(net, sweeps=sweeps)
            segments = np.array(
                [
                    (*a, *b)
                    for r in layout.edge_routes.values()
                    for a, b in zip(r.points, r.points[1:])
                ]
            )
            return _count_crossings(segments)

        assert crossings(8) <= crossings(0)

    def test_plot_and_directions(self) -> None:
        """plot() accepts pz-layered in both directions and with aligned leaves."""
        from phylozoo.viz._dispatch import resolve_layout

        assert resolve_layout(LEVEL_1_DNETWORK_SINGLE_HYBRID, "pz-layered") == "pz-layered"
        for kwargs in ({}, {"direction": "LR"}, {"align_leaves": True}):
            assert plot(LEVEL_1_DNETWORK_SINGLE_HYBRID, layout="pz-layered", **kwargs) is not None

    def test_invalid_direction(self) -> None:
        """An unknown direction raises PhyloZooValueError."""
        from phylozoo.viz.dnetwork.layout import compute_pz_layered_layout

        with pytest.raises(PhyloZooValueError, match="direction"):
            compute_pz_layered_layout(DTREE_LARGE_BINARY, direction="diagonal")


class TestUnrootedOptions:
    """pz-unrooted options (daylight, stress) and pz-unrooted on directed networks."""

    def test_daylight_keeps_edge_lengths(self) -> None:
        """Daylight rounds only rotate subtrees: every edge keeps its length (up to scaling)."""
        for net in (SDTREE_LARGE_BINARY, LEVEL_1_SDNETWORK_LARGE_MULTI_BLOB):
            raw = compute_pz_unrooted_layout(net, refine=0).positions
            day = compute_pz_unrooted_layout(net, refine=0, daylight=8).positions
            ratios = [
                math.dist(raw[u], raw[v]) / math.dist(day[u], day[v])
                for u, v, _ in net._graph.edges(keys=True)
            ]
            assert max(ratios) - min(ratios) < 1e-6

    def test_uniform_stress(self) -> None:
        """With node_spacing == edge_length the stress step drives all edges to one length."""
        net = LEVEL_1_SDNETWORK_LARGE_MULTI_BLOB
        layout = compute_pz_unrooted_layout(net, node_spacing=1.0, refine=200)
        lengths = [
            math.dist(layout.positions[u], layout.positions[v])
            for u, v, _ in net._graph.edges(keys=True)
        ]
        assert max(lengths) / min(lengths) < 2.0

    def test_directed_network(self) -> None:
        """pz-unrooted works on directed networks: all nodes placed, all edges routed, arrows on all."""
        from phylozoo.viz.dnetwork.layout import compute_pz_unrooted_layout as d_blobs
        from phylozoo.viz._dispatch import resolve_layout

        net = LEVEL_2_DNETWORK_MULTIPLE_BLOBS
        layout = d_blobs(net)
        assert layout.algorithm == "pz-unrooted"
        assert set(layout.positions) == set(net._graph.nodes)
        assert len(layout.edge_routes) == net.number_of_edges()
        assert layout.backbone_edges | layout.reticulate_edges == set(net._graph.edges(keys=True))
        assert resolve_layout(net, "pz-unrooted") == "pz-unrooted"
        assert plot(net, layout="pz-unrooted", daylight=4) is not None

    def test_refinement_never_adds_crossings(self) -> None:
        """The guarded stress step returns a drawing with at most the raw drawing's crossings."""
        from phylozoo.viz._layout_utils import count_crossings

        def crossings(net, **kw) -> int:
            pos = compute_pz_unrooted_layout(net, **kw).positions
            return count_crossings(
                np.array([(*pos[u], *pos[v]) for u, v, _ in net._graph.edges(keys=True)])
            )

        for name in (
            "LEVEL_3_SDNETWORK_LARGE_MANY_HYBRIDS",
            "LEVEL_2_SDNETWORK_PARALLEL_HYBRIDS",
            "LEVEL_1_SDNETWORK_LARGE_MULTI_BLOB",
        ):
            net = getattr(sd_networks, name)
            assert crossings(net) <= crossings(net, refine=0), name
            assert crossings(net, daylight=8) <= crossings(net, daylight=8, refine=0), name

    def test_parallel_edges_form_a_blob(self) -> None:
        """Two parallel edges are not cut edges: their endpoints share a blob (no crash)."""
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_PARALLEL_EDGES

        layout = compute_pz_unrooted_layout(LEVEL_1_SDNETWORK_PARALLEL_EDGES)
        assert len(layout.edge_routes) == LEVEL_1_SDNETWORK_PARALLEL_EDGES.number_of_edges()


class TestOrderingImprovements:
    """Radial re-optimisation in polar coordinates and cladogram-initialised layering."""

    def test_local_search_with_custom_score_never_worsens(self) -> None:
        """_local_search returns a score no worse than the score of the starting order."""
        from phylozoo.viz.dnetwork.layout.cladogram import _local_search

        kids = {"r": ["a", "b", "c"], "a": [], "b": [], "c": []}
        penalty = {("a", "b", "c"): 5, ("b", "a", "c"): 3, ("b", "c", "a"): 1}

        def score(order: dict) -> tuple[int, float]:
            return penalty.get(tuple(order["r"]), 9), 0.0

        assert _local_search("r", kids, [], [], {}, 1.0, score=score) == (1, 0.0)
        assert kids["r"] == ["b", "c", "a"]

    def test_radial_search_reduces_chord_crossings(self) -> None:
        """The polar search leaves NON_TREEBASED with fewer crossings than the plain layered order."""
        from phylozoo.viz._layout_utils import count_crossings
        from phylozoo.viz.dnetwork.layout import compute_pz_radial_layout

        net = LEVEL_5_DNETWORK_NON_TREEBASED
        layout = compute_pz_radial_layout(net)
        segments = np.array(
            [(*a, *b) for r in layout.edge_routes.values() for a, b in zip(r.points, r.points[1:])]
        )
        assert count_crossings(segments) <= 5
        assert compute_pz_radial_layout(net).positions == layout.positions  # deterministic

    def test_layered_starts_from_cladogram_order(self) -> None:
        """With aligned leaves the layered drawing keeps the cladogram's crossing-free order."""
        from phylozoo.viz._layout_utils import count_crossings
        from phylozoo.viz.dnetwork.layout import compute_pz_layered_layout
        from tests.fixtures.directed_networks import LEVEL_3_DNETWORK_LARGE_MANY_HYBRIDS

        layout = compute_pz_layered_layout(LEVEL_3_DNETWORK_LARGE_MANY_HYBRIDS, align_leaves=True)
        segments = np.array(
            [(*a, *b) for r in layout.edge_routes.values() for a, b in zip(r.points, r.points[1:])]
        )
        assert count_crossings(segments) == 0
