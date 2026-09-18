"""
Tests for styling in viz.

This test suite covers the DNetStyle class and default styles.
"""

from phylozoo.viz.dnetwork import DNetStyle, default_style


class TestDNetStyle:
    """Test DNetStyle class."""

    def test_default_values(self) -> None:
        """Test default style values."""
        style = DNetStyle()

        # Defaults match SDNetStyle so directed and semi-directed plots look alike.
        assert style.node_color == "white"
        assert style.leaf_color == "#0a0a0a"
        assert style.hybrid_color == "#fcc0bc"
        assert style.leaf_size == 100.0
        assert style.edge_color == "gray"
        assert style.hybrid_edge_color == "red"
        assert style.arrow_head_size == 12.0
        assert style.with_labels is True
        assert style.label_rotation is None
        assert style.arrows is None

    def test_matches_sdnet_style(self) -> None:
        """Directed and semi-directed default styles share their colours."""
        from phylozoo.viz.sdnetwork import SDNetStyle

        d, sd = DNetStyle(), SDNetStyle()
        for attr in ("node_color", "leaf_color", "hybrid_color", "hybrid_edge_color", "edge_color"):
            assert getattr(d, attr) == getattr(sd, attr)

    def test_copy_keeps_arrows(self) -> None:
        """copy() carries the arrows and label_rotation options."""
        style = DNetStyle(arrows="all", label_rotation=0.0).copy()
        assert style.arrows == "all"
        assert style.label_rotation == 0.0

    def test_leaf_size_default_uses_node_size(self) -> None:
        """Test that leaf_size=None defaults to node_size for leaves."""
        from phylozoo.viz import plot
        from phylozoo.core.network.dnetwork import DirectedPhyNetwork

        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"})],
        )
        style = DNetStyle(node_size=1000, leaf_size=None)
        ax = plot(net, style=style, layout="pz-cladogram", trials=100)
        assert ax is not None

    def test_custom_values(self) -> None:
        """Test custom style values."""
        style = DNetStyle(
            node_color="blue",
            leaf_color="green",
            edge_width=3.0,
            arrow_head_size=24.0,
        )

        assert style.node_color == "blue"
        assert style.leaf_color == "green"
        assert style.edge_width == 3.0
        assert style.arrow_head_size == 24.0

    def test_copy(self) -> None:
        """Test style copying."""
        style1 = DNetStyle(node_color="blue")
        style2 = style1.copy()

        assert style2.node_color == "blue"
        assert style2.arrow_head_size == style1.arrow_head_size
        assert style2 is not style1


class TestDefaultStyles:
    """Test default style functions."""

    def test_default_style(self) -> None:
        """Test default_style function."""
        style = default_style()

        assert isinstance(style, DNetStyle)
        assert style.node_color == "white"


class TestLabelPlacement:
    """Fixed-rotation labels sit beside the node on its outward side."""

    def test_fixed_rotation_alignment(self) -> None:
        """Rightward direction -> left-aligned, vertically centred; downward -> centred, top."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        from phylozoo.viz._render import draw_label

        _, ax = plt.subplots()
        style = DNetStyle(label_rotation=0.0)
        right = draw_label(ax, (1.0, 0.0), "x", style, anchor=(0.0, 0.0), rotation=0.0)
        assert (right.get_ha(), right.get_va()) == ("left", "center")
        down = draw_label(ax, (0.0, -1.0), "x", style, anchor=(0.0, 0.0), rotation=0.0)
        assert (down.get_ha(), down.get_va()) == ("center", "top")
        auto = draw_label(ax, (1.0, 0.0), "x", style, anchor=(0.0, 0.0), rotation=None)
        assert (auto.get_ha(), auto.get_va()) == ("left", "center")
        plt.close("all")


class TestRootStyle:
    """root_color / root_size options."""

    def test_defaults_and_copy(self) -> None:
        """Unset by default and carried by copy()."""
        style = DNetStyle()
        assert style.root_color is None and style.root_size is None
        copied = DNetStyle(root_color="gold", root_size=300.0).copy()
        assert (copied.root_color, copied.root_size) == ("gold", 300.0)

    def test_render_uses_root_style(self) -> None:
        """The root marker takes root_color/root_size when set, node values otherwise."""
        from phylozoo.viz._render import _get_node_color, _get_node_size

        plain = DNetStyle()
        assert _get_node_color("root", plain) == plain.node_color
        assert _get_node_size("root", plain) == plain.node_size
        styled = DNetStyle(root_color="gold", root_size=300.0)
        assert _get_node_color("root", styled) == "gold"
        assert _get_node_size("root", styled) == 300.0

    def test_unrooted_layout_marks_root(self) -> None:
        """pz-unrooted draws the root larger than the other internal nodes; other layouts do not."""
        from phylozoo.viz import plot
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID as net

        def radii(layout: str) -> dict[float, int]:
            ax = plot(net, layout=layout)
            counts: dict[float, int] = {}
            for patch in ax.patches:
                if hasattr(patch, "get_radius"):
                    counts[round(patch.get_radius(), 6)] = (
                        counts.get(round(patch.get_radius(), 6), 0) + 1
                    )
            return counts

        unrooted, cladogram = radii("pz-unrooted"), radii("pz-cladogram")
        assert (
            len(unrooted) == len(cladogram) + 1
        )  # one extra (root) radius in the unrooted drawing
        assert unrooted[max(unrooted)] == 1
