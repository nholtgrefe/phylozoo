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

        assert style.node_color == 'lightblue'
        assert style.leaf_color == 'lightblue'
        assert style.hybrid_color == 'lightblue'
        assert style.leaf_size is None
        assert style.edge_color == 'gray'
        assert style.hybrid_edge_color == 'red'
        assert style.with_labels is True

    def test_leaf_size_default_uses_node_size(self) -> None:
        """Test that leaf_size=None defaults to node_size for leaves."""
        from phylozoo.viz import plot
        from phylozoo.core.network.dnetwork import DirectedPhyNetwork

        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})],
        )
        style = DNetStyle(node_size=1000, leaf_size=None)
        ax = plot(net, style=style, layout='pz-dag', trials=100)
        assert ax is not None

    def test_custom_values(self) -> None:
        """Test custom style values."""
        style = DNetStyle(
            node_color='blue',
            leaf_color='green',
            edge_width=3.0,
        )

        assert style.node_color == 'blue'
        assert style.leaf_color == 'green'
        assert style.edge_width == 3.0

    def test_copy(self) -> None:
        """Test style copying."""
        style1 = DNetStyle(node_color='blue')
        style2 = style1.copy()

        assert style2.node_color == 'blue'
        assert style2 is not style1


class TestDefaultStyles:
    """Test default style functions."""

    def test_default_style(self) -> None:
        """Test default_style function."""
        style = default_style()

        assert isinstance(style, DNetStyle)
        assert style.node_color == 'lightblue'
