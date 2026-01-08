"""
Tests for styling in visualize2.

This test suite covers the NetworkStyle class and default styles.
"""

from phylozoo.visualize.dnetwork.styling import NetworkStyle, default_style


class TestNetworkStyle:
    """Test NetworkStyle class."""

    def test_default_values(self) -> None:
        """Test default style values."""
        style = NetworkStyle()

        assert style.node_color == 'lightblue'
        assert style.leaf_color == 'lightgreen'
        assert style.hybrid_color == 'salmon'
        assert style.edge_color == 'gray'
        assert style.hybrid_edge_color == 'red'
        assert style.with_labels is True

    def test_custom_values(self) -> None:
        """Test custom style values."""
        style = NetworkStyle(
            node_color='blue',
            leaf_color='green',
            edge_width=3.0,
        )

        assert style.node_color == 'blue'
        assert style.leaf_color == 'green'
        assert style.edge_width == 3.0

    def test_copy(self) -> None:
        """Test style copying."""
        style1 = NetworkStyle(node_color='blue')
        style2 = style1.copy()

        assert style2.node_color == 'blue'
        assert style2 is not style1


class TestDefaultStyles:
    """Test default style functions."""

    def test_default_style(self) -> None:
        """Test default_style function."""
        style = default_style()

        assert isinstance(style, NetworkStyle)
        assert style.node_color == 'lightblue'
