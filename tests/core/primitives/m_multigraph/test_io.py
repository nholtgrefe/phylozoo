"""
Tests for mixed multi-graph I/O functions.

This test suite covers the IOMixin integration for MixedMultiGraph, including
save, load, convert, and format conversion methods for phylozoo-dot format.
"""

import os
import tempfile

import pytest

from phylozoo.core.primitives.m_multigraph import MixedMultiGraph


class TestMixedMultiGraphPhyloZooDOTIO:
    """Test IOMixin methods for MixedMultiGraph with phylozoo-dot format."""

    def test_to_string(self) -> None:
        """Test converting graph to phylozoo-dot string."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2, weight=1.0)
        G.add_undirected_edge(2, 3, weight=2.0)

        pzdot_str = G.to_string()

        assert "graph" in pzdot_str
        assert "1 -> 2" in pzdot_str
        assert "2 -- 3" in pzdot_str

    def test_to_string_with_attributes(self) -> None:
        """Test phylozoo-dot string with node and edge attributes."""
        G = MixedMultiGraph(attributes={"name": "test"})
        G.add_node(1, label="Node1", shape="box")
        G.add_node(2, label="Node2")
        G.add_directed_edge(1, 2, weight=1.0, color="blue")
        G.add_undirected_edge(2, 3, weight=2.0, color="red")

        pzdot_str = G.to_string()

        assert "graph" in pzdot_str
        assert "name=test" in pzdot_str
        assert "label=Node1" in pzdot_str or 'label="Node1"' in pzdot_str
        assert "weight=1.0" in pzdot_str
        assert "weight=2.0" in pzdot_str

    def test_to_string_parallel_edges(self) -> None:
        """Test phylozoo-dot string with parallel edges."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2, weight=1.0)
        G.add_directed_edge(1, 2, weight=2.0)  # Parallel directed edge
        G.add_undirected_edge(2, 3, weight=3.0)
        G.add_undirected_edge(2, 3, weight=4.0)  # Parallel undirected edge

        pzdot_str = G.to_string()

        assert "graph" in pzdot_str
        assert pzdot_str.count("1 -> 2") == 2  # Two directed edges
        assert pzdot_str.count("2 -- 3") == 2  # Two undirected edges

    def test_save_basic(self) -> None:
        """Test basic phylozoo-dot file saving."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)
        G.add_undirected_edge(2, 3)

        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.pzdot")
            G.save(file_path, overwrite=True)

            assert os.path.exists(file_path)

            # Read and verify content
            with open(file_path, "r") as f:
                content = f.read()

            assert "graph" in content
            assert "1 -> 2" in content
            assert "2 -- 3" in content

    def test_from_string(self) -> None:
        """Test parsing phylozoo-dot string to graph."""
        pzdot_str = """graph {
    1 [label="Node1"];
    2 [label="Node2"];
    1 -> 2 [weight=1.0];
    2 -- 3 [weight=2.0];
}"""

        G = MixedMultiGraph.from_string(pzdot_str, format="phylozoo-dot")

        assert G.number_of_nodes() == 3
        assert G.number_of_edges() == 2
        assert 1 in G
        assert 2 in G
        assert 3 in G

    def test_from_string_with_attributes(self) -> None:
        """Test parsing phylozoo-dot string with attributes."""
        pzdot_str = """graph {
    name=test_graph;
    1 [label="Node1", shape=box];
    2 [label="Node2"];
    1 -> 2 [weight=1.0, color=blue];
    2 -- 3 [weight=2.0, color=red];
}"""

        G = MixedMultiGraph.from_string(pzdot_str, format="phylozoo-dot")

        assert G.number_of_nodes() == 3
        assert G.number_of_edges() == 2
        # Check graph attributes
        assert G._directed.graph.get("name") == "test_graph"
        # Check node attributes
        assert G._directed.nodes[1].get("label") == "Node1"
        assert G._directed.nodes[1].get("shape") == "box"

    def test_from_string_parallel_edges(self) -> None:
        """Test parsing phylozoo-dot string with parallel edges."""
        pzdot_str = """graph {
    1 -> 2 [weight=1.0, key=0];
    1 -> 2 [weight=2.0, key=1];
    2 -- 3 [weight=3.0, key=0];
    2 -- 3 [weight=4.0, key=1];
}"""

        G = MixedMultiGraph.from_string(pzdot_str, format="phylozoo-dot")

        assert G.number_of_edges() == 4
        assert G._directed.number_of_edges(1, 2) == 2
        assert G._undirected.number_of_edges(2, 3) == 2

    def test_load_basic(self) -> None:
        """Test loading graph from phylozoo-dot file."""
        pzdot_content = """graph {
    1 -> 2 [weight=1.0];
    2 -- 3 [weight=2.0];
}"""

        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.pzdot")
            with open(file_path, "w") as f:
                f.write(pzdot_content)

            G = MixedMultiGraph.load(file_path)

            assert G.number_of_nodes() == 3
            assert G.number_of_edges() == 2

    def test_round_trip(self) -> None:
        """Test round-trip conversion (to_string -> from_string)."""
        G = MixedMultiGraph(attributes={"name": "test"})
        G.add_node(1, label="Node1", shape="box")
        G.add_node(2, label="Node2")
        G.add_directed_edge(1, 2, weight=1.0, color="blue")
        G.add_undirected_edge(2, 3, weight=2.0, color="red")
        G.add_directed_edge(3, 4, weight=3.0)

        pzdot_str = G.to_string()
        G2 = MixedMultiGraph.from_string(pzdot_str, format="phylozoo-dot")

        assert G2.number_of_nodes() == G.number_of_nodes()
        assert G2.number_of_edges() == G.number_of_edges()
        assert G2._directed.number_of_edges(1, 2) == G._directed.number_of_edges(1, 2)
        assert G2._undirected.number_of_edges(2, 3) == G._undirected.number_of_edges(2, 3)

    def test_round_trip_parallel_edges(self) -> None:
        """Test round-trip with parallel edges."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2, weight=1.0)
        G.add_directed_edge(1, 2, weight=2.0)  # Parallel directed edge
        G.add_undirected_edge(2, 3, weight=3.0)
        G.add_undirected_edge(2, 3, weight=4.0)  # Parallel undirected edge

        pzdot_str = G.to_string()
        G2 = MixedMultiGraph.from_string(pzdot_str, format="phylozoo-dot")

        assert G2.number_of_edges() == G.number_of_edges()
        assert G2._directed.number_of_edges(1, 2) == 2
        assert G2._undirected.number_of_edges(2, 3) == 2

    def test_round_trip_file(self) -> None:
        """Test round-trip conversion via file."""
        G = MixedMultiGraph(attributes={"name": "test"})
        G.add_directed_edge(1, 2, weight=1.0)
        G.add_undirected_edge(2, 3, weight=2.0)

        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.pzdot")
            G.save(file_path, overwrite=True)

            G2 = MixedMultiGraph.load(file_path)

            assert G2.number_of_nodes() == G.number_of_nodes()
            assert G2.number_of_edges() == G.number_of_edges()

    def test_empty_graph(self) -> None:
        """Test empty graph I/O."""
        G = MixedMultiGraph()

        pzdot_str = G.to_string()
        G2 = MixedMultiGraph.from_string(pzdot_str, format="phylozoo-dot")

        assert G2.number_of_nodes() == 0
        assert G2.number_of_edges() == 0

    def test_single_node(self) -> None:
        """Test graph with single node."""
        G = MixedMultiGraph()
        G.add_node(1, label="Single")

        pzdot_str = G.to_string()
        G2 = MixedMultiGraph.from_string(pzdot_str, format="phylozoo-dot")

        assert G2.number_of_nodes() == 1
        assert G2.number_of_edges() == 0
        assert 1 in G2

    def test_string_node_ids(self) -> None:
        """Test graph with string node IDs."""
        G = MixedMultiGraph()
        G.add_directed_edge("A", "B", weight=1.0)
        G.add_undirected_edge("B", "C", weight=2.0)

        pzdot_str = G.to_string()
        G2 = MixedMultiGraph.from_string(pzdot_str, format="phylozoo-dot")

        assert G2.number_of_nodes() == 3
        assert G2.number_of_edges() == 2
        assert "A" in G2
        assert "B" in G2
        assert "C" in G2

    def test_mixed_node_types(self) -> None:
        """Test graph with mixed node ID types."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, "A", weight=1.0)
        G.add_undirected_edge("A", 2.5, weight=2.0)

        pzdot_str = G.to_string()
        G2 = MixedMultiGraph.from_string(pzdot_str, format="phylozoo-dot")

        assert G2.number_of_nodes() == 3
        assert G2.number_of_edges() == 2

    def test_graph_attributes(self) -> None:
        """Test graph-level attributes."""
        G = MixedMultiGraph(attributes={"name": "test_graph", "version": "1.0"})
        G.add_directed_edge(1, 2)

        pzdot_str = G.to_string()
        G2 = MixedMultiGraph.from_string(pzdot_str, format="phylozoo-dot")

        assert G2._directed.graph.get("name") == "test_graph"
        assert G2._directed.graph.get("version") == "1.0"

    def test_node_attributes(self) -> None:
        """Test node attributes preservation."""
        G = MixedMultiGraph()
        G.add_node(1, label="Node1", shape="box", color="red")
        G.add_node(2, label="Node2", shape="circle")
        G.add_directed_edge(1, 2)

        pzdot_str = G.to_string()
        G2 = MixedMultiGraph.from_string(pzdot_str, format="phylozoo-dot")

        assert G2._directed.nodes[1].get("label") == "Node1"
        assert G2._directed.nodes[1].get("shape") == "box"
        assert G2._directed.nodes[1].get("color") == "red"
        assert G2._directed.nodes[2].get("label") == "Node2"
        assert G2._directed.nodes[2].get("shape") == "circle"

    def test_edge_attributes(self) -> None:
        """Test edge attributes preservation."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2, weight=1.0, label="edge1", color="blue")
        G.add_undirected_edge(2, 3, weight=2.0, label="edge2", color="red")

        pzdot_str = G.to_string()
        G2 = MixedMultiGraph.from_string(pzdot_str, format="phylozoo-dot")

        # Check directed edge attributes
        edges_1_2 = [
            (u, v, k, d)
            for u, v, k, d in G2._directed.edges(keys=True, data=True)
            if u == 1 and v == 2
        ]
        assert len(edges_1_2) == 1
        assert edges_1_2[0][3].get("weight") == 1.0
        assert edges_1_2[0][3].get("label") == "edge1"
        assert edges_1_2[0][3].get("color") == "blue"

        # Check undirected edge attributes
        edges_2_3 = [
            (u, v, k, d)
            for u, v, k, d in G2._undirected.edges(keys=True, data=True)
            if (u == 2 and v == 3) or (u == 3 and v == 2)
        ]
        assert len(edges_2_3) == 1
        assert edges_2_3[0][3].get("weight") == 2.0
        assert edges_2_3[0][3].get("label") == "edge2"
        assert edges_2_3[0][3].get("color") == "red"

    def test_invalid_format(self) -> None:
        """Test error handling for invalid format."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)

        from phylozoo.utils.exceptions import PhyloZooFormatError

        with pytest.raises(PhyloZooFormatError, match="Format 'invalid' not supported"):
            G.to_string(format="invalid")

    def test_malformed_string(self) -> None:
        """Test error handling for malformed phylozoo-dot string."""
        from phylozoo.utils.exceptions import PhyloZooParseError

        malformed = "not a valid graph"

        with pytest.raises(PhyloZooParseError, match="Could not find graph declaration"):
            MixedMultiGraph.from_string(malformed, format="phylozoo-dot")

    def test_unmatched_braces(self) -> None:
        """Test error handling for unmatched braces."""
        from phylozoo.utils.exceptions import PhyloZooParseError

        malformed = "graph { 1 -> 2"  # Missing closing brace

        with pytest.raises(PhyloZooParseError, match="Unmatched braces"):
            MixedMultiGraph.from_string(malformed, format="phylozoo-dot")

    def test_convert(self) -> None:
        """Test format conversion between files."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2, weight=1.0)
        G.add_undirected_edge(2, 3, weight=2.0)

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "input.pzdot")
            output_path = os.path.join(tmpdir, "output.pzdot")

            G.save(input_path, overwrite=True)
            MixedMultiGraph.convert(input_path, output_path, overwrite=True)

            G2 = MixedMultiGraph.load(output_path)
            assert G2.number_of_nodes() == G.number_of_nodes()
            assert G2.number_of_edges() == G.number_of_edges()

    def test_convert_string(self) -> None:
        """Test format conversion between strings."""
        pzdot_str = """graph {
    1 -> 2 [weight=1.0];
    2 -- 3 [weight=2.0];
}"""

        # Convert to same format (should be identity)
        result = MixedMultiGraph.convert_string(pzdot_str, "phylozoo-dot", "phylozoo-dot")

        G1 = MixedMultiGraph.from_string(pzdot_str, format="phylozoo-dot")
        G2 = MixedMultiGraph.from_string(result, format="phylozoo-dot")

        assert G2.number_of_nodes() == G1.number_of_nodes()
        assert G2.number_of_edges() == G1.number_of_edges()

    def test_auto_detect_format(self) -> None:
        """Test automatic format detection from file extension."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)

        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.pzdot")
            G.save(file_path, overwrite=True)

            # Load without specifying format (should auto-detect)
            G2 = MixedMultiGraph.load(file_path)
            assert G2.number_of_edges() == 1

    def test_explicit_format_override(self) -> None:
        """Test explicit format override even with extension."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)

        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.txt")
            G.save(file_path, format="phylozoo-dot", overwrite=True)

            # Load with explicit format
            G2 = MixedMultiGraph.load(file_path, format="phylozoo-dot")
            assert G2.number_of_edges() == 1


class TestPhyloZooDOTRoundTrip:
    """Round-trip fidelity for phylozoo-dot (no spurious node labels)."""

    @staticmethod
    def _dsig(g: MixedMultiGraph) -> list:
        return sorted(g.directed_edges_iter(keys=True))

    @staticmethod
    def _usig(g: MixedMultiGraph) -> list:
        return sorted((tuple(sorted((a, b))), k) for a, b, k in g.undirected_edges_iter(keys=True))

    def test_labelless_nodes_have_no_spurious_attribute(self) -> None:
        """Nodes without a label must not gain a (non-string) ``label`` on round-trip.

        The writer used to inject ``label=<node id>``; the reader parsed it back as
        an int, leaving every node with a non-string label that broke downstream
        label validation (e.g. ``attach_leaves_to_generator``).
        """
        G = MixedMultiGraph(directed_edges=[(1, 2), (3, 1)], undirected_edges=[(2, 3)])
        G2 = MixedMultiGraph.from_string(G.to_string("phylozoo-dot"), "phylozoo-dot")

        assert "label=" not in G.to_string("phylozoo-dot")
        assert all(data == {} for _, data in G2.nodes(data=True))

    def test_roundtrip_preserves_edges_including_parallel(self) -> None:
        """Directed and undirected edges, with parallel multiplicities, round-trip."""
        G = MixedMultiGraph(
            directed_edges=[(1, 2), (1, 2)],  # parallel directed
            undirected_edges=[(3, 4), (3, 4), (2, 3)],  # parallel undirected
        )
        G2 = MixedMultiGraph.from_string(G.to_string("phylozoo-dot"), "phylozoo-dot")

        assert self._dsig(G) == self._dsig(G2)
        assert self._usig(G) == self._usig(G2)
        assert {type(n).__name__ for n in G2.nodes()} == {"int"}

    def test_real_string_label_is_preserved(self) -> None:
        """A genuine string label survives the round-trip unchanged."""
        G = MixedMultiGraph(directed_edges=[(1, 2)])
        G.add_node(1, label="taxonA")
        G2 = MixedMultiGraph.from_string(G.to_string("phylozoo-dot"), "phylozoo-dot")
        assert dict(G2.nodes(data=True))[1].get("label") == "taxonA"

    def test_numeric_label_is_read_as_string(self) -> None:
        """An (already-saved) ``label=<int>`` is read back as a string, not an int."""
        pzdot = "graph {\n  1 [label=1];\n  2 [label=2];\n  1 -> 2;\n}"
        G = MixedMultiGraph.from_string(pzdot, "phylozoo-dot")
        labels = {data["label"] for _, data in G.nodes(data=True) if "label" in data}
        assert labels == {"1", "2"}
        assert all(isinstance(lbl, str) for lbl in labels)


class TestMixedMultiGraphStandardDOTIO:
    """The standard ``dot`` format for MixedMultiGraph (digraph + dir=none)."""

    @staticmethod
    def _dsig(g: MixedMultiGraph) -> list:
        return sorted(g.directed_edges_iter(keys=True))

    @staticmethod
    def _usig(g: MixedMultiGraph) -> list:
        return sorted((tuple(sorted((a, b))), k) for a, b, k in g.undirected_edges_iter(keys=True))

    def test_dot_is_a_supported_format(self) -> None:
        assert "dot" in MixedMultiGraph._supported_formats
        assert MixedMultiGraph._default_format == "phylozoo-dot"

    def test_to_dot_uses_digraph_with_dir_none(self) -> None:
        """Undirected edges become ``-> [dir=none]`` inside a ``digraph``."""
        G = MixedMultiGraph(directed_edges=[(1, 2)], undirected_edges=[(2, 3)])
        dot_str = G.to_string("dot")
        assert dot_str.lstrip().startswith("digraph")
        assert "1 -> 2" in dot_str
        assert "2 -> 3 [dir=none]" in dot_str
        assert "--" not in dot_str  # undirected is encoded via dir=none, not '--'

    def test_round_trip_directed_and_undirected_with_parallel(self) -> None:
        G = MixedMultiGraph(
            directed_edges=[(1, 2), (1, 2)],  # parallel directed
            undirected_edges=[(3, 4), (3, 4), (2, 3)],  # parallel undirected
        )
        G2 = MixedMultiGraph.from_string(G.to_string("dot"), "dot")
        assert self._dsig(G) == self._dsig(G2)
        assert self._usig(G) == self._usig(G2)

    def test_no_spurious_node_labels(self) -> None:
        G = MixedMultiGraph(directed_edges=[(1, 2)], undirected_edges=[(2, 3)])
        G2 = MixedMultiGraph.from_string(G.to_string("dot"), "dot")
        assert all(data == {} for _, data in G2.nodes(data=True))

    def test_pzdot_and_dot_describe_the_same_graph(self) -> None:
        """Loading a graph via ``phylozoo-dot`` and via ``dot`` yields the same edges."""
        G = MixedMultiGraph(directed_edges=[(1, 2), (4, 5)], undirected_edges=[(2, 3), (3, 1)])
        via_pz = MixedMultiGraph.from_string(G.to_string("phylozoo-dot"), "phylozoo-dot")
        via_dot = MixedMultiGraph.from_string(G.to_string("dot"), "dot")
        assert self._dsig(via_pz) == self._dsig(via_dot)
        assert self._usig(via_pz) == self._usig(via_dot)

    def test_dot_output_is_valid_graphviz(self) -> None:
        """The ``dot`` output parses with Graphviz (unlike ``phylozoo-dot``)."""
        pgv = pytest.importorskip("pygraphviz")
        G = MixedMultiGraph(directed_edges=[(1, 2), (1, 2)], undirected_edges=[(2, 3)])
        pgv.AGraph(string=G.to_string("dot"))  # must not raise
        with pytest.raises(Exception):
            pgv.AGraph(string=G.to_string("phylozoo-dot"))  # mixed -> / -- is not valid DOT
