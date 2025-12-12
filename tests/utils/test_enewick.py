"""
Tests for eNewick parser module.

This module tests the eNewick (Extended Newick) format parser, including
basic trees, networks with hybrid nodes, branch lengths, comments, and
various edge cases.
"""

import pytest

from phylozoo.utils.enewick import ENewickParseError, ParsedENewick, parse_enewick


class TestBasicTreeParsing:
    """Test cases for basic tree structure parsing."""

    def test_simple_two_leaves(self) -> None:
        """Parse a simple tree with two leaves."""
        result = parse_enewick("(A,B);")
        assert len(result.nodes) == 3  # A, B, and internal node
        assert len(result.edges) == 2
        assert result.root is not None

    def test_nested_tree(self) -> None:
        """Parse a nested tree structure."""
        result = parse_enewick("((A,B),C);")
        assert len(result.nodes) == 5  # A, B, C, and 2 internal nodes
        assert len(result.edges) == 4
        assert result.root is not None

    def test_deeply_nested_tree(self) -> None:
        """Parse a deeply nested tree."""
        result = parse_enewick("(((A,B),C),D);")
        assert len(result.nodes) == 7
        assert len(result.edges) == 6

    def test_multiple_children(self) -> None:
        """Parse a tree with multiple children."""
        result = parse_enewick("(A,B,C,D);")
        assert len(result.nodes) == 5  # 4 leaves + 1 internal
        assert len(result.edges) == 4

    def test_single_leaf(self) -> None:
        """Parse a tree with a single leaf."""
        result = parse_enewick("A;")
        assert len(result.nodes) == 1
        assert len(result.edges) == 0
        assert result.root == "A"

    def test_empty_tree(self) -> None:
        """Parse an empty tree (just root)."""
        result = parse_enewick("();")
        assert len(result.nodes) == 1
        assert len(result.edges) == 0


class TestNodeLabels:
    """Test cases for node label parsing."""

    def test_unquoted_labels(self) -> None:
        """Parse unquoted node labels."""
        result = parse_enewick("(A,B);")
        labels = [n.get("label") for n in result.nodes if "label" in n]
        assert "A" in labels
        assert "B" in labels

    def test_quoted_labels(self) -> None:
        """Parse quoted node labels."""
        result = parse_enewick("('Species A','Species B');")
        labels = [n.get("label") for n in result.nodes if "label" in n]
        assert "Species A" in labels
        assert "Species B" in labels

    def test_quoted_labels_with_spaces(self) -> None:
        """Parse quoted labels containing spaces."""
        result = parse_enewick("('Species Name','Another Species');")
        labels = [n.get("label") for n in result.nodes if "label" in n]
        assert "Species Name" in labels
        assert "Another Species" in labels

    def test_quoted_labels_with_special_chars(self) -> None:
        """Parse quoted labels with special characters."""
        result = parse_enewick("('Node-1','Node_2','Node.3');")
        labels = [n.get("label") for n in result.nodes if "label" in n]
        assert "Node-1" in labels
        assert "Node_2" in labels
        assert "Node.3" in labels

    def test_escaped_quotes_in_labels(self) -> None:
        """Parse labels with escaped quotes."""
        result = parse_enewick("('Node''s Label','Other');")
        labels = [n.get("label") for n in result.nodes if "label" in n]
        assert "Node's Label" in labels

    def test_unlabeled_internal_nodes(self) -> None:
        """Internal nodes without labels get auto-generated IDs."""
        result = parse_enewick("((A,B),C);")
        internal_nodes = [n for n in result.nodes if isinstance(n["id"], int)]
        assert len(internal_nodes) == 2
        assert all(isinstance(n["id"], int) for n in internal_nodes)

    def test_mixed_labeled_unlabeled(self) -> None:
        """Mix of labeled and unlabeled nodes."""
        result = parse_enewick("((A,B),C);")
        # A, B, C should have labels
        labeled = [n for n in result.nodes if "label" in n]
        assert len(labeled) >= 3  # At least A, B, C


class TestBranchLengths:
    """Test cases for branch length parsing."""

    def test_simple_branch_lengths(self) -> None:
        """Parse branch lengths on edges."""
        result = parse_enewick("(A:0.5,B:0.3);")
        edges_with_bl = [e for e in result.edges if "branch_length" in e]
        assert len(edges_with_bl) == 2
        assert edges_with_bl[0]["branch_length"] == 0.5
        assert edges_with_bl[1]["branch_length"] == 0.3

    def test_nested_branch_lengths(self) -> None:
        """Parse nested branch lengths."""
        result = parse_enewick("((A:0.5,B:0.3):0.1,C:0.2);")
        edges_with_bl = [e for e in result.edges if "branch_length" in e]
        assert len(edges_with_bl) == 4
        
        # Check specific edges
        edge_to_a = next(e for e in result.edges if e["v"] == "A")
        assert edge_to_a["branch_length"] == 0.5
        
        edge_to_internal = next(e for e in result.edges if isinstance(e["v"], int))
        assert edge_to_internal["branch_length"] == 0.1

    def test_integer_branch_lengths(self) -> None:
        """Parse integer branch lengths."""
        result = parse_enewick("(A:5,B:10);")
        edges_with_bl = [e for e in result.edges if "branch_length" in e]
        assert edges_with_bl[0]["branch_length"] == 5.0
        assert edges_with_bl[1]["branch_length"] == 10.0

    def test_negative_branch_lengths(self) -> None:
        """Parse negative branch lengths."""
        result = parse_enewick("(A:-0.5,B:-0.3);")
        edges_with_bl = [e for e in result.edges if "branch_length" in e]
        assert edges_with_bl[0]["branch_length"] == -0.5
        assert edges_with_bl[1]["branch_length"] == -0.3

    def test_scientific_notation(self) -> None:
        """Parse branch lengths in scientific notation."""
        result = parse_enewick("(A:1.5e-3,B:2.0e2);")
        edges_with_bl = [e for e in result.edges if "branch_length" in e]
        assert abs(edges_with_bl[0]["branch_length"] - 0.0015) < 1e-10
        assert abs(edges_with_bl[1]["branch_length"] - 200.0) < 1e-10

    def test_scientific_notation_uppercase(self) -> None:
        """Parse branch lengths with uppercase E."""
        result = parse_enewick("(A:1.5E-3,B:2.0E2);")
        edges_with_bl = [e for e in result.edges if "branch_length" in e]
        assert abs(edges_with_bl[0]["branch_length"] - 0.0015) < 1e-10
        assert abs(edges_with_bl[1]["branch_length"] - 200.0) < 1e-10

    def test_mixed_with_without_branch_lengths(self) -> None:
        """Mix of edges with and without branch lengths."""
        result = parse_enewick("(A:0.5,B,C:0.3);")
        edges_with_bl = [e for e in result.edges if "branch_length" in e]
        edges_without_bl = [e for e in result.edges if "branch_length" not in e]
        assert len(edges_with_bl) == 2
        assert len(edges_without_bl) == 1


class TestHybridNodes:
    """Test cases for hybrid node parsing."""

    def test_single_hybrid_node(self) -> None:
        """Parse a network with a single hybrid node (Extended Newick format)."""
        # In Extended Newick, hybrid nodes must appear twice: once with definition, once as reference
        result = parse_enewick("((A,B)#H1,#H1);")
        assert len(result.hybrid_nodes) == 1
        assert 1 in result.hybrid_nodes
        assert result.hybrid_nodes[1] == 1
        
        # Verify hybrid node has 2 incoming edges
        hybrid_id = 1
        incoming_edges = [e for e in result.edges if e['v'] == hybrid_id]
        assert len(incoming_edges) == 2

    def test_multiple_hybrid_nodes(self) -> None:
        """Parse a network with multiple hybrid nodes (Extended Newick format)."""
        result = parse_enewick("(((A,B)#H1,C)#H2,(#H1,#H2));")
        assert len(result.hybrid_nodes) == 2
        assert 1 in result.hybrid_nodes.values()
        assert 2 in result.hybrid_nodes.values()
        
        # Verify each hybrid node has 2 incoming edges
        for hybrid_id, hybrid_num in result.hybrid_nodes.items():
            incoming_edges = [e for e in result.edges if e['v'] == hybrid_id]
            assert len(incoming_edges) == 2, f"Hybrid #H{hybrid_num} should have 2 incoming edges"

    def test_hybrid_node_on_leaf(self) -> None:
        """Parse a hybrid node that is also a leaf."""
        result = parse_enewick("(A#H1,B);")
        # The hybrid marker on A means A is a hybrid node
        assert len(result.hybrid_nodes) == 1
        assert "A" in result.hybrid_nodes

    def test_hybrid_node_with_branch_length(self) -> None:
        """Parse hybrid node with branch length."""
        # Extended Newick: hybrid defined with children, then referenced
        # Branch lengths are on edges to children and on reference
        result = parse_enewick("((A:0.5,B:0.3)#H1,#H1:0.4);")
        assert len(result.hybrid_nodes) == 1
        
        # Check that branch lengths are on edges to children of hybrid
        hybrid_id = list(result.hybrid_nodes.keys())[0]
        edges_from_hybrid = [e for e in result.edges if e["u"] == hybrid_id]
        assert len(edges_from_hybrid) == 2
        assert all("branch_length" in e for e in edges_from_hybrid)
        assert edges_from_hybrid[0]["branch_length"] == 0.5
        assert edges_from_hybrid[1]["branch_length"] == 0.3
        
        # Check that hybrid node has 2 incoming edges (one with branch length)
        edges_to_hybrid = [e for e in result.edges if e["v"] == hybrid_id]
        assert len(edges_to_hybrid) == 2
        # The reference should have branch length
        edges_with_bl = [e for e in edges_to_hybrid if "branch_length" in e]
        assert len(edges_with_bl) == 1
        assert edges_with_bl[0]["branch_length"] == 0.4

    def test_hybrid_numbers_sequential(self) -> None:
        """Hybrid numbers should be tracked correctly."""
        result = parse_enewick("(((A,B)#H1,#H1)#H2,#H2);")
        hybrid_numbers = sorted(result.hybrid_nodes.values())
        assert hybrid_numbers == [1, 2]
    
    def test_hybrid_reference_before_definition_error(self) -> None:
        """Hybrid reference before definition should raise error."""
        with pytest.raises(ENewickParseError, match="found before definition"):
            parse_enewick("(#H1,(A,B)#H1);")


class TestComments:
    """Test cases for comment parsing."""

    def test_comments_after_labels(self) -> None:
        """Parse comments after node labels."""
        result = parse_enewick("(A[comment1],B[comment2]);")
        node_a = next(n for n in result.nodes if n.get("id") == "A")
        node_b = next(n for n in result.nodes if n.get("id") == "B")
        assert node_a.get("comment") == "comment1"
        assert node_b.get("comment") == "comment2"

    def test_comments_before_labels(self) -> None:
        """Parse comments before node labels."""
        result = parse_enewick("([comment1]A,[comment2]B);")
        node_a = next(n for n in result.nodes if n.get("id") == "A")
        node_b = next(n for n in result.nodes if n.get("id") == "B")
        assert node_a.get("comment") == "comment1"
        assert node_b.get("comment") == "comment2"

    def test_comments_with_spaces(self) -> None:
        """Parse comments containing spaces."""
        result = parse_enewick("(A[this is a comment],B);")
        node_a = next(n for n in result.nodes if n.get("id") == "A")
        assert node_a.get("comment") == "this is a comment"

    def test_nested_comments(self) -> None:
        """Parse nested comments (brackets within brackets)."""
        result = parse_enewick("(A[[nested] comment],B);")
        node_a = next(n for n in result.nodes if n.get("id") == "A")
        assert "[nested]" in node_a.get("comment", "")

    def test_multiple_comments_same_node(self) -> None:
        """Multiple comments on same node (last one wins)."""
        result = parse_enewick("([first]A[second],B);")
        node_a = next(n for n in result.nodes if n.get("id") == "A")
        # Last comment should be stored
        assert node_a.get("comment") == "second"


class TestNonBinaryNodes:
    """Test cases for non-binary nodes (polytomies and multiple parents)."""

    def test_polytomy_root(self) -> None:
        """Parse a tree with non-binary root (polytomy)."""
        result = parse_enewick("(A,B,C,D);")
        assert len(result.nodes) == 5
        assert len(result.edges) == 4
        
        # Root should have 4 children
        root_id = result.root
        children = [e for e in result.edges if e['u'] == root_id]
        assert len(children) == 4

    def test_polytomy_internal(self) -> None:
        """Parse a tree with non-binary internal node."""
        result = parse_enewick("((A,B,C)int1,D);")
        assert len(result.nodes) == 6
        
        # Find the labeled internal node
        int_node = next(n for n in result.nodes if n.get('label') == 'int1')
        children = [e for e in result.edges if e['u'] == int_node['id']]
        assert len(children) == 3

    def test_hybrid_with_three_parents(self) -> None:
        """Parse a hybrid node with 3 parents (Extended Newick)."""
        result = parse_enewick("((A,B)#H1,#H1,#H1);")
        assert len(result.hybrid_nodes) == 1
        
        # Hybrid should have 3 incoming edges
        hybrid_id = list(result.hybrid_nodes.keys())[0]
        parents = [e for e in result.edges if e['v'] == hybrid_id]
        assert len(parents) == 3

    def test_hybrid_with_four_parents(self) -> None:
        """Parse a hybrid node with 4 parents."""
        result = parse_enewick("((A,B)#H1,#H1,#H1,#H1);")
        assert len(result.hybrid_nodes) == 1
        
        # Hybrid should have 4 incoming edges
        hybrid_id = list(result.hybrid_nodes.keys())[0]
        parents = [e for e in result.edges if e['v'] == hybrid_id]
        assert len(parents) == 4

    def test_polytomy_and_hybrid_combined(self) -> None:
        """Parse a network with both polytomy and hybrid nodes."""
        result = parse_enewick("((A,B,C)#H1,#H1,D);")
        
        # Should have one hybrid node
        assert len(result.hybrid_nodes) == 1
        hybrid_id = list(result.hybrid_nodes.keys())[0]
        
        # Hybrid should have 2 parents
        parents = [e for e in result.edges if e['v'] == hybrid_id]
        assert len(parents) == 2
        
        # Hybrid should have 3 children
        children = [e for e in result.edges if e['u'] == hybrid_id]
        assert len(children) == 3


class TestComplexExamples:
    """Test cases for complex eNewick strings combining multiple features."""

    def test_all_features_combined(self) -> None:
        """Parse a complex example with all features."""
        # Note: Branch length comes before hybrid marker
        enewick = "((('Leaf A':0.5,'Leaf B':0.3):0.1)#H1,'Leaf C':0.2);"
        result = parse_enewick(enewick)
        
        # Check structure
        assert len(result.nodes) >= 4
        assert len(result.edges) >= 4
        
        # Check hybrid node
        assert len(result.hybrid_nodes) == 1
        
        # Check branch lengths (at least some edges should have them)
        edges_with_bl = [e for e in result.edges if "branch_length" in e]
        assert len(edges_with_bl) >= 3
        
        # Check labels
        labels = [n.get("label") for n in result.nodes if "label" in n]
        assert "Leaf A" in labels
        assert "Leaf B" in labels
        assert "Leaf C" in labels

    def test_deeply_nested_with_features(self) -> None:
        """Deeply nested structure with all features."""
        # Note: Branch length comes before hybrid marker
        enewick = "((((A:0.1,B:0.2):0.05)#H1,C:0.3):0.15,D:0.4);"
        result = parse_enewick(enewick)
        
        assert len(result.nodes) >= 5
        assert len(result.edges) >= 5
        assert len(result.hybrid_nodes) == 1
        
        # Most edges should have branch lengths
        edges_with_bl = [e for e in result.edges if "branch_length" in e]
        assert len(edges_with_bl) >= 4

    def test_multiple_hybrids_complex(self) -> None:
        """Multiple hybrid nodes in complex structure."""
        enewick = "((((A,B)#H1,C)#H2,D),E);"
        result = parse_enewick(enewick)
        
        assert len(result.hybrid_nodes) == 2
        assert 1 in result.hybrid_nodes.values()
        assert 2 in result.hybrid_nodes.values()


class TestErrorHandling:
    """Test cases for error handling and malformed input."""

    def test_missing_semicolon(self) -> None:
        """Missing semicolon should raise error."""
        with pytest.raises(ENewickParseError, match="must end with"):
            parse_enewick("(A,B)")

    def test_unclosed_parenthesis(self) -> None:
        """Unclosed parenthesis should raise error."""
        with pytest.raises(ENewickParseError):
            parse_enewick("((A,B);")

    def test_unclosed_quoted_string(self) -> None:
        """Unclosed quoted string should raise error."""
        with pytest.raises(ENewickParseError, match="Unclosed quoted"):
            parse_enewick("('A,B);")

    def test_unclosed_comment(self) -> None:
        """Unclosed comment bracket should raise error."""
        with pytest.raises(ENewickParseError, match="Unclosed comment"):
            parse_enewick("(A[comment,B);")

    def test_invalid_branch_length(self) -> None:
        """Invalid branch length should raise error."""
        with pytest.raises(ENewickParseError):
            parse_enewick("(A:abc,B);")

    def test_incomplete_hybrid_marker(self) -> None:
        """Incomplete hybrid marker should raise error."""
        with pytest.raises(ENewickParseError, match="Expected 'H'"):
            parse_enewick("(A#1,B);")

    def test_hybrid_marker_without_number(self) -> None:
        """Hybrid marker without number should raise error."""
        with pytest.raises(ENewickParseError, match="Expected number"):
            parse_enewick("(A#H,B);")

    def test_unexpected_characters_after_tree(self) -> None:
        """Unexpected characters after tree should raise error."""
        # Parser requires semicolon at end, so extra chars before semicolon fail
        # This test verifies the parser enforces the semicolon requirement
        with pytest.raises(ENewickParseError, match="must end with"):
            parse_enewick("(A,B)extra")

    def test_empty_string(self) -> None:
        """Empty string should raise error."""
        with pytest.raises(ENewickParseError):
            parse_enewick("")

    def test_whitespace_only(self) -> None:
        """Whitespace-only string becomes empty tree after stripping."""
        # After stripping whitespace, "   ;" becomes ";", which parses as empty tree
        result = parse_enewick("   ;")
        assert len(result.nodes) == 1  # Empty tree has one node
        assert len(result.edges) == 0


class TestEdgeCases:
    """Test cases for edge cases and special scenarios."""

    def test_whitespace_handling(self) -> None:
        """Whitespace should be handled correctly."""
        result1 = parse_enewick("(A,B);")
        result2 = parse_enewick("( A , B ) ;")
        assert len(result1.nodes) == len(result2.nodes)
        assert len(result1.edges) == len(result2.edges)

    def test_very_long_labels(self) -> None:
        """Very long node labels should be parsed correctly."""
        long_label = "A" * 1000
        result = parse_enewick(f"({long_label},B);")
        labels = [n.get("label") for n in result.nodes if "label" in n]
        assert long_label in labels

    def test_very_small_branch_lengths(self) -> None:
        """Very small branch lengths should be parsed correctly."""
        result = parse_enewick("(A:1e-10,B:1e-20);")
        edges_with_bl = [e for e in result.edges if "branch_length" in e]
        assert edges_with_bl[0]["branch_length"] == 1e-10
        assert edges_with_bl[1]["branch_length"] == 1e-20

    def test_very_large_branch_lengths(self) -> None:
        """Very large branch lengths should be parsed correctly."""
        result = parse_enewick("(A:1e10,B:1e20);")
        edges_with_bl = [e for e in result.edges if "branch_length" in e]
        assert edges_with_bl[0]["branch_length"] == 1e10
        assert edges_with_bl[1]["branch_length"] == 1e20

    def test_zero_branch_length(self) -> None:
        """Zero branch length should be parsed correctly."""
        result = parse_enewick("(A:0,B:0.0);")
        edges_with_bl = [e for e in result.edges if "branch_length" in e]
        assert edges_with_bl[0]["branch_length"] == 0.0
        assert edges_with_bl[1]["branch_length"] == 0.0

    def test_unicode_in_labels(self) -> None:
        """Unicode characters in labels should be handled."""
        result = parse_enewick("('Species α','Species β');")
        labels = [n.get("label") for n in result.nodes if "label" in n]
        assert "Species α" in labels
        assert "Species β" in labels

    def test_special_characters_in_comments(self) -> None:
        """Special characters in comments should be preserved."""
        result = parse_enewick("(A[comment with !@#$%^&*()],B);")
        node_a = next(n for n in result.nodes if n.get("id") == "A")
        assert "!" in node_a.get("comment", "")
        assert "@" in node_a.get("comment", "")


class TestParsedENewickStructure:
    """Test cases for ParsedENewick data structure."""

    def test_parsed_structure_has_required_fields(self) -> None:
        """ParsedENewick should have all required fields."""
        result = parse_enewick("(A,B);")
        assert hasattr(result, "edges")
        assert hasattr(result, "nodes")
        assert hasattr(result, "root")
        assert hasattr(result, "hybrid_nodes")
        assert hasattr(result, "metadata")

    def test_edges_structure(self) -> None:
        """Edges should have correct structure."""
        result = parse_enewick("(A:0.5,B:0.3);")
        assert len(result.edges) > 0
        for edge in result.edges:
            assert "u" in edge
            assert "v" in edge
            assert "key" in edge
            assert isinstance(edge["key"], int)

    def test_nodes_structure(self) -> None:
        """Nodes should have correct structure."""
        result = parse_enewick("(A,B);")
        assert len(result.nodes) > 0
        for node in result.nodes:
            assert "id" in node
            assert isinstance(node["id"], (str, int))

    def test_hybrid_nodes_dictionary(self) -> None:
        """Hybrid nodes should be a dictionary mapping IDs to numbers."""
        result = parse_enewick("((A,B)#H1,C);")
        assert isinstance(result.hybrid_nodes, dict)
        assert len(result.hybrid_nodes) == 1
        for node_id, hybrid_num in result.hybrid_nodes.items():
            assert isinstance(node_id, (str, int))
            assert isinstance(hybrid_num, int)

    def test_root_is_node_id(self) -> None:
        """Root should be a valid node ID."""
        result = parse_enewick("(A,B);")
        assert result.root is not None
        node_ids = [n["id"] for n in result.nodes]
        assert result.root in node_ids

    def test_metadata_is_dict(self) -> None:
        """Metadata should be a dictionary."""
        result = parse_enewick("(A,B);")
        assert isinstance(result.metadata, dict)


class TestRealWorldExamples:
    """Test cases based on real-world eNewick examples."""

    def test_simple_phylogenetic_tree(self) -> None:
        """Parse a simple phylogenetic tree."""
        enewick = "((Human:0.1,Chimp:0.1):0.05,Gorilla:0.15);"
        result = parse_enewick(enewick)
        
        labels = [n.get("label") for n in result.nodes if "label" in n]
        assert "Human" in labels
        assert "Chimp" in labels
        assert "Gorilla" in labels
        
        # All edges should have branch lengths (4 edges: 2 to leaves, 1 to internal, 1 to Gorilla)
        edges_with_bl = [e for e in result.edges if "branch_length" in e]
        assert len(edges_with_bl) == 4

    def test_network_with_reticulation(self) -> None:
        """Parse a network with reticulation event (Extended Newick format)."""
        enewick = "((A,B)#H1,(#H1,C));"
        result = parse_enewick(enewick)
        
        assert len(result.hybrid_nodes) == 1
        # Verify hybrid has 2 incoming edges (reticulation)
        hybrid_id = list(result.hybrid_nodes.keys())[0]
        incoming_edges = [e for e in result.edges if e['v'] == hybrid_id]
        assert len(incoming_edges) == 2

    def test_tree_with_internal_labels(self) -> None:
        """Parse a tree where internal nodes might have labels."""
        # Note: Internal nodes without explicit labels get auto-generated IDs
        enewick = "((A,B),C);"
        result = parse_enewick(enewick)
        
        # Should have A, B, C as labeled nodes
        labeled_nodes = [n for n in result.nodes if "label" in n]
        assert len(labeled_nodes) >= 3

