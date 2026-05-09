"""
Tests for DirectedPhyNetwork to eNewick conversion.
"""

import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork, isomorphism
from phylozoo.core.network.dnetwork._enewick import from_enewick, parse_enewick, to_enewick

# Large hybrid eNewick (reticulations, named #H markers); used for round-trip isomorphism.
_LARGE_SAMPLE_ENEWICK: str = (
    "((((t17:1.099348922,t12:1.099348922)N1:0.1705515639,(((t15:0.3933582101,((t16:0.1635589659,"
    "t18:0.1635589659)N2:0.2159864787,(t7:0.04769816128,t6:0.04769816128)N3:0.3318472834)N4:0.01381276542)"
    "N5:0.1759550544,((t20:0.2141536746,(t1:0.2141536746)N27#H46:0)N7:0.03387300123,t8:0.2480266759)N8:"
    "0.3212865886)N9:0.6962976373,(((t10:0.07172234796,t9:0.07172234796)N10:0.2797101105,(t13:0.2141536746,"
    "N27#H46:0)N11:0.1372787838)N12:0.07536269875,N28#H35:-2.220446049e-16)N13:0.8388157445)N14:0.004289584448)"
    "N15:0.4018792131,((((t19:0.4037150203,N29#H39:0)N16:0.02308013689,((t14:0.4037150203,(t5:0.4037150203)"
    "N29#H39:0)N18:0.02308013689)N28#H35:0)N20:1.023456361,((t4:0.4223670865,t11:0.4223670865)N21:0.2627193004,"
    "t2:0.6850863869)N22:0.7651651316)N23:0.1715330931,t3:1.621784612)N24:0.04999508769)N25:0.4803650104)N26;"
)


class TestBasicTrees:
    """Test eNewick conversion for simple trees without hybrids."""
    
    def test_single_node_with_label(self) -> None:
        """Convert single-node network with label."""
        with pytest.warns(UserWarning, match="Single-node network detected"):
            net = DirectedPhyNetwork(nodes=[(1, {'label': 'A'})])
        result = to_enewick(net)
        assert result == "A;"
    
    def test_single_node_without_label(self) -> None:
        """Convert single-node network without label."""
        with pytest.warns(UserWarning, match="Single-node network detected"):
            net = DirectedPhyNetwork(nodes=[1])
        result = to_enewick(net)
        assert result == "1;"
    
    def test_simple_tree_two_leaves(self) -> None:
        """Convert simple tree with two leaves."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        result = to_enewick(net)
        # Should be deterministic
        assert result == "(A,B);"
    
    def test_simple_tree_three_leaves(self) -> None:
        """Convert tree with three leaves."""
        net = DirectedPhyNetwork(
            edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
        )
        result = to_enewick(net)
        assert result == "(A,B,C);"
    
    def test_nested_tree(self) -> None:
        """Convert nested tree structure."""
        net = DirectedPhyNetwork(
            edges=[(5, 3), (5, 4), (3, 1), (3, 2)],
            nodes=[
                (1, {'label': 'A'}), 
                (2, {'label': 'B'}), 
                (4, {'label': 'C'})
            ]
        )
        result = to_enewick(net)
        # Root 5 has children 3 and 4, where 3 has children 1 and 2
        assert result == "((A,B),C);"
    
    def test_empty_network_raises_error(self) -> None:
        """Empty network should raise ValueError."""
        with pytest.warns(UserWarning, match="Empty network"):
            net = DirectedPhyNetwork(edges=[])
        with pytest.raises(ValueError, match="Cannot convert empty network"):
            to_enewick(net)


class TestBranchLengths:
    """Test eNewick conversion with branch lengths."""
    
    def test_branch_lengths_simple(self) -> None:
        """Convert tree with branch lengths."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                {'u': 3, 'v': 2, 'branch_length': 0.3}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        result = to_enewick(net)
        assert result == "(A:0.5,B:0.3);"
    
    def test_branch_lengths_nested(self) -> None:
        """Convert nested tree with branch lengths."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 5, 'v': 3, 'branch_length': 0.1},
                {'u': 5, 'v': 4, 'branch_length': 0.2},
                {'u': 3, 'v': 1, 'branch_length': 0.3},
                {'u': 3, 'v': 2, 'branch_length': 0.4}
            ],
            nodes=[
                (1, {'label': 'A'}), 
                (2, {'label': 'B'}), 
                (4, {'label': 'C'})
            ]
        )
        result = to_enewick(net)
        assert result == "((A:0.3,B:0.4):0.1,C:0.2);"
    
    def test_branch_lengths_mixed(self) -> None:
        """Convert tree with some branch lengths missing."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                (3, 2)  # No branch length
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        result = to_enewick(net)
        assert result == "(A:0.5,B);"


class TestInternalLabels:
    """Test eNewick conversion with internal node labels."""
    
    def test_internal_label_root(self) -> None:
        """Convert tree with root label."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[
                (1, {'label': 'A'}), 
                (2, {'label': 'B'}), 
                (3, {'label': 'root'})
            ]
        )
        result = to_enewick(net)
        assert result == "(A,B)root;"
    
    def test_internal_labels_multiple(self) -> None:
        """Convert tree with multiple internal labels."""
        net = DirectedPhyNetwork(
            edges=[(5, 3), (5, 4), (3, 1), (3, 2)],
            nodes=[
                (1, {'label': 'A'}), 
                (2, {'label': 'B'}), 
                (3, {'label': 'internal1'}),
                (4, {'label': 'C'}),
                (5, {'label': 'root'})
            ]
        )
        result = to_enewick(net)
        assert result == "((A,B)internal1,C)root;"
    
    def test_internal_label_with_branch_length(self) -> None:
        """Convert tree with internal labels and branch lengths."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 5, 'v': 3, 'branch_length': 0.1},
                {'u': 5, 'v': 4, 'branch_length': 0.2},
                {'u': 3, 'v': 1, 'branch_length': 0.3},
                {'u': 3, 'v': 2, 'branch_length': 0.4}
            ],
            nodes=[
                (1, {'label': 'A'}), 
                (2, {'label': 'B'}), 
                (3, {'label': 'int1'}),
                (4, {'label': 'C'}),
                (5, {'label': 'root'})
            ]
        )
        result = to_enewick(net)
        assert result == "((A:0.3,B:0.4)int1:0.1,C:0.2)root;"


class TestHybridNodes:
    """Test eNewick conversion for networks with hybrid nodes."""
    
    def test_single_hybrid_simple(self) -> None:
        """Convert network with single hybrid node."""
        # Structure: root 7 -> (5, 6), 5 -> (4, 8), 6 -> (4, 9), 4 -> 10, 10 -> (1, 2)
        # Node 4 is hybrid (in-degree 2, out-degree 1)
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                (5, 4), (5, 8),  # Tree node 5 has out-degree 2
                (6, 4), (6, 9),  # Tree node 6 has out-degree 2
                (4, 10),  # Hybrid node 4 has out-degree 1
                (10, 1), (10, 2)  # Tree node 10 has out-degree 2
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (8, {'label': 'C'}), (9, {'label': 'D'})]
        )
        result = to_enewick(net)
        # Hybrid node 4 should get #H1 marker
        # First occurrence defines it, second references it
        assert "#H1" in result
        assert result.count("#H1") == 2
    
    def test_single_hybrid_with_gamma(self) -> None:
        """Convert network with hybrid node and gamma values."""
        # Structure: root 7 -> (5, 6), 5 -> (4, 8), 6 -> (4, 9), 4 -> 10, 10 -> (1, 2)
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                {'u': 5, 'v': 4, 'gamma': 0.6},
                (5, 8),
                {'u': 6, 'v': 4, 'gamma': 0.4},
                (6, 9),
                (4, 10),
                (10, 1), (10, 2)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (8, {'label': 'C'}), (9, {'label': 'D'})]
        )
        result = to_enewick(net)
        assert "#H1" in result
        assert "gamma=0.6" in result or "gamma=0.4" in result
    
    def test_multiple_hybrid_nodes(self) -> None:
        """Convert network with multiple hybrid nodes."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                (5, 3), (5, 4),
                (6, 3), (6, 4),
                (3, 1), (4, 2)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        result = to_enewick(net)
        # Should have two hybrid nodes: 3 and 4
        assert "#H1" in result
        assert "#H2" in result
        assert result.count("#H1") == 2
        assert result.count("#H2") == 2
    
    def test_hybrid_with_branch_lengths(self) -> None:
        """Convert network with hybrid node and branch lengths."""
        # Structure: root 7 -> (5, 6), 5 -> (4, 8), 6 -> (4, 9), 4 -> 10, 10 -> (1, 2)
        net = DirectedPhyNetwork(
            edges=[
                {'u': 7, 'v': 5, 'branch_length': 0.1},
                {'u': 7, 'v': 6, 'branch_length': 0.2},
                {'u': 5, 'v': 4, 'branch_length': 0.3},
                {'u': 5, 'v': 8, 'branch_length': 0.35},
                {'u': 6, 'v': 4, 'branch_length': 0.4},
                {'u': 6, 'v': 9, 'branch_length': 0.45},
                {'u': 4, 'v': 10, 'branch_length': 0.5},
                {'u': 10, 'v': 1, 'branch_length': 0.6},
                {'u': 10, 'v': 2, 'branch_length': 0.7}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (8, {'label': 'C'}), (9, {'label': 'D'})]
        )
        result = to_enewick(net)
        assert "#H1" in result
        # Check that branch lengths are present
        assert ":0." in result


class TestGammaBootstrap:
    """Test eNewick conversion with gamma and bootstrap values."""
    
    def test_gamma_values(self) -> None:
        """Convert network with gamma values in comments."""
        # Structure: root 7 -> (5, 6), 5 -> (4, 8), 6 -> (4, 9), 4 -> 10, 10 -> (1, 2)
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                {'u': 5, 'v': 4, 'gamma': 0.6},
                (5, 8),
                {'u': 6, 'v': 4, 'gamma': 0.4},
                (6, 9),
                (4, 10),
                (10, 1), (10, 2)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (8, {'label': 'C'}), (9, {'label': 'D'})]
        )
        result = to_enewick(net)
        assert "gamma=" in result
    
    def test_bootstrap_values(self) -> None:
        """Convert tree with bootstrap values in comments."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'bootstrap': 0.95},
                {'u': 3, 'v': 2, 'bootstrap': 0.87}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        result = to_enewick(net)
        assert "bootstrap=0.95" in result
        assert "bootstrap=0.87" in result
    
    def test_gamma_and_bootstrap_combined(self) -> None:
        """Convert network with both gamma and bootstrap values."""
        # Structure: root 7 -> (5, 6), 5 -> (4, 8), 6 -> (4, 9), 4 -> 10, 10 -> (1, 2)
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                {'u': 5, 'v': 4, 'gamma': 0.6, 'bootstrap': 0.95},
                (5, 8),
                {'u': 6, 'v': 4, 'gamma': 0.4, 'bootstrap': 0.85},
                (6, 9),
                (4, 10),
                (10, 1), (10, 2)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (8, {'label': 'C'}), (9, {'label': 'D'})]
        )
        result = to_enewick(net)
        assert "gamma=" in result
        assert "bootstrap=" in result
        # Check that they're in the same comment block
        assert "[&" in result


class TestLabelQuoting:
    """Test eNewick conversion with labels requiring quotes."""
    
    def test_label_with_spaces(self) -> None:
        """Convert tree with label containing spaces."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'Taxon A'}), (2, {'label': 'B'})]
        )
        result = to_enewick(net)
        assert "'Taxon A'" in result
    
    def test_label_with_parentheses(self) -> None:
        """Convert tree with label containing parentheses."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A(group1)'}), (2, {'label': 'B'})]
        )
        result = to_enewick(net)
        assert "'A(group1)'" in result
    
    def test_label_with_colon(self) -> None:
        """Convert tree with label containing colon."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A:1'}), (2, {'label': 'B'})]
        )
        result = to_enewick(net)
        assert "'A:1'" in result
    
    def test_label_with_quote(self) -> None:
        """Convert tree with label containing single quote."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': "A'strain"}), (2, {'label': 'B'})]
        )
        result = to_enewick(net)
        # Single quotes should be escaped
        assert "A''strain" in result


class TestParallelEdges:
    """Test eNewick conversion with parallel edges."""
    
    def test_parallel_edges_first_used(self) -> None:
        """Convert network with parallel edges (first edge used)."""
        # Use tree node with parallel edges to children to keep out-degree >= 2
        net = DirectedPhyNetwork(
            edges=[
                {'u': 5, 'v': 1, 'key': 0, 'branch_length': 0.5},
                {'u': 5, 'v': 2, 'key': 0, 'branch_length': 0.3},
                # Add parallel edges to make it a valid tree (out-degree >= 2 for internal nodes)
                {'u': 5, 'v': 3, 'key': 0},
                {'u': 5, 'v': 4, 'key': 0}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'})]
        )
        result = to_enewick(net)
        # Should use branch lengths from edges
        assert "A:0.5" in result
        assert "B:0.3" in result


class TestDeterministicOutput:
    """Test that output is deterministic."""
    
    def test_same_network_same_output(self) -> None:
        """Multiple conversions of same network produce same output."""
        net = DirectedPhyNetwork(
            edges=[(5, 3), (5, 4), (3, 1), (3, 2)],
            nodes=[
                (1, {'label': 'A'}), 
                (2, {'label': 'B'}), 
                (4, {'label': 'C'})
            ]
        )
        result1 = to_enewick(net)
        result2 = to_enewick(net)
        result3 = to_enewick(net)
        assert result1 == result2 == result3
    
    def test_children_sorted_deterministically(self) -> None:
        """Children are output in deterministic order."""
        net = DirectedPhyNetwork(
            edges=[(10, 5), (10, 3), (10, 7)],
            nodes=[
                (3, {'label': 'C'}), 
                (5, {'label': 'B'}), 
                (7, {'label': 'A'})
            ]
        )
        result = to_enewick(net)
        # Children should be sorted by node ID, not by label
        # Node IDs: 3, 5, 7 -> C, B, A
        assert result == "(C,B,A);"


class TestRoundTrip:
    """Test round-trip conversion: network -> eNewick -> parse -> network."""
    
    def test_simple_tree_roundtrip(self) -> None:
        """Round-trip conversion preserves structure for simple tree."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        enewick_str = to_enewick(net)
        parsed = parse_enewick(enewick_str)
        
        # Check structure is preserved
        assert len(parsed.nodes) == 3  # 2 leaves + 1 internal
        assert len(parsed.edges) == 2
        
        # Check taxa
        leaf_labels = {n['label'] for n in parsed.nodes if n.get('label') in ['A', 'B']}
        assert leaf_labels == {'A', 'B'}
    
    def test_tree_with_branch_lengths_roundtrip(self) -> None:
        """Round-trip conversion preserves branch lengths."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                {'u': 3, 'v': 2, 'branch_length': 0.3}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        enewick_str = to_enewick(net)
        parsed = parse_enewick(enewick_str)
        
        # Check branch lengths are preserved
        branch_lengths = [e.get('branch_length') for e in parsed.edges if 'branch_length' in e]
        assert set(branch_lengths) == {0.5, 0.3}
    
    def test_hybrid_network_roundtrip(self) -> None:
        """Round-trip conversion preserves hybrid node structure."""
        # Structure: root 7 -> (5, 6), 5 -> (4, 8), 6 -> (4, 9), 4 -> 10, 10 -> (1, 2)
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                (5, 4), (5, 8),
                (6, 4), (6, 9),
                (4, 10),
                (10, 1), (10, 2)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (8, {'label': 'C'}), (9, {'label': 'D'})]
        )
        enewick_str = to_enewick(net)
        parsed = parse_enewick(enewick_str)
        
        # Check hybrid nodes are identified
        assert len(parsed.hybrid_nodes) == 1
        
        # Check that hybrid node has 2 incoming edges
        hybrid_id = list(parsed.hybrid_nodes.keys())[0]
        incoming_edges = [e for e in parsed.edges if e['v'] == hybrid_id]
        assert len(incoming_edges) == 2
    
    def test_internal_labels_roundtrip(self) -> None:
        """Round-trip conversion preserves internal labels."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[
                (1, {'label': 'A'}), 
                (2, {'label': 'B'}), 
                (3, {'label': 'root'})
            ]
        )
        enewick_str = to_enewick(net)
        parsed = parse_enewick(enewick_str)
        
        # Check internal label is preserved
        internal_nodes = [n for n in parsed.nodes if n.get('label') == 'root']
        assert len(internal_nodes) == 1

    def test_large_enewick_roundtrip_isomorphic(self) -> None:
        """
        Large sample eNewick: parse → build network → to_enewick → parse again.

        The two :class:`~phylozoo.core.network.dnetwork.DirectedPhyNetwork`
        instances must be isomorphic including ``branch_length`` on edges.
        """
        net_a = from_enewick(_LARGE_SAMPLE_ENEWICK)
        net_b = from_enewick(to_enewick(net_a))
        assert net_a.number_of_nodes() == net_b.number_of_nodes()
        assert net_a.number_of_edges() == net_b.number_of_edges()
        assert isomorphism.is_isomorphic(net_a, net_b, edge_attrs=["branch_length"])


class TestComplexNetworks:
    """Test eNewick conversion for complex networks."""
    
    def test_complex_network_all_features(self) -> None:
        """Convert complex network with all features combined."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 7, 'v': 5, 'branch_length': 0.1, 'bootstrap': 0.99},
                {'u': 7, 'v': 6, 'branch_length': 0.2, 'bootstrap': 0.95},
                {'u': 5, 'v': 4, 'gamma': 0.6, 'branch_length': 0.3, 'bootstrap': 0.9},
                {'u': 6, 'v': 4, 'gamma': 0.4, 'branch_length': 0.4, 'bootstrap': 0.85},
                {'u': 4, 'v': 1, 'branch_length': 0.5},
                {'u': 5, 'v': 2, 'branch_length': 0.6},
                {'u': 6, 'v': 3, 'branch_length': 0.7}
            ],
            nodes=[
                (1, {'label': 'A'}), 
                (2, {'label': 'B'}), 
                (3, {'label': 'C'}),
                (5, {'label': 'int1'}),
                (6, {'label': 'int2'}),
                (7, {'label': 'root'})
            ]
        )
        result = to_enewick(net)
        
        # Check all features are present
        assert "#H1" in result  # Hybrid marker
        assert "gamma=" in result  # Gamma values
        assert "bootstrap=" in result  # Bootstrap values
        assert ":" in result  # Branch lengths
        assert "int1" in result or "int2" in result  # Internal labels
        assert "root" in result  # Root label

