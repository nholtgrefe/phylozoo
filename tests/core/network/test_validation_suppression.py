"""
Tests for validation suppression on network classes.

This module tests that the validation_aware decorator works correctly
on DirectedPhyNetwork, MixedPhyNetwork, and SemiDirectedPhyNetwork.
"""

import warnings

import pytest

from phylozoo.core.network import DirectedPhyNetwork, SemiDirectedPhyNetwork
from phylozoo.core.network.sdnetwork.base import MixedPhyNetwork
from phylozoo.utils.validation import no_validation


class TestDirectedPhyNetworkValidationSuppression:
    """Test validation suppression for DirectedPhyNetwork."""

    def test_validate_runs_by_default(self) -> None:
        """Validation runs during initialization by default."""
        # This should raise ValueError because network is disconnected
        with pytest.raises(ValueError, match="not connected"):
            DirectedPhyNetwork(
                edges=[(1, 2), (3, 4)],
                nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"}), (4, {"label": "D"})],
            )

    def test_validate_suppressed_during_init(self) -> None:
        """Validation is suppressed during initialization in no_validation context."""
        # This should succeed because validation is suppressed
        with no_validation():
            net = DirectedPhyNetwork(
                edges=[(1, 2), (3, 4)],  # Disconnected, would normally fail
                nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"}), (4, {"label": "D"})],
            )
        # Network was created despite being invalid
        assert net.number_of_nodes() == 4

    def test_validate_suppressed_explicitly(self) -> None:
        """Validation can be suppressed explicitly for validate method."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})])
        # Modify to make invalid (add disconnected component)
        net._graph.add_edge(5, 6)
        net._graph.add_node(5)
        net._graph.add_node(6)
        
        # Should raise error normally
        with pytest.raises(ValueError, match="not connected"):
            net.validate()
        
        # Should not raise when suppressed
        with no_validation(only=["validate"]):
            net.validate()  # Should return None silently

    def test_helper_methods_can_be_suppressed(self) -> None:
        """Helper validation methods (_validate_*) can be suppressed independently."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})])
        # Add disconnected component - this would normally be caught by _validate_structural_constraints
        net._graph.add_edge(5, 6)
        net._graph.add_node(5)
        net._graph.add_node(6)
        # Add labels for the disconnected nodes so they're valid leaves
        net._node_to_label[5] = "E"
        net._node_to_label[6] = "F"
        net._label_to_node["E"] = 5
        net._label_to_node["F"] = 6
        
        # When not suppressed, should fail on connectivity
        with pytest.raises(ValueError, match="not connected"):
            net.validate()
        
        # When structural constraints are suppressed, connectivity check is skipped
        # The network structure is invalid but validation passes because structural checks are suppressed
        with no_validation(only=["_validate_structural_constraints"]):
            net.validate()  # Should succeed because connectivity check is suppressed


class TestMixedPhyNetworkValidationSuppression:
    """Test validation suppression for MixedPhyNetwork."""

    def test_validate_runs_by_default(self) -> None:
        """Validation runs during initialization by default."""
        # This should raise ValueError because network is disconnected
        with pytest.raises(ValueError, match="not connected"):
            MixedPhyNetwork(
                directed_edges=[],
                undirected_edges=[(1, 2), (3, 4)],
                nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"}), (4, {"label": "D"})],
            )

    def test_validate_suppressed_during_init(self) -> None:
        """Validation is suppressed during initialization in no_validation context."""
        # This should succeed because validation is suppressed
        with no_validation():
            net = MixedPhyNetwork(
                directed_edges=[],
                undirected_edges=[(1, 2), (3, 4)],  # Disconnected, would normally fail
                nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"}), (4, {"label": "D"})],
            )
        # Network was created despite being invalid
        assert net.number_of_nodes() == 4

    @pytest.mark.filterwarnings("ignore:Validity is not fully checked.*MixedPhyNetworks:UserWarning")
    def test_validate_suppressed_explicitly(self) -> None:
        """Validation can be suppressed explicitly for validate method."""
        # Create valid network first (internal node needs degree >= 3)
        net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (4, {"label": "C"})],
        )
        # Modify to make invalid (add self-loop which violates structural constraints)
        net._graph._undirected.add_edge(1, 1)
        
        # Should raise error normally (before reaching mixed-network constraint check)
        with pytest.raises(ValueError, match="Self-loops"):
            net.validate()
        
        # Should not raise when suppressed
        with no_validation(only=["validate"]):
            net.validate()  # Should return None silently


class TestSemiDirectedPhyNetworkValidationSuppression:
    """Test validation suppression for SemiDirectedPhyNetwork."""

    def test_validate_runs_by_default(self) -> None:
        """Validation runs during initialization by default."""
        # This should raise ValueError because network is disconnected
        with pytest.raises(ValueError, match="not connected"):
            SemiDirectedPhyNetwork(
                directed_edges=[],
                undirected_edges=[(1, 2), (3, 4)],
                nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"}), (4, {"label": "D"})],
            )

    def test_validate_suppressed_during_init(self) -> None:
        """Validation is suppressed during initialization in no_validation context."""
        # This should succeed because validation is suppressed
        with no_validation():
            net = SemiDirectedPhyNetwork(
                directed_edges=[],
                undirected_edges=[(1, 2), (3, 4)],  # Disconnected, would normally fail
                nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"}), (4, {"label": "D"})],
            )
        # Network was created despite being invalid
        assert net.number_of_nodes() == 4

    def test_validate_suppressed_explicitly(self) -> None:
        """Validation can be suppressed explicitly for validate method."""
        # Create valid network first (internal node needs degree >= 3)
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (4, {"label": "C"})],
        )
        # Modify to make invalid (add self-loop which violates structural constraints)
        net._graph._undirected.add_edge(1, 1)
        
        # Should raise error normally
        with pytest.raises(ValueError, match="Self-loops"):
            net.validate()
        
        # Should not raise when suppressed
        with no_validation(only=["validate"]):
            net.validate()  # Should return None silently


class TestNestedValidationSuppression:
    """Test nested validation suppression contexts."""

    def test_nested_contexts_work(self) -> None:
        """Nested no_validation contexts work correctly."""
        # Outer context suppresses validation
        with no_validation():
            net1 = DirectedPhyNetwork(
                edges=[(1, 2), (3, 4)],  # Disconnected
                nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"}), (4, {"label": "D"})],
            )
            # Inner context also suppresses
            with no_validation():
                net2 = DirectedPhyNetwork(
                    edges=[(5, 6), (7, 8)],  # Also disconnected
                    nodes=[(5, {"label": "E"}), (6, {"label": "F"}), (7, {"label": "G"}), (8, {"label": "H"})],
                )
            # Still in outer context
            net3 = DirectedPhyNetwork(
                edges=[(9, 10)],  # Single edge, valid
                nodes=[(9, {"label": "I"}), (10, {"label": "J"})],
            )
        # All networks created successfully
        assert net1.number_of_nodes() == 4
        assert net2.number_of_nodes() == 4
        assert net3.number_of_nodes() == 2

    def test_validation_restored_after_context(self) -> None:
        """Validation is restored after exiting context."""
        # Create invalid network with suppression
        with no_validation():
            net = DirectedPhyNetwork(
                edges=[(1, 2), (3, 4)],  # Disconnected
                nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"}), (4, {"label": "D"})],
            )
        
        # After context, validation should work normally
        with pytest.raises(ValueError, match="not connected"):
            net.validate()


class TestWildcardPatterns:
    """Test that wildcard patterns work for helper methods."""

    def test_wildcard_suppresses_helpers(self) -> None:
        """Wildcard pattern _validate_struct* suppresses structural constraint methods."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})])
        # Add disconnected component - this would normally be caught by _validate_structural_constraints
        net._graph.add_edge(5, 6)
        net._graph.add_node(5)
        net._graph.add_node(6)
        # Add labels for the disconnected nodes so they're valid leaves
        net._node_to_label[5] = "E"
        net._node_to_label[6] = "F"
        net._label_to_node["E"] = 5
        net._label_to_node["F"] = 6
        
        # When not suppressed, should fail on connectivity
        with pytest.raises(ValueError, match="not connected"):
            net.validate()
        
        # When suppressed using wildcard pattern _validate_struct*, connectivity check is skipped
        # The disconnected nodes have valid degrees (5 has in-degree 1, out-degree 0; 6 has in-degree 1, out-degree 0)
        # So validate() should succeed when structural constraints are suppressed
        with no_validation(only=["_validate_struct*"]):
            net.validate()  # Should succeed because connectivity check is suppressed via wildcard pattern


class TestDefaultSuppression:
    """Test that default suppression works correctly."""

    def test_default_suppresses_validate(self) -> None:
        """Default suppression pattern suppresses validate method."""
        # Using no_validation() without arguments should use default (validate)
        with no_validation():
            net = DirectedPhyNetwork(
                edges=[(1, 2), (3, 4)],  # Disconnected
                nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"}), (4, {"label": "D"})],
            )
        # Should succeed
        assert net.number_of_nodes() == 4
        
        # But validate should still work when called explicitly outside context
        with pytest.raises(ValueError, match="not connected"):
            net.validate()

