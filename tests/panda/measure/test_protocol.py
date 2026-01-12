"""
Tests for the DiversityMeasure protocol.

This test suite verifies that the protocol is correctly defined and that
implementations satisfy the protocol interface.
"""

from typing import Any, Set

import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.panda.measure.protocol import DiversityMeasure
from phylozoo.panda.measure.all_paths import AllPathsDiversity


class TestDiversityMeasureProtocol:
    """Test cases for the DiversityMeasure protocol."""

    def test_all_paths_satisfies_protocol(self) -> None:
        """Test that AllPathsDiversity satisfies the DiversityMeasure protocol."""
        measure = AllPathsDiversity()
        
        # Check that it has the required methods
        assert hasattr(measure, 'compute_diversity')
        assert hasattr(measure, 'solve_maximization')
        
        # Check method signatures
        import inspect
        compute_sig = inspect.signature(measure.compute_diversity)
        assert 'network' in compute_sig.parameters
        assert 'taxa' in compute_sig.parameters
        
        solve_sig = inspect.signature(measure.solve_maximization)
        assert 'network' in solve_sig.parameters
        assert 'k' in solve_sig.parameters

    def test_protocol_is_structural(self) -> None:
        """
        Test that Protocol works via structural typing.
        
        Protocols don't require explicit inheritance - any class that
        implements the required methods satisfies the protocol.
        """
        class CustomMeasure:
            """A custom measure that implements the protocol without explicit inheritance."""
            
            def compute_diversity(
                self,
                network: DirectedPhyNetwork,
                taxa: Set[str],
                **kwargs: Any,
            ) -> float:
                """Compute diversity."""
                return 1.0
            
            def solve_maximization(
                self,
                network: DirectedPhyNetwork,
                k: int,
                **kwargs: Any,
            ) -> tuple[float, Set[str]]:
                """Solve maximization."""
                return (1.0, set())
        
        # This should work - the class satisfies the protocol structurally
        measure: DiversityMeasure = CustomMeasure()
        assert measure.compute_diversity(None, set()) == 1.0  # type: ignore

