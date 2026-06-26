"""Tests for the semi-directed-generator construction step (``semidirect_generators``)."""

import pytest

from phylozoo.core.network.dnetwork.generator.construction import (
    all_level_k_generators as d_all_level_k_generators,
)
from phylozoo.core.network.sdnetwork.generator.construction import (
    all_level_k_generators,
    semidirect_generators,
)


class TestSemidirectGenerators:
    @pytest.mark.parametrize("k", [1, 2, 3])
    def test_reproduces_all_level_k(self, k: int) -> None:
        """Semi-directing the complete level-k d-generators yields the complete sd set."""
        result = semidirect_generators(d_all_level_k_generators(k))
        assert len(result) == len(all_level_k_generators(k))

    def test_mixed_levels_are_allowed(self) -> None:
        mixed = d_all_level_k_generators(2) | d_all_level_k_generators(3)
        result = semidirect_generators(mixed)
        assert len(result) == len(all_level_k_generators(2)) + len(all_level_k_generators(3))
