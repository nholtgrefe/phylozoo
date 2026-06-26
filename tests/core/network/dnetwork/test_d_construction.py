"""Tests for the directed-generator construction step (``gambette_step``)."""

import pytest

from phylozoo.core.network.dnetwork.generator.construction import (
    all_level_k_generators,
    gambette_step,
)


class TestGambetteStep:
    @pytest.mark.parametrize("k", [1, 2, 3])
    def test_step_from_complete_prev_level_reproduces_all_level_k(self, k: int) -> None:
        """One step on the complete level-(k-1) set yields the complete level-k set.

        ``k == 1`` exercises the level-0 input (the R1 rule bootstraps the single
        node into a level-1 generator).
        """
        stepped = gambette_step(all_level_k_generators(k - 1))
        assert len(stepped) == len(all_level_k_generators(k))

    def test_mixed_levels_are_allowed(self) -> None:
        """Inputs of different levels are each advanced independently."""
        mixed = all_level_k_generators(1) | all_level_k_generators(2)
        result = gambette_step(mixed)
        # level-1 -> level-2 and level-2 -> level-3 (disjoint).
        assert len(result) == len(all_level_k_generators(2)) + len(all_level_k_generators(3))
