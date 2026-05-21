"""
Round-trip tests connecting DistanceMatrix ↔ WeightedSplitSystem.

These tests verify that the conversions::

    WeightedSplitSystem  →  distances_from_splitsystem  →  DistanceMatrix
    DistanceMatrix       →  split_decomposition          →  WeightedSplitSystem

are inverses of each other on the appropriate subsets, and that the
SplitSystem → tree → induced_splits round-trip is also lossless.
"""

import numpy as np
import pytest

from phylozoo.core.distance import DistanceMatrix, split_decomposition
from phylozoo.core.split import Split
from phylozoo.core.split.algorithms import distances_from_splitsystem, tree_from_splitsystem
from phylozoo.core.split.classifications import is_tree_compatible
from phylozoo.core.split.splitsystem import SplitSystem
from phylozoo.core.split.weighted_splitsystem import WeightedSplitSystem

# ---------------------------------------------------------------------------
# Helpers / shared data
# ---------------------------------------------------------------------------

TABLE1_MATRIX = np.array(
    [
        [0, 4, 5, 7, 13, 8, 6],
        [4, 0, 1, 3, 9, 12, 10],
        [5, 1, 0, 2, 8, 13, 11],
        [7, 3, 2, 0, 6, 11, 13],
        [13, 9, 8, 6, 0, 5, 7],
        [8, 12, 13, 11, 5, 0, 2],
        [6, 10, 11, 13, 7, 2, 0],
    ],
    dtype=np.float64,
)
TABLE1_LABELS = list("ABCDEFG")


def _cherries_system() -> WeightedSplitSystem:
    """Caterpillar tree on {A,B,C,D}: ((A,B),(C,D)) with all edge lengths positive.

    Splits and weights:
      {A}|{B,C,D} = 1,  {B}|{A,C,D} = 1
      {C}|{A,B,D} = 1,  {D}|{A,B,C} = 1
      {A,B}|{C,D} = 2   (internal split)
    """
    return WeightedSplitSystem(
        {
            Split({"A"}, {"B", "C", "D"}): 1.0,
            Split({"B"}, {"A", "C", "D"}): 1.0,
            Split({"C"}, {"A", "B", "D"}): 1.0,
            Split({"D"}, {"A", "B", "C"}): 1.0,
            Split({"A", "B"}, {"C", "D"}): 2.0,
        }
    )


# ---------------------------------------------------------------------------
# WeightedSplitSystem → DistanceMatrix → split_decomposition → same system
# ---------------------------------------------------------------------------


class TestSplitSystemToDistancesAndBack:
    """Split system → distances_from_splitsystem → split_decomposition recovers the system."""

    def test_cherries_tree_splits_recovered(self) -> None:
        """All five splits of a 4-taxon caterpillar tree are recovered."""
        system = _cherries_system()
        dm = distances_from_splitsystem(system)
        recovered, residual = split_decomposition(dm)

        assert recovered.splits == system.splits

    def test_cherries_tree_weights_recovered(self) -> None:
        """Isolation indices match the original edge weights after round-trip."""
        system = _cherries_system()
        dm = distances_from_splitsystem(system)
        recovered, _ = split_decomposition(dm)

        for sp in system.splits:
            assert recovered.get_weight(sp) == pytest.approx(system.get_weight(sp))

    def test_cherries_tree_residual_zero(self) -> None:
        """A tree split system has zero residual after round-trip."""
        system = _cherries_system()
        dm = distances_from_splitsystem(system)
        _, residual = split_decomposition(dm)

        assert np.allclose(residual.np_array, 0.0, atol=1e-9)

    def test_path4_system_recovered(self) -> None:
        """Path-4 d-splits {1}|rest, {1,2}|{3,4}, rest|{4} are recovered.

        Nodes 2 and 3 are internal in the caterpillar, so their trivial splits
        have zero pendant-edge length and do NOT appear as d-splits.
        """
        path4_splits = WeightedSplitSystem(
            {
                Split({1}, {2, 3, 4}): 1.0,
                Split({1, 2}, {3, 4}): 1.0,
                Split({1, 2, 3}, {4}): 1.0,
            }
        )
        dm = distances_from_splitsystem(path4_splits)
        recovered, _ = split_decomposition(dm)

        assert recovered.splits == path4_splits.splits
        for sp in path4_splits.splits:
            assert recovered.get_weight(sp) == pytest.approx(path4_splits.get_weight(sp))

    def test_single_split_two_taxa(self) -> None:
        """For n=2 the unique split is recovered with the correct weight."""
        system = WeightedSplitSystem({Split({"X"}, {"Y"}): 7.0})
        dm = distances_from_splitsystem(system)
        recovered, _ = split_decomposition(dm)

        sp = Split({"X"}, {"Y"})
        assert sp in recovered.splits
        assert recovered.get_weight(sp) == pytest.approx(7.0)


# ---------------------------------------------------------------------------
# DistanceMatrix → split_decomposition → distances_from_splitsystem → same d
# ---------------------------------------------------------------------------


class TestDistancesToSplitsAndBack:
    """split_decomposition → distances_from_splitsystem reconstructs the input."""

    def test_table1_totally_decomposable(self) -> None:
        """Table 1 (Bandelt & Dress 1992) has zero residual; d^1 == d."""
        dm = DistanceMatrix(TABLE1_MATRIX, labels=TABLE1_LABELS)
        system, residual = split_decomposition(dm)

        assert np.allclose(residual.np_array, 0.0, atol=1e-9), "residual should be zero"

        d1_dm = distances_from_splitsystem(system)
        d1 = np.array(
            [
                [d1_dm.get_distance(TABLE1_LABELS[r], TABLE1_LABELS[c]) for c in range(7)]
                for r in range(7)
            ]
        )
        assert np.allclose(d1, TABLE1_MATRIX, atol=1e-9)

    def test_path4_tree_metric(self) -> None:
        """Path-tree metric on {1,2,3,4} is reconstructed exactly via d-splits."""
        matrix = np.array(
            [[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]], dtype=np.float64
        )
        dm = DistanceMatrix(matrix, labels=[1, 2, 3, 4])
        system, residual = split_decomposition(dm)

        d1_dm = distances_from_splitsystem(system)
        d1 = np.array(
            [[d1_dm.get_distance(r + 1, c + 1) for c in range(4)] for r in range(4)],
            dtype=np.float64,
        )
        reconstructed = d1 + residual.np_array
        assert np.allclose(reconstructed, matrix, atol=1e-9)

    def test_distances_plus_residual_equals_input(self) -> None:
        """d^1 + d^0 == d for every distance matrix (decomposition identity)."""
        rng = np.random.default_rng(0)
        n = 5
        raw = rng.uniform(1, 10, (n, n))
        raw = (raw + raw.T) / 2
        np.fill_diagonal(raw, 0)
        dm = DistanceMatrix(raw, labels=list("ABCDE"))

        system, residual = split_decomposition(dm)

        labels = list(dm.labels)
        if len(system.splits) > 0:
            d1_dm = distances_from_splitsystem(system)
            d1 = np.array(
                [[d1_dm.get_distance(labels[r], labels[c]) for c in range(n)] for r in range(n)]
            )
        else:
            d1 = np.zeros((n, n))

        assert np.allclose(d1 + residual.np_array, raw, atol=1e-9)

    def test_cherries_tree_distances_exact(self) -> None:
        """Distances computed from a tree system decompose back to the same distances."""
        system = _cherries_system()
        dm = distances_from_splitsystem(system)
        recovered_system, residual = split_decomposition(dm)

        labels = list(dm.labels)
        n = len(labels)
        d1_dm = distances_from_splitsystem(recovered_system)
        d1 = np.array(
            [[d1_dm.get_distance(labels[r], labels[c]) for c in range(n)] for r in range(n)]
        )
        assert np.allclose(d1 + residual.np_array, dm.np_array, atol=1e-9)


# ---------------------------------------------------------------------------
# SplitSystem → tree_from_splitsystem → induced_splits → same splits
# ---------------------------------------------------------------------------


class TestSplitSystemToTreeAndBack:
    """Compatible SplitSystem → tree_from_splitsystem → induced_splits is lossless."""

    def _make_compatible_system(self) -> SplitSystem:
        """4-taxon compatible system: 4 trivials + 1 internal split."""
        return SplitSystem(
            [
                Split({"A"}, {"B", "C", "D"}),
                Split({"B"}, {"A", "C", "D"}),
                Split({"C"}, {"A", "B", "D"}),
                Split({"D"}, {"A", "B", "C"}),
                Split({"A", "B"}, {"C", "D"}),
            ]
        )

    def test_tree_is_tree_compatible(self) -> None:
        system = self._make_compatible_system()
        assert is_tree_compatible(system)

    def test_induced_splits_match(self) -> None:
        from phylozoo.core.network.sdnetwork.derivations import induced_splits

        system = self._make_compatible_system()
        tree = tree_from_splitsystem(system)
        result = induced_splits(tree)

        assert result.splits == system.splits

    def test_three_taxa_trivial_system(self) -> None:
        """Three trivial splits on a star topology round-trip correctly."""
        from phylozoo.core.network.sdnetwork.derivations import induced_splits

        system = SplitSystem(
            [
                Split({"X"}, {"Y", "Z"}),
                Split({"Y"}, {"X", "Z"}),
                Split({"Z"}, {"X", "Y"}),
            ]
        )
        tree = tree_from_splitsystem(system)
        result = induced_splits(tree)

        assert result.splits == system.splits

    def test_five_taxa_caterpillar(self) -> None:
        """Five-taxon caterpillar round-trips all splits."""
        from phylozoo.core.network.sdnetwork.derivations import induced_splits

        # Caterpillar: ((A,B),(C,(D,E)))
        system = SplitSystem(
            [
                Split({"A"}, {"B", "C", "D", "E"}),
                Split({"B"}, {"A", "C", "D", "E"}),
                Split({"C"}, {"A", "B", "D", "E"}),
                Split({"D"}, {"A", "B", "C", "E"}),
                Split({"E"}, {"A", "B", "C", "D"}),
                Split({"A", "B"}, {"C", "D", "E"}),
                Split({"D", "E"}, {"A", "B", "C"}),
            ]
        )
        tree = tree_from_splitsystem(system)
        result = induced_splits(tree)

        assert result.splits == system.splits
