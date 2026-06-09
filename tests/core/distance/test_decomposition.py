"""
Tests for split decomposition of distance matrices.

Reference data is taken from Table 1 of Bandelt & Dress (1992):
  "Split Decomposition: A New and Useful Approach to Phylogenetic Analysis of
   Distance Data", Mol. Phylogenet. Evol. 1(3):242-252.

Seven taxa A–G with distance matrix d (see TABLE1_MATRIX below).
The four d-splits (coded by their minority parts) and isolation indices are:

    EFG  →  alpha = 6   ({A,B,C,D} | {E,F,G})
    AFG  →  alpha = 4   ({B,C,D,E} | {A,F,G})
    DEF  →  alpha = 2   ({A,B,C,G} | {D,E,F})
    CDE  →  alpha = 1   ({A,B,F,G} | {C,D,E})

The residual is zero (d = d^1 exactly), so Table 1 is totally decomposable.
"""

import numpy as np
import pytest

from phylozoo.core.distance import DistanceMatrix, isolation_index, split_decomposition
from phylozoo.core.distance.classifications import is_totally_decomposable, is_tree_metric
from phylozoo.core.split import Split
from phylozoo.utils.exceptions import PhyloZooValueError

# ---------------------------------------------------------------------------
# Shared fixtures / constants
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


@pytest.fixture()
def table1_dm() -> DistanceMatrix:
    return DistanceMatrix(TABLE1_MATRIX, labels=TABLE1_LABELS)


@pytest.fixture()
def path4_dm() -> DistanceMatrix:
    """Path-tree metric on four taxa: 1 — 2 — 3 — 4."""
    matrix = np.array([[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]], dtype=np.float64)
    return DistanceMatrix(matrix, labels=[1, 2, 3, 4])


# ---------------------------------------------------------------------------
# isolation_index
# ---------------------------------------------------------------------------


class TestIsolationIndex:
    """Tests for the isolation_index helper."""

    def test_efg_abcd(self, table1_dm: DistanceMatrix) -> None:
        s = Split({"E", "F", "G"}, {"A", "B", "C", "D"})
        assert isolation_index(table1_dm, s) == pytest.approx(6.0)

    def test_afg_bcde(self, table1_dm: DistanceMatrix) -> None:
        s = Split({"A", "F", "G"}, {"B", "C", "D", "E"})
        assert isolation_index(table1_dm, s) == pytest.approx(4.0)

    def test_def_abcg(self, table1_dm: DistanceMatrix) -> None:
        s = Split({"D", "E", "F"}, {"A", "B", "C", "G"})
        assert isolation_index(table1_dm, s) == pytest.approx(2.0)

    def test_cde_abfg(self, table1_dm: DistanceMatrix) -> None:
        s = Split({"C", "D", "E"}, {"A", "B", "F", "G"})
        assert isolation_index(table1_dm, s) == pytest.approx(1.0)

    def test_non_dsplit_returns_zero(self, table1_dm: DistanceMatrix) -> None:
        """A partition that is not a d-split has isolation index 0."""
        # {A,E} | {B,C,D,F,G} — not a d-split for this matrix
        s = Split({"A", "E"}, {"B", "C", "D", "F", "G"})
        assert isolation_index(table1_dm, s) == pytest.approx(0.0, abs=1e-10)

    def test_isolation_index_nonnegative(self, table1_dm: DistanceMatrix) -> None:
        """Isolation index is always non-negative."""
        labels = TABLE1_LABELS
        for k in range(1, len(labels)):
            set1 = set(labels[:k])
            set2 = set(labels[k:])
            s = Split(set1, set2)
            assert isolation_index(table1_dm, s) >= 0.0

    def test_tree_metric_trivial_split(self, path4_dm: DistanceMatrix) -> None:
        """On a tree metric the trivial split has positive isolation index."""
        s = Split({1}, {2, 3, 4})
        assert isolation_index(path4_dm, s) > 0.0

    def test_unknown_element_raises(self, table1_dm: DistanceMatrix) -> None:
        s = Split({"A", "Z"}, {"B", "C", "D", "E", "F", "G"})
        with pytest.raises(PhyloZooValueError, match="not present in the distance matrix"):
            isolation_index(table1_dm, s)

    def test_two_taxa(self) -> None:
        """For n=2 the single split has isolation index equal to the distance."""
        matrix = np.array([[0, 5], [5, 0]], dtype=np.float64)
        dm = DistanceMatrix(matrix, labels=["X", "Y"])
        s = Split({"X"}, {"Y"})
        assert isolation_index(dm, s) == pytest.approx(5.0)


# ---------------------------------------------------------------------------
# split_decomposition — split discovery
# ---------------------------------------------------------------------------


class TestSplitDecompositionSplits:
    """Tests that split_decomposition finds the correct d-splits."""

    def test_table1_split_count(self, table1_dm: DistanceMatrix) -> None:
        system, _ = split_decomposition(table1_dm)
        assert len(system.splits) == 4

    def test_table1_all_expected_splits_present(self, table1_dm: DistanceMatrix) -> None:
        system, _ = split_decomposition(table1_dm)
        expected = [
            Split({"E", "F", "G"}, {"A", "B", "C", "D"}),
            Split({"A", "F", "G"}, {"B", "C", "D", "E"}),
            Split({"D", "E", "F"}, {"A", "B", "C", "G"}),
            Split({"C", "D", "E"}, {"A", "B", "F", "G"}),
        ]
        for sp in expected:
            assert sp in system.splits, f"Expected split {sp} not found"

    def test_table1_isolation_indices(self, table1_dm: DistanceMatrix) -> None:
        system, _ = split_decomposition(table1_dm)
        expected_weights = {
            Split({"E", "F", "G"}, {"A", "B", "C", "D"}): 6.0,
            Split({"A", "F", "G"}, {"B", "C", "D", "E"}): 4.0,
            Split({"D", "E", "F"}, {"A", "B", "C", "G"}): 2.0,
            Split({"C", "D", "E"}, {"A", "B", "F", "G"}): 1.0,
        }
        for sp, expected_alpha in expected_weights.items():
            assert system.get_weight(sp) == pytest.approx(expected_alpha), (
                f"Split {sp}: expected alpha={expected_alpha}, " f"got {system.get_weight(sp)}"
            )

    def test_tree_metric_leaf_trivial_splits(self, path4_dm: DistanceMatrix) -> None:
        """On the path tree 1-2-3-4 only the leaf trivial splits are d-splits.
        Internal nodes (2 and 3) have isolation index 0 because the path
        triangle equality d_2_1 + d_2_3 = d_1_3 holds, and similarly for 3."""
        system, _ = split_decomposition(path4_dm)
        trivial_splits = {sp for sp in system.splits if sp.is_trivial}
        # Only leaves 1 and 4 yield trivial d-splits
        assert len(trivial_splits) == 2
        leaf_elements = {frozenset({1}), frozenset({4})}
        for sp in trivial_splits:
            singleton = frozenset(sp.set1) if len(sp.set1) == 1 else frozenset(sp.set2)
            assert singleton in leaf_elements

    def test_two_taxa(self) -> None:
        matrix = np.array([[0, 7], [7, 0]], dtype=np.float64)
        dm = DistanceMatrix(matrix, labels=["P", "Q"])
        system, residual = split_decomposition(dm)
        assert len(system.splits) == 1
        sp = next(iter(system.splits))
        assert sp.elements == {"P", "Q"}
        assert system.get_weight(sp) == pytest.approx(7.0)
        assert np.allclose(residual.np_array, 0.0)

    def test_single_taxon(self) -> None:
        matrix = np.array([[0]], dtype=np.float64)
        dm = DistanceMatrix(matrix, labels=["A"])
        system, residual = split_decomposition(dm)
        assert len(system.splits) == 0
        assert np.allclose(residual.np_array, 0.0)


# ---------------------------------------------------------------------------
# split_decomposition — reconstruction identity
# ---------------------------------------------------------------------------


class TestSplitDecompositionReconstruction:
    """d^1 + d^0 == d (reconstruction identity)."""

    def _d1_matrix(self, system, labels: list, n: int) -> np.ndarray:
        from phylozoo.core.split.algorithms import distances_from_splitsystem

        if len(system.splits) == 0:
            return np.zeros((n, n), dtype=np.float64)
        d1_raw = distances_from_splitsystem(system)
        return np.array(
            [[d1_raw.get_distance(labels[r], labels[c]) for c in range(n)] for r in range(n)],
            dtype=np.float64,
        )

    def test_table1_reconstruction(self, table1_dm: DistanceMatrix) -> None:
        system, residual = split_decomposition(table1_dm)
        n = len(table1_dm)
        d1 = self._d1_matrix(system, TABLE1_LABELS, n)
        reconstructed = d1 + residual.np_array
        assert np.allclose(reconstructed, TABLE1_MATRIX, atol=1e-9)

    def test_path4_reconstruction(self, path4_dm: DistanceMatrix) -> None:
        labels = list(path4_dm.labels)
        n = len(path4_dm)
        system, residual = split_decomposition(path4_dm)
        d1 = self._d1_matrix(system, labels, n)
        reconstructed = d1 + residual.np_array
        assert np.allclose(reconstructed, path4_dm.np_array, atol=1e-9)

    def test_residual_is_nonnegative(self, table1_dm: DistanceMatrix) -> None:
        """Split decomposition approximates d from below: d^0 >= 0."""
        _, residual = split_decomposition(table1_dm)
        assert np.all(residual.np_array >= -1e-10)

    def test_residual_symmetric(self, table1_dm: DistanceMatrix) -> None:
        _, residual = split_decomposition(table1_dm)
        assert np.allclose(residual.np_array, residual.np_array.T)

    def test_residual_zero_diagonal(self, table1_dm: DistanceMatrix) -> None:
        _, residual = split_decomposition(table1_dm)
        assert np.allclose(np.diag(residual.np_array), 0.0)


# ---------------------------------------------------------------------------
# is_tree_metric
# ---------------------------------------------------------------------------


class TestIsTreeMetric:
    def test_path_tree_is_tree_metric(self, path4_dm: DistanceMatrix) -> None:
        assert is_tree_metric(path4_dm) is True

    def test_table1_is_not_tree_metric(self, table1_dm: DistanceMatrix) -> None:
        assert is_tree_metric(table1_dm) is False

    def test_two_taxa_is_tree_metric(self) -> None:
        matrix = np.array([[0, 3], [3, 0]], dtype=np.float64)
        dm = DistanceMatrix(matrix)
        assert is_tree_metric(dm) is True

    def test_three_taxa_is_tree_metric(self) -> None:
        # Any metric on 3 taxa is a tree metric (star tree)
        matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]], dtype=np.float64)
        dm = DistanceMatrix(matrix)
        assert is_tree_metric(dm) is True

    def test_star_tree_is_tree_metric(self) -> None:
        """Star tree: all taxa at distance 1 from center (equidistant)."""
        # 4 taxa all at distance 2 from each other: star with edge length 1
        matrix = np.array(
            [[0, 2, 2, 2], [2, 0, 2, 2], [2, 2, 0, 2], [2, 2, 2, 0]], dtype=np.float64
        )
        dm = DistanceMatrix(matrix)
        assert is_tree_metric(dm) is True

    def test_non_tree_metric_four_taxa(self) -> None:
        """A submatrix of Table 1 (taxa A,B,E,G) fails the 4-point condition.
        s1=d_AB+d_EG=4+7=11, s2=d_AE+d_BG=13+10=23, s3=d_AG+d_BE=6+9=15;
        max=23 is achieved only by s2 so count=1 < 2."""
        matrix = np.array(
            [[0, 4, 13, 6], [4, 0, 9, 10], [13, 9, 0, 7], [6, 10, 7, 0]],
            dtype=np.float64,
        )
        dm = DistanceMatrix(matrix, labels=["A", "B", "E", "G"])
        assert is_tree_metric(dm) is False


# ---------------------------------------------------------------------------
# is_totally_decomposable
# ---------------------------------------------------------------------------


class TestIsTotallyDecomposable:
    def test_table1_totally_decomposable(self, table1_dm: DistanceMatrix) -> None:
        assert is_totally_decomposable(table1_dm) is True

    def test_path4_totally_decomposable(self, path4_dm: DistanceMatrix) -> None:
        # Tree metrics are totally decomposable
        assert is_totally_decomposable(path4_dm) is True

    def test_totally_decomposable_with_residual(self) -> None:
        """A metric with a non-zero residual is not totally decomposable."""
        # Random symmetric positive matrix that is NOT totally decomposable
        # Use the ribosomal RNA example insight: randomly chosen metrics
        # tend to have large residuals.
        rng = np.random.default_rng(42)
        n = 5
        raw = rng.uniform(1, 10, (n, n))
        raw = (raw + raw.T) / 2
        np.fill_diagonal(raw, 0)
        # Force triangle inequality by adding enough to off-diagonals... or just
        # use a known non-totally-decomposable example:
        # Take Table 1 and add an indecomposable perturbation.
        perturbed = TABLE1_MATRIX.copy()
        # Add a small symmetric perturbation to a single off-diagonal entry
        # that breaks total decomposability
        perturbed[0, 4] += 0.5
        perturbed[4, 0] += 0.5
        dm = DistanceMatrix(perturbed, labels=TABLE1_LABELS)
        assert is_totally_decomposable(dm) is False

    def test_two_taxa_totally_decomposable(self) -> None:
        matrix = np.array([[0, 4], [4, 0]], dtype=np.float64)
        dm = DistanceMatrix(matrix)
        assert is_totally_decomposable(dm) is True
