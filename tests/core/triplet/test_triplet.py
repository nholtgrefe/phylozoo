"""
Tests for Triplet class.
"""

import pytest

from phylozoo.core.triplet import Triplet
from phylozoo.core.split.base import Split


class TestTripletInit:
    """Tests for Triplet initialization."""

    def test_init_from_split(self) -> None:
        """Test creating a triplet from a split."""
        split = Split({1}, {2, 3})
        triplet = Triplet(split)

        assert triplet.taxa == frozenset({1, 2, 3})
        assert triplet.split == split
        assert triplet.is_resolved()
        assert not triplet.is_star()

    @pytest.mark.parametrize("taxa_input", [{1, 2, 3}, frozenset({1, 2, 3})])
    def test_init_from_set_like_star(self, taxa_input: set[int] | frozenset[int]) -> None:
        """Test creating a star triplet from set or frozenset of taxa."""
        triplet = Triplet(taxa_input)
        assert triplet.taxa == frozenset({1, 2, 3})
        assert triplet.split is None
        assert not triplet.is_resolved()
        assert triplet.is_star()

    def test_init_non_trivial_split_error(self) -> None:
        """Test that non-trivial splits raise ValueError."""
        # A 2|2 split is non-trivial; not a valid triplet split.
        non_trivial_split = Split({1, 2}, {3, 4})

        with pytest.raises(ValueError, match="exactly 3 elements"):
            Triplet(non_trivial_split)

    def test_init_wrong_number_taxa_split(self) -> None:
        """Test that splits with wrong number of taxa raise ValueError."""
        wrong_split = Split({1}, {2, 3, 4})

        with pytest.raises(ValueError, match="Split must have exactly 3 elements"):
            Triplet(wrong_split)

    @pytest.mark.parametrize("taxa", [{1, 2}, {1, 2, 3, 4}])
    def test_init_wrong_number_taxa_star(self, taxa: set[int]) -> None:
        """Test that star trees with wrong number of taxa raise ValueError."""
        with pytest.raises(ValueError, match="Taxa must have exactly 3 elements"):
            Triplet(taxa)


class TestTripletProperties:
    """Tests for Triplet properties."""

    def test_taxa_property(self) -> None:
        """Test taxa property."""
        triplet = Triplet(Split({1}, {2, 3}))
        assert triplet.taxa == frozenset({1, 2, 3})
        assert isinstance(triplet.taxa, frozenset)

    def test_split_property(self) -> None:
        """Test split property."""
        split = Split({1}, {2, 3})
        triplet = Triplet(split)
        assert triplet.split == split

        star_triplet = Triplet({1, 2, 3})
        assert star_triplet.split is None

    def test_immutability(self) -> None:
        """Test that triplet is immutable."""
        triplet = Triplet(Split({1}, {2, 3}))

        with pytest.raises(AttributeError, match="Cannot modify attribute"):
            triplet._taxa = frozenset({5, 6, 7})

        with pytest.raises(AttributeError, match="Cannot modify attribute"):
            triplet._split = None


class TestTripletMethods:
    """Tests for Triplet methods."""

    @pytest.mark.parametrize(
        "triplet_input,expected_resolved",
        [
            (Split({1}, {2, 3}), True),
            ({1, 2, 3}, False),
        ],
    )
    def test_is_resolved(self, triplet_input: Split | set[int], expected_resolved: bool) -> None:
        """Test is_resolved method."""
        triplet = Triplet(triplet_input)
        assert triplet.is_resolved() == expected_resolved

    @pytest.mark.parametrize(
        "triplet_input,expected_star",
        [
            (Split({1}, {2, 3}), False),
            ({1, 2, 3}, True),
        ],
    )
    def test_is_star(self, triplet_input: Split | set[int], expected_star: bool) -> None:
        """Test is_star method."""
        triplet = Triplet(triplet_input)
        assert triplet.is_star() == expected_star

    @pytest.mark.parametrize(
        "triplet_input",
        [Split({1}, {2, 3}), {1, 2, 3}],
    )
    def test_copy(self, triplet_input: Split | set[int]) -> None:
        """Test copy method for resolved and star triplets."""
        original = Triplet(triplet_input)
        copied = original.copy()
        assert copied is not original
        assert copied.taxa == original.taxa
        assert copied.split == original.split
        assert copied == original

    def test_to_network_resolved(self) -> None:
        """Test to_network for resolved triplet."""
        triplet = Triplet(Split({1}, {2, 3}))
        network = triplet.to_network()

        assert network.number_of_nodes() == 5  # 3 leaves + 1 internal + 1 root
        assert network.number_of_edges() == 4  # 1 outgroup + 1 root-to-internal + 2 cherry
        assert network.taxa == {"1", "2", "3"}
        assert len(network.leaves) == 3
        # The root has the outgroup as a direct leaf child and an internal node.
        root = network.root_node
        children = list(network.children(root))
        assert len(children) == 2
        leaf_children = [c for c in children if c in network.leaves]
        assert len(leaf_children) == 1
        assert network.get_label(leaf_children[0]) == "1"

    def test_to_network_star(self) -> None:
        """Test to_network for star tree."""
        triplet = Triplet({1, 2, 3})
        network = triplet.to_network()

        assert network.number_of_nodes() == 4  # 3 leaves + 1 root
        assert network.number_of_edges() == 3  # 3 leaf edges
        assert network.taxa == {"1", "2", "3"}
        assert len(network.leaves) == 3
        # All three leaves are direct children of the root.
        root = network.root_node
        children = set(network.children(root))
        assert children == network.leaves


class TestTripletOutgroupCherry:
    """Tests for Triplet outgroup and cherry properties."""

    def test_outgroup_resolved(self) -> None:
        """Test outgroup for resolved triplet."""
        triplet = Triplet(Split({1}, {2, 3}))
        assert triplet.outgroup == frozenset({1})

    def test_cherry_resolved(self) -> None:
        """Test cherry for resolved triplet."""
        triplet = Triplet(Split({1}, {2, 3}))
        assert triplet.cherry == frozenset({2, 3})

    def test_outgroup_star(self) -> None:
        """Test outgroup is None for star triplet."""
        triplet = Triplet({1, 2, 3})
        assert triplet.outgroup is None

    def test_cherry_star(self) -> None:
        """Test cherry is None for star triplet."""
        triplet = Triplet({1, 2, 3})
        assert triplet.cherry is None

    @pytest.mark.parametrize(
        "split,expected_outgroup,expected_cherry",
        [
            (Split({1}, {2, 3}), frozenset({1}), frozenset({2, 3})),
            (Split({2, 3}, {1}), frozenset({1}), frozenset({2, 3})),
            (Split({5}, {6, 7}), frozenset({5}), frozenset({6, 7})),
        ],
    )
    def test_outgroup_cherry_independent_of_split_argument_order(
        self,
        split: Split,
        expected_outgroup: frozenset[int],
        expected_cherry: frozenset[int],
    ) -> None:
        """Test that outgroup and cherry are correct regardless of split argument order."""
        triplet = Triplet(split)
        assert triplet.outgroup == expected_outgroup
        assert triplet.cherry == expected_cherry


class TestTripletEquality:
    """Tests for Triplet equality and hashing."""

    def test_eq_same_split(self) -> None:
        """Test equality with same split."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({1}, {2, 3}))

        assert t1 == t2
        assert hash(t1) == hash(t2)

    def test_eq_different_split_same_taxa(self) -> None:
        """Test equality with different split but same taxa."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))

        assert t1 != t2

    def test_eq_same_star(self) -> None:
        """Test equality with same star tree."""
        t1 = Triplet({1, 2, 3})
        t2 = Triplet({1, 2, 3})

        assert t1 == t2
        assert hash(t1) == hash(t2)

    def test_eq_different_taxa(self) -> None:
        """Test equality with different taxa."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({4}, {5, 6}))

        assert t1 != t2

    def test_eq_resolved_vs_star(self) -> None:
        """Test equality between resolved and star on same taxa."""
        resolved = Triplet(Split({1}, {2, 3}))
        star = Triplet({1, 2, 3})

        assert resolved != star

    def test_eq_different_type(self) -> None:
        """Test equality with different type."""
        triplet = Triplet(Split({1}, {2, 3}))

        assert triplet != "not a triplet"
        assert triplet != 42
        assert triplet != Split({1}, {2, 3})


class TestTripletRepr:
    """Tests for Triplet string representation."""

    @pytest.mark.parametrize(
        "triplet_input,expect_split_in_repr,check_taxa_in_repr",
        [
            (Split({1}, {2, 3}), True, False),
            ({1, 2, 3}, False, True),
        ],
    )
    def test_repr(
        self,
        triplet_input: Split | set[int],
        expect_split_in_repr: bool,
        check_taxa_in_repr: bool,
    ) -> None:
        """Test string representation of resolved and star triplets."""
        triplet = Triplet(triplet_input)
        repr_str = repr(triplet)
        assert repr_str.startswith("Triplet(")
        if expect_split_in_repr:
            assert "Split" in repr_str
        else:
            assert "Split" not in repr_str
        if check_taxa_in_repr:
            assert "1" in repr_str and "2" in repr_str and "3" in repr_str


class TestTripletHashability:
    """Tests for Triplet hashability."""

    def test_hashable(self) -> None:
        """Test that triplets are hashable."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet({1, 2, 3})

        # Should not raise
        hash(t1)
        hash(t2)

    def test_hash_consistency(self) -> None:
        """Test that hash is consistent."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({1}, {2, 3}))

        assert hash(t1) == hash(t2)

    def test_can_use_in_set(self) -> None:
        """Test that triplets can be used in sets."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        t3 = Triplet({1, 2, 3})

        triplet_set = {t1, t2, t3}
        assert len(triplet_set) == 3
        assert t1 in triplet_set
        assert t2 in triplet_set
        assert t3 in triplet_set

    def test_can_use_in_dict(self) -> None:
        """Test that triplets can be used as dictionary keys."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))

        triplet_dict = {t1: "first", t2: "second"}
        assert triplet_dict[t1] == "first"
        assert triplet_dict[t2] == "second"
