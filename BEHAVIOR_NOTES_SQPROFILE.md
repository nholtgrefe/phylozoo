# Behavior Notes for SqQuartetProfile and SqQuartetProfileSet

## Summary

Comprehensive tests have been added for `SqQuartetProfile` and `SqQuartetProfileSet`. All 38 tests pass. Below are behaviors that might be unexpected or that you might want to change.

## Unexpected Behaviors

### 1. **SqQuartetProfileSet: Multiple Quartets with Same Taxa**

**Issue**: When multiple `Quartet` objects with the same taxa are passed to `SqQuartetProfileSet`, only the **last one is kept** instead of being merged into a single profile with multiple quartets.

**Example**:
```python
q1 = Quartet(Split({1, 2}, {3, 4}))
q2 = Quartet(Split({1, 3}, {2, 4}))  # Same taxa, different quartet
profileset = SqQuartetProfileSet(profiles=[q1, q2])
# Result: Only one profile with q2 (the last one), not a merged profile with both quartets
```

**Root Cause**: The implementation converts each `Quartet` to a `SqQuartetProfile` separately before passing to the parent class. The parent class (`QuartetProfileSet`) groups profiles by taxa, but when multiple profiles have the same taxa, it **overwrites** (keeps the last one) rather than merging.

**Expected Behavior**: Ideally, quartets with the same taxa should be merged into a single `SqQuartetProfile` containing both quartets (if they form a valid 2-quartet profile).

**Workaround**: Create the `SqQuartetProfile` explicitly:
```python
sq_profile = SqQuartetProfile([q1, q2], reticulation_leaf='A')  # if needed
profileset = SqQuartetProfileSet(profiles=[sq_profile])
```

### 2. **SqQuartetProfileSet: Quartet Weights Lost During Auto-Conversion**

**Issue**: When a `Quartet` is auto-converted to `SqQuartetProfile`, the quartet weight becomes the profile weight, but within the profile, the quartet gets weight 1.0 (default when creating from a list).

**Example**:
```python
q1 = Quartet(Split({1, 2}, {3, 4}))
profileset = SqQuartetProfileSet(profiles=[(q1, 0.6)])
profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
# Profile weight: 0.6 (correct)
# Quartet weight within profile: 1.0 (not 0.6)
```

**Root Cause**: The conversion uses `SqQuartetProfile([obj])` which assigns default weight 1.0 to the quartet.

**Expected Behavior**: The quartet weight should be preserved within the profile, or at least the profile weight should reflect the quartet weight.

**Note**: This is less critical since the profile weight is preserved, but the internal quartet weight doesn't match.

### 3. **SqQuartetProfile: Two Quartets Always Form Circular Ordering**

**Observation**: Any two resolved quartets on the same 4 taxa will always have at least one common circular ordering. This means the validation that checks for circular ordering will always pass for valid resolved quartets.

**Implication**: The test case for "two quartets without circular ordering" is not achievable with resolved quartets. The validation is still useful as a sanity check, but it will never actually fail for valid input.

## Design Decisions

### 1. **Reticulation Leaf Optional for 2-Quartet Profiles**

The `reticulation_leaf` parameter is optional even when a profile has 2 quartets. This allows creating 2-quartet profiles without specifying which leaf is the reticulation leaf.

**Question**: Should `reticulation_leaf` be required for 2-quartet profiles, or is optional acceptable?

### 2. **No Automatic Merging in SqQuartetProfileSet**

The current implementation doesn't automatically merge quartets with the same taxa when they're provided as separate `Quartet` objects. This is consistent with the parent class behavior, but might be unexpected.

**Question**: Should `SqQuartetProfileSet` automatically merge quartets with the same taxa into a single profile?

## Recommendations

1. **Fix Quartet Merging**: Modify `SqQuartetProfileSet.__init__` to detect when multiple quartets have the same taxa and merge them into a single `SqQuartetProfile` before passing to the parent class.

2. **Preserve Quartet Weights**: When auto-converting `Quartet` to `SqQuartetProfile`, preserve the quartet weight by using `SqQuartetProfile({quartet: weight})` instead of `SqQuartetProfile([quartet])`.

3. **Consider Making Reticulation Leaf Required**: For 2-quartet profiles, consider making `reticulation_leaf` a required parameter to ensure all 2-quartet profiles have this information.

## Test Coverage

All tests pass (38 tests total):
- **SqQuartetProfile**: 20 tests covering initialization, properties, edge cases
- **SqQuartetProfileSet**: 18 tests covering initialization, properties, edge cases

The tests document the current behavior, including the unexpected behaviors mentioned above.

