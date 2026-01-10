#!/usr/bin/env python3
"""Test sqprofileset_from_network on LEVEL_1_SDNETWORK_TWO_HYBRIDS_SEPARATE."""

from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_TWO_HYBRIDS_SEPARATE
from phylozoo.inference.squirrel.qsimilarity import sqprofileset_from_network

# Get the profileset
network = LEVEL_1_SDNETWORK_TWO_HYBRIDS_SEPARATE
profileset = sqprofileset_from_network(network)

# Print results
print("=" * 80)
print("SqQuartetProfileSet from LEVEL_1_SDNETWORK_TWO_HYBRIDS_SEPARATE")
print("=" * 80)
print()
print(f"Repr: {repr(profileset)}")
print()
print(f"Number of profiles: {len(profileset)}")
print(f"Taxa: {sorted(profileset.taxa)}")
print()
print("Profiles (sorted by taxa set):")
print("-" * 80)
for taxa_set, (profile, weight) in sorted(profileset.profiles.items()):
    taxa_list = sorted(taxa_set)
    quartets_dict = dict(profile.quartets)
    ret_leaf = profile.reticulation_leaf
    print(f"Taxa: {taxa_list}")
    print(f"  Weight: {weight}")
    print(f"  Quartets: {quartets_dict}")
    print(f"  Reticulation leaf: {ret_leaf}")
    print()

