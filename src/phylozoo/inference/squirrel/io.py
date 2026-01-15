"""
I/O module for SqQuartetProfileSet.

This module provides serialization/deserialization functions for SqQuartetProfileSet
using the PhyloZoo format (.pz extension).

The PhyloZoo format is a JSON-based format that serializes:
- Taxa: List of all taxa in the profileset
- Profiles: List of profiles, each containing:
  - taxa: The four taxa for this profile
  - quartets: List of quartets with their weights
  - reticulation_leaf: Optional reticulation leaf identifier
  - profile_weight: Weight of the profile

Example Format
--------------
A minimal example of the PhyloZoo format:

.. code-block:: json

    {
        "taxa": ["A", "B", "C", "D"],
        "profiles": [
            {
                "taxa": ["A", "B", "C", "D"],
                "quartets": [
                    {
                        "quartet": {
                            "type": "resolved",
                            "split": {
                                "set1": ["A", "B"],
                                "set2": ["C", "D"]
                            }
                        },
                        "weight": 1.0
                    }
                ],
                "reticulation_leaf": null,
                "profile_weight": 1.0
            }
        ]
    }

For star quartets, the quartet structure is:

.. code-block:: json

    {
        "quartet": {
            "type": "star",
            "taxa": ["A", "B", "C", "D"]
        },
        "weight": 0.5
    }
"""

from __future__ import annotations

import json
from typing import Any

from typing import TYPE_CHECKING

from ...core.quartet.base import Quartet
from ...core.split.base import Split
from ...utils.exceptions import PhyloZooFormatError, PhyloZooParseError
from .sqprofile import SqQuartetProfile

if TYPE_CHECKING:
    from .sqprofileset import SqQuartetProfileSet


def to_pz(profileset: 'SqQuartetProfileSet', **kwargs: Any) -> str:
    """
    Convert SqQuartetProfileSet to PhyloZoo format string.
    
    Serializes the profileset to a JSON string containing all profiles,
    quartets, weights, and taxa information.
    
    Parameters
    ----------
    profileset : SqQuartetProfileSet
        The profile set to serialize.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    str
        JSON string representation of the profileset.
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> from phylozoo.core.quartet.base import Quartet
    >>> q1 = Quartet(Split({1, 2}, {3, 4}))
    >>> profileset = SqQuartetProfileSet(profiles=[q1])
    >>> pz_string = to_pz(profileset)
    >>> isinstance(pz_string, str)
    True
    """
    # Serialize profiles: extract quartets as basic data
    profiles_data = []
    for four_taxa, (profile, weight) in profileset.profiles.items():
        # Serialize quartets in this profile
        quartets_data = []
        for q, q_weight in profile.quartets.items():
            if q.is_resolved():
                # Resolved quartet: serialize split as (set1, set2)
                split = q.split
                if split is None:
                    raise PhyloZooFormatError("Resolved quartet should have a split")
                quartet_data = {
                    'type': 'resolved',
                    'split': {
                        'set1': sorted(split.set1),  # Convert to list for JSON
                        'set2': sorted(split.set2),
                    },
                }
            else:
                # Star quartet: serialize as taxa list
                quartet_data = {
                    'type': 'star',
                    'taxa': sorted(q.taxa),  # Convert to list for JSON
                }
            quartets_data.append({
                'quartet': quartet_data,
                'weight': q_weight,
            })
        
        profile_data = {
            'taxa': sorted(four_taxa),  # Convert to list for JSON
            'quartets': quartets_data,
            'reticulation_leaf': profile.reticulation_leaf,
            'profile_weight': weight,
        }
        profiles_data.append(profile_data)
    
    # Serialize taxa
    taxa_data = sorted(profileset.taxa)  # Convert to list for JSON
    
    # Create JSON structure
    data = {
        'profiles': profiles_data,
        'taxa': taxa_data,
    }
    
    return json.dumps(data, indent=None)


def from_pz(pz_string: str, **kwargs: Any) -> 'SqQuartetProfileSet':
    """
    Parse PhyloZoo format string and create SqQuartetProfileSet.
    
    Deserializes a JSON string back into a SqQuartetProfileSet object.
    
    Parameters
    ----------
    pz_string : str
        JSON string containing profileset data.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    SqQuartetProfileSet
        Parsed squirrel quartet profile set.
    
    Raises
    ------
    PhyloZooParseError
        If the JSON string is malformed or cannot be parsed.
    PhyloZooFormatError
        If the parsed data structure is invalid.
    
    Examples
    --------
    >>> pz_string = '{"profiles": [...], "taxa": [...]}'
    >>> profileset = from_pz(pz_string)
    >>> isinstance(profileset, SqQuartetProfileSet)
    True
    """
    try:
        data = json.loads(pz_string)
    except json.JSONDecodeError as e:
        raise PhyloZooParseError(f"Invalid JSON in PhyloZoo format: {e}") from e
    
    # Validate structure
    if not isinstance(data, dict):
        raise PhyloZooFormatError("PhyloZoo format must be a JSON object")
    
    if 'profiles' not in data or 'taxa' not in data:
        raise PhyloZooFormatError("PhyloZoo format must contain 'profiles' and 'taxa' keys")
    
    # Reconstruct profiles
    reconstructed_profiles = []
    for profile_data in data['profiles']:
        # Validate profile structure
        if 'taxa' not in profile_data or 'quartets' not in profile_data:
            raise PhyloZooFormatError("Profile must contain 'taxa' and 'quartets' keys")
        
        # Reconstruct quartets
        quartets_list = []
        for q_data in profile_data['quartets']:
            if 'quartet' not in q_data or 'weight' not in q_data:
                raise PhyloZooFormatError("Quartet data must contain 'quartet' and 'weight' keys")
            
            quartet_info = q_data['quartet']
            q_weight = q_data['weight']
            
            if quartet_info['type'] == 'resolved':
                # Reconstruct resolved quartet from split
                split_data = quartet_info['split']
                split = Split(set(split_data['set1']), set(split_data['set2']))
                quartet = Quartet(split)
            elif quartet_info['type'] == 'star':
                # Reconstruct star quartet from taxa
                quartet = Quartet(frozenset(quartet_info['taxa']))
            else:
                raise PhyloZooFormatError(f"Unknown quartet type: {quartet_info['type']}")
            
            quartets_list.append((quartet, q_weight))
        
        # Reconstruct SqQuartetProfile
        sq_profile = SqQuartetProfile(
            quartets=quartets_list,
            reticulation_leaf=profile_data.get('reticulation_leaf'),
        )
        
        reconstructed_profiles.append((sq_profile, profile_data['profile_weight']))
    
    # Reconstruct profileset
    # Import here to avoid circular imports and ensure it's available in worker process
    from .sqprofileset import SqQuartetProfileSet
    
    taxa = frozenset(data['taxa'])
    
    return SqQuartetProfileSet(
        profiles=reconstructed_profiles,
        taxa=taxa,
    )


def _register_formats() -> None:
    """Register format handlers with FormatRegistry."""
    from ...utils.io import FormatRegistry
    from .sqprofileset import SqQuartetProfileSet
    
    FormatRegistry.register(
        SqQuartetProfileSet,
        'pz',
        reader=from_pz,
        writer=to_pz,
        extensions=['.pz'],
        default=True,
    )


# Register formats when module is imported
_register_formats()
