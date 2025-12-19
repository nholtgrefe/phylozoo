"""
Fixture networks for testing SemiDirectedPhyNetwork.

Each network is defined as a module constant and includes metadata
documenting its properties. Networks are organized into categories:
- Trees: Networks with no hybrid nodes
- Simple Hybrids: Networks with 1-2 hybrid nodes
- Multiple Blobs: Networks with multiple bi-edge connected components
- High Level: Networks with level >= 2
- Large Networks: Networks with 30+ nodes
- Parallel Edges: Networks with multiple edges between same node pairs
- Special Cases: Deep, wide, and other edge cases

Note: These networks are converted from DirectedPhyNetwork fixtures
using to_lsa_network and to_sd_network transformations.
"""

from typing import Any

from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork

# ============================================================================
# TREES
# ============================================================================

SDTREE_EMPTY = SemiDirectedPhyNetwork(
)
"""Converted from DTREE_EMPTY.

Properties:
- Nodes: 0, Edges: 0
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 0
"""

SDTREE_LARGE_BINARY = SemiDirectedPhyNetwork(
    undirected_edges=[
        (21, 1),
        (21, 2),
        (21, 31),
        (31, 22),
        (31, 36),
        (22, 3),
        (22, 4),
        (23, 5),
        (23, 6),
        (23, 32),
        (32, 24),
        (32, 36),
        (24, 7),
        (24, 8),
        (25, 9),
        (25, 10),
        (25, 33),
        (33, 26),
        (33, 37),
        (26, 11),
        (26, 12),
        (27, 13),
        (27, 14),
        (27, 34),
        (34, 28),
        (34, 37),
        (28, 15),
        (28, 16),
        (29, 17),
        (29, 18),
        (29, 35),
        (35, 30),
        (35, 38),
        (30, 19),
        (30, 20),
        (36, 39),
        (37, 39),
        (38, 1038),
        (38, 40),
        (39, 40),
        (40, 1040)
    ],
    nodes=[
        (1, {'label': 'L1'}),
        (2, {'label': 'L2'}),
        (3, {'label': 'L3'}),
        (4, {'label': 'L4'}),
        (5, {'label': 'L5'}),
        (6, {'label': 'L6'}),
        (7, {'label': 'L7'}),
        (8, {'label': 'L8'}),
        (9, {'label': 'L9'}),
        (10, {'label': 'L10'}),
        (11, {'label': 'L11'}),
        (12, {'label': 'L12'}),
        (13, {'label': 'L13'}),
        (14, {'label': 'L14'}),
        (15, {'label': 'L15'}),
        (16, {'label': 'L16'}),
        (17, {'label': 'L17'}),
        (18, {'label': 'L18'}),
        (19, {'label': 'L19'}),
        (20, {'label': 'L20'}),
        (1038, {'label': 'Dummy1038'}),
        (1040, {'label': 'Dummy1040'})
    ]
)
"""Converted from DTREE_LARGE_BINARY.

Properties:
- Nodes: 42, Edges: 41
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 0
"""

SDTREE_MEDIUM_BINARY = SemiDirectedPhyNetwork(
    undirected_edges=[
        (11, 1),
        (11, 2),
        (11, 16),
        (16, 12),
        (16, 19),
        (12, 3),
        (12, 4),
        (13, 5),
        (13, 6),
        (13, 17),
        (17, 14),
        (17, 19),
        (14, 7),
        (14, 8),
        (15, 9),
        (15, 10),
        (15, 18),
        (18, 1018),
        (18, 20),
        (19, 20),
        (20, 1020)
    ],
    nodes=[
        (1, {'label': 'L1'}),
        (2, {'label': 'L2'}),
        (3, {'label': 'L3'}),
        (4, {'label': 'L4'}),
        (5, {'label': 'L5'}),
        (6, {'label': 'L6'}),
        (7, {'label': 'L7'}),
        (8, {'label': 'L8'}),
        (9, {'label': 'L9'}),
        (10, {'label': 'L10'}),
        (1018, {'label': 'Dummy1018'}),
        (1020, {'label': 'Dummy1020'})
    ]
)
"""Converted from DTREE_MEDIUM_BINARY.

Properties:
- Nodes: 22, Edges: 21
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 0
"""

SDTREE_NON_BINARY_LARGE = SemiDirectedPhyNetwork(
    undirected_edges=[
        (100, 20),
        (100, 21),
        (100, 22),
        (100, 23),
        (100, 24),
        (20, 10),
        (20, 11),
        (20, 12),
        (21, 13),
        (21, 14),
        (22, 15),
        (22, 16),
        (23, 17),
        (23, 18),
        (24, 19),
        (24, 27),
        (10, 1),
        (10, 2),
        (11, 3),
        (11, 25),
        (12, 4),
        (12, 26),
        (13, 5),
        (13, 28),
        (14, 6),
        (14, 29),
        (15, 7),
        (15, 30),
        (16, 8),
        (16, 31),
        (17, 9),
        (17, 32),
        (18, 33),
        (18, 34),
        (19, 35),
        (19, 36)
    ],
    nodes=[
        (27, {'label': 'L12'}),
        (1, {'label': 'L1'}),
        (2, {'label': 'L2'}),
        (3, {'label': 'L3'}),
        (25, {'label': 'L10'}),
        (4, {'label': 'L4'}),
        (26, {'label': 'L11'}),
        (5, {'label': 'L5'}),
        (28, {'label': 'L13'}),
        (6, {'label': 'L6'}),
        (29, {'label': 'L14'}),
        (7, {'label': 'L7'}),
        (30, {'label': 'L15'}),
        (8, {'label': 'L8'}),
        (31, {'label': 'L16'}),
        (9, {'label': 'L9'}),
        (32, {'label': 'L17'}),
        (33, {'label': 'L18'}),
        (34, {'label': 'L19'}),
        (35, {'label': 'L20'}),
        (36, {'label': 'L21'})
    ]
)
"""Converted from DTREE_NON_BINARY_LARGE.

Properties:
- Nodes: 37, Edges: 36
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 0
"""

SDTREE_NON_BINARY_SMALL = SemiDirectedPhyNetwork(
    undirected_edges=[
        (10, 5),
        (10, 6),
        (10, 7),
        (10, 8),
        (5, 1),
        (5, 2),
        (6, 3),
        (6, 9),
        (7, 4),
        (7, 11)
    ],
    nodes=[
        (8, {'label': 'E'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'}),
        (3, {'label': 'C'}),
        (9, {'label': 'F'}),
        (4, {'label': 'D'}),
        (11, {'label': 'G'})
    ]
)
"""Converted from DTREE_NON_BINARY_SMALL.

Properties:
- Nodes: 11, Edges: 10
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 0
"""

SDTREE_SINGLE_NODE = SemiDirectedPhyNetwork(
    nodes=[
        (1, {'label': 'A'})
    ]
)
"""Converted from DTREE_SINGLE_NODE.

Properties:
- Nodes: 1, Edges: 0
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 0
"""

SDTREE_SMALL_BINARY = SemiDirectedPhyNetwork(
    undirected_edges=[
        (3, 1),
        (3, 2),
        (3, 4)
    ],
    nodes=[
        (1, {'label': 'A'}),
        (2, {'label': 'B'}),
        (4, {'label': 'C'})
    ]
)
"""Converted from DTREE_SMALL_BINARY.

Properties:
- Nodes: 4, Edges: 3
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 0
"""

# ============================================================================
# SIMPLE HYBRIDS
# ============================================================================

LEVEL_1_SDNETWORK_SINGLE_HYBRID = SemiDirectedPhyNetwork(
    directed_edges=[
        (5, 4),
        (6, 4)
    ],
    undirected_edges=[
        (5, 3),
        (5, 6),
        (6, 7),
        (4, 8),
        (8, 1),
        (8, 2)
    ],
    nodes=[
        (3, {'label': 'C'}),
        (7, {'label': 'D'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'})
    ]
)
"""Converted from LEVEL_1_DNETWORK_SINGLE_HYBRID.

Properties:
- Nodes: 8, Edges: 8
- Level: 1, Vertex level: 1, Reticulation number: 1
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 1
"""

LEVEL_1_SDNETWORK_SINGLE_HYBRID_BINARY = SemiDirectedPhyNetwork(
    directed_edges=[
        (5, 4),
        (6, 4)
    ],
    undirected_edges=[
        (5, 1),
        (5, 6),
        (6, 2),
        (4, 7),
        (7, 3),
        (7, 9)
    ],
    nodes=[
        (1, {'label': 'A'}),
        (2, {'label': 'B'}),
        (3, {'label': 'C'}),
        (9, {'label': 'D'})
    ]
)
"""Converted from LEVEL_1_DNETWORK_SINGLE_HYBRID_BINARY.

Properties:
- Nodes: 8, Edges: 8
- Level: 1, Vertex level: 1, Reticulation number: 1
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 1
"""

LEVEL_1_SDNETWORK_TWO_HYBRIDS_SEPARATE = SemiDirectedPhyNetwork(
    directed_edges=[
        (5, 4),
        (6, 4),
        (7, 9),
        (8, 9)
    ],
    undirected_edges=[
        (5, 11),
        (5, 6),
        (11, 7),
        (11, 8),
        (6, 3),
        (7, 12),
        (8, 13),
        (4, 1),
        (9, 2)
    ],
    nodes=[
        (3, {'label': 'C'}),
        (1, {'label': 'A'}),
        (12, {'label': 'D'}),
        (13, {'label': 'E'}),
        (2, {'label': 'B'})
    ]
)
"""Converted from LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE.

Properties:
- Nodes: 12, Edges: 13
- Level: 1, Vertex level: 1, Reticulation number: 2
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 2
"""

LEVEL_2_SDNETWORK_DIAMOND_HYBRID = SemiDirectedPhyNetwork(
    directed_edges=[
        (7, 5),
        (7, 6),
        (8, 5),
        (8, 6),
        (11, 6)
    ],
    undirected_edges=[
        (10, 7),
        (10, 8),
        (10, 11),
        (11, 14),
        (11, 15),
        (5, 12),
        (12, 1),
        (12, 2),
        (6, 13),
        (13, 3),
        (13, 4)
    ],
    nodes=[
        (14, {'label': 'E'}),
        (15, {'label': 'F'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'}),
        (3, {'label': 'C'}),
        (4, {'label': 'D'})
    ]
)
"""Converted from LEVEL_2_DNETWORK_DIAMOND_HYBRID.

Properties:
- Nodes: 14, Edges: 16
- Level: 3, Vertex level: 2, Reticulation number: 3
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 2
"""

LEVEL_2_SDNETWORK_NESTED_HYBRIDS = SemiDirectedPhyNetwork(
    directed_edges=[
        (15, 12),
        (16, 12),
        (10, 8),
        (10, 7),
        (11, 8),
        (5, 7)
    ],
    undirected_edges=[
        (15, 21),
        (15, 16),
        (21, 22),
        (21, 23),
        (16, 4),
        (12, 13),
        (13, 10),
        (13, 11),
        (11, 17),
        (8, 9),
        (9, 5),
        (9, 6),
        (5, 18),
        (6, 3),
        (6, 19),
        (7, 14),
        (14, 1),
        (14, 2)
    ],
    nodes=[
        (4, {'label': 'D'}),
        (22, {'label': 'H'}),
        (23, {'label': 'I'}),
        (17, {'label': 'E'}),
        (18, {'label': 'F'}),
        (3, {'label': 'C'}),
        (19, {'label': 'G'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'})
    ]
)
"""Converted from LEVEL_2_DNETWORK_NESTED_HYBRIDS.

Properties:
- Nodes: 22, Edges: 24
- Level: 2, Vertex level: 2, Reticulation number: 3
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 3
"""

LEVEL_2_SDNETWORK_TRIANGLE_HYBRID = SemiDirectedPhyNetwork(
    directed_edges=[
        (7, 6),
        (5, 4),
        (5, 6),
        (6, 4)
    ],
    undirected_edges=[
        (7, 5),
        (7, 3),
        (4, 9),
        (9, 1),
        (9, 2)
    ],
    nodes=[
        (3, {'label': 'C'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'})
    ]
)
"""Converted from LEVEL_2_DNETWORK_TRIANGLE_HYBRID.

Properties:
- Nodes: 8, Edges: 9
- Level: 2, Vertex level: 2, Reticulation number: 2
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 2
"""

LEVEL_2_SDNETWORK_TWO_HYBRIDS_SAME_BLOB = SemiDirectedPhyNetwork(
    directed_edges=[
        (10, 8),
        (10, 9),
        (11, 8),
        (8, 9)
    ],
    undirected_edges=[
        (15, 10),
        (15, 11),
        (15, 4),
        (11, 3),
        (9, 12),
        (12, 1),
        (12, 2)
    ],
    nodes=[
        (4, {'label': 'D'}),
        (3, {'label': 'C'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'})
    ]
)
"""Converted from LEVEL_2_DNETWORK_TWO_HYBRIDS_SAME_BLOB.

Properties:
- Nodes: 10, Edges: 11
- Level: 2, Vertex level: 2, Reticulation number: 2
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 2
"""

# ============================================================================
# MULTIPLE BLOBS
# ============================================================================

LEVEL_1_SDNETWORK_CHAIN_OF_BLOBS = SemiDirectedPhyNetwork(
    directed_edges=[
        (20, 10),
        (21, 10),
        (22, 11),
        (23, 11),
        (24, 12),
        (25, 12),
        (27, 13),
        (28, 13)
    ],
    undirected_edges=[
        (20, 2),
        (20, 21),
        (21, 35),
        (10, 14),
        (14, 1),
        (14, 18),
        (18, 30),
        (18, 32),
        (30, 22),
        (30, 23),
        (22, 36),
        (23, 37),
        (11, 15),
        (15, 3),
        (15, 19),
        (19, 40),
        (19, 33),
        (40, 24),
        (40, 25),
        (24, 38),
        (25, 39),
        (12, 16),
        (16, 4),
        (16, 29),
        (29, 46),
        (29, 34),
        (46, 27),
        (46, 28),
        (27, 6),
        (28, 41),
        (13, 17),
        (17, 5),
        (17, 31)
    ],
    nodes=[
        (2, {'label': 'L2'}),
        (35, {'label': 'X35'}),
        (1, {'label': 'L1'}),
        (32, {'label': '32'}),
        (36, {'label': 'X36'}),
        (37, {'label': 'X37'}),
        (3, {'label': 'L3'}),
        (40, {'label': 'X40'}),
        (33, {'label': '33'}),
        (38, {'label': 'X38'}),
        (39, {'label': 'X39'}),
        (4, {'label': 'L4'}),
        (34, {'label': '34'}),
        (6, {'label': 'L6'}),
        (41, {'label': 'X41'}),
        (5, {'label': 'L5'}),
        (31, {'label': '31'})
    ]
)
"""Converted from LEVEL_1_DNETWORK_CHAIN_OF_BLOBS.

Properties:
- Nodes: 38, Edges: 41
- Level: 1, Vertex level: 1, Reticulation number: 4
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 4
"""

LEVEL_1_SDNETWORK_FIVE_BLOBS = SemiDirectedPhyNetwork(
    directed_edges=[
        (50, 40),
        (51, 40),
        (52, 41),
        (53, 41),
        (54, 42),
        (55, 42),
        (56, 43),
        (57, 43),
        (58, 44),
        (59, 44)
    ],
    undirected_edges=[
        (50, 6),
        (50, 51),
        (51, 60),
        (60, 52),
        (60, 53),
        (52, 7),
        (53, 70),
        (40, 45),
        (45, 1),
        (45, 11),
        (70, 54),
        (70, 55),
        (54, 8),
        (55, 80),
        (41, 46),
        (46, 2),
        (46, 12),
        (80, 56),
        (80, 57),
        (56, 9),
        (57, 90),
        (42, 47),
        (47, 3),
        (47, 13),
        (90, 58),
        (90, 59),
        (58, 10),
        (59, 16),
        (43, 48),
        (48, 4),
        (48, 14),
        (44, 49),
        (49, 5),
        (49, 15)
    ],
    nodes=[
        (6, {'label': 'L6'}),
        (1, {'label': 'L1'}),
        (11, {'label': 'L11'}),
        (7, {'label': 'L7'}),
        (2, {'label': 'L2'}),
        (12, {'label': 'L12'}),
        (8, {'label': 'L8'}),
        (3, {'label': 'L3'}),
        (13, {'label': 'L13'}),
        (9, {'label': 'L9'}),
        (4, {'label': 'L4'}),
        (14, {'label': 'L14'}),
        (10, {'label': 'L10'}),
        (16, {'label': 'L16'}),
        (5, {'label': 'L5'}),
        (15, {'label': 'L15'})
    ]
)
"""Converted from LEVEL_1_DNETWORK_FIVE_BLOBS.

Properties:
- Nodes: 40, Edges: 44
- Level: 1, Vertex level: 1, Reticulation number: 5
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 5
"""

LEVEL_1_SDNETWORK_STAR_BLOBS = SemiDirectedPhyNetwork(
    directed_edges=[
        (51, 40),
        (52, 40),
        (61, 41),
        (62, 41),
        (71, 42),
        (72, 42),
        (81, 43),
        (82, 43)
    ],
    undirected_edges=[
        (100, 50),
        (100, 60),
        (100, 70),
        (100, 80),
        (50, 51),
        (50, 52),
        (60, 61),
        (60, 62),
        (70, 71),
        (70, 72),
        (80, 81),
        (80, 82),
        (51, 5),
        (52, 15),
        (61, 6),
        (62, 16),
        (71, 7),
        (72, 9),
        (81, 8),
        (82, 14),
        (40, 44),
        (44, 1),
        (44, 10),
        (41, 45),
        (45, 2),
        (45, 11),
        (42, 46),
        (46, 3),
        (46, 12),
        (43, 47),
        (47, 4),
        (47, 13)
    ],
    nodes=[
        (5, {'label': 'L5'}),
        (15, {'label': 'L15'}),
        (6, {'label': 'L6'}),
        (16, {'label': 'L16'}),
        (7, {'label': 'L7'}),
        (9, {'label': 'L9'}),
        (8, {'label': 'L8'}),
        (14, {'label': 'L14'}),
        (1, {'label': 'L1'}),
        (10, {'label': 'L10'}),
        (2, {'label': 'L2'}),
        (11, {'label': 'L11'}),
        (3, {'label': 'L3'}),
        (12, {'label': 'L12'}),
        (4, {'label': 'L4'}),
        (13, {'label': 'L13'})
    ]
)
"""Converted from LEVEL_1_DNETWORK_STAR_BLOBS.

Properties:
- Nodes: 37, Edges: 40
- Level: 1, Vertex level: 1, Reticulation number: 4
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 4
"""

LEVEL_1_SDNETWORK_TWO_BLOBS = SemiDirectedPhyNetwork(
    directed_edges=[
        (5, 4),
        (6, 4),
        (7, 9),
        (8, 9)
    ],
    undirected_edges=[
        (5, 3),
        (5, 6),
        (6, 11),
        (11, 7),
        (11, 8),
        (7, 12),
        (8, 13),
        (4, 14),
        (14, 1),
        (14, 16),
        (9, 15),
        (15, 2),
        (15, 17)
    ],
    nodes=[
        (3, {'label': 'C'}),
        (1, {'label': 'A'}),
        (16, {'label': 'F'}),
        (12, {'label': 'D'}),
        (13, {'label': 'E'}),
        (2, {'label': 'B'}),
        (17, {'label': 'G'})
    ]
)
"""Converted from LEVEL_1_DNETWORK_TWO_BLOBS.

Properties:
- Nodes: 16, Edges: 17
- Level: 1, Vertex level: 1, Reticulation number: 2
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 2
"""

LEVEL_2_SDNETWORK_THREE_BLOBS = SemiDirectedPhyNetwork(
    directed_edges=[
        (20, 15),
        (21, 15),
        (22, 16),
        (22, 17),
        (23, 16),
        (29, 17)
    ],
    undirected_edges=[
        (30, 20),
        (30, 21),
        (30, 6),
        (20, 5),
        (21, 25),
        (25, 22),
        (25, 23),
        (23, 26),
        (15, 28),
        (28, 1),
        (28, 32),
        (26, 27),
        (26, 7),
        (27, 3),
        (27, 4),
        (16, 29),
        (29, 34),
        (17, 31),
        (31, 2),
        (31, 33)
    ],
    nodes=[
        (6, {'label': 'F'}),
        (5, {'label': 'E'}),
        (1, {'label': 'A'}),
        (32, {'label': 'H'}),
        (7, {'label': 'G'}),
        (34, {'label': 'J'}),
        (2, {'label': 'B'}),
        (33, {'label': 'I'}),
        (3, {'label': 'C'}),
        (4, {'label': 'D'})
    ]
)
"""Converted from LEVEL_2_DNETWORK_THREE_BLOBS.

Properties:
- Nodes: 24, Edges: 26
- Level: 2, Vertex level: 2, Reticulation number: 3
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 3
"""

SDTREE_BALANCED = SemiDirectedPhyNetwork(
    undirected_edges=[
        (20, 10),
        (20, 11),
        (20, 12),
        (20, 13),
        (20, 14),
        (10, 1),
        (10, 6),
        (11, 2),
        (11, 7),
        (12, 3),
        (12, 8),
        (13, 4),
        (13, 9),
        (14, 5),
        (14, 15)
    ],
    nodes=[
        (1, {'label': 'L1'}),
        (6, {'label': 'L6'}),
        (2, {'label': 'L2'}),
        (7, {'label': 'L7'}),
        (3, {'label': 'L3'}),
        (8, {'label': 'L8'}),
        (4, {'label': 'L4'}),
        (9, {'label': 'L9'}),
        (5, {'label': 'L5'}),
        (15, {'label': 'L10'})
    ]
)
"""Converted from DTREE_BALANCED.

Properties:
- Nodes: 16, Edges: 15
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 0
"""

# ============================================================================
# HIGH LEVEL
# ============================================================================

LEVEL_2_SDNETWORK_MULTIPLE_BLOBS = SemiDirectedPhyNetwork(
    directed_edges=[
        (40, 30),
        (41, 30),
        (31, 20),
        (31, 21),
        (32, 20),
        (42, 33),
        (43, 33),
        (61, 21),
        (34, 22),
        (34, 23),
        (35, 22),
        (64, 23)
    ],
    undirected_edges=[
        (40, 3),
        (40, 41),
        (41, 12),
        (41, 9),
        (30, 60),
        (60, 31),
        (60, 32),
        (31, 11),
        (32, 45),
        (45, 42),
        (45, 43),
        (42, 4),
        (43, 13),
        (43, 10),
        (20, 61),
        (61, 16),
        (21, 66),
        (66, 62),
        (66, 7),
        (62, 1),
        (62, 17),
        (33, 63),
        (63, 34),
        (63, 35),
        (34, 68),
        (35, 6),
        (22, 64),
        (64, 14),
        (23, 67),
        (67, 65),
        (67, 8),
        (65, 2),
        (65, 15)
    ],
    nodes=[
        (3, {'label': 'C'}),
        (12, {'label': 'L'}),
        (9, {'label': 'I'}),
        (11, {'label': 'K'}),
        (16, {'label': 'P'}),
        (7, {'label': 'G'}),
        (1, {'label': 'A'}),
        (17, {'label': 'Q'}),
        (4, {'label': 'D'}),
        (13, {'label': 'M'}),
        (10, {'label': 'J'}),
        (68, {'label': 'E'}),
        (6, {'label': 'F'}),
        (14, {'label': 'N'}),
        (8, {'label': 'H'}),
        (2, {'label': 'B'}),
        (15, {'label': 'O'})
    ]
)
"""Converted from LEVEL_2_DNETWORK_MULTIPLE_BLOBS.

Properties:
- Nodes: 40, Edges: 45
- Level: 2, Vertex level: 2, Reticulation number: 6
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 6
"""

LEVEL_2_SDNETWORK_PARALLEL_HYBRIDS = SemiDirectedPhyNetwork(
    directed_edges=[
        (20, 15),
        (20, 10),
        (21, 15),
        (22, 15),
        (16, 10),
        (17, 10)
    ],
    undirected_edges=[
        (30, 20),
        (30, 21),
        (30, 22),
        (21, 3),
        (22, 4),
        (15, 50),
        (50, 16),
        (50, 17),
        (16, 5),
        (17, 6),
        (10, 51),
        (51, 1),
        (51, 2)
    ],
    nodes=[
        (3, {'label': 'C'}),
        (4, {'label': 'D'}),
        (5, {'label': 'E'}),
        (6, {'label': 'F'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'})
    ]
)
"""Converted from LEVEL_2_DNETWORK_PARALLEL_HYBRIDS.

Properties:
- Nodes: 16, Edges: 19
- Level: 4, Vertex level: 2, Reticulation number: 4
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 2
"""

LEVEL_2_SDNETWORK_SINGLE_BLOB = SemiDirectedPhyNetwork(
    directed_edges=[
        (25, 20),
        (26, 20),
        (18, 15),
        (19, 15),
        (12, 10),
        (12, 8),
        (13, 10),
        (29, 8)
    ],
    undirected_edges=[
        (25, 22),
        (25, 26),
        (26, 5),
        (20, 27),
        (27, 18),
        (27, 19),
        (18, 21),
        (19, 4),
        (15, 28),
        (28, 12),
        (28, 13),
        (13, 3),
        (10, 29),
        (29, 23),
        (8, 31),
        (31, 1),
        (31, 2)
    ],
    nodes=[
        (22, {'label': 'L7'}),
        (5, {'label': 'L5'}),
        (21, {'label': 'L6'}),
        (4, {'label': 'L4'}),
        (3, {'label': 'L3'}),
        (23, {'label': 'L8'}),
        (1, {'label': 'L1'}),
        (2, {'label': 'L2'})
    ]
)
"""Converted from LEVEL_2_DNETWORK_SINGLE_BLOB.

Properties:
- Nodes: 22, Edges: 25
- Level: 2, Vertex level: 2, Reticulation number: 4
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 4
"""

LEVEL_3_DNETWORK = SemiDirectedPhyNetwork(
    directed_edges=[
        (25, 20),
        (26, 20),
        (18, 15),
        (19, 15),
        (232, 234),
        (233, 234),
        (12, 10),
        (12, 8),
        (13, 10),
        (29, 8)
    ],
    undirected_edges=[
        (25, 22),
        (25, 26),
        (26, 5),
        (20, 27),
        (27, 18),
        (27, 19),
        (18, 21),
        (19, 4),
        (15, 28),
        (28, 232),
        (28, 233),
        (232, 12),
        (233, 13),
        (13, 3),
        (234, 235),
        (235, 236),
        (235, 237),
        (10, 29),
        (29, 23),
        (8, 31),
        (31, 1),
        (31, 2)
    ],
    nodes=[
        (22, {'label': 'L7'}),
        (5, {'label': 'L5'}),
        (21, {'label': 'L6'}),
        (4, {'label': 'L4'}),
        (3, {'label': 'L3'}),
        (236, {'label': 'L236'}),
        (237, {'label': 'L237'}),
        (23, {'label': 'L8'}),
        (1, {'label': 'L1'}),
        (2, {'label': 'L2'})
    ]
)
"""Converted from LEVEL_3_DNETWORK.

Properties:
- Nodes: 28, Edges: 32
- Level: 3, Vertex level: 3, Reticulation number: 5
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 5
"""

LEVEL_3_SDNETWORK_HIGH_VERTEX_LEVEL = SemiDirectedPhyNetwork(
    directed_edges=[
        (40, 30),
        (40, 20),
        (41, 30),
        (42, 30),
        (31, 20),
        (31, 10),
        (32, 20),
        (21, 10),
        (22, 10)
    ],
    undirected_edges=[
        (50, 40),
        (50, 41),
        (50, 42),
        (41, 7),
        (42, 4),
        (30, 60),
        (60, 31),
        (60, 32),
        (32, 3),
        (20, 61),
        (61, 21),
        (61, 22),
        (21, 5),
        (22, 6),
        (10, 62),
        (62, 1),
        (62, 2)
    ],
    nodes=[
        (7, {'label': '7'}),
        (4, {'label': 'D'}),
        (3, {'label': 'C'}),
        (5, {'label': 'E'}),
        (6, {'label': 'F'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'})
    ]
)
"""Converted from LEVEL_3_DNETWORK_HIGH_VERTEX_LEVEL.

Properties:
- Nodes: 21, Edges: 26
- Level: 6, Vertex level: 3, Reticulation number: 6
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 3
"""

# ============================================================================
# LARGE NETWORKS
# ============================================================================

LEVEL_1_SDNETWORK_LARGE_DEEP_HYBRID_CHAIN = SemiDirectedPhyNetwork(
    directed_edges=[
        (90, 80),
        (91, 80),
        (81, 70),
        (82, 70),
        (71, 60),
        (72, 60),
        (61, 50),
        (62, 50),
        (51, 40),
        (52, 40),
        (41, 30),
        (42, 30)
    ],
    undirected_edges=[
        (90, 3),
        (90, 91),
        (91, 24),
        (80, 101),
        (101, 81),
        (101, 82),
        (101, 9),
        (81, 22),
        (81, 4),
        (82, 23),
        (70, 102),
        (102, 71),
        (102, 72),
        (102, 10),
        (71, 20),
        (71, 5),
        (72, 21),
        (60, 103),
        (103, 61),
        (103, 62),
        (103, 11),
        (61, 18),
        (61, 6),
        (62, 19),
        (50, 104),
        (104, 51),
        (104, 52),
        (104, 12),
        (51, 16),
        (51, 7),
        (52, 17),
        (40, 105),
        (105, 41),
        (105, 42),
        (105, 13),
        (41, 14),
        (41, 8),
        (42, 15),
        (30, 106),
        (106, 1),
        (106, 2)
    ],
    nodes=[
        (3, {'label': 'L3'}),
        (24, {'label': 'L24'}),
        (9, {'label': 'L9'}),
        (22, {'label': 'L22'}),
        (4, {'label': 'L4'}),
        (23, {'label': 'L23'}),
        (10, {'label': 'L10'}),
        (20, {'label': 'L20'}),
        (5, {'label': 'L5'}),
        (21, {'label': 'L21'}),
        (11, {'label': 'L11'}),
        (18, {'label': 'L18'}),
        (6, {'label': 'L6'}),
        (19, {'label': 'L19'}),
        (12, {'label': 'L12'}),
        (16, {'label': 'L16'}),
        (7, {'label': 'L7'}),
        (17, {'label': 'L17'}),
        (13, {'label': 'L13'}),
        (14, {'label': 'L14'}),
        (8, {'label': 'L8'}),
        (15, {'label': 'L15'}),
        (1, {'label': 'L1'}),
        (2, {'label': 'L2'})
    ]
)
"""Converted from LEVEL_1_DNETWORK_LARGE_DEEP_HYBRID_CHAIN.

Properties:
- Nodes: 48, Edges: 53
- Level: 1, Vertex level: 1, Reticulation number: 6
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 6
"""

LEVEL_1_SDNETWORK_LARGE_FEW_HYBRIDS = SemiDirectedPhyNetwork(
    directed_edges=[
        (31, 1054),
        (32, 1054),
        (34, 1055),
        (35, 1055)
    ],
    undirected_edges=[
        (31, 1),
        (31, 2),
        (31, 46),
        (46, 32),
        (46, 54),
        (32, 3),
        (32, 4),
        (33, 5),
        (33, 6),
        (33, 47),
        (47, 34),
        (47, 54),
        (34, 7),
        (34, 8),
        (35, 9),
        (35, 10),
        (35, 48),
        (48, 36),
        (48, 55),
        (36, 11),
        (36, 12),
        (37, 13),
        (37, 14),
        (37, 49),
        (49, 38),
        (49, 55),
        (38, 15),
        (38, 16),
        (39, 17),
        (39, 18),
        (39, 50),
        (50, 40),
        (50, 56),
        (40, 19),
        (40, 20),
        (41, 21),
        (41, 22),
        (41, 51),
        (51, 42),
        (51, 56),
        (42, 23),
        (42, 24),
        (43, 25),
        (43, 26),
        (43, 52),
        (52, 44),
        (52, 57),
        (44, 27),
        (44, 28),
        (45, 29),
        (45, 30),
        (45, 53),
        (53, 1053),
        (53, 57),
        (54, 58),
        (55, 58),
        (56, 59),
        (57, 59),
        (58, 59),
        (1054, 1056),
        (1056, 1058),
        (1056, 1059),
        (1055, 1057),
        (1057, 1060),
        (1057, 1061)
    ],
    nodes=[
        (1, {'label': 'L1'}),
        (2, {'label': 'L2'}),
        (3, {'label': 'L3'}),
        (4, {'label': 'L4'}),
        (5, {'label': 'L5'}),
        (6, {'label': 'L6'}),
        (7, {'label': 'L7'}),
        (8, {'label': 'L8'}),
        (9, {'label': 'L9'}),
        (10, {'label': 'L10'}),
        (11, {'label': 'L11'}),
        (12, {'label': 'L12'}),
        (13, {'label': 'L13'}),
        (14, {'label': 'L14'}),
        (15, {'label': 'L15'}),
        (16, {'label': 'L16'}),
        (17, {'label': 'L17'}),
        (18, {'label': 'L18'}),
        (19, {'label': 'L19'}),
        (20, {'label': 'L20'}),
        (21, {'label': 'L21'}),
        (22, {'label': 'L22'}),
        (23, {'label': 'L23'}),
        (24, {'label': 'L24'}),
        (25, {'label': 'L25'}),
        (26, {'label': 'L26'}),
        (27, {'label': 'L27'}),
        (28, {'label': 'L28'}),
        (29, {'label': 'L29'}),
        (30, {'label': 'L30'}),
        (1053, {'label': 'Dummy1053'}),
        (1058, {'label': 'HL0'}),
        (1059, {'label': 'HL1'}),
        (1060, {'label': 'HL2'}),
        (1061, {'label': 'HL3'})
    ]
)
"""Converted from LEVEL_1_DNETWORK_LARGE_FEW_HYBRIDS.

Properties:
- Nodes: 68, Edges: 69
- Level: 1, Vertex level: 1, Reticulation number: 2
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 2
"""

LEVEL_1_SDNETWORK_LARGE_GLUED = SemiDirectedPhyNetwork(
    directed_edges=[
        (100, 90),
        (101, 90),
        (104, 95),
        (105, 95)
    ],
    undirected_edges=[
        (100, 2),
        (100, 101),
        (101, 110),
        (110, 102),
        (110, 103),
        (102, 91),
        (102, 92),
        (103, 93),
        (103, 94),
        (103, 120),
        (90, 96),
        (96, 1),
        (96, 18),
        (91, 3),
        (91, 4),
        (92, 5),
        (92, 6),
        (93, 7),
        (93, 8),
        (94, 9),
        (94, 10),
        (120, 104),
        (120, 105),
        (104, 12),
        (105, 130),
        (130, 106),
        (130, 107),
        (130, 108),
        (106, 13),
        (106, 14),
        (107, 15),
        (107, 20),
        (108, 16),
        (108, 17),
        (95, 97),
        (97, 11),
        (97, 19)
    ],
    nodes=[
        (2, {'label': 'L2'}),
        (1, {'label': 'L1'}),
        (18, {'label': 'L18'}),
        (3, {'label': 'L3'}),
        (4, {'label': 'L4'}),
        (5, {'label': 'L5'}),
        (6, {'label': 'L6'}),
        (7, {'label': 'L7'}),
        (8, {'label': 'L8'}),
        (9, {'label': 'L9'}),
        (10, {'label': 'L10'}),
        (12, {'label': 'L12'}),
        (11, {'label': 'L11'}),
        (19, {'label': 'L19'}),
        (13, {'label': 'L13'}),
        (14, {'label': 'L14'}),
        (15, {'label': 'L15'}),
        (20, {'label': 'L20'}),
        (16, {'label': 'L16'}),
        (17, {'label': 'L17'})
    ]
)
"""Converted from LEVEL_1_DNETWORK_LARGE_GLUED.

Properties:
- Nodes: 40, Edges: 41
- Level: 1, Vertex level: 1, Reticulation number: 2
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 2
"""

LEVEL_1_SDNETWORK_LARGE_MULTI_BLOB = SemiDirectedPhyNetwork(
    directed_edges=[
        (10, 8),
        (11, 8),
        (30, 28),
        (31, 28),
        (50, 48),
        (51, 48),
        (70, 68),
        (71, 68)
    ],
    undirected_edges=[
        (15, 10),
        (15, 11),
        (15, 35),
        (10, 2),
        (11, 3),
        (35, 30),
        (35, 31),
        (35, 55),
        (30, 22),
        (31, 23),
        (55, 50),
        (55, 51),
        (55, 75),
        (8, 12),
        (12, 1),
        (12, 4),
        (50, 42),
        (51, 43),
        (75, 70),
        (75, 71),
        (28, 32),
        (32, 21),
        (32, 24),
        (70, 62),
        (71, 63),
        (48, 52),
        (52, 41),
        (52, 44),
        (68, 72),
        (72, 61),
        (72, 64)
    ],
    nodes=[
        (2, {'label': 'B0L2'}),
        (3, {'label': 'B0L3'}),
        (1, {'label': 'B0L1'}),
        (4, {'label': 'B0L4'}),
        (22, {'label': 'B1L2'}),
        (23, {'label': 'B1L3'}),
        (21, {'label': 'B1L1'}),
        (24, {'label': 'B1L4'}),
        (42, {'label': 'B2L2'}),
        (43, {'label': 'B2L3'}),
        (41, {'label': 'B2L1'}),
        (44, {'label': 'B2L4'}),
        (62, {'label': 'B3L2'}),
        (63, {'label': 'B3L3'}),
        (61, {'label': 'B3L1'}),
        (64, {'label': 'B3L4'})
    ]
)
"""Converted from LEVEL_1_DNETWORK_LARGE_MULTI_BLOB.

Properties:
- Nodes: 36, Edges: 39
- Level: 1, Vertex level: 1, Reticulation number: 4
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 4
"""

LEVEL_3_SDNETWORK_LARGE_MANY_HYBRIDS = SemiDirectedPhyNetwork(
    directed_edges=[
        (26, 1050),
        (39, 1054),
        (27, 1050),
        (29, 1051),
        (30, 1051),
        (32, 1052),
        (33, 1052),
        (35, 1053),
        (36, 1053),
        (38, 1054)
    ],
    undirected_edges=[
        (26, 1),
        (26, 2),
        (26, 39),
        (39, 27),
        (39, 46),
        (27, 3),
        (27, 4),
        (28, 5),
        (28, 6),
        (28, 40),
        (40, 29),
        (40, 46),
        (29, 7),
        (29, 8),
        (30, 9),
        (30, 10),
        (30, 41),
        (41, 31),
        (41, 47),
        (31, 11),
        (31, 12),
        (32, 13),
        (32, 14),
        (32, 42),
        (42, 33),
        (42, 47),
        (33, 15),
        (33, 16),
        (34, 17),
        (34, 18),
        (34, 43),
        (43, 35),
        (43, 48),
        (35, 19),
        (35, 20),
        (36, 21),
        (36, 22),
        (36, 44),
        (44, 37),
        (44, 48),
        (37, 23),
        (37, 24),
        (38, 25),
        (38, 1038),
        (38, 45),
        (45, 1045),
        (45, 49),
        (46, 50),
        (47, 50),
        (48, 51),
        (49, 1049),
        (49, 51),
        (50, 51),
        (1050, 1055),
        (1055, 1060),
        (1055, 1061),
        (1051, 1056),
        (1056, 1062),
        (1056, 1063),
        (1052, 1057),
        (1057, 1064),
        (1057, 1065),
        (1053, 1058),
        (1058, 1066),
        (1058, 1067),
        (1054, 1059),
        (1059, 1068),
        (1059, 1069)
    ],
    nodes=[
        (1, {'label': 'L1'}),
        (2, {'label': 'L2'}),
        (3, {'label': 'L3'}),
        (4, {'label': 'L4'}),
        (5, {'label': 'L5'}),
        (6, {'label': 'L6'}),
        (7, {'label': 'L7'}),
        (8, {'label': 'L8'}),
        (9, {'label': 'L9'}),
        (10, {'label': 'L10'}),
        (11, {'label': 'L11'}),
        (12, {'label': 'L12'}),
        (13, {'label': 'L13'}),
        (14, {'label': 'L14'}),
        (15, {'label': 'L15'}),
        (16, {'label': 'L16'}),
        (17, {'label': 'L17'}),
        (18, {'label': 'L18'}),
        (19, {'label': 'L19'}),
        (20, {'label': 'L20'}),
        (21, {'label': 'L21'}),
        (22, {'label': 'L22'}),
        (23, {'label': 'L23'}),
        (24, {'label': 'L24'}),
        (25, {'label': 'L25'}),
        (1038, {'label': 'Dummy1038'}),
        (1045, {'label': 'Dummy1045'}),
        (1049, {'label': 'Dummy1049'}),
        (1060, {'label': 'HL0'}),
        (1061, {'label': 'HL1'}),
        (1062, {'label': 'HL2'}),
        (1063, {'label': 'HL3'}),
        (1064, {'label': 'HL4'}),
        (1065, {'label': 'HL5'}),
        (1066, {'label': 'HL6'}),
        (1067, {'label': 'HL7'}),
        (1068, {'label': 'HL8'}),
        (1069, {'label': 'HL9'})
    ]
)
"""Converted from LEVEL_3_DNETWORK_LARGE_MANY_HYBRIDS.

Properties:
- Nodes: 74, Edges: 78
- Level: 3, Vertex level: 3, Reticulation number: 5
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 5
"""

SDTREE_VERY_LARGE = SemiDirectedPhyNetwork(
    undirected_edges=[
        (51, 1),
        (51, 2),
        (51, 76),
        (76, 52),
        (76, 89),
        (52, 3),
        (52, 4),
        (53, 5),
        (53, 6),
        (53, 77),
        (77, 54),
        (77, 89),
        (54, 7),
        (54, 8),
        (55, 9),
        (55, 10),
        (55, 78),
        (78, 56),
        (78, 90),
        (56, 11),
        (56, 12),
        (57, 13),
        (57, 14),
        (57, 79),
        (79, 58),
        (79, 90),
        (58, 15),
        (58, 16),
        (59, 17),
        (59, 18),
        (59, 80),
        (80, 60),
        (80, 91),
        (60, 19),
        (60, 20),
        (61, 21),
        (61, 22),
        (61, 81),
        (81, 62),
        (81, 91),
        (62, 23),
        (62, 24),
        (63, 25),
        (63, 26),
        (63, 82),
        (82, 64),
        (82, 92),
        (64, 27),
        (64, 28),
        (65, 29),
        (65, 30),
        (65, 83),
        (83, 66),
        (83, 92),
        (66, 31),
        (66, 32),
        (67, 33),
        (67, 34),
        (67, 84),
        (84, 68),
        (84, 93),
        (68, 35),
        (68, 36),
        (69, 37),
        (69, 38),
        (69, 85),
        (85, 70),
        (85, 93),
        (70, 39),
        (70, 40),
        (71, 41),
        (71, 42),
        (71, 86),
        (86, 72),
        (86, 94),
        (72, 43),
        (72, 44),
        (73, 45),
        (73, 46),
        (73, 87),
        (87, 74),
        (87, 94),
        (74, 47),
        (74, 48),
        (75, 49),
        (75, 50),
        (75, 88),
        (88, 1088),
        (88, 95),
        (89, 96),
        (90, 96),
        (91, 97),
        (92, 97),
        (93, 98),
        (94, 98),
        (95, 1095),
        (95, 99),
        (96, 100),
        (97, 100),
        (98, 101),
        (99, 1099),
        (99, 101),
        (100, 101)
    ],
    nodes=[
        (1, {'label': 'L1'}),
        (2, {'label': 'L2'}),
        (3, {'label': 'L3'}),
        (4, {'label': 'L4'}),
        (5, {'label': 'L5'}),
        (6, {'label': 'L6'}),
        (7, {'label': 'L7'}),
        (8, {'label': 'L8'}),
        (9, {'label': 'L9'}),
        (10, {'label': 'L10'}),
        (11, {'label': 'L11'}),
        (12, {'label': 'L12'}),
        (13, {'label': 'L13'}),
        (14, {'label': 'L14'}),
        (15, {'label': 'L15'}),
        (16, {'label': 'L16'}),
        (17, {'label': 'L17'}),
        (18, {'label': 'L18'}),
        (19, {'label': 'L19'}),
        (20, {'label': 'L20'}),
        (21, {'label': 'L21'}),
        (22, {'label': 'L22'}),
        (23, {'label': 'L23'}),
        (24, {'label': 'L24'}),
        (25, {'label': 'L25'}),
        (26, {'label': 'L26'}),
        (27, {'label': 'L27'}),
        (28, {'label': 'L28'}),
        (29, {'label': 'L29'}),
        (30, {'label': 'L30'}),
        (31, {'label': 'L31'}),
        (32, {'label': 'L32'}),
        (33, {'label': 'L33'}),
        (34, {'label': 'L34'}),
        (35, {'label': 'L35'}),
        (36, {'label': 'L36'}),
        (37, {'label': 'L37'}),
        (38, {'label': 'L38'}),
        (39, {'label': 'L39'}),
        (40, {'label': 'L40'}),
        (41, {'label': 'L41'}),
        (42, {'label': 'L42'}),
        (43, {'label': 'L43'}),
        (44, {'label': 'L44'}),
        (45, {'label': 'L45'}),
        (46, {'label': 'L46'}),
        (47, {'label': 'L47'}),
        (48, {'label': 'L48'}),
        (49, {'label': 'L49'}),
        (50, {'label': 'L50'}),
        (1088, {'label': 'Dummy1088'}),
        (1095, {'label': 'Dummy1095'}),
        (1099, {'label': 'Dummy1099'})
    ]
)
"""Converted from DTREE_VERY_LARGE.

Properties:
- Nodes: 104, Edges: 103
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 0
"""

# ============================================================================
# PARALLEL EDGES
# ============================================================================

LEVEL_1_SDNETWORK_PARALLEL_EDGES = SemiDirectedPhyNetwork(
    directed_edges=[
        (10, 5),
        (10, 5, 1)
    ],
    undirected_edges=[
        (10, 6),
        (6, 3),
        (6, 4),
        (5, 11),
        (11, 1),
        (11, 2)
    ],
    nodes=[
        (3, {'label': 'C'}),
        (4, {'label': 'D'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'})
    ]
)
"""Converted from LEVEL_1_DNETWORK_PARALLEL_EDGES.

Properties:
- Nodes: 8, Edges: 8
- Level: 1, Vertex level: 1, Reticulation number: 1
- Is tree: False, Is binary: True
- Has parallel edges: True
- Number of hybrid nodes: 1
"""

LEVEL_1_SDNETWORK_PARALLEL_EDGES_HYBRID = SemiDirectedPhyNetwork(
    directed_edges=[
        (10, 8),
        (10, 8, 1),
        (11, 8),
        (11, 8, 1)
    ],
    undirected_edges=[
        (10, 3),
        (10, 11),
        (11, 4),
        (8, 9),
        (9, 5),
        (9, 6),
        (5, 1),
        (5, 2)
    ],
    nodes=[
        (3, {'label': 'C'}),
        (4, {'label': 'D'}),
        (6, {'label': 'E'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'})
    ]
)
"""Converted from LEVEL_1_DNETWORK_PARALLEL_EDGES_HYBRID.

Properties:
- Nodes: 10, Edges: 12
- Level: 3, Vertex level: 1, Reticulation number: 3
- Is tree: False, Is binary: False
- Has parallel edges: True
- Number of hybrid nodes: 1
"""

LEVEL_1_SDNETWORK_PARALLEL_IN_BLOB = SemiDirectedPhyNetwork(
    directed_edges=[
        (15, 12),
        (15, 12, 1),
        (16, 12),
        (16, 12, 1),
        (13, 10),
        (13, 10, 1)
    ],
    undirected_edges=[
        (15, 3),
        (15, 16),
        (16, 4),
        (12, 13),
        (10, 11),
        (11, 1),
        (11, 2)
    ],
    nodes=[
        (3, {'label': 'C'}),
        (4, {'label': 'D'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'})
    ]
)
"""Converted from LEVEL_1_DNETWORK_PARALLEL_IN_BLOB.

Properties:
- Nodes: 10, Edges: 13
- Level: 3, Vertex level: 1, Reticulation number: 4
- Is tree: False, Is binary: False
- Has parallel edges: True
- Number of hybrid nodes: 2
"""

LEVEL_2_SDNETWORK_MANY_PARALLEL_EDGES = SemiDirectedPhyNetwork(
    directed_edges=[
        (20, 15),
        (20, 15, 1),
        (20, 15, 2),
        (16, 10),
        (16, 10, 1),
        (25, 10),
        (25, 10, 1)
    ],
    undirected_edges=[
        (20, 16),
        (16, 4),
        (15, 25),
        (25, 3),
        (10, 11),
        (11, 5),
        (11, 7),
        (5, 1),
        (5, 2)
    ],
    nodes=[
        (4, {'label': 'D'}),
        (3, {'label': 'C'}),
        (7, {'label': 'E'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'})
    ]
)
"""Converted from LEVEL_2_DNETWORK_MANY_PARALLEL_EDGES.

Properties:
- Nodes: 12, Edges: 16
- Level: 5, Vertex level: 2, Reticulation number: 5
- Is tree: False, Is binary: False
- Has parallel edges: True
- Number of hybrid nodes: 2
"""

# ============================================================================
# SPECIAL CASES
# ============================================================================

LEVEL_1_SDNETWORK_MIXED_TOPOLOGY = SemiDirectedPhyNetwork(
    directed_edges=[
        (81, 71),
        (82, 71),
        (51, 40),
        (52, 40)
    ],
    undirected_edges=[
        (100, 80),
        (100, 81),
        (100, 82),
        (100, 83),
        (80, 50),
        (80, 10),
        (81, 12),
        (82, 13),
        (83, 72),
        (83, 73),
        (83, 74),
        (50, 1),
        (50, 11),
        (72, 4),
        (72, 7),
        (72, 25),
        (73, 5),
        (73, 8),
        (74, 6),
        (74, 9),
        (71, 90),
        (90, 51),
        (90, 52),
        (51, 3),
        (52, 14),
        (40, 91),
        (91, 2),
        (91, 15)
    ],
    nodes=[
        (10, {'label': 'L10'}),
        (12, {'label': 'L12'}),
        (13, {'label': 'L13'}),
        (1, {'label': 'L1'}),
        (11, {'label': 'L11'}),
        (4, {'label': 'L4'}),
        (7, {'label': 'L7'}),
        (25, {'label': 'L16'}),
        (5, {'label': 'L5'}),
        (8, {'label': 'L8'}),
        (6, {'label': 'L6'}),
        (9, {'label': 'L9'}),
        (3, {'label': 'L3'}),
        (14, {'label': 'L14'}),
        (2, {'label': 'L2'}),
        (15, {'label': 'L15'})
    ]
)
"""Converted from LEVEL_1_DNETWORK_MIXED_TOPOLOGY.

Properties:
- Nodes: 31, Edges: 32
- Level: 1, Vertex level: 1, Reticulation number: 2
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 2
"""

LEVEL_1_SDNETWORK_NON_LSA = SemiDirectedPhyNetwork(
    directed_edges=[
        (5, 4),
        (6, 4)
    ],
    undirected_edges=[
        (10, 5),
        (10, 6),
        (10, 136),
        (5, 3),
        (6, 7),
        (4, 8),
        (8, 1),
        (8, 2)
    ],
    nodes=[
        (136, {'label': 'E'}),
        (3, {'label': 'C'}),
        (7, {'label': 'D'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'})
    ]
)
"""Converted from LEVEL_1_DNETWORK_NON_LSA.

Properties:
- Nodes: 10, Edges: 10
- Level: 1, Vertex level: 1, Reticulation number: 1
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 1
"""

LEVEL_2_SDNETWORK_NON_LSA = SemiDirectedPhyNetwork(
    directed_edges=[
        (25, 20),
        (26, 20),
        (142, 143),
        (142, 143, 1),
        (18, 15),
        (19, 15),
        (12, 10),
        (12, 8),
        (13, 10),
        (29, 8)
    ],
    undirected_edges=[
        (30, 25),
        (30, 26),
        (30, 142),
        (25, 22),
        (26, 5),
        (20, 27),
        (27, 18),
        (27, 19),
        (18, 21),
        (19, 4),
        (15, 28),
        (28, 12),
        (28, 13),
        (13, 3),
        (10, 29),
        (29, 23),
        (8, 31),
        (31, 1),
        (31, 2),
        (143, 144)
    ],
    nodes=[
        (22, {'label': 'L7'}),
        (5, {'label': 'L5'}),
        (21, {'label': 'L6'}),
        (4, {'label': 'L4'}),
        (3, {'label': 'L3'}),
        (23, {'label': 'L8'}),
        (1, {'label': 'L1'}),
        (2, {'label': 'L2'}),
        (144, {'label': 'L144'})
    ]
)
"""Converted from LEVEL_2_DNETWORK_NON_LSA.

Properties:
- Nodes: 26, Edges: 30
- Level: 2, Vertex level: 2, Reticulation number: 5
- Is tree: False, Is binary: True
- Has parallel edges: True
- Number of hybrid nodes: 5
"""

LEVEL_2_SDNETWORK_NON_LSA_NESTED = SemiDirectedPhyNetwork(
    directed_edges=[
        (15, 12),
        (16, 12),
        (10, 8),
        (10, 7),
        (11, 8),
        (5, 7)
    ],
    undirected_edges=[
        (20, 15),
        (20, 16),
        (20, 155),
        (15, 21),
        (16, 4),
        (21, 22),
        (21, 23),
        (12, 13),
        (13, 10),
        (13, 11),
        (11, 17),
        (8, 9),
        (9, 5),
        (9, 6),
        (5, 18),
        (6, 3),
        (6, 19),
        (7, 14),
        (14, 1),
        (14, 2)
    ],
    nodes=[
        (155, {'label': 'J'}),
        (4, {'label': 'D'}),
        (22, {'label': 'H'}),
        (23, {'label': 'I'}),
        (17, {'label': 'E'}),
        (1, {'label': 'A'}),
        (2, {'label': 'B'}),
        (18, {'label': 'F'}),
        (3, {'label': 'C'}),
        (19, {'label': 'G'})
    ]
)
"""Converted from LEVEL_2_DNETWORK_NON_LSA_NESTED.

Properties:
- Nodes: 24, Edges: 26
- Level: 2, Vertex level: 2, Reticulation number: 3
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 3
"""

LEVEL_3_SDNETWORK_CHAIN_HYBRIDS = SemiDirectedPhyNetwork(
    directed_edges=[
        (40, 30),
        (40, 20),
        (41, 30),
        (31, 20),
        (31, 10),
        (21, 10)
    ],
    undirected_edges=[
        (40, 41),
        (41, 3),
        (30, 60),
        (60, 31),
        (60, 4),
        (20, 61),
        (61, 21),
        (61, 5),
        (21, 7),
        (10, 62),
        (62, 1),
        (62, 2)
    ],
    nodes=[
        (3, {'label': 'L3'}),
        (4, {'label': 'L4'}),
        (5, {'label': 'L5'}),
        (7, {'label': 'L7'}),
        (1, {'label': 'L1'}),
        (2, {'label': 'L2'})
    ]
)
"""Converted from LEVEL_3_DNETWORK_CHAIN_HYBRIDS.

Properties:
- Nodes: 16, Edges: 18
- Level: 3, Vertex level: 3, Reticulation number: 3
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 3
"""

SDTREE_SIMPLE = SemiDirectedPhyNetwork(
    undirected_edges=[
        (15, 1),
        (15, 2),
        (15, 3)
    ],
    nodes=[
        (1, {'label': 'A'}),
        (2, {'label': 'B'}),
        (3, {'label': 'C'})
    ]
)
"""Converted from DTREE_SIMPLE.

Properties:
- Nodes: 4, Edges: 3
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 0
"""

# ============================================================================
# METADATA REGISTRY
# ============================================================================

NETWORK_METADATA: dict[str, dict[str, Any]] = {
    'LEVEL_1_SDNETWORK_CHAIN_OF_BLOBS': {
        'category': 'multiple_blobs',
        'nodes': 38,
        'edges': 41,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 4,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 4,
    },
    'LEVEL_1_SDNETWORK_FIVE_BLOBS': {
        'category': 'multiple_blobs',
        'nodes': 40,
        'edges': 44,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 5,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 5,
    },
    'LEVEL_1_SDNETWORK_LARGE_DEEP_HYBRID_CHAIN': {
        'category': 'large_networks',
        'nodes': 48,
        'edges': 53,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 6,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 6,
    },
    'LEVEL_1_SDNETWORK_LARGE_FEW_HYBRIDS': {
        'category': 'large_networks',
        'nodes': 68,
        'edges': 69,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 2,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 2,
    },
    'LEVEL_1_SDNETWORK_LARGE_GLUED': {
        'category': 'large_networks',
        'nodes': 40,
        'edges': 41,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 2,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 2,
    },
    'LEVEL_1_SDNETWORK_LARGE_MULTI_BLOB': {
        'category': 'large_networks',
        'nodes': 36,
        'edges': 39,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 4,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 4,
    },
    'LEVEL_1_SDNETWORK_MIXED_TOPOLOGY': {
        'category': 'special_cases',
        'nodes': 31,
        'edges': 32,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 2,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 2,
    },
    'LEVEL_1_SDNETWORK_NON_LSA': {
        'category': 'non_lsa',
        'nodes': 10,
        'edges': 10,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 1,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 1,
    },
    'LEVEL_1_SDNETWORK_PARALLEL_EDGES': {
        'category': 'parallel_edges',
        'nodes': 8,
        'edges': 8,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 1,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': True,
        'num_hybrids': 1,
    },
    'LEVEL_1_SDNETWORK_PARALLEL_EDGES_HYBRID': {
        'category': 'parallel_edges',
        'nodes': 10,
        'edges': 12,
        'level': 3,
        'vertex_level': 1,
        'reticulation_number': 3,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': True,
        'num_hybrids': 1,
    },
    'LEVEL_1_SDNETWORK_PARALLEL_IN_BLOB': {
        'category': 'parallel_edges',
        'nodes': 10,
        'edges': 13,
        'level': 3,
        'vertex_level': 1,
        'reticulation_number': 4,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': True,
        'num_hybrids': 2,
    },
    'LEVEL_1_SDNETWORK_SINGLE_HYBRID': {
        'category': 'simple_hybrids',
        'nodes': 8,
        'edges': 8,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 1,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 1,
    },
    'LEVEL_1_SDNETWORK_SINGLE_HYBRID_BINARY': {
        'category': 'simple_hybrids',
        'nodes': 8,
        'edges': 8,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 1,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 1,
    },
    'LEVEL_1_SDNETWORK_STAR_BLOBS': {
        'category': 'multiple_blobs',
        'nodes': 37,
        'edges': 40,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 4,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 4,
    },
    'LEVEL_1_SDNETWORK_TWO_BLOBS': {
        'category': 'multiple_blobs',
        'nodes': 16,
        'edges': 17,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 2,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 2,
    },
    'LEVEL_1_SDNETWORK_TWO_HYBRIDS_SEPARATE': {
        'category': 'simple_hybrids',
        'nodes': 12,
        'edges': 13,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 2,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 2,
    },
    'LEVEL_2_SDNETWORK_DIAMOND_HYBRID': {
        'category': 'simple_hybrids',
        'nodes': 14,
        'edges': 16,
        'level': 3,
        'vertex_level': 2,
        'reticulation_number': 3,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 2,
    },
    'LEVEL_2_SDNETWORK_MANY_PARALLEL_EDGES': {
        'category': 'parallel_edges',
        'nodes': 12,
        'edges': 16,
        'level': 5,
        'vertex_level': 2,
        'reticulation_number': 5,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': True,
        'num_hybrids': 2,
    },
    'LEVEL_2_SDNETWORK_MULTIPLE_BLOBS': {
        'category': 'high_level',
        'nodes': 40,
        'edges': 45,
        'level': 2,
        'vertex_level': 2,
        'reticulation_number': 6,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 6,
    },
    'LEVEL_2_SDNETWORK_NESTED_HYBRIDS': {
        'category': 'simple_hybrids',
        'nodes': 22,
        'edges': 24,
        'level': 2,
        'vertex_level': 2,
        'reticulation_number': 3,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 3,
    },
    'LEVEL_2_SDNETWORK_NON_LSA': {
        'category': 'non_lsa',
        'nodes': 26,
        'edges': 30,
        'level': 2,
        'vertex_level': 2,
        'reticulation_number': 5,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': True,
        'num_hybrids': 5,
    },
    'LEVEL_2_SDNETWORK_NON_LSA_NESTED': {
        'category': 'non_lsa',
        'nodes': 24,
        'edges': 26,
        'level': 2,
        'vertex_level': 2,
        'reticulation_number': 3,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 3,
    },
    'LEVEL_2_SDNETWORK_PARALLEL_HYBRIDS': {
        'category': 'high_level',
        'nodes': 16,
        'edges': 19,
        'level': 4,
        'vertex_level': 2,
        'reticulation_number': 4,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 2,
    },
    'LEVEL_2_SDNETWORK_SINGLE_BLOB': {
        'category': 'high_level',
        'nodes': 22,
        'edges': 25,
        'level': 2,
        'vertex_level': 2,
        'reticulation_number': 4,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 4,
    },
    'LEVEL_2_SDNETWORK_THREE_BLOBS': {
        'category': 'multiple_blobs',
        'nodes': 24,
        'edges': 26,
        'level': 2,
        'vertex_level': 2,
        'reticulation_number': 3,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 3,
    },
    'LEVEL_2_SDNETWORK_TRIANGLE_HYBRID': {
        'category': 'simple_hybrids',
        'nodes': 8,
        'edges': 9,
        'level': 2,
        'vertex_level': 2,
        'reticulation_number': 2,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 2,
    },
    'LEVEL_2_SDNETWORK_TWO_HYBRIDS_SAME_BLOB': {
        'category': 'simple_hybrids',
        'nodes': 10,
        'edges': 11,
        'level': 2,
        'vertex_level': 2,
        'reticulation_number': 2,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 2,
    },
    'LEVEL_3_DNETWORK': {
        'category': 'high_level',
        'nodes': 28,
        'edges': 32,
        'level': 3,
        'vertex_level': 3,
        'reticulation_number': 5,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 5,
    },
    'LEVEL_3_SDNETWORK_CHAIN_HYBRIDS': {
        'category': 'special_cases',
        'nodes': 16,
        'edges': 18,
        'level': 3,
        'vertex_level': 3,
        'reticulation_number': 3,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 3,
    },
    'LEVEL_3_SDNETWORK_HIGH_VERTEX_LEVEL': {
        'category': 'high_level',
        'nodes': 21,
        'edges': 26,
        'level': 6,
        'vertex_level': 3,
        'reticulation_number': 6,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 3,
    },
    'LEVEL_3_SDNETWORK_LARGE_MANY_HYBRIDS': {
        'category': 'large_networks',
        'nodes': 74,
        'edges': 78,
        'level': 3,
        'vertex_level': 3,
        'reticulation_number': 5,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 5,
    },
    'SDTREE_BALANCED': {
        'category': 'multiple_blobs',
        'nodes': 16,
        'edges': 15,
        'level': 0,
        'vertex_level': 0,
        'reticulation_number': 0,
        'is_tree': True,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 0,
    },
    'SDTREE_EMPTY': {
        'category': 'trees',
        'nodes': 0,
        'edges': 0,
        'level': 0,
        'vertex_level': 0,
        'reticulation_number': 0,
        'is_tree': True,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 0,
    },
    'SDTREE_LARGE_BINARY': {
        'category': 'trees',
        'nodes': 42,
        'edges': 41,
        'level': 0,
        'vertex_level': 0,
        'reticulation_number': 0,
        'is_tree': True,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 0,
    },
    'SDTREE_MEDIUM_BINARY': {
        'category': 'trees',
        'nodes': 22,
        'edges': 21,
        'level': 0,
        'vertex_level': 0,
        'reticulation_number': 0,
        'is_tree': True,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 0,
    },
    'SDTREE_NON_BINARY_LARGE': {
        'category': 'trees',
        'nodes': 37,
        'edges': 36,
        'level': 0,
        'vertex_level': 0,
        'reticulation_number': 0,
        'is_tree': True,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 0,
    },
    'SDTREE_NON_BINARY_SMALL': {
        'category': 'trees',
        'nodes': 11,
        'edges': 10,
        'level': 0,
        'vertex_level': 0,
        'reticulation_number': 0,
        'is_tree': True,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 0,
    },
    'SDTREE_SIMPLE': {
        'category': 'special_cases',
        'nodes': 4,
        'edges': 3,
        'level': 0,
        'vertex_level': 0,
        'reticulation_number': 0,
        'is_tree': True,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 0,
    },
    'SDTREE_SINGLE_NODE': {
        'category': 'trees',
        'nodes': 1,
        'edges': 0,
        'level': 0,
        'vertex_level': 0,
        'reticulation_number': 0,
        'is_tree': True,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 0,
    },
    'SDTREE_SMALL_BINARY': {
        'category': 'trees',
        'nodes': 4,
        'edges': 3,
        'level': 0,
        'vertex_level': 0,
        'reticulation_number': 0,
        'is_tree': True,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 0,
    },
    'SDTREE_VERY_LARGE': {
        'category': 'large_networks',
        'nodes': 104,
        'edges': 103,
        'level': 0,
        'vertex_level': 0,
        'reticulation_number': 0,
        'is_tree': True,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 0,
    },
}


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_networks_by_category(category: str) -> list[SemiDirectedPhyNetwork]:
    """
    Get all networks in a specific category.
    
    Parameters
    ----------
    category : str
        Category name (e.g., 'trees', 'simple_hybrids', 'high_level').
    
    Returns
    -------
    list[SemiDirectedPhyNetwork]
        List of networks in the specified category.
    
    Examples
    --------
    >>> trees = get_networks_by_category('trees')
    >>> len(trees) >= 5
    True
    """
    return [globals()[k] for k, v in NETWORK_METADATA.items() 
            if v['category'] == category]


def get_networks_with_property(**kwargs: Any) -> list[SemiDirectedPhyNetwork]:
    """
    Get networks matching all specified properties.
    
    Parameters
    ----------
    **kwargs : Any
        Property-value pairs to match (e.g., level=2, is_tree=True).
    
    Returns
    -------
    list[SemiDirectedPhyNetwork]
        List of networks matching all specified properties.
    
    Examples
    --------
    >>> level_2_nets = get_networks_with_property(level=2)
    >>> all(net.number_of_nodes() > 0 for net in level_2_nets)
    True
    
    >>> trees = get_networks_with_property(is_tree=True, is_binary=True)
    >>> len(trees) >= 3
    True
    """
    matching = []
    for name, meta in NETWORK_METADATA.items():
        if all(meta.get(k) == v for k, v in kwargs.items()):
            matching.append(globals()[name])
    return matching


def get_all_networks() -> list[SemiDirectedPhyNetwork]:
    """
    Get all fixture networks.
    
    Returns
    -------
    list[SemiDirectedPhyNetwork]
        List of all networks.
    
    Examples
    --------
    >>> all_nets = get_all_networks()
    >>> len(all_nets) >= 40
    True
    """
    return [globals()[k] for k in NETWORK_METADATA.keys()]


def get_network_metadata(network_name: str) -> dict[str, Any]:
    """
    Get metadata for a specific network.
    
    Parameters
    ----------
    network_name : str
        Name of the network (e.g., 'SDTREE_SMALL_BINARY').
    
    Returns
    -------
    dict[str, Any]
        Metadata dictionary for the network.
    
    Examples
    --------
    >>> meta = get_network_metadata('SDTREE_SMALL_BINARY')
    >>> meta['category']
    'trees'
    """
    return NETWORK_METADATA.get(network_name, {})