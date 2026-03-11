"""
Fixture networks for testing DirectedPhyNetwork.

Each network is defined as a module constant and includes metadata
documenting its properties. Networks are organized into categories:
- Trees: Networks with no hybrid nodes
- Simple Hybrids: Networks with 1-2 hybrid nodes
- Multiple Blobs: Networks with multiple bi-edge connected components
- High Level: Networks with level >= 2
- Large Networks: Networks with 30+ nodes
- Parallel Edges: Networks with multiple edges between same node pairs
- Special Cases: Deep, wide, and other edge cases
"""

from typing import Any

from phylozoo.core.network import DirectedPhyNetwork

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def _build_balanced_binary_tree(num_leaves: int) -> tuple[list[tuple[int, int]], list[tuple[int, dict]]]:
    """
    Build a balanced binary tree with specified number of leaves.
    
    Parameters
    ----------
    num_leaves : int
        Number of leaf nodes to create.
    
    Returns
    -------
    tuple[list[tuple[int, int]], list[tuple[int, dict]]]
        Tuple of (edges, nodes) for the tree.
    """
    edges = []
    nodes = [(i, {'label': f'L{i}'}) for i in range(1, num_leaves + 1)]
    
    # Build tree bottom-up
    current_level = list(range(1, num_leaves + 1))
    next_node_id = num_leaves + 1
    
    while len(current_level) > 1:
        next_level = []
        for i in range(0, len(current_level), 2):
            if i + 1 < len(current_level):
                # Pair exists
                parent = next_node_id
                edges.append((parent, current_level[i]))
                edges.append((parent, current_level[i + 1]))
                next_level.append(parent)
                next_node_id += 1
            else:
                # Odd one out - connect to a new parent with a sibling leaf to avoid pass-through
                parent = next_node_id
                edges.append((parent, current_level[i]))
                # Add a dummy leaf to ensure parent has out-degree >= 2
                dummy_leaf = next_node_id + 1000  # Use high ID to avoid conflicts
                edges.append((parent, dummy_leaf))
                nodes.append((dummy_leaf, {'label': f'Dummy{dummy_leaf}'}))
                next_level.append(parent)
                next_node_id += 1
        current_level = next_level
    
    return edges, nodes


def _create_simple_hybrid(root_id: int, hybrid_id: int, parent1_id: int, parent2_id: int,
                         tree_node_id: int, leaf_ids: list[int], extra_leaf_ids: list[int]) -> list[tuple[int, int]]:
    """
    Create a simple hybrid structure: root -> parents -> hybrid -> tree_node -> leaves.
    
    This pattern ensures the hybrid node (in-degree 2) has out-degree 1, which is valid.
    
    Parameters
    ----------
    root_id : int
        ID of root node.
    hybrid_id : int
        ID of the hybrid node.
    parent1_id, parent2_id : int
        IDs of the two parent nodes leading to the hybrid.
    tree_node_id : int
        ID of tree node that receives from hybrid.
    leaf_ids : list[int]
        IDs of leaves that come from tree_node.
    extra_leaf_ids : list[int]
        IDs of extra leaves from parents (to satisfy degree constraints).
    
    Returns
    -------
    list[tuple[int, int]]
        List of edges forming the hybrid structure.
    """
    edges = [
        (root_id, parent1_id), (root_id, parent2_id),  # Root to parents
        (parent1_id, hybrid_id), (parent2_id, hybrid_id),  # Parents to hybrid
        (hybrid_id, tree_node_id),  # Hybrid to tree node (out-degree 1)
    ]
    
    # Tree node to leaves
    for leaf in leaf_ids:
        edges.append((tree_node_id, leaf))
    
    # Extra leaves from parents (to ensure out-degree >= 2)
    if len(extra_leaf_ids) >= 2:
        edges.append((parent1_id, extra_leaf_ids[0]))
        edges.append((parent2_id, extra_leaf_ids[1]))
    
    return edges


# ============================================================================
# TREES
# ============================================================================

DTREE_EMPTY = DirectedPhyNetwork(edges=[], nodes=[])
"""Empty network with no nodes or edges.

Properties:
- Nodes: 0, Edges: 0
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 0"""

DTREE_SINGLE_NODE = DirectedPhyNetwork(
    nodes=[(1, {'label': 'A'})]
)
"""Single node network (root is also a leaf).

Properties:
- Nodes: 1, Edges: 0
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 0"""

DTREE_SMALL_BINARY = DirectedPhyNetwork(
    edges=[(5, 3), (5, 4), (3, 1), (3, 2)],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
)
"""Small binary tree with 3 leaves.

Properties:
- Nodes: 5, Edges: 4
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 0"""

# Build medium binary tree with 10 leaves using helper
_medium_tree_edges, _medium_tree_nodes = _build_balanced_binary_tree(10)
DTREE_MEDIUM_BINARY = DirectedPhyNetwork(edges=_medium_tree_edges, nodes=_medium_tree_nodes)
"""Medium binary tree with 10 leaves.

Properties:
- Nodes: 23, Edges: 22
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 0"""

_large_tree_edges, _large_tree_nodes = _build_balanced_binary_tree(20)
DTREE_LARGE_BINARY = DirectedPhyNetwork(edges=_large_tree_edges, nodes=_large_tree_nodes)
"""Large binary tree with 20 leaves.

Properties:
- Nodes: 43, Edges: 42
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 0"""

DTREE_NON_BINARY_SMALL = DirectedPhyNetwork(
    edges=[(10, 5), (10, 6), (10, 7), (10, 8), (5, 1), (5, 2), (6, 3), (6, 9), (7, 4), (7, 11)],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), 
           (4, {'label': 'D'}), (8, {'label': 'E'}), (9, {'label': 'F'}), (11, {'label': 'G'})]
)
"""Small non-binary tree (root has out-degree 4).

Properties:
- Nodes: 11, Edges: 10
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 0"""

DTREE_NON_BINARY_LARGE = DirectedPhyNetwork(
    edges=[
        (100, 20), (100, 21), (100, 22), (100, 23), (100, 24),  # Root with out-degree 5
        (20, 10), (20, 11), (20, 12),  # Out-degree 3
        (21, 13), (21, 14),
        (22, 15), (22, 16),
        (23, 17), (23, 18),
        (24, 19), (24, 27),  # Node 24 now has out-degree 2
        (10, 1), (10, 2), (11, 3), (11, 25), (12, 4), (12, 26), 
        (13, 5), (13, 28), (14, 6), (14, 29),  # Nodes 13, 14 now have out-degree 2
        (15, 7), (15, 30), (16, 8), (16, 31),  # Nodes 15, 16 now have out-degree 2
        (17, 9), (17, 32), (18, 33), (18, 34), (19, 35), (19, 36)  # Nodes 17, 18, 19 have out-degree 2
    ],
    nodes=[(i, {'label': f'L{i}'}) for i in range(1, 10)] + 
          [(i, {'label': f'L{i-15}'}) for i in range(25, 37)]
)
"""Large non-binary tree with varied out-degrees.

Properties:
- Nodes: 37, Edges: 36
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 0"""

# ============================================================================
# SIMPLE HYBRIDS
# ============================================================================

LEVEL_1_DNETWORK_SINGLE_HYBRID = DirectedPhyNetwork(
    edges=[
        (10, 5), (10, 6),  # Root to tree nodes
        (5, 4), (6, 4),    # Both lead to hybrid 4 (in-degree 2)
        (4, 8),            # Hybrid to tree node
        (8, 1), (8, 2),    # Tree node to leaves
        (5, 3), (6, 7)     # Additional leaves to satisfy degree constraints
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
)
"""Network with single hybrid node, level 1.

Properties:
- Nodes: 9, Edges: 9
- Level: 1, Vertex level: 1, Reticulation number: 1
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 1"""

LEVEL_1_DNETWORK_SINGLE_HYBRID_BINARY = DirectedPhyNetwork(
    edges=[
        (8, 5), (8, 6),
        (5, 4), (6, 4),
        (4, 7),
        (5, 1), (6, 2), (7, 3), (7, 9)
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (9, {'label': 'D'})]
)
"""Binary network with single hybrid node.

Properties:
- Nodes: 9, Edges: 9
- Level: 1, Vertex level: 1, Reticulation number: 1
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 1"""

LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE = DirectedPhyNetwork(
    edges=[
        # First blob with hybrid 4
        (10, 5), (10, 6), (5, 4), (6, 4), (4, 1),
        # Cut-edge connecting
        (5, 11),
        # Second blob with hybrid 9
        (11, 7), (11, 8), (7, 9), (8, 9), (9, 2),
        # Additional leaves for degree constraints
        (6, 3), (7, 12), (8, 13)
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (12, {'label': 'D'}), (13, {'label': 'E'})]
)
"""Network with two hybrids in separate blobs.

Properties:
- Nodes: 13, Edges: 14
- Level: 1, Vertex level: 1, Reticulation number: 2
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 2"""

LEVEL_2_DNETWORK_TWO_HYBRIDS_SAME_BLOB = DirectedPhyNetwork(
    edges=[
        (15, 10), (15, 11),  # Root
        (10, 8), (11, 8),    # First hybrid
        (8, 9), (10, 9),     # Second hybrid (receives from 8 and 10)
        (9, 12),             # Tree node
        (12, 1), (12, 2),    # Leaves
        (11, 3), (15, 4)     # Additional leaves
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'})]
)
"""Network with two hybrids in same blob.

Properties:
- Nodes: 10, Edges: 11
- Level: 2, Vertex level: 2, Reticulation number: 2
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 2"""

LEVEL_2_DNETWORK_TRIANGLE_HYBRID = DirectedPhyNetwork(
    edges=[
        (7, 5), (7, 6),   # Root to two nodes
        (5, 4), (6, 4),   # Both to hybrid (forms triangle with edge below)
        (5, 6),           # Edge making triangle
        (4, 9),           # Hybrid to tree node
        (9, 1), (9, 2),   # Tree node to leaves
        (7, 3)            # Extra leaf
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
)
"""Simple triangle reticulation (cycle of length 3).

Properties:
- Nodes: 8, Edges: 9
- Level: 2, Vertex level: 2, Reticulation number: 2
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 2"""

LEVEL_2_DNETWORK_NESTED_HYBRIDS = DirectedPhyNetwork(
    edges=[
        (20, 15), (20, 16),      # Root
        (15, 12), (16, 12),      # Hybrid 12
        (15, 21),                # Node 15 extra edge
        (12, 13),                # Hybrid to tree node
        (13, 10), (13, 11),      # Tree node splits
        (10, 8), (11, 8),        # Hybrid 8
        (8, 9),                  # Hybrid to tree node
        (9, 5), (9, 6),          # Tree node splits
        (10, 7),                 # Node 10 to create hybrid 7
        (5, 7),                  # Node 5 to create hybrid 7
        (7, 14),                 # Hybrid to tree node
        (14, 1), (14, 2),        # Leaves
        (5, 18), (6, 3), (6, 19), (16, 4), (11, 17), (21, 22), (21, 23)  # More leaves
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'}), (17, {'label': 'E'}), (18, {'label': 'F'}), (19, {'label': 'G'}), (22, {'label': 'H'}), (23, {'label': 'I'})]
)
"""Nested hybrid structure.

Properties:
- Nodes: 23, Edges: 25
- Level: 2, Vertex level: 2, Reticulation number: 3
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 3"""

LEVEL_2_DNETWORK_DIAMOND_HYBRID = DirectedPhyNetwork(
    edges=[
        (10, 7), (10, 8), (10, 11),  # Root with out-degree 3
        (7, 5), (8, 5),              # Hybrid 5 from 7,8
        (5, 12),                     # Hybrid 5 to tree node 12
        (12, 1), (12, 2),            # Tree node 12 to leaves
        (7, 6), (8, 6), (11, 6),     # Hybrid 6 from 7,8,11
        (6, 13),                     # Hybrid 6 to tree node 13
        (13, 3), (13, 4),            # Tree node 13 to leaves
        (11, 14), (11, 15)           # Node 11 to additional leaves (now out-degree 3)
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'}),
           (14, {'label': 'E'}), (15, {'label': 'F'})]
)
"""Diamond-shaped hybrid structure.

Properties:
- Nodes: 14, Edges: 16
- Level: 3, Vertex level: 2, Reticulation number: 3
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 2"""

# ============================================================================
# NON-TREE-BASED
# ============================================================================

LEVEL_5_DNETWORK_NON_TREEBASED = DirectedPhyNetwork(
    edges=[
        (100, 50),
        (100, 51),
        (50, 11),
        (50, 52),
        (51, 12),
        (51, 53),
        (52, 13),
        (52, 14),
        (53, 14),
        (53, 13),
        (11, 20),
        (11, 21),
        (12, 21),
        (12, 22),
        (13, 20),
        (14, 22),
        (20, 1),
        (21, 2),
        (22, 3),
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})],
)
"""Binary network that is NOT tree-based (violates omnian condition).

6 omnians (11, 12, 13, 14, 52, 53) with only 3 distinct children (20, 21, 22),
so |children(S)| < |S| for S = all omnians.

Properties:
- Nodes: 15, Edges: 19
- Level: 5, Vertex level: 5, Reticulation number: 5
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 5
- Is tree-based: False"""

# ============================================================================
# MULTIPLE BLOBS
# ============================================================================

LEVEL_1_DNETWORK_TWO_BLOBS = DirectedPhyNetwork(
    edges=[
        # Blob 1: hybrid node 4
        (10, 5), (10, 6), (5, 4), (6, 4), (4, 14), (14, 1), (14, 16),  # Hybrid to tree node with 2 children
        # Bridge
        (6, 11),
        # Blob 2: hybrid node 9
        (11, 7), (11, 8), (7, 9), (8, 9), (9, 15), (15, 2), (15, 17),  # Hybrid to tree node with 2 children
        # Extra leaves
        (5, 3), (7, 12), (8, 13)
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (12, {'label': 'D'}), (13, {'label': 'E'}),
           (16, {'label': 'F'}), (17, {'label': 'G'})]
)
"""Two separate blobs, each with level 1.

Properties:
- Nodes: 17, Edges: 18
- Level: 1, Vertex level: 1, Reticulation number: 2
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 2"""

LEVEL_2_DNETWORK_THREE_BLOBS = DirectedPhyNetwork(
    edges=[
        # Blob 1: level 1 hybrid 15
        (30, 20), (30, 21), (20, 15), (21, 15), (15, 28), (28, 1), (28, 32),  # Hybrid to tree with 2 children
        # Bridge to blob 2
        (21, 25),
        # Blob 2: level 2 (two hybrids in same blob)
        (25, 22), (25, 23), (22, 16), (23, 16),  # Hybrid 16
        (16, 29),  # Hybrid 16 to tree node 29
        (29, 17), (29, 34), (22, 17),  # Tree node (with 2 children) and parent to hybrid 17
        (17, 31), (31, 2), (31, 33),  # Hybrid 17 to tree with 2 children
        # Bridge to blob 3
        (23, 26),
        # Blob 3: tree (no hybrids)
        (26, 27), (27, 3), (27, 4),
        # Extra leaves
        (20, 5), (30, 6), (26, 7)  # Added leaf to node 26 (now out-degree 2)
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), 
           (4, {'label': 'D'}), (5, {'label': 'E'}), (6, {'label': 'F'}), (7, {'label': 'G'}),
           (32, {'label': 'H'}), (33, {'label': 'I'}), (34, {'label': 'J'})]
)
"""Three blobs with varying levels (1, 2, 0).

Properties:
- Nodes: 24, Edges: 26
- Level: 2, Vertex level: 2, Reticulation number: 3
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 3"""

LEVEL_1_DNETWORK_FIVE_BLOBS = DirectedPhyNetwork(
    edges=[
        # Blob 1: hybrid 40 -> tree 45 -> leaves
        (100, 50), (100, 51), (50, 40), (51, 40), (40, 45), (45, 1), (45, 11),
        # Bridge
        (51, 60),
        # Blob 2: hybrid 41 -> tree 46 -> leaves
        (60, 52), (60, 53), (52, 41), (53, 41), (41, 46), (46, 2), (46, 12),
        # Bridge
        (53, 70),
        # Blob 3: hybrid 42 -> tree 47 -> leaves
        (70, 54), (70, 55), (54, 42), (55, 42), (42, 47), (47, 3), (47, 13),
        # Bridge
        (55, 80),
        # Blob 4: hybrid 43 -> tree 48 -> leaves
        (80, 56), (80, 57), (56, 43), (57, 43), (43, 48), (48, 4), (48, 14),
        # Bridge
        (57, 90),
        # Blob 5: hybrid 44 -> tree 49 -> leaves
        (90, 58), (90, 59), (58, 44), (59, 44), (44, 49), (49, 5), (49, 15),
        # Extra leaves
        (50, 6), (52, 7), (54, 8), (56, 9), (58, 10), (59, 16)
    ],
    nodes=[(i, {'label': f'L{i}'}) for i in range(1, 17)]
)
"""Large network with 5 separate blobs.

Properties:
- Nodes: 41, Edges: 45
- Level: 1, Vertex level: 1, Reticulation number: 5
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 5"""

DTREE_BALANCED = DirectedPhyNetwork(
    edges=[
        (20, 10), (20, 11), (20, 12), (20, 13), (20, 14),
        (10, 1), (10, 6), (11, 2), (11, 7), (12, 3), (12, 8), (13, 4), (13, 9), (14, 5), (14, 15)
    ],
    nodes=[(i, {'label': f'L{i}'}) for i in range(1, 10)] + [(15, {'label': 'L10'})]
)
"""Network where most blobs are trivial 1-blobs (leaves).

Properties:
- Nodes: 16, Edges: 15
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 0"""

LEVEL_1_DNETWORK_CHAIN_OF_BLOBS = DirectedPhyNetwork(
    edges=[
        # Chain of 4 blobs - each hybrid connects to tree node first
        (50, 20), (50, 21), (20, 10), (21, 10),  # Hybrid 10
        (10, 14), (14, 1), (14, 18),  # Hybrid to tree (with 2 children)
        (18, 30), (18, 32),  # Tree continues to bridge and has 2 children
        (30, 22), (30, 23), (22, 11), (23, 11),  # Hybrid 11
        (11, 15), (15, 3), (15, 19),  # Hybrid to tree (with 2 children)
        (19, 40), (19, 33),  # Tree continues to bridge and has 2 children
        (40, 24), (40, 25), (24, 12), (25, 12),  # Hybrid 12
        (12, 16), (16, 4), (16, 29),  # Hybrid to tree (with 2 children)
        (29, 46), (29, 34),  # Tree continues to bridge and has 2 children
        (46, 27), (46, 28), (27, 13), (28, 13),  # Hybrid 13
        (13, 17), (17, 5), (17, 31),  # Hybrid to tree (with 2 children)
        # Extra leaves to fix pass-through nodes
        (20, 2), (21, 35), (22, 36), (23, 37), (24, 38), (25, 39), (27, 6), (28, 41)
    ],
    nodes=[(i, {'label': f'L{i}'}) for i in range(1, 7)] + [(i, {'label': f'X{i}'}) for i in range(35, 42)]
)
"""Chain of 4 blobs connected sequentially.

Properties:
- Nodes: 39, Edges: 42
- Level: 1, Vertex level: 1, Reticulation number: 4
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 4"""

LEVEL_1_DNETWORK_STAR_BLOBS = DirectedPhyNetwork(
    edges=[
        # Central hub node
        (100, 50), (100, 60), (100, 70), (100, 80),
        # Blob 1: hybrid 40 -> tree 44 -> leaves
        (50, 51), (50, 52), (51, 40), (52, 40), (40, 44), (44, 1), (44, 10), (52, 15),
        # Blob 2: hybrid 41 -> tree 45 -> leaves
        (60, 61), (60, 62), (61, 41), (62, 41), (41, 45), (45, 2), (45, 11), (62, 16),
        # Blob 3: hybrid 42 -> tree 46 -> leaves
        (70, 71), (70, 72), (71, 42), (72, 42), (42, 46), (46, 3), (46, 12), (72, 9),
        # Blob 4: hybrid 43 -> tree 47 -> leaves
        (80, 81), (80, 82), (81, 43), (82, 43), (43, 47), (47, 4), (47, 13), (82, 14),
        # Extra leaves from spokes
        (51, 5), (61, 6), (71, 7), (81, 8)
    ],
    nodes=[(i, {'label': f'L{i}'}) for i in range(1, 17)]
)
"""Star topology with central hub and 4 blob branches.

Properties:
- Nodes: 37, Edges: 40
- Level: 1, Vertex level: 1, Reticulation number: 4
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 4"""

# ============================================================================
# HIGH LEVEL NETWORKS
# ============================================================================

LEVEL_2_DNETWORK_SINGLE_BLOB = DirectedPhyNetwork(
    edges=[
        (30, 25), (30, 26),              # Root
        (25, 20), (25, 22), (26, 20),   # Node 25 branches (to hybrid 20 and leaf 22)
        (20, 27),                        # Hybrid 20 to tree node 27
        (27, 18), (27, 19),              # Tree node 27 splits
        (18, 15), (18, 21), (19, 15),   # Tree node 18 branches (to hybrid 15 and leaf 21)
        (15, 28),                        # Hybrid 15 to tree node 28
        (28, 12), (28, 13),              # Tree node 28 splits
        (12, 10), (13, 10),              # Hybrid 10
        (10, 29),                        # Hybrid 10 to tree node 29
        (29, 8), (29, 23), (12, 8),     # Tree node 29 branches (to hybrid 8 and leaf 23)
        (8, 31),                         # Hybrid 8 to tree node 31
        (31, 1), (31, 2),                # Tree node 31 to leaves
        (13, 3), (19, 4), (26, 5)        # More leaves
    ],
    nodes=[(i, {'label': f'L{i}'}) for i in range(1, 6)] + [(21, {'label': 'L6'}), (22, {'label': 'L7'}), (23, {'label': 'L8'})]
)
"""Single blob with level 3.

Properties:
- Nodes: 23, Edges: 26
- Level: 2, Vertex level: 2, Reticulation number: 4
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 4"""

LEVEL_2_DNETWORK_MULTIPLE_BLOBS = DirectedPhyNetwork(
    edges=[
        # Blob 1: level 2 - fixing all hybrids
        (50, 40), (50, 41), (40, 30), (41, 30),  # Hybrid 30
        (30, 60),                                  # Hybrid 30 to tree node 60
        (60, 31), (60, 32),                        # Tree node 60 splits
        (31, 20), (31, 11), (32, 20),              # Node 31 branches
        (20, 61),                                  # Hybrid 20 to tree node 61
        (61, 21), (61, 16), (31, 21),             # Tree node 61 branches, hybrid 21
        (21, 66),                                  # Hybrid 21 to tree node 66
        (66, 62), (66, 7),                        # Tree node 66 branches
        (62, 1), (62, 17),                        # Tree node 62 branches
        # Bridge
        (32, 45),
        # Blob 2: level 2 - fixing all hybrids
        (45, 42), (45, 43), (42, 33), (43, 33),  # Hybrid 33
        (33, 63),                                  # Hybrid 33 to tree node 63
        (63, 34), (63, 35),                        # Tree node 63 splits
        (41, 12), (43, 13),                        # Nodes 41, 43 branch
        (34, 22), (34, 68), (35, 22), (35, 6),    # Nodes 34, 35 branch (use 68 instead of 5)
        (22, 64),                                  # Hybrid 22 to tree node 64
        (64, 23), (64, 14), (34, 23),             # Tree node 64 branches, hybrid 23
        (23, 67),                                  # Hybrid 23 to tree node 67
        (67, 65), (67, 8),                        # Tree node 67 branches
        (65, 2), (65, 15),                        # Tree node 65 branches
        # Extra leaves
        (40, 3), (41, 9), (42, 4), (43, 10)
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'}), (6, {'label': 'F'}), (7, {'label': 'G'}), (8, {'label': 'H'}), (9, {'label': 'I'}), (10, {'label': 'J'}), (11, {'label': 'K'}), (12, {'label': 'L'}), (13, {'label': 'M'}), (14, {'label': 'N'}), (15, {'label': 'O'}), (16, {'label': 'P'}), (17, {'label': 'Q'}), (68, {'label': 'E'})]
)
"""Two blobs, each with level 2.

Properties:
- Nodes: 41, Edges: 46
- Level: 2, Vertex level: 2, Reticulation number: 6
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 6"""

LEVEL_3_DNETWORK_HIGH_VERTEX_LEVEL = DirectedPhyNetwork(
    edges=[
        (50, 40), (50, 41), (50, 42),    # Root with out-degree 3
        (40, 30), (41, 30), (41, 7), (42, 30),  # Node 41 branches
        (30, 60),                        # Hybrid 30 to tree node 60
        (60, 31), (60, 32),              # Tree node 60 splits
        (31, 20), (32, 20), (40, 20),    # Hybrid 20 with in-degree 3
        (20, 61),                        # Hybrid 20 to tree node 61
        (61, 21), (61, 22),              # Tree node 61 splits
        (21, 10), (21, 5), (22, 10), (22, 6), (31, 10),  # Nodes 21, 22 branch (to hybrid 10 and leaves)
        (10, 62),                        # Hybrid 10 to tree node 62
        (62, 1), (62, 2),                # Tree node 62 to leaves
        (32, 3), (42, 4)                 # More leaves
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'}), (5, {'label': 'E'}), (6, {'label': 'F'})]
)
"""Network with high vertex level (3 hybrids in same blob).

Properties:
- Nodes: 21, Edges: 26
- Level: 6, Vertex level: 3, Reticulation number: 6
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 3"""

LEVEL_2_DNETWORK_PARALLEL_HYBRIDS = DirectedPhyNetwork(
    edges=[
        (30, 20), (30, 21), (30, 22),    # Root
        (20, 15), (21, 15), (22, 15),    # Hybrid 15 (in-degree 3)
        (15, 50),                        # Hybrid 15 to tree node 50
        (50, 16), (50, 17),              # Tree node 50 splits
        (16, 10), (16, 5), (17, 10), (17, 6), (20, 10),  # Nodes 16, 17 branch (to hybrid 10 and leaves)
        (10, 51),                        # Hybrid 10 to tree node 51
        (51, 1), (51, 2),                # Tree node 51 to leaves
        (21, 3), (22, 4)                 # More leaves
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'}), (5, {'label': 'E'}), (6, {'label': 'F'})]
)
"""Level 2 network with parallel hybrid pattern.

Properties:
- Nodes: 16, Edges: 19
- Level: 4, Vertex level: 2, Reticulation number: 4
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 2"""

# ============================================================================
# LARGE NETWORKS
# ============================================================================

def _build_large_tree_with_hybrids(num_leaves: int, num_hybrids: int) -> tuple[list, list]:
    """Helper to build large network with specified hybrids.
    
    Creates a base tree and adds hybrid nodes by:
    1. Identifying internal nodes to become hybrid parents
    2. Creating new hybrid nodes
    3. Ensuring each hybrid has in-degree >= 2 and out-degree = 1
    """
    # Start with balanced tree
    edges, nodes = _build_balanced_binary_tree(num_leaves)
    
    # Find max node ID
    max_node = max([n for n, _ in nodes] + [u for u, v in edges] + [v for u, v in edges])
    
    # Find internal nodes to use as hybrid parents
    all_nodes = set([u for u, v in edges] + [v for u, v in edges])
    leaf_nodes = set([n for n, _ in nodes])
    internal_nodes = list(all_nodes - leaf_nodes)
    
    # Add hybrids by creating new hybrid nodes
    hybrid_node_start = max_node + 1
    tree_node_start = hybrid_node_start + num_hybrids
    
    for i in range(min(num_hybrids, len(internal_nodes) // 3)):
        # Get two parent nodes for this hybrid
        parent1_idx = i * 3
        parent2_idx = i * 3 + 1
        if parent2_idx < len(internal_nodes):
            parent1 = internal_nodes[parent1_idx]
            parent2 = internal_nodes[parent2_idx]
            hybrid = hybrid_node_start + i
            tree_node = tree_node_start + i
            
            # Add edges: parent1 -> hybrid, parent2 -> hybrid
            edges.append((parent1, hybrid))
            edges.append((parent2, hybrid))
            # Hybrid -> tree_node (out-degree 1)
            edges.append((hybrid, tree_node))
            # Tree_node -> two leaves (create new leaves)
            leaf1 = tree_node_start + num_hybrids + i * 2
            leaf2 = tree_node_start + num_hybrids + i * 2 + 1
            edges.append((tree_node, leaf1))
            edges.append((tree_node, leaf2))
            nodes.append((leaf1, {'label': f'HL{i*2}'}))
            nodes.append((leaf2, {'label': f'HL{i*2+1}'}))
    
    return edges, nodes

_large_many_hybrids_edges, _large_many_hybrids_nodes = _build_large_tree_with_hybrids(25, 5)
LEVEL_3_DNETWORK_LARGE_MANY_HYBRIDS = DirectedPhyNetwork(
    edges=_large_many_hybrids_edges,
    nodes=_large_many_hybrids_nodes
)
"""Large network with 25+ leaves and 5 hybrid nodes.

Properties:
- Nodes: 75, Edges: 79
- Level: 3, Vertex level: 3, Reticulation number: 5
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 5"""

_large_few_hybrids_edges, _large_few_hybrids_nodes = _build_large_tree_with_hybrids(30, 2)
LEVEL_1_DNETWORK_LARGE_FEW_HYBRIDS = DirectedPhyNetwork(
    edges=_large_few_hybrids_edges,
    nodes=_large_few_hybrids_nodes
)
"""Large network with 30+ leaves and only 2 hybrid nodes.

Properties:
- Nodes: 69, Edges: 70
- Level: 1, Vertex level: 1, Reticulation number: 2
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 2"""

# Build large multi-blob network
_blob_edges = []
_blob_nodes = []
for blob_idx in range(4):
    base = blob_idx * 20
    # Each blob has a hybrid with proper structure: hybrid -> tree_node -> leaf
    _blob_edges.extend([
        (base + 15, base + 10), (base + 15, base + 11),  # Parents to hybrid 8
        (base + 10, base + 8), (base + 11, base + 8),    # Hybrid 8
        (base + 8, base + 12),                           # Hybrid to tree node 12
        (base + 12, base + 1), (base + 12, base + 4),  # Tree node 12 branches (to 2 leaves)
        (base + 10, base + 2), (base + 11, base + 3)     # Parent nodes to extra leaves
    ])
    _blob_nodes.extend([
        (base + 1, {'label': f'B{blob_idx}L1'}),
        (base + 2, {'label': f'B{blob_idx}L2'}),
        (base + 3, {'label': f'B{blob_idx}L3'}),
        (base + 4, {'label': f'B{blob_idx}L4'})
    ])
    # Connect blobs with bridges
    if blob_idx > 0:
        prev_base = (blob_idx - 1) * 20
        _blob_edges.append((prev_base + 15, base + 15))

LEVEL_1_DNETWORK_LARGE_MULTI_BLOB = DirectedPhyNetwork(edges=_blob_edges, nodes=_blob_nodes)
"""Large network with 40+ nodes and 4 blobs.

Properties:
- Nodes: 36, Edges: 39
- Level: 1, Vertex level: 1, Reticulation number: 4
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 4"""

_very_large_edges, _very_large_nodes = _build_balanced_binary_tree(50)
DTREE_VERY_LARGE = DirectedPhyNetwork(edges=_very_large_edges, nodes=_very_large_nodes)
"""Very large tree with 50+ leaves.

Properties:
- Nodes: 105, Edges: 104
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 0"""

# Large network built by gluing smaller networks
LEVEL_1_DNETWORK_LARGE_GLUED = DirectedPhyNetwork(
    edges=[
        # Section 1: Small hybrid network - hybrid 90 -> tree 96 -> leaves
        (200, 100), (200, 101), (100, 90), (101, 90), (90, 96), (96, 1), (96, 18), (100, 2),
        # Bridge
        (101, 110),
        # Section 2: Binary tree
        (110, 102), (110, 103), (102, 91), (102, 92), (103, 93), (103, 94),
        (91, 3), (91, 4), (92, 5), (92, 6), (93, 7), (93, 8), (94, 9), (94, 10),
        # Bridge
        (103, 120),
        # Section 3: Another hybrid - hybrid 95 -> tree 97 -> leaves
        (120, 104), (120, 105), (104, 95), (105, 95), (95, 97), (97, 11), (97, 19), (104, 12),
        # Bridge
        (105, 130),
        # Section 4: Non-binary tree
        (130, 106), (130, 107), (130, 108), (106, 13), (106, 14), 
        (107, 15), (107, 20), (108, 16), (108, 17)
    ],
    nodes=[(i, {'label': f'L{i}'}) for i in range(1, 21)]
)
"""Large network built by gluing smaller networks together.

Properties:
- Nodes: 41, Edges: 42
- Level: 1, Vertex level: 1, Reticulation number: 2
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 2"""

LEVEL_1_DNETWORK_LARGE_DEEP_HYBRID_CHAIN = DirectedPhyNetwork(
    edges=[
        # Chain of 6 hybrid nodes - each hybrid -> tree node
        (100, 90), (100, 91), (90, 80), (91, 80), (91, 24),  # Node 91 branches
        (80, 101), (101, 81), (101, 82), (101, 9), # Hybrid 80 -> tree 101 -> nodes (add leaf)
        (81, 70), (81, 22), (82, 70), (82, 23),    # Nodes 81, 82 branch
        (70, 102), (102, 71), (102, 72), (102, 10), # Hybrid 70 -> tree 102 -> nodes (add leaf)
        (71, 60), (71, 20), (72, 60), (72, 21),    # Nodes 71, 72 branch
        (60, 103), (103, 61), (103, 62), (103, 11), # Hybrid 60 -> tree 103 -> nodes (add leaf)
        (61, 50), (61, 18), (62, 50), (62, 19),    # Nodes 61, 62 branch
        (50, 104), (104, 51), (104, 52), (104, 12), # Hybrid 50 -> tree 104 -> nodes (add leaf)
        (51, 40), (51, 16), (52, 40), (52, 17),    # Nodes 51, 52 branch
        (40, 105), (105, 41), (105, 42), (105, 13), # Hybrid 40 -> tree 105 -> nodes (add leaf)
        (41, 30), (41, 14), (42, 30), (42, 15),    # Nodes 41, 42 branch
        (30, 106), (106, 1), (106, 2),             # Hybrid 30 -> tree 106 -> leaves
        (90, 3), (81, 4), (71, 5),
        (61, 6), (51, 7), (41, 8)
    ],
    nodes=[(i, {'label': f'L{i}'}) for i in range(1, 25)]
)
"""Large network with chain of 6 hybrid nodes.

Properties:
- Nodes: 49, Edges: 54
- Level: 1, Vertex level: 1, Reticulation number: 6
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 6"""

# ============================================================================
# PARALLEL EDGES
# ============================================================================

LEVEL_1_DNETWORK_PARALLEL_EDGES = DirectedPhyNetwork(
    edges=[
        (10, 5, 0), (10, 5, 1),  # Parallel edges to hybrid 5
        (10, 6, 0),               # Edge to node 6 (with explicit key)
        (5, 11),                  # Hybrid 5 to tree node 11
        (11, 1), (11, 2),         # Tree node 11 branches
        (6, 3), (6, 4)
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'})]
)
"""Tree with parallel edges between root and child.

Properties:
- Nodes: 8, Edges: 8
- Level: 1, Vertex level: 1, Reticulation number: 1
- Is tree: False, Is binary: False
- Has parallel edges: True
- Number of hybrid nodes: 1"""

LEVEL_1_DNETWORK_PARALLEL_EDGES_HYBRID = DirectedPhyNetwork(
    edges=[
        (15, 10), (15, 11),
        (10, 8, 0), (10, 8, 1),  # Parallel edges to hybrid
        (11, 8, 0), (11, 8, 1),  # More parallel edges to same hybrid
        (8, 9),                  # Hybrid to tree node
        (9, 5), (9, 6),          # Tree node branches
        (5, 1), (5, 2),
        (10, 3), (11, 4)
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'}), (6, {'label': 'E'})]
)
"""Hybrid network with parallel edges to hybrid node.

Properties:
- Nodes: 11, Edges: 13
- Level: 3, Vertex level: 1, Reticulation number: 3
- Is tree: False, Is binary: False
- Has parallel edges: True
- Number of hybrid nodes: 1"""

LEVEL_2_DNETWORK_MANY_PARALLEL_EDGES = DirectedPhyNetwork(
    edges=[
        (20, 15, 0), (20, 15, 1), (20, 15, 2),  # Triple parallel to hybrid 15
        (20, 16),
        (15, 25),                                 # Hybrid 15 to tree node 25
        (25, 10, 0), (25, 10, 1),                # Tree node 25 to hybrid 10 (double parallel)
        (16, 10, 0), (16, 10, 1),                # Double parallel to hybrid
        (10, 11),                                 # Hybrid 10 to tree node
        (11, 5), (11, 7),                         # Tree node branches
        (5, 1), (5, 2),
        (25, 3), (16, 4)                         # Tree node 25 and node 16 to leaves
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'}), (7, {'label': 'E'})]
)
"""Network with multiple sets of parallel edges.

Properties:
- Nodes: 12, Edges: 16
- Level: 5, Vertex level: 2, Reticulation number: 5
- Is tree: False, Is binary: False
- Has parallel edges: True
- Number of hybrid nodes: 2"""

LEVEL_1_DNETWORK_PARALLEL_IN_BLOB = DirectedPhyNetwork(
    edges=[
        (20, 15), (20, 16),
        (15, 12, 0), (15, 12, 1),  # Parallel within blob
        (16, 12, 0), (16, 12, 1),  # More parallel to hybrid
        (12, 13),                  # Hybrid to tree node
        (13, 10, 0), (13, 10, 1),  # Parallel from tree node to hybrid 10
        (10, 11),                  # Hybrid 10 to tree node 11
        (11, 1), (11, 2),          # Tree node 11 to leaves
        (15, 3), (16, 4)
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'})]
)
"""Network with parallel edges concentrated in one blob.

Properties:
- Nodes: 11, Edges: 14
- Level: 3, Vertex level: 1, Reticulation number: 4
- Is tree: False, Is binary: False
- Has parallel edges: True
- Number of hybrid nodes: 2"""

# ============================================================================
# SPECIAL CASES
# ============================================================================

DTREE_SIMPLE = DirectedPhyNetwork(
    edges=[
        # Very deep chain: 3 levels (simplified to avoid pass-through nodes)
        (30, 29), (29, 15), (29, 3),
        (15, 1), (15, 2)
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
)
"""Very deep network (3 levels deep).

Properties:
- Nodes: 6, Edges: 5
- Level: 0, Vertex level: 0, Reticulation number: 0
- Is tree: True, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 0"""

# DTREE_WIDE removed as requested

LEVEL_3_DNETWORK_CHAIN_HYBRIDS = DirectedPhyNetwork(
    edges=[
        (50, 40), (50, 41),
        (40, 30), (41, 30),  # Hybrid 30
        (30, 60),            # Hybrid 30 to tree node 60
        (60, 31), (60, 4),   # Tree node 60 branches (not pass-through)
        (31, 20), (40, 20),  # Hybrid 20 (31 and 40)
        (20, 61),            # Hybrid 20 to tree node 61
        (61, 21), (61, 5),   # Tree node 61 branches (not pass-through)
        (21, 10), (21, 7), (31, 10),  # Node 21 branches (to hybrid 10 and leaf 7)
        (10, 62),            # Hybrid 10 to tree node 62
        (62, 1), (62, 2),    # Tree node 62 to leaves
        (41, 3)
    ],
    nodes=[(i, {'label': f'L{i}'}) for i in [1, 2, 3, 4, 5, 7]]  # Leaves: 1, 2, 3, 4, 5, 7 (6 not used)
)
"""Chain of sequential hybrid nodes.

Properties:
- Nodes: 17, Edges: 19
- Level: 3, Vertex level: 3, Reticulation number: 3
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 3"""

LEVEL_1_DNETWORK_MIXED_TOPOLOGY = DirectedPhyNetwork(
    edges=[
        # Mix of deep, wide, and hybrid structures
        (100, 80), (100, 81), (100, 82), (100, 83),  # Wide root
        # Deep branch 1 (simplified to avoid pass-through)
        (80, 50), (50, 1), (50, 11), (80, 10),  # Node 50 branches
        # Hybrid branch 2 - hybrid 71 and hybrid 40 both need tree nodes
        (81, 71), (81, 12), (82, 71), (82, 13),  # Nodes 81, 82 branch
        (71, 90),            # Hybrid 71 to tree node 90
        (90, 51), (90, 52),  # Tree node 90 splits
        (51, 40), (52, 40), (52, 14),  # Node 52 branches
        (40, 91),            # Hybrid 40 to tree node 91
        (91, 2), (91, 15),   # Tree node 91 branches
        (51, 3),
        # Wide branch 3
        (83, 72), (83, 73), (83, 74), (72, 4), (72, 7), (72, 25), (73, 5), (73, 8), (74, 6), (74, 9)
    ],
    nodes=[(i, {'label': f'L{i}'}) for i in range(1, 16)] + [(25, {'label': 'L16'})]
)
"""Complex topology mixing deep, wide, and hybrid structures.

Properties:
- Nodes: 31, Edges: 32
- Level: 1, Vertex level: 1, Reticulation number: 2
- Is tree: False, Is binary: False
- Has parallel edges: False
- Number of hybrid nodes: 2"""

# ============================================================================
# NON-LSA NETWORKS
# ============================================================================

LEVEL_1_DNETWORK_NON_LSA = DirectedPhyNetwork(
    edges=[
        (10, 5),
        (10, 6),
        (5, 4),
        (5, 3),
        (6, 4),
        (6, 7),
        (4, 8),
        (8, 1),
        (8, 2),
        (132, 133),
        (133, 134),
        {'u': 133, 'v': 134, 'key': 1},
        (134, 135),
        (135, 10),
        (135, 136)
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'}), (136, {'label': 'E'})]
)
"""Non-LSA network based on LEVEL_1_DNETWORK_SINGLE_HYBRID.

Created by adding a new root (132) with structure: 132 -> 133 -> 134 (parallel edges, twice), 134 -> 135 -> original root (10).
This creates a non-LSA network with a 2-blob structure above the original network.

Properties:
- Nodes: 14, Edges: 15
- Level: 1, Vertex level: 1, Reticulation number: 2
- Is tree: False, Is binary: False
- Has parallel edges: True
- Number of hybrid nodes: 2"""

LEVEL_2_DNETWORK_NON_LSA = DirectedPhyNetwork(
    edges=[
        (30, 25),
        (30, 26),
        (25, 20),
        (25, 22),
        (26, 20),
        (26, 5),
        (20, 27),
        (27, 18),
        (27, 19),
        (18, 15),
        (18, 21),
        (19, 15),
        (19, 4),
        (15, 28),
        (28, 12),
        (28, 13),
        (12, 10),
        (12, 8),
        (13, 10),
        (13, 3),
        (10, 29),
        (8, 31),
        (29, 8),
        (29, 23),
        (31, 1),
        (31, 2),
        (141, 142),
        (142, 143),
        {'u': 142, 'v': 143, 'key': 1},
        (142, 30),
        (143, 144)
    ],
    nodes=[(1, {'label': 'L1'}), (2, {'label': 'L2'}), (3, {'label': 'L3'}), (4, {'label': 'L4'}), (5, {'label': 'L5'}), (21, {'label': 'L6'}), (22, {'label': 'L7'}), (23, {'label': 'L8'}), (144, {'label': 'L144'})]
)
"""Non-LSA network based on LEVEL_2_DNETWORK_SINGLE_BLOB.

Created by adding a new root (141) with structure: 141 -> 142 -> 143 (parallel edges), 142 -> old_root (30).
This creates a non-LSA network with a 2-blob structure above the original network.

Properties:
- Nodes: 27, Edges: 31
- Level: 2, Vertex level: 2, Reticulation number: 5
- Is tree: False, Is binary: False
- Has parallel edges: True
- Number of hybrid nodes: 5"""

LEVEL_2_DNETWORK_NON_LSA_NESTED = DirectedPhyNetwork(
    edges=[
        (20, 15),
        (20, 16),
        (15, 12),
        (15, 21),
        (16, 12),
        (16, 4),
        (12, 13),
        (21, 22),
        (21, 23),
        (13, 10),
        (13, 11),
        (10, 8),
        (10, 7),
        (11, 8),
        (11, 17),
        (8, 9),
        (7, 14),
        (9, 5),
        (9, 6),
        (5, 7),
        (5, 18),
        (6, 3),
        (6, 19),
        (14, 1),
        (14, 2),
        (151, 152),
        (152, 153),
        {'u': 152, 'v': 153, 'key': 1},
        (153, 154),
        (154, 20),
        (154, 155)
    ],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'}), (17, {'label': 'E'}), (18, {'label': 'F'}), (19, {'label': 'G'}), (22, {'label': 'H'}), (23, {'label': 'I'}), (155, {'label': 'J'})]
)
"""Non-LSA network based on LEVEL_2_DNETWORK_NESTED_HYBRIDS.

Created by adding a new root (151) with structure: 151 -> 152 -> 153 (parallel edges, twice), 153 -> 154 -> original root (20).
This creates a non-LSA network with a 2-blob structure above the original network.

Properties:
- Nodes: 28, Edges: 31
- Level: 2, Vertex level: 2, Reticulation number: 4
- Is tree: False, Is binary: False
- Has parallel edges: True
- Number of hybrid nodes: 4"""

LEVEL_3_DNETWORK = DirectedPhyNetwork(
    edges=[
        (30, 25),
        (30, 26),
        (25, 20),
        (25, 22),
        (26, 20),
        (26, 5),
        (20, 27),
        (27, 18),
        (27, 19),
        (18, 15),
        (18, 21),
        (19, 15),
        (19, 4),
        (15, 28),
        (28, 232),
        (28, 233),
        (232, 12),
        (232, 234),
        (12, 10),
        (12, 8),
        (233, 13),
        (233, 234),
        (13, 10),
        (13, 3),
        (10, 29),
        (8, 31),
        (29, 8),
        (29, 23),
        (31, 1),
        (31, 2),
        (234, 235),
        (235, 236),
        (235, 237)
    ],
    nodes=[(1, {'label': 'L1'}), (2, {'label': 'L2'}), (3, {'label': 'L3'}), (4, {'label': 'L4'}), (5, {'label': 'L5'}), (21, {'label': 'L6'}), (22, {'label': 'L7'}), (23, {'label': 'L8'}), (236, {'label': 'L236'}), (237, {'label': 'L237'})]
)
"""Level-3 network created from LEVEL_2_DNETWORK_SINGLE_BLOB.

Created by subdividing two edges within the largest blob and adding a new hybrid node.
The new hybrid connects the subdivided edges, increasing the blob level from 2 to 3.

Properties:
- Nodes: 29, Edges: 33
- Level: 3, Vertex level: 3, Reticulation number: 5
- Is tree: False, Is binary: True
- Has parallel edges: False
- Number of hybrid nodes: 5"""

# ============================================================================
# METADATA REGISTRY
# ============================================================================

NETWORK_METADATA: dict[str, dict[str, Any]] = {
    # Trees
    'DTREE_EMPTY': {
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
    'DTREE_SINGLE_NODE': {
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
    'DTREE_SMALL_BINARY': {
        'category': 'trees',
        'nodes': 5,
        'edges': 4,
        'level': 0,
        'vertex_level': 0,
        'reticulation_number': 0,
        'is_tree': True,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 0,
    },
    'DTREE_MEDIUM_BINARY': {
        'category': 'trees',
        'nodes': 23,
        'edges': 22,
        'level': 0,
        'vertex_level': 0,
        'reticulation_number': 0,
        'is_tree': True,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 0,
    },
    'DTREE_LARGE_BINARY': {
        'category': 'trees',
        'nodes': 43,
        'edges': 42,
        'level': 0,
        'vertex_level': 0,
        'reticulation_number': 0,
        'is_tree': True,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 0,
    },
    'DTREE_NON_BINARY_SMALL': {
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
    'DTREE_NON_BINARY_LARGE': {
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
    # Simple Hybrids
    'LEVEL_1_DNETWORK_SINGLE_HYBRID': {
        'category': 'simple_hybrids',
        'nodes': 9,
        'edges': 9,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 1,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 1,
    },
    'LEVEL_1_DNETWORK_SINGLE_HYBRID_BINARY': {
        'category': 'simple_hybrids',
        'nodes': 9,
        'edges': 9,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 1,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 1,
    },
    'LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE': {
        'category': 'simple_hybrids',
        'nodes': 13,
        'edges': 14,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 2,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 2,
    },
    'LEVEL_2_DNETWORK_TWO_HYBRIDS_SAME_BLOB': {
        'category': 'simple_hybrids',
        'nodes': 10,
        'edges': 11,
        'level': 2,
        'vertex_level': 2,
        'reticulation_number': 2,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 2,
    },
    'LEVEL_2_DNETWORK_TRIANGLE_HYBRID': {
        'category': 'simple_hybrids',
        'nodes': 8,
        'edges': 9,
        'level': 2,
        'vertex_level': 2,
        'reticulation_number': 2,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 2,
    },
    'LEVEL_2_DNETWORK_NESTED_HYBRIDS': {
        'category': 'simple_hybrids',
        'nodes': 23,
        'edges': 25,
        'level': 2,
        'vertex_level': 2,
        'reticulation_number': 3,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 3,
    },
    'LEVEL_2_DNETWORK_DIAMOND_HYBRID': {
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
    'LEVEL_5_DNETWORK_NON_TREEBASED': {
        'category': 'non_treebased',
        'nodes': 15,
        'edges': 19,
        'level': 5,
        'vertex_level': 5,
        'reticulation_number': 5,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 5,
    },
    # Multiple Blobs
    'LEVEL_1_DNETWORK_TWO_BLOBS': {
        'category': 'multiple_blobs',
        'nodes': 17,
        'edges': 18,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 2,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 2,
    },
    'LEVEL_2_DNETWORK_THREE_BLOBS': {
        'category': 'multiple_blobs',
        'nodes': 24,
        'edges': 26,
        'level': 2,
        'vertex_level': 2,
        'reticulation_number': 3,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 3,
    },
    'LEVEL_1_DNETWORK_FIVE_BLOBS': {
        'category': 'multiple_blobs',
        'nodes': 41,
        'edges': 45,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 5,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 5,
    },
    'DTREE_BALANCED': {
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
    'LEVEL_1_DNETWORK_CHAIN_OF_BLOBS': {
        'category': 'multiple_blobs',
        'nodes': 39,
        'edges': 42,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 4,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 4,
    },
    'LEVEL_1_DNETWORK_STAR_BLOBS': {
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
    # High Level
    'LEVEL_2_DNETWORK_SINGLE_BLOB': {
        'category': 'high_level',
        'nodes': 23,
        'edges': 26,
        'level': 2,
        'vertex_level': 2,
        'reticulation_number': 4,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 4,
    },
    'LEVEL_2_DNETWORK_MULTIPLE_BLOBS': {
        'category': 'high_level',
        'nodes': 41,
        'edges': 46,
        'level': 2,
        'vertex_level': 2,
        'reticulation_number': 6,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 6,
    },
    'LEVEL_3_DNETWORK_HIGH_VERTEX_LEVEL': {
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
    'LEVEL_2_DNETWORK_PARALLEL_HYBRIDS': {
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
    # Large Networks
    'LEVEL_3_DNETWORK_LARGE_MANY_HYBRIDS': {
        'category': 'large_networks',
        'nodes': 75,
        'edges': 79,
        'level': 3,
        'vertex_level': 3,
        'reticulation_number': 5,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 5,
    },
    'LEVEL_1_DNETWORK_LARGE_FEW_HYBRIDS': {
        'category': 'large_networks',
        'nodes': 69,
        'edges': 70,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 2,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 2,
    },
    'LEVEL_1_DNETWORK_LARGE_MULTI_BLOB': {
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
    'DTREE_VERY_LARGE': {
        'category': 'large_networks',
        'nodes': 105,
        'edges': 104,
        'level': 0,
        'vertex_level': 0,
        'reticulation_number': 0,
        'is_tree': True,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 0,
    },
    'LEVEL_1_DNETWORK_LARGE_GLUED': {
        'category': 'large_networks',
        'nodes': 41,
        'edges': 42,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 2,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 2,
    },
    'LEVEL_1_DNETWORK_LARGE_DEEP_HYBRID_CHAIN': {
        'category': 'large_networks',
        'nodes': 49,
        'edges': 54,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 6,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 6,
    },
    # Parallel Edges
    'LEVEL_1_DNETWORK_PARALLEL_EDGES': {
        'category': 'parallel_edges',
        'nodes': 8,
        'edges': 8,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 1,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': True,
        'num_hybrids': 1,
    },
    'LEVEL_1_DNETWORK_PARALLEL_EDGES_HYBRID': {
        'category': 'parallel_edges',
        'nodes': 11,
        'edges': 13,
        'level': 3,
        'vertex_level': 1,
        'reticulation_number': 3,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': True,
        'num_hybrids': 1,
    },
    'LEVEL_2_DNETWORK_MANY_PARALLEL_EDGES': {
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
    'LEVEL_1_DNETWORK_PARALLEL_IN_BLOB': {
        'category': 'parallel_edges',
        'nodes': 11,
        'edges': 14,
        'level': 3,
        'vertex_level': 1,
        'reticulation_number': 4,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': True,
        'num_hybrids': 2,
    },
    # Special Cases
    'DTREE_SIMPLE': {
        'category': 'special_cases',
        'nodes': 6,
        'edges': 5,
        'level': 0,
        'vertex_level': 0,
        'reticulation_number': 0,
        'is_tree': True,
        'is_binary': False,
        'has_parallel_edges': False,
        'num_hybrids': 0,
    },
    'LEVEL_3_DNETWORK_CHAIN_HYBRIDS': {
        'category': 'special_cases',
        'nodes': 17,
        'edges': 19,
        'level': 3,
        'vertex_level': 3,
        'reticulation_number': 3,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 3,
    },
    'LEVEL_1_DNETWORK_MIXED_TOPOLOGY': {
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
    # Non-LSA Networks
    'LEVEL_1_DNETWORK_NON_LSA': {
        'category': 'non_lsa',
        'nodes': 14,
        'edges': 15,
        'level': 1,
        'vertex_level': 1,
        'reticulation_number': 2,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': True,
        'num_hybrids': 2,
    },
    'LEVEL_2_DNETWORK_NON_LSA': {
        'category': 'non_lsa',
        'nodes': 27,
        'edges': 31,
        'level': 2,
        'vertex_level': 2,
        'reticulation_number': 5,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': True,
        'num_hybrids': 5,
    },
    'LEVEL_2_DNETWORK_NON_LSA_NESTED': {
        'category': 'non_lsa',
        'nodes': 28,
        'edges': 31,
        'level': 2,
        'vertex_level': 2,
        'reticulation_number': 4,
        'is_tree': False,
        'is_binary': False,
        'has_parallel_edges': True,
        'num_hybrids': 4,
    },
    'LEVEL_3_DNETWORK': {
        'category': 'high_level',
        'nodes': 29,
        'edges': 33,
        'level': 3,
        'vertex_level': 3,
        'reticulation_number': 5,
        'is_tree': False,
        'is_binary': True,
        'has_parallel_edges': False,
        'num_hybrids': 5,
    },
}


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_networks_by_category(category: str) -> list[DirectedPhyNetwork]:
    """
    Get all networks in a specific category.
    
    Parameters
    ----------
    category : str
        Category name (e.g., 'trees', 'simple_hybrids', 'high_level').
    
    Returns
    -------
    list[DirectedPhyNetwork]
        List of networks in the specified category.
    
    Examples
    --------
    >>> trees = get_networks_by_category('trees')
    >>> len(trees) >= 5
    True
    """
    return [globals()[k] for k, v in NETWORK_METADATA.items() 
            if v['category'] == category]


def get_networks_with_property(**kwargs: Any) -> list[DirectedPhyNetwork]:
    """
    Get networks matching all specified properties.
    
    Parameters
    ----------
    **kwargs : Any
        Property-value pairs to match (e.g., level=2, is_tree=True).
    
    Returns
    -------
    list[DirectedPhyNetwork]
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


def get_all_networks() -> list[DirectedPhyNetwork]:
    """
    Get all fixture networks.
    
    Returns
    -------
    list[DirectedPhyNetwork]
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
        Name of the network (e.g., 'DTREE_SMALL_BINARY').
    
    Returns
    -------
    dict[str, Any]
        Metadata dictionary for the network.
    
    Examples
    --------
    >>> meta = get_network_metadata('DTREE_SMALL_BINARY')
    >>> meta['category']
    'trees'
    """
    return NETWORK_METADATA.get(network_name, {})

