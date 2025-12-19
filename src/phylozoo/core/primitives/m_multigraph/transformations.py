"""
Graph transformations module.

This module provides functions to transform MixedMultiGraph instances
(e.g., identify nodes, orient edges, suppress degree-2 nodes, etc.).
"""

from typing import TYPE_CHECKING, Any, TypeVar
from collections import deque

import networkx as nx

if TYPE_CHECKING:
    from ..d_multigraph import DirectedMultiGraph
    from . import MixedMultiGraph
else:
    from . import MixedMultiGraph
    from ..d_multigraph import DirectedMultiGraph

T = TypeVar('T')


def identify_two_nodes(graph: 'MixedMultiGraph', u: T, v: T) -> None:
    """
    Identify two nodes u and v by keeping node u.
    
    Modifies the graph in place. All edges incident to v are moved to u,
    and v is removed.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to modify.
    u : T
        Node to keep.
    v : T
        Node to identify with u (will be removed).
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> G.add_directed_edge(1, 2)
    0
    >>> G.add_undirected_edge(2, 3)
    0
    >>> identify_two_nodes(G, 1, 2)
    >>> list(G.nodes())
    [1, 3]
    """
    # Use NetworkX's contracted_nodes on each graph
    nx.contracted_nodes(graph._undirected, u, v, self_loops=False, copy=False)
    nx.contracted_nodes(graph._directed, u, v, self_loops=False, copy=False)
    nx.contracted_nodes(graph._combined, u, v, self_loops=False, copy=False)
    
    # Clean up any self-loops that might have been created
    # (contracted_nodes with self_loops=False should handle this, but be safe)
    if graph._undirected.has_edge(u, u):
        for k in list(graph._undirected[u][u].keys()):
            graph._undirected.remove_edge(u, u, key=k)
            graph._combined.remove_edge(u, u, key=k)
    if graph._directed.has_edge(u, u):
        for k in list(graph._directed[u][u].keys()):
            graph._directed.remove_edge(u, u, key=k)
            graph._combined.remove_edge(u, u, key=k)


def identify_node_set(graph: 'MixedMultiGraph', nodes: list[T] | set[T]) -> None:
    """
    Identify all nodes in the set by keeping the first node.
    
    Modifies the graph in place.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to modify.
    nodes : list[T] | set[T]
        Iterable of nodes to identify. The first node will be kept.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> G.add_undirected_edge(1, 2)
    0
    >>> G.add_undirected_edge(2, 3)
    0
    >>> identify_node_set(G, [1, 2, 3])
    >>> len(G.nodes()) <= 3
    True
    """
    nodes_list = list(nodes)
    if len(nodes_list) < 2:
        return
    
    for i in range(1, len(nodes_list)):
        identify_two_nodes(graph, nodes_list[0], nodes_list[i])


def orient_away_from_vertex(graph: 'MixedMultiGraph', root: T) -> 'DirectedMultiGraph':
    """
    Orient all edges in a mixed multigraph away from a root vertex using BFS.
    
    This function performs a BFS from the given root vertex, orienting all undirected edges
    away from the current vertex. When a directed edge (u, v) is encountered, v is stored
    in a separate queue. After processing all undirected edges in the BFS, the function
    continues BFS from vertices in the directed queue, orienting undirected edges away
    from those vertices. Directed edges cannot be reoriented.
    
    The result is a DirectedMultiGraph where all edges are directed away from the root.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The mixed multigraph to orient.
    root : T
        The root vertex from which to orient all edges.
    
    Returns
    -------
    DirectedMultiGraph
        A new directed multigraph with all edges oriented away from the root.
    
    Raises
    ------
    ValueError
        If the root vertex is not in the graph, or if the graph contains cycles that
        prevent a valid orientation.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> G.add_undirected_edge(1, 2)
    0
    >>> G.add_undirected_edge(2, 3)
    0
    >>> dm = orient_away_from_vertex(G, 1)
    >>> list(dm.edges())
    [(1, 2), (2, 3)]
    """
    if root not in graph.nodes():
        raise ValueError(f"Root vertex {root} not found in graph")
    
    # Create a new DirectedMultiGraph
    from ..d_multigraph import DirectedMultiGraph
    dm = DirectedMultiGraph()
    
    # Add all nodes with attributes
    for node in graph.nodes():
        # Get node attributes from either graph (they should be the same)
        node_attrs = {}
        if node in graph._undirected:
            node_attrs.update(graph._undirected.nodes[node])
        if node in graph._directed:
            node_attrs.update(graph._directed.nodes[node])
        dm.add_node(node, **node_attrs)
    
    # BFS queues: one for undirected exploration, one for directed exploration
    undirected_queue = deque([root])
    directed_queue = deque()
    visited = {root}
    processed_undirected_edges = set()  # Track processed undirected edges (u, v, key)
    
    # Process undirected edges first
    while undirected_queue or directed_queue:
        # Process undirected queue first
        while undirected_queue:
            current = undirected_queue.popleft()
            
            # Process all undirected edges incident to current
            for u, v, key, data in graph._undirected.edges(current, keys=True, data=True):
                # Normalize edge representation to avoid processing twice
                edge_key = (min(u, v), max(u, v), key)
                if edge_key in processed_undirected_edges:
                    continue
                processed_undirected_edges.add(edge_key)
                
                neighbor = v if u == current else u
                
                if neighbor not in visited:
                    visited.add(neighbor)
                    # Orient edge away from current: current -> neighbor
                    dm.add_edge(current, neighbor, key=key, **data)
                    # Add neighbor to undirected queue to continue BFS
                    undirected_queue.append(neighbor)
                else:
                    # Neighbor already visited - orient away from current
                    if not dm.has_edge(current, neighbor, key=key):
                        dm.add_edge(current, neighbor, key=key, **data)
            
            # Process all directed edges incident to current
            # Outgoing edges: current -> v
            for u, v, key, data in graph._directed.edges(current, keys=True, data=True):
                if u == current:  # Outgoing edge
                    if v not in visited:
                        visited.add(v)
                        dm.add_edge(current, v, key=key, **data)
                        # Add to directed queue (not undirected queue)
                        directed_queue.append(v)
                    else:
                        # Already visited, but add edge if not present
                        if not dm.has_edge(current, v, key=key):
                            dm.add_edge(current, v, key=key, **data)
            
            # Incoming edges: u -> current
            for u, v, key, data in graph._directed.edges(keys=True, data=True):
                if v == current:  # Incoming edge
                    # Check if this edge points towards a vertex already visited via undirected BFS
                    # This would create an invalid orientation (edge pointing towards BFS path)
                    if current in visited and u not in visited:
                        # Check if current was reached via undirected edges
                        current_reached_via_undirected = any(
                            (min(current, n), max(current, n), k) in processed_undirected_edges
                            for n in visited
                            for k in range(graph._undirected.number_of_edges(current, n) if graph._undirected.has_edge(current, n) else 0)
                        )
                        if current_reached_via_undirected:
                            raise ValueError(
                                f"Directed edge ({u}, {current}, key={key}) points towards vertex "
                                f"{current} which is already in the BFS path from root {root}."
                            )
                    if u not in visited:
                        visited.add(u)
                        dm.add_edge(u, current, key=key, **data)
                        # Add to directed queue
                        directed_queue.append(u)
                    else:
                        # Already visited, but add edge if not present
                        if not dm.has_edge(u, current, key=key):
                            dm.add_edge(u, current, key=key, **data)
        
        # Process directed queue (continue BFS from vertices reached via directed edges)
        while directed_queue:
            current = directed_queue.popleft()
            
            # Process undirected edges incident to current (orient away from current)
            for u, v, key, data in graph._undirected.edges(current, keys=True, data=True):
                # Normalize edge representation to avoid processing twice
                edge_key = (min(u, v), max(u, v), key)
                if edge_key in processed_undirected_edges:
                    continue
                processed_undirected_edges.add(edge_key)
                
                neighbor = v if u == current else u
                
                if neighbor not in visited:
                    visited.add(neighbor)
                    dm.add_edge(current, neighbor, key=key, **data)
                    directed_queue.append(neighbor)
                else:
                    # Already visited, but add edge if not present
                    if not dm.has_edge(current, neighbor, key=key):
                        dm.add_edge(current, neighbor, key=key, **data)
            
            # Process directed edges from current
            for u, v, key, data in graph._directed.edges(current, keys=True, data=True):
                if u == current:  # Outgoing edge
                    if v not in visited:
                        visited.add(v)
                        dm.add_edge(current, v, key=key, **data)
                        directed_queue.append(v)
                    else:
                        if not dm.has_edge(current, v, key=key):
                            dm.add_edge(current, v, key=key, **data)
    
    # Add any remaining directed edges that weren't reached in BFS
    # (shouldn't happen in a connected graph, but handle it)
    for u, v, key, data in graph._directed.edges(keys=True, data=True):
        if not dm.has_edge(u, v, key=key):
            dm.add_edge(u, v, key=key, **data)
    
    # Check for unreachable undirected edges
    for u, v, key in graph._undirected.edges(keys=True):
        edge_key = (min(u, v), max(u, v), key)
        if edge_key not in processed_undirected_edges:
            raise ValueError(
                f"Undirected edge ({u}, {v}, key={key}) is not reachable from root {root}. "
                f"This may indicate disconnected components."
            )
    
    # Verify all nodes are included
    if dm.number_of_nodes() != graph.number_of_nodes():
        raise ValueError(
            f"Orientation failed: resulting graph has {dm.number_of_nodes()} nodes, "
            f"expected {graph.number_of_nodes()}. This may indicate disconnected components "
            f"or cycles that prevent a valid orientation."
        )
    
    return dm


def suppress_degree2_node(graph: 'MixedMultiGraph', node: T, merged_attrs: dict[str, Any] | None = None) -> None:
    """
    Suppress a single degree-2 node in a mixed multigraph in place.
    
    A degree-2 node has exactly two incident edges. Suppression connects the two
    neighbors directly, preserving edge types according to these rules:
    - undirected + undirected -> undirected edge
    - directed_in + directed_out (u->x, x->v) -> directed edge (u->v)
    - directed_out + directed_in (x->u, v->x) -> directed edge (v->u)
    - directed_in + undirected (u->x, x—v) -> directed edge (u->v)
    - undirected + directed_in (u—x, v->x) -> directed edge (u->v)
    - directed_out + undirected (x->u, x—v) -> directed edge (u->v)
    - undirected + directed_out (u—x, x->v) -> directed edge (u->v)
    
    Invalid combinations that raise ValueError:
    - directed_in + directed_in: Multiple incoming directed edges
    - directed_out + directed_out: Multiple outgoing directed edges
    
    This operation modifies the graph in place. Suppression may create parallel edges.
    
    Edge attributes are handled as follows:
    - If `merged_attrs` is provided: these attributes are used directly for the new edge.
      This allows the caller to apply special merging logic before suppression.
    - If `merged_attrs` is None: attributes are merged by taking the first edge's data,
      then the second edge's data overriding. The order of edges is determined by edge
      type priority: directed_in edges first, then directed_out edges, then undirected
      edges. If both edges are of the same type, the order may be non-deterministic
      (depends on graph iteration order). For attributes present in both edges, the
      second edge's value overrides the first.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The mixed multigraph to modify. **This function modifies the graph in place.**
    node : T
        The degree-2 node to suppress.
    merged_attrs : dict[str, Any] | None, optional
        Pre-merged attributes to use for the resulting edge. If None, attributes are merged
        by taking edge1_data first, then edge2_data overriding. 
        
        When provided, these attributes will be used directly for the new edge created
        during suppression. This is useful when special attribute handling is needed.
    
    Raises
    ------
    ValueError
        If the node is not degree-2, or has an invalid edge configuration (e.g., multiple
        directed edges in the same direction, or more than 2 incident edges).
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> G.add_undirected_edge(1, 2)
    0
    >>> G.add_undirected_edge(2, 3)
    0
    >>> suppress_degree2_node(G, 2)
    >>> list(G.edges())
    [(1, 3)]
    """
    # Check that node exists
    if node not in graph.nodes():
        raise ValueError(f"Node {node} not found in graph")
    
    # Verify node is degree-2 using the public API
    if graph.degree(node) != 2:
        raise ValueError(
            f"Node {node} has degree {graph.degree(node)}, expected degree 2"
        )
    
    # Collect incident edges using the public API
    # Note: We already verified degree == 2, so we know there are exactly 2 incident edges
    undirected_edges = list(graph.incident_undirected_edges(node, keys=True, data=True))
    directed_in = list(graph.incident_parent_edges(node, keys=True, data=True))
    directed_out = list(graph.incident_child_edges(node, keys=True, data=True))
    
    # Determine neighbors and edge types
    neighbors = []
    edge_types = []  # 'undirected', 'directed_in', 'directed_out'
    
    for u, v, key, data in directed_in:
        neighbors.append((u, key, data))
        edge_types.append('directed_in')
    
    for u, v, key, data in directed_out:
        neighbors.append((v, key, data))
        edge_types.append('directed_out')
    
    for u, v, key, data in undirected_edges:
        neighbor = v if u == node else u
        neighbors.append((neighbor, key, data))
        edge_types.append('undirected')
    
    # Defensive check: we already verified degree == 2, so this should always be true
    # This check catches potential inconsistencies in the graph implementation
    if len(neighbors) != 2:
        raise ValueError(
            f"Node {node} has {len(neighbors)} incident edges, expected exactly 2. "
            f"This may indicate an inconsistency in the graph implementation."
        )
    
    (n1, k1, d1), (n2, k2, d2) = neighbors
    type1, type2 = edge_types
    
    # Remove the node and its incident edges
    graph.remove_node(node)
    
    # Determine resulting edge type and direction
    # Use provided merged_attrs if given, otherwise merge attributes (prefer d1, then d2)
    if merged_attrs is None:
        merged_attrs = {}
        if d1:
            merged_attrs.update(d1)
        if d2:
            merged_attrs.update(d2)
    
    if type1 == 'undirected' and type2 == 'undirected':
        # undirected + undirected -> undirected
        graph.add_undirected_edge(n1, n2, key=k1 if k1 == k2 else None, **merged_attrs)
    elif type1 == 'directed_in' and type2 == 'directed_out':
        # u->x, x->v -> u->v
        graph.add_directed_edge(n1, n2, key=k1 if k1 == k2 else None, **merged_attrs)
    elif type1 == 'directed_out' and type2 == 'directed_in':
        # x->u, v->x -> v->u (n1 is target of x->n1, n2 is source of n2->x)
        graph.add_directed_edge(n2, n1, key=k1 if k1 == k2 else None, **merged_attrs)
    elif type1 == 'directed_in' and type2 == 'undirected':
        # u->x, x—v -> u->v
        graph.add_directed_edge(n1, n2, key=k1 if k1 == k2 else None, **merged_attrs)
    elif type1 == 'undirected' and type2 == 'directed_in':
        # u—x, v->x -> u->v
        graph.add_directed_edge(n1, n2, key=k1 if k1 == k2 else None, **merged_attrs)
    elif type1 == 'directed_out' and type2 == 'undirected':
        # x->u, x—v -> u->v (from target of directed edge to neighbor of undirected edge)
        # n1 is the target of x->n1, n2 is neighbor of x—n2
        graph.add_directed_edge(n1, n2, key=k1 if k1 == k2 else None, **merged_attrs)
    elif type1 == 'undirected' and type2 == 'directed_out':
        # u—x, x->v -> u->v
        graph.add_directed_edge(n1, n2, key=k1 if k1 == k2 else None, **merged_attrs)
    elif type1 == 'directed_in' and type2 == 'directed_in':
        # Multiple incoming directed edges - invalid for degree-2 node
        raise ValueError(
            f"Node {node} has multiple incoming directed edges, "
            f"cannot determine suppression direction"
        )
    elif type1 == 'directed_out' and type2 == 'directed_out':
        # Multiple outgoing directed edges - invalid for degree-2 node
        raise ValueError(
            f"Node {node} has multiple outgoing directed edges, "
            f"cannot determine suppression direction"
        )
    else:
        raise ValueError(
            f"Unexpected edge type combination for node {node}: "
            f"{type1} and {type2}"
        )


def identify_parallel_edge(graph: 'MixedMultiGraph', u: T, v: T, merged_attrs: dict[str, Any] | None = None) -> None:
    """
    Identify all parallel edges between two nodes by keeping one edge.
    
    This function removes all parallel edges between u and v except one,
    effectively merging them into a single edge. The first edge (lowest key)
    is kept, and all others are removed.
    
    For mixed multigraphs, parallel edges can be either:
    - Multiple directed edges from u to v (if any directed edges exist)
    - Multiple undirected edges between u and v (if any undirected edges exist)
    
    Note: Directed and undirected edges between the same nodes are mutually exclusive
    in MixedMultiGraph, so this function handles only one type at a time.
    
    Edge attributes are handled as follows:
    - If `merged_attrs` is provided: these attributes are used directly for the kept edge.
      This allows the caller to apply special merging logic before identification.
    - If `merged_attrs` is None: attributes are merged by taking the first edge's data,
      then subsequent edges' data overriding. For attributes present in multiple edges,
      the last edge's value overrides earlier values.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The mixed multigraph to modify. **This function modifies the graph in place.**
    u : T
        First node.
    v : T
        Second node.
    merged_attrs : dict[str, Any] | None, optional
        Pre-merged attributes to use for the kept edge. If None, attributes are merged
        by taking the first edge's data first, then subsequent edges' data overriding.
        
        When provided, these attributes will be used directly for the kept edge.
        This is useful when special attribute handling is needed.
    
    Raises
    ------
    ValueError
        If either node is not in the graph, or if no edges exist between u and v.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> G.add_directed_edge(1, 2, weight=1.0)
    0
    >>> G.add_directed_edge(1, 2, weight=2.0)
    1
    >>> G.add_directed_edge(1, 2, label='test')
    2
    >>> identify_parallel_edge(G, 1, 2)
    >>> G.number_of_edges(1, 2)
    1
    >>> # With undirected edges
    >>> G2 = MixedMultiGraph()
    >>> G2.add_undirected_edge(1, 2, weight=1.0)
    0
    >>> G2.add_undirected_edge(1, 2, weight=2.0)
    1
    >>> identify_parallel_edge(G2, 1, 2, merged_attrs={'weight': 3.0, 'merged': True})
    >>> edge_data = G2._undirected[1][2][0]
    >>> edge_data['weight']
    3.0
    >>> edge_data['merged']
    True
    """
    # Check that nodes exist
    if u not in graph.nodes():
        raise ValueError(f"Node {u} not found in graph")
    if v not in graph.nodes():
        raise ValueError(f"Node {v} not found in graph")
    
    # Check if there are any edges between u and v
    if not graph.has_edge(u, v):
        raise ValueError(f"No edges exist between nodes {u} and {v}")
    
    # Check if edges are directed or undirected (mutually exclusive)
    has_directed = graph._directed.has_edge(u, v)
    has_undirected = graph._undirected.has_edge(u, v) or graph._undirected.has_edge(v, u)
    
    if has_directed:
        # Handle parallel directed edges
        edges_dict = graph._directed[u].get(v, {})
        if not edges_dict:
            raise ValueError(f"No directed edges exist between nodes {u} and {v}")
        
        num_edges = len(edges_dict)
        if num_edges <= 1:
            # No parallel edges, nothing to do
            return
        
        # Collect all edge keys and data
        edge_keys = sorted(edges_dict.keys())
        first_key = edge_keys[0]
        first_data = edges_dict[first_key]
        
        # Determine merged attributes
        if merged_attrs is None:
            merged_attrs = {}
            # Start with first edge's data
            if first_data:
                merged_attrs.update(first_data)
            # Then override with subsequent edges' data
            for key in edge_keys[1:]:
                edge_data = edges_dict[key]
                if edge_data:
                    merged_attrs.update(edge_data)
        
        # Remove all edges between u and v
        for key in edge_keys:
            graph.remove_directed_edge(u, v, key=key)
        
        # Add back a single edge with merged attributes
        graph.add_directed_edge(u, v, key=first_key, **merged_attrs)
    
    elif has_undirected:
        # Handle parallel undirected edges
        # Check both directions (u, v) and (v, u) since undirected edges are symmetric
        edges_dict = graph._undirected[u].get(v, {})
        if not edges_dict:
            edges_dict = graph._undirected[v].get(u, {})
        
        if not edges_dict:
            raise ValueError(f"No undirected edges exist between nodes {u} and {v}")
        
        num_edges = len(edges_dict)
        if num_edges <= 1:
            # No parallel edges, nothing to do
            return
        
        # Collect all edge keys and data
        edge_keys = sorted(edges_dict.keys())
        first_key = edge_keys[0]
        first_data = edges_dict[first_key]
        
        # Determine merged attributes
        if merged_attrs is None:
            merged_attrs = {}
            # Start with first edge's data
            if first_data:
                merged_attrs.update(first_data)
            # Then override with subsequent edges' data
            for key in edge_keys[1:]:
                edge_data = edges_dict[key]
                if edge_data:
                    merged_attrs.update(edge_data)
        
        # Remove all edges between u and v (handle both directions)
        for key in edge_keys:
            # Try removing as (u, v) first, then (v, u)
            try:
                graph.remove_edge(u, v, key=key)
            except (ValueError, KeyError):
                graph.remove_edge(v, u, key=key)
        
        # Add back a single edge with merged attributes
        graph.add_undirected_edge(u, v, key=first_key, **merged_attrs)
    
    else:
        raise ValueError(f"No edges exist between nodes {u} and {v}")
