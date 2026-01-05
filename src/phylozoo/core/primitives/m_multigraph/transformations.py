"""
Graph transformations module.

This module provides functions to transform MixedMultiGraph instances
(e.g., identify nodes, orient edges, suppress degree-2 nodes, etc.).
"""

from typing import TYPE_CHECKING, Any, TypeVar, Iterable
from collections import deque

import networkx as nx

if TYPE_CHECKING:
    from ..d_multigraph import DirectedMultiGraph
    from . import MixedMultiGraph
else:
    from . import MixedMultiGraph
    from ..d_multigraph import DirectedMultiGraph

T = TypeVar('T')


def identify_vertices(graph: 'MixedMultiGraph', vertices: list[T], merged_attrs: dict[str, Any] | None = None) -> None:
    """
    Identify multiple vertices by keeping the first vertex.
    
    This function identifies all vertices in the list with the first vertex.
    All edges incident to the other vertices are moved to the first vertex,
    and the other vertices are removed. The first vertex's attributes are
    preserved (or replaced with merged_attrs if provided).
    
    This operation modifies the graph in place. Identification may create
    new parallel edges. Self-loops are not created (edges from a vertex
    to itself are removed).
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The mixed multigraph to modify. **This function modifies the graph in place.**
    vertices : list[T]
        List of vertices to identify. The first vertex will be kept, and all
        others will be merged into it.
    merged_attrs : dict[str, Any] | None, optional
        Attributes to use for the kept vertex. If None, the first vertex's
        attributes are preserved. If provided, these attributes replace the
        first vertex's attributes.
    
    Raises
    ------
    ValueError
        If the vertices list is empty, if any vertex is not in the graph, if
        identification would create both directed and undirected edges between
        the same pair of nodes, or if identification would create edges in both
        directions between the same pair of nodes (u->v and v->u), which is not allowed.
    
    Notes
    -----
    This function does not create self-loops. Any edges that would become
    self-loops after identification are removed.
    
    Identification may create parallel edges. However, if identification would
    result in both directed and undirected edges between the same pair of nodes
    (e.g., both a directed edge u->v and an undirected edge u-v), or edges in
    both directions (e.g., both u->v and v->u), a ValueError is raised as this
    violates the graph's constraints.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> G.add_directed_edge(1, 2)
    0
    >>> G.add_undirected_edge(2, 3)
    0
    >>> G.add_undirected_edge(3, 4)
    0
    >>> identify_vertices(G, [1, 2, 3])
    >>> list(G.nodes())
    [1, 4]
    >>> G.has_edge(1, 4)
    True
    """
    if not vertices:
        raise ValueError("Vertices list cannot be empty")
    
    vertices_list = list(vertices)
    if len(vertices_list) < 2:
        return  # Nothing to identify
    
    first_vertex = vertices_list[0]
    other_vertices = vertices_list[1:]
    
    # Check that all vertices exist
    for v in vertices_list:
        if v not in graph.nodes():
            raise ValueError(f"Vertex {v} not found in graph")
    
    # Before merging, check if identification would create conflicting edge types
    # Collect all edges that would be created, tracking their types
    edges_to_create: dict[tuple[T, T], set[str]] = {}  # (u, v) -> set of edge types ('directed' or 'undirected')
    
    for other_v in other_vertices:
        # Check outgoing directed edges from other_v
        for u, v, key, data in graph.incident_child_edges(other_v, keys=True, data=True):
            if v == first_vertex or v in other_vertices:
                continue  # Skip self-loops and edges to other vertices being merged
            edge_key = (first_vertex, v)
            if edge_key not in edges_to_create:
                edges_to_create[edge_key] = set()
            edges_to_create[edge_key].add('directed')
        
        # Check incoming directed edges to other_v
        for u, v, key, data in graph.incident_parent_edges(other_v, keys=True, data=True):
            if u == first_vertex or u in other_vertices:
                continue  # Skip self-loops and edges from other vertices being merged
            edge_key = (u, first_vertex)
            if edge_key not in edges_to_create:
                edges_to_create[edge_key] = set()
            edges_to_create[edge_key].add('directed')
        
        # Check undirected edges incident to other_v
        for u, v, key, data in graph.incident_undirected_edges(other_v, keys=True, data=True):
            neighbor = v if u == other_v else u
            if neighbor == first_vertex or neighbor in other_vertices:
                continue  # Skip self-loops and edges to other vertices being merged
            # Normalize edge key (undirected edges are symmetric)
            edge_key = (min(first_vertex, neighbor), max(first_vertex, neighbor))
            if edge_key not in edges_to_create:
                edges_to_create[edge_key] = set()
            edges_to_create[edge_key].add('undirected')
    
    # Check for conflicts: both directed and undirected edges between same nodes
    # Also check for bidirectional directed edges (u->v and v->u)
    # Also check against existing edges from first_vertex
    for (u, v), edge_types in edges_to_create.items():
        # Check if edges_to_create has both types (directed and undirected)
        if 'directed' in edge_types and 'undirected' in edge_types:
            raise ValueError(
                f"Identification would create both directed and undirected edges between {u} and {v}, "
                f"which violates mutual exclusivity."
            )
        
        # Check for bidirectional directed edges (u->v and v->u)
        reverse_key = (v, u)
        if reverse_key in edges_to_create:
            reverse_types = edges_to_create[reverse_key]
            if 'directed' in edge_types and 'directed' in reverse_types:
                raise ValueError(
                    f"Identification would create edges in both directions between {u} and {v}, "
                    f"which is not allowed."
                )
        
        # Check against existing edges from first_vertex
        # Check if first_vertex already has a directed edge to/from this node
        if u == first_vertex:
            if graph._directed.has_edge(first_vertex, v):
                if 'undirected' in edge_types:
                    raise ValueError(
                        f"Identification would create both directed and undirected edges between {first_vertex} and {v}, "
                        f"which violates mutual exclusivity."
                    )
                # Check for bidirectional: first_vertex->v exists, and we'd create v->first_vertex
                if 'directed' in edge_types:
                    # Check if reverse direction would be created
                    if (v, first_vertex) in edges_to_create and 'directed' in edges_to_create[(v, first_vertex)]:
                        raise ValueError(
                            f"Identification would create edges in both directions between {first_vertex} and {v}, "
                            f"which is not allowed."
                        )
            if graph._directed.has_edge(v, first_vertex):
                if 'undirected' in edge_types:
                    raise ValueError(
                        f"Identification would create both directed and undirected edges between {first_vertex} and {v}, "
                        f"which violates mutual exclusivity."
                    )
                # Check for bidirectional: v->first_vertex exists, and we'd create first_vertex->v
                if 'directed' in edge_types:
                    raise ValueError(
                        f"Identification would create edges in both directions between {first_vertex} and {v}, "
                        f"which is not allowed."
                    )
        elif v == first_vertex:
            if graph._directed.has_edge(first_vertex, u):
                if 'undirected' in edge_types:
                    raise ValueError(
                        f"Identification would create both directed and undirected edges between {first_vertex} and {u}, "
                        f"which violates mutual exclusivity."
                    )
                # Check for bidirectional: first_vertex->u exists, and we'd create u->first_vertex
                if 'directed' in edge_types:
                    raise ValueError(
                        f"Identification would create edges in both directions between {first_vertex} and {u}, "
                        f"which is not allowed."
                    )
            if graph._directed.has_edge(u, first_vertex):
                if 'undirected' in edge_types:
                    raise ValueError(
                        f"Identification would create both directed and undirected edges between {first_vertex} and {u}, "
                        f"which violates mutual exclusivity."
                    )
                # Check for bidirectional: u->first_vertex exists, and we'd create first_vertex->u
                if 'directed' in edge_types:
                    # Check if reverse direction would be created
                    if (first_vertex, u) in edges_to_create and 'directed' in edges_to_create[(first_vertex, u)]:
                        raise ValueError(
                            f"Identification would create edges in both directions between {first_vertex} and {u}, "
                            f"which is not allowed."
                        )
        
        # Check if first_vertex already has an undirected edge to this node
        if u == first_vertex:
            if graph._undirected.has_edge(first_vertex, v):
                if 'directed' in edge_types:
                    raise ValueError(
                        f"Identification would create both directed and undirected edges between {first_vertex} and {v}, "
                        f"which violates mutual exclusivity."
                    )
        elif v == first_vertex:
            if graph._undirected.has_edge(first_vertex, u):
                if 'directed' in edge_types:
                    raise ValueError(
                        f"Identification would create both directed and undirected edges between {first_vertex} and {u}, "
                        f"which violates mutual exclusivity."
                    )
        
        # Also check reverse direction for undirected edges
        reverse_key_undir = (v, u) if u < v else (u, v)
        if reverse_key_undir in edges_to_create:
            reverse_types = edges_to_create[reverse_key_undir]
            if 'directed' in edge_types and 'undirected' in reverse_types:
                raise ValueError(
                    f"Identification would create both directed and undirected edges between {u} and {v}, "
                    f"which violates mutual exclusivity."
                )
            if 'undirected' in edge_types and 'directed' in reverse_types:
                raise ValueError(
                    f"Identification would create both directed and undirected edges between {u} and {v}, "
                    f"which violates mutual exclusivity."
                )
    
    # Collect node attributes from first vertex
    first_vertex_attrs = {}
    if first_vertex in graph._undirected.nodes():
        first_vertex_attrs.update(graph._undirected.nodes[first_vertex])
    if first_vertex in graph._directed.nodes():
        first_vertex_attrs.update(graph._directed.nodes[first_vertex])
    
    # Merge all other vertices into first vertex
    for other_v in other_vertices:
        # Collect all edges incident to other_v before removal
        directed_edges_to_add: list[tuple[T, T, int, dict[str, Any]]] = []
        undirected_edges_to_add: list[tuple[T, T, int, dict[str, Any]]] = []
        
        # Outgoing directed edges from other_v
        for u, v, key, data in graph.incident_child_edges(other_v, keys=True, data=True):
            if v == first_vertex or v in other_vertices:
                continue  # Skip self-loops and edges to other vertices being merged
            directed_edges_to_add.append((first_vertex, v, key, data or {}))
        
        # Incoming directed edges to other_v
        for u, v, key, data in graph.incident_parent_edges(other_v, keys=True, data=True):
            if u == first_vertex or u in other_vertices:
                continue  # Skip self-loops and edges from other vertices being merged
            directed_edges_to_add.append((u, first_vertex, key, data or {}))
        
        # Undirected edges incident to other_v
        for u, v, key, data in graph.incident_undirected_edges(other_v, keys=True, data=True):
            neighbor = v if u == other_v else u
            if neighbor == first_vertex or neighbor in other_vertices:
                continue  # Skip self-loops and edges to other vertices being merged
            undirected_edges_to_add.append((first_vertex, neighbor, key, data or {}))
        
        # Remove the vertex (this removes all its incident edges)
        graph.remove_node(other_v)
        
        # Add edges to first_vertex
        for u, v, key, data in directed_edges_to_add:
            # Check if edge with this key already exists - if so, let it auto-generate
            # This allows parallel edges to be created
            if graph.has_edge(u, v, key=key):
                # Key conflict - let add_directed_edge auto-generate a new key
                graph.add_directed_edge(u, v, key=None, **data)
            else:
                graph.add_directed_edge(u, v, key=key, **data)
        
        for u, v, key, data in undirected_edges_to_add:
            # Check if edge with this key already exists - if so, let it auto-generate
            # This allows parallel edges to be created
            if graph.has_edge(u, v, key=key):
                # Key conflict - let add_undirected_edge auto-generate a new key
                graph.add_undirected_edge(u, v, key=None, **data)
            else:
                graph.add_undirected_edge(u, v, key=key, **data)
    
    # Update first vertex attributes
    if merged_attrs is not None:
        # Replace attributes with merged_attrs
        for attr_name, attr_value in merged_attrs.items():
            if first_vertex in graph._undirected.nodes():
                graph._undirected.nodes[first_vertex][attr_name] = attr_value
            if first_vertex in graph._directed.nodes():
                graph._directed.nodes[first_vertex][attr_name] = attr_value
            if first_vertex in graph._combined.nodes():
                graph._combined.nodes[first_vertex][attr_name] = attr_value
        # Remove attributes not in merged_attrs
        attrs_to_remove = set(first_vertex_attrs.keys()) - set(merged_attrs.keys())
        for attr_name in attrs_to_remove:
            if first_vertex in graph._undirected.nodes() and attr_name in graph._undirected.nodes[first_vertex]:
                del graph._undirected.nodes[first_vertex][attr_name]
            if first_vertex in graph._directed.nodes() and attr_name in graph._directed.nodes[first_vertex]:
                del graph._directed.nodes[first_vertex][attr_name]
            if first_vertex in graph._combined.nodes() and attr_name in graph._combined.nodes[first_vertex]:
                del graph._combined.nodes[first_vertex][attr_name]
    # Otherwise, first vertex's attributes are already preserved


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
                # Use string comparison to handle mixed types (e.g., int and str node IDs)
                if str(u) <= str(v):
                    edge_key = (u, v, key)
                else:
                    edge_key = (v, u, key)
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
                # Use string comparison to handle mixed types (e.g., int and str node IDs)
                if str(u) <= str(v):
                    edge_key = (u, v, key)
                else:
                    edge_key = (v, u, key)
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
        # Use string comparison to handle mixed types (e.g., int and str node IDs)
        if str(u) <= str(v):
            edge_key = (u, v, key)
        else:
            edge_key = (v, u, key)
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
    - directed_in + undirected (u->x, x—v) -> undirected edge (u—v)
    - undirected + directed_out (u—x, x->v) -> directed edge (u->v)
    
    Invalid combinations that raise ValueError:
    - directed_in + directed_in: Multiple incoming directed edges
    - directed_out + directed_out: Multiple outgoing directed edges
    
    Note: Edge types are processed in order (directed_in first, then undirected, then
    directed_out), so the order in the rules above is deterministic.
    
    This operation modifies the graph in place. Suppression may create parallel edges.
    
    Edge attributes are handled as follows:
    - If `merged_attrs` is provided: these attributes are used directly for the new edge.
      This allows the caller to apply special merging logic before suppression.
    - If `merged_attrs` is None: attributes are merged by taking the first edge's data,
      then the second edge's data overriding. The order of edges is determined by edge
      type priority: directed_in edges first, then undirected edges, then directed_out
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
    # Order: directed_in first, then undirected, then directed_out
    neighbors = []
    edge_types = []  # 'undirected', 'directed_in', 'directed_out'
    
    for u, v, key, data in directed_in:
        neighbors.append((u, key, data))
        edge_types.append('directed_in')
    
    for u, v, key, data in undirected_edges:
        neighbor = v if u == node else u
        neighbors.append((neighbor, key, data))
        edge_types.append('undirected')
    
    for u, v, key, data in directed_out:
        neighbors.append((v, key, data))
        edge_types.append('directed_out')
    
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
        # Always use key=None to allow parallel edges to be created if the edge already exists
        graph.add_undirected_edge(n1, n2, key=None, **merged_attrs)
    elif type1 == 'directed_in' and type2 == 'directed_out':
        # u->x, x->v -> u->v
        # n1 is source of incoming (u), n2 is target of outgoing (v) -> u->v
        # Always use key=None to allow parallel edges to be created if the edge already exists
        graph.add_directed_edge(n1, n2, key=None, **merged_attrs)
    elif type1 == 'directed_in' and type2 == 'undirected':
        # u->x, x—v -> undirected edge
        # n1 is source of incoming (u), n2 is neighbor of undirected (v) -> undirected u—v
        # Always use key=None to allow parallel edges to be created if the edge already exists
        graph.add_undirected_edge(n1, n2, key=None, **merged_attrs)
    elif type1 == 'undirected' and type2 == 'directed_out':
        # u—x, x->v -> directed edge
        # n1 is neighbor of undirected (u), n2 is target of outgoing (v) -> u->v
        # Always use key=None to allow parallel edges to be created if the edge already exists
        graph.add_directed_edge(n1, n2, key=None, **merged_attrs)
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


def subgraph(graph: 'MixedMultiGraph', nodes: Iterable[T]) -> 'MixedMultiGraph':
    """
    Return the induced subgraph of `graph` on the given `nodes`.

    The returned object is a new MixedMultiGraph instance containing only
    the specified nodes and any edges (both undirected and directed,
    including parallel edges and their keys/attributes) whose endpoints are
    both in `nodes`. Node and edge attributes are preserved.

    Parameters
    ----------
    graph : MixedMultiGraph
        Source mixed multigraph.
    nodes : Iterable[T]
        Iterable of node identifiers to include in the induced subgraph.

    Returns
    -------
    MixedMultiGraph
        A new MixedMultiGraph containing the induced subgraph.

    Raises
    ------
    ValueError
        If any node in `nodes` is not present in `graph`.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> from phylozoo.core.primitives.m_multigraph.transformations import subgraph
    >>> G = MixedMultiGraph()
    >>> G.add_node(1, label='A')
    >>> G.add_node(2, label='B')
    >>> G.add_undirected_edge(1, 2, weight=1.0)
    0
    >>> H = subgraph(G, [1, 2])
    >>> list(H.nodes())
    [1, 2]
    >>> list(H.undirected_edges_iter(keys=True, data=True))
    [(1, 2, 0, {'weight': 1.0})]
    """
    nodes_set = set(nodes)

    # Empty nodes -> return empty graph
    if not nodes_set:
        return MixedMultiGraph()

    # Validate nodes exist in source graph
    for n in nodes_set:
        if n not in graph.nodes():
            raise ValueError(f"Node {n} not found in graph")

    new_graph = MixedMultiGraph()

    # Preserve node attributes
    for n in nodes_set:
        node_attrs: dict[str, Any] = {}
        if n in graph._undirected.nodes():
            node_attrs.update(dict(graph._undirected.nodes[n]))
        if n in graph._directed.nodes():
            node_attrs.update(dict(graph._directed.nodes[n]))
        if node_attrs:
            new_graph.add_node(n, **node_attrs)
        else:
            new_graph.add_node(n)

    # Preserve undirected edges (keys and data) where both endpoints are in nodes_set
    for u, v, key, data in graph._undirected.edges(keys=True, data=True):
        if u in nodes_set and v in nodes_set:
            edge_data = dict(data) if data else {}
            new_graph.add_undirected_edge(u, v, key=key, **edge_data)

    # Preserve directed edges (keys and data) where both endpoints are in nodes_set
    for u, v, key, data in graph._directed.edges(keys=True, data=True):
        if u in nodes_set and v in nodes_set:
            edge_data = dict(data) if data else {}
            new_graph.add_directed_edge(u, v, key=key, **edge_data)

    return new_graph
