"""
Graph transformations module.

This module provides functions to transform DirectedMultiGraph instances
(e.g., identify nodes, suppress degree-2 nodes, etc.).
"""

from typing import Any, TypeVar, Iterable

import networkx as nx

from . import DirectedMultiGraph

T = TypeVar('T')


def identify_vertices(graph: 'DirectedMultiGraph', vertices: list[T], merged_attrs: dict[str, Any] | None = None) -> None:
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
    graph : DirectedMultiGraph
        The directed multigraph to modify. **This function modifies the graph in place.**
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
        If the vertices list is empty, if any vertex is not in the graph, or
        if identification would create edges in both directions between the
        same pair of nodes (u->v and v->u), which is not allowed.
    
    Notes
    -----
    This function does not create self-loops. Any edges that would become
    self-loops after identification are removed.
    
    Identification may create parallel edges. However, if identification would
    result in edges in both directions between the same pair of nodes (e.g.,
    both u->v and v->u), a ValueError is raised as this violates the graph's
    constraints.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2)
    0
    >>> G.add_edge(2, 3)
    0
    >>> G.add_edge(3, 4)
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
    
    # Before merging, check if identification would create bidirectional edges
    # After merging, edges from other_vertices will be moved to first_vertex
    # We need to check if this would create both u->v and v->u for any external node pair
    
    # Collect all external nodes that would have edges to/from first_vertex after merging
    external_outgoing_targets: set[T] = set()  # Nodes that would receive edges FROM first_vertex
    external_incoming_sources: set[T] = set()  # Nodes that would have edges TO first_vertex
    
    for other_v in other_vertices:
        # Check outgoing edges from other_v (will become first_vertex -> target)
        for u, v, key, data in graph.incident_child_edges(other_v, keys=True, data=True):
            if v == first_vertex or v in other_vertices:
                continue  # Skip self-loops and edges to vertices being merged
            external_outgoing_targets.add(v)
        
        # Check incoming edges to other_v (will become source -> first_vertex)
        for u, v, key, data in graph.incident_parent_edges(other_v, keys=True, data=True):
            if u == first_vertex or u in other_vertices:
                continue  # Skip self-loops and edges from vertices being merged
            external_incoming_sources.add(u)
    
    # Check if first_vertex already has edges that would conflict
    # Check outgoing edges from first_vertex
    for u, v, key, data in graph.incident_child_edges(first_vertex, keys=True, data=True):
        if v not in vertices_list:  # External target
            # After merge: first_vertex -> v (already exists)
            # Check if merging would also create v -> first_vertex
            if v in external_incoming_sources:
                raise ValueError(
                    f"Identification would create edges in both directions between {first_vertex} and {v}, "
                    f"which is not allowed."
                )
    
    # Check incoming edges to first_vertex
    for u, v, key, data in graph.incident_parent_edges(first_vertex, keys=True, data=True):
        if u not in vertices_list:  # External source
            # After merge: u -> first_vertex (already exists)
            # Check if merging would also create first_vertex -> u
            if u in external_outgoing_targets:
                raise ValueError(
                    f"Identification would create edges in both directions between {first_vertex} and {u}, "
                    f"which is not allowed."
                )
    
    # Also check if edges being merged would create bidirectional edges between external nodes
    # This handles the case where other_vertex has u->v and first_vertex has v->u
    for target in external_outgoing_targets:
        # After merge: first_vertex -> target
        # Check if first_vertex already has target -> first_vertex
        if graph.has_edge(target, first_vertex) and target not in vertices_list:
            raise ValueError(
                f"Identification would create edges in both directions between {first_vertex} and {target}, "
                f"which is not allowed."
            )
    
    for source in external_incoming_sources:
        # After merge: source -> first_vertex
        # Check if first_vertex already has first_vertex -> source
        if graph.has_edge(first_vertex, source) and source not in vertices_list:
            raise ValueError(
                f"Identification would create edges in both directions between {first_vertex} and {source}, "
                f"which is not allowed."
            )
    
    # Collect node attributes from first vertex
    first_vertex_attrs = {}
    if first_vertex in graph._graph.nodes():
        first_vertex_attrs.update(graph._graph.nodes[first_vertex])
    
    # Merge all other vertices into first vertex
    for other_v in other_vertices:
        # Collect all edges incident to other_v before removal
        edges_to_add: list[tuple[T, T, int, dict[str, Any]]] = []
        
        # Outgoing edges from other_v
        for u, v, key, data in graph.incident_child_edges(other_v, keys=True, data=True):
            if v == first_vertex or v in other_vertices:
                continue  # Skip self-loops and edges to other vertices being merged
            edges_to_add.append((first_vertex, v, key, data or {}))
        
        # Incoming edges to other_v
        for u, v, key, data in graph.incident_parent_edges(other_v, keys=True, data=True):
            if u == first_vertex or u in other_vertices:
                continue  # Skip self-loops and edges from other vertices being merged
            edges_to_add.append((u, first_vertex, key, data or {}))
        
        # Remove the vertex (this removes all its incident edges)
        graph.remove_node(other_v)
        
        # Add edges to first_vertex
        for u, v, key, data in edges_to_add:
            # Check if edge with this key already exists - if so, let it auto-generate
            # This allows parallel edges to be created
            if graph.has_edge(u, v, key=key):
                # Key conflict - let add_edge auto-generate a new key
                graph.add_edge(u, v, key=None, **data)
            else:
                graph.add_edge(u, v, key=key, **data)
    
    # Update first vertex attributes
    if merged_attrs is not None:
        # Replace attributes with merged_attrs
        for attr_name, attr_value in merged_attrs.items():
            graph._graph.nodes[first_vertex][attr_name] = attr_value
            graph._combined.nodes[first_vertex][attr_name] = attr_value
        # Remove attributes not in merged_attrs
        attrs_to_remove = set(first_vertex_attrs.keys()) - set(merged_attrs.keys())
        for attr_name in attrs_to_remove:
            if attr_name in graph._graph.nodes[first_vertex]:
                del graph._graph.nodes[first_vertex][attr_name]
            if attr_name in graph._combined.nodes[first_vertex]:
                del graph._combined.nodes[first_vertex][attr_name]
    # Otherwise, first vertex's attributes are already preserved


def suppress_degree2_node(graph: 'DirectedMultiGraph', node: T, merged_attrs: dict[str, Any] | None = None) -> None:
    """
    Suppress a single degree-2 node in a directed multigraph in place.
    
    A degree-2 node has exactly two incident edges. For a directed multigraph, the only
    valid configuration is indegree=1 and outdegree=1 (u->x->v). Suppression connects
    the two neighbors directly: u->x->v becomes u->v.
    
    Invalid configurations that raise ValueError:
    - indegree=2, outdegree=0: Multiple incoming directed edges
    - indegree=0, outdegree=2: Multiple outgoing directed edges
    
    This operation modifies the graph in place. Suppression may create parallel edges.
    
    Edge attributes are handled as follows:
    - If `merged_attrs` is provided: these attributes are used directly for the new edge.
      This allows the caller to apply special merging logic before suppression.
    - If `merged_attrs` is None: attributes are merged by taking the incoming edge's data
      first, then the outgoing edge's data overriding. For attributes present in both edges,
      the outgoing edge's value overrides the incoming edge's value.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The directed multigraph to modify. **This function modifies the graph in place.**
    node : T
        The degree-2 node to suppress.
    merged_attrs : dict[str, Any] | None, optional
        Pre-merged attributes to use for the resulting edge. If None, attributes are merged
        by taking incoming edge data first, then outgoing edge data overriding.
        
        When provided, these attributes will be used directly for the new edge created
        during suppression. This is useful when special attribute handling is needed.
    
    Raises
    ------
    ValueError
        If the node is not degree-2, or has an invalid edge configuration (indegree=2
        with outdegree=0, or indegree=0 with outdegree=2).
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2)
    0
    >>> G.add_edge(2, 3)
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
    
    # Check indegree and outdegree
    in_deg = graph.indegree(node)
    out_deg = graph.outdegree(node)
    
    # Invalid configurations
    if in_deg == 2 and out_deg == 0:
        raise ValueError(
            f"Node {node} has indegree 2 and outdegree 0. "
            f"Cannot suppress node with multiple incoming directed edges."
        )
    if in_deg == 0 and out_deg == 2:
        raise ValueError(
            f"Node {node} has indegree 0 and outdegree 2. "
            f"Cannot suppress node with multiple outgoing directed edges."
        )
    
    # Valid configuration: indegree=1, outdegree=1
    if in_deg != 1 or out_deg != 1:
        raise ValueError(
            f"Node {node} has indegree {in_deg} and outdegree {out_deg}. "
            f"Expected indegree 1 and outdegree 1 for suppression."
        )
    
    # Collect incident edges using the public API
    directed_in = list(graph.incident_parent_edges(node, keys=True, data=True))
    directed_out = list(graph.incident_child_edges(node, keys=True, data=True))
    
    # We should have exactly one incoming and one outgoing edge
    if len(directed_in) != 1 or len(directed_out) != 1:
        raise ValueError(
            f"Node {node} has {len(directed_in)} incoming and {len(directed_out)} outgoing edges. "
            f"Expected exactly 1 of each."
        )
    
    # Get the neighbors and edge data
    (u, _, k1, d1) = directed_in[0]  # u->node
    (_, v, k2, d2) = directed_out[0]  # node->v
    
    # Remove the node and its incident edges
    graph.remove_node(node)
    
    # Determine resulting edge attributes
    # Use provided merged_attrs if given, otherwise merge attributes (incoming first, then outgoing)
    if merged_attrs is None:
        merged_attrs = {}
        if d1:
            merged_attrs.update(d1)
        if d2:
            merged_attrs.update(d2)
    
    # Add the new edge u->v
    graph.add_edge(u, v, key=k1 if k1 == k2 else None, **merged_attrs)


def identify_parallel_edge(graph: 'DirectedMultiGraph', u: T, v: T, merged_attrs: dict[str, Any] | None = None) -> None:
    """
    Identify all parallel edges between two nodes by keeping one edge.
    
    This function removes all parallel edges between u and v except one,
    effectively merging them into a single edge. The first edge (lowest key)
    is kept, and all others are removed.
    
    Edge attributes are handled as follows:
    - If `merged_attrs` is provided: these attributes are used directly for the kept edge.
      This allows the caller to apply special merging logic before identification.
    - If `merged_attrs` is None: attributes are merged by taking the first edge's data,
      then subsequent edges' data overriding. For attributes present in multiple edges,
      the last edge's value overrides earlier values.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The directed multigraph to modify. **This function modifies the graph in place.**
    u : T
        Source node.
    v : T
        Target node.
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
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2, weight=1.0)
    0
    >>> G.add_edge(1, 2, weight=2.0)
    1
    >>> G.add_edge(1, 2, label='test')
    2
    >>> identify_parallel_edge(G, 1, 2)
    >>> G.number_of_edges(1, 2)
    1
    >>> # With custom merged attributes
    >>> G2 = DirectedMultiGraph()
    >>> G2.add_edge(1, 2, weight=1.0)
    0
    >>> G2.add_edge(1, 2, weight=2.0)
    1
    >>> identify_parallel_edge(G2, 1, 2, merged_attrs={'weight': 3.0, 'merged': True})
    >>> edge_data = G2._graph[1][2][0]
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
    
    # Get all parallel edges between u and v
    edges_dict = graph._graph[u].get(v, {})
    if not edges_dict:
        raise ValueError(f"No edges exist between nodes {u} and {v}")
    
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
    else:
        # Use provided merged_attrs directly
        pass
    
    # Remove all edges between u and v
    for key in edge_keys:
        graph.remove_edge(u, v, key=key)
    
    # Add back a single edge with merged attributes
    graph.add_edge(u, v, key=first_key, **merged_attrs)


def subgraph(graph: 'DirectedMultiGraph', nodes: Iterable[T]) -> 'DirectedMultiGraph':
    """
    Return the induced subgraph of `graph` on the given `nodes`.

    The returned object is a new DirectedMultiGraph instance containing only
    the specified nodes and any edges (including parallel edges and their
    keys/attributes) whose endpoints are both in `nodes`. Node and edge
    attributes are preserved.

    Parameters
    ----------
    graph : DirectedMultiGraph
        Source directed multigraph.
    nodes : Iterable[T]
        Iterable of node identifiers to include in the induced subgraph.

    Returns
    -------
    DirectedMultiGraph
        A new DirectedMultiGraph containing the induced subgraph.

    Raises
    ------
    ValueError
        If any node in `nodes` is not present in `graph`.

    Examples
    --------
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2, weight=1.0)
    0
    >>> G.add_edge(2, 3, weight=2.0)
    0
    >>> H = subgraph(G, [1, 2])
    >>> list(H.nodes())
    [1, 2]
    >>> list(H.edges_iter(keys=True, data=True))
    [(1, 2, 0, {'weight': 1.0})]
    """
    nodes_set = set(nodes)

    # Empty nodes -> return empty graph
    if not nodes_set:
        return DirectedMultiGraph()

    # Validate nodes exist in source graph
    for n in nodes_set:
        if n not in graph.nodes():
            raise ValueError(f"Node {n} not found in graph")

    new_graph = DirectedMultiGraph()

    # Preserve node attributes
    for n in nodes_set:
        attrs = dict(graph._graph.nodes[n]) if n in graph._graph.nodes else {}
        new_graph.add_node(n, **attrs)

    # Preserve edges (including keys and data) where both endpoints are in nodes_set
    for u, v, key, data in graph._graph.edges(keys=True, data=True):
        if u in nodes_set and v in nodes_set:
            # data can be None; ensure dict
            edge_data = dict(data) if data else {}
            # Preserve the same key where possible
            new_graph.add_edge(u, v, key=key, **edge_data)

    return new_graph
