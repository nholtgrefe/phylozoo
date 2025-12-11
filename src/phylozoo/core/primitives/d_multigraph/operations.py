"""
Operations for DirectedMultiGraph.

This module provides functions for working with DirectedMultiGraph instances,
following NetworkX-style function-based API.
"""

from typing import Any, Dict, Iterator, List, Set, TypeVar

import networkx as nx

from . import DirectedMultiGraph

T = TypeVar('T')


def number_of_connected_components(graph: 'DirectedMultiGraph') -> int:
    """
    Return the number of weakly connected components.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to analyze.
    
    Returns
    -------
    int
        Number of connected components.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2)
    0
    >>> G.add_edge(3, 4)
    0
    >>> number_of_connected_components(G)
    2
    """
    return nx.number_connected_components(graph._combined)


def is_connected(graph: 'DirectedMultiGraph') -> bool:
    """
    Check if graph is weakly connected.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to check.
    
    Returns
    -------
    bool
        True if graph is connected, False otherwise.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2)
    0
    >>> G.add_edge(2, 3)
    0
    >>> is_connected(G)
    True
    """
    return nx.is_connected(graph._combined)


def connected_components(graph: 'DirectedMultiGraph') -> Iterator[set[T]]:
    """
    Get weakly connected components.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to analyze.
    
    Returns
    -------
    Iterator[set[T]]
        Iterator over sets of nodes in each component.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2)
    0
    >>> G.add_edge(3, 4)
    0
    >>> list(connected_components(G))
    [{1, 2}, {3, 4}]
    """
    return nx.connected_components(graph._combined)


def has_self_loops(graph: 'DirectedMultiGraph') -> bool:
    """
    Check whether the directed multigraph contains any self-loops.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to inspect.
    
    Returns
    -------
    bool
        True if at least one self-loop exists, False otherwise.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> from phylozoo.core.primitives.d_multigraph.operations import has_self_loops
    >>> G = DirectedMultiGraph()
    >>> has_self_loops(G)
    False
    >>> G.add_edge(1, 1)
    0
    >>> has_self_loops(G)
    True
    """
    return nx.number_of_selfloops(graph._graph) > 0


def identify_two_nodes(graph: 'DirectedMultiGraph', u: T, v: T) -> None:
    """
    Identify two nodes u and v by keeping node u.
    
    Modifies the graph in place. All edges incident to v are moved to u,
    and v is removed.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to modify.
    u : T
        Node to keep.
    v : T
        Node to identify with u (will be removed).
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2)
    0
    >>> G.add_edge(2, 3)
    0
    >>> identify_two_nodes(G, 1, 2)
    >>> list(G.nodes())
    [1, 3]
    """
    # Use NetworkX's contracted_nodes on each graph
    nx.contracted_nodes(graph._graph, u, v, self_loops=False, copy=False)
    nx.contracted_nodes(graph._combined, u, v, self_loops=False, copy=False)
    
    # Clean up any self-loops that might have been created
    if graph._graph.has_edge(u, u):
        for k in list(graph._graph[u][u].keys()):
            graph._graph.remove_edge(u, u, key=k)
            graph._combined.remove_edge(u, u, key=k)


def identify_node_set(graph: 'DirectedMultiGraph', nodes: list[T] | set[T]) -> None:
    """
    Identify all nodes in the set by keeping the first node.
    
    Modifies the graph in place.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to modify.
    nodes : list[T] | set[T]
        Iterable of nodes to identify. The first node will be kept.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2)
    0
    >>> G.add_edge(2, 3)
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



