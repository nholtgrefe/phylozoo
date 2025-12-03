"""
Directed multi-graph module.

This module provides the DirectedMultiGraph class for working with directed multi-graphs.
"""

from typing import Any, List, Optional, Set, Tuple, TypeVar

import networkx as nx

T = TypeVar('T')


class DirectedMultiGraph(nx.MultiDiGraph):
    """
    Directed multi-graph where all edges are directed.

    This class allows parallel directed edges (multiple directed edges between
    the same nodes), and each parallel edge can have different parameters
    (weights, attributes, etc.) via edge keys.

    Since this inherits from nx.MultiDiGraph, most NetworkX methods work
    directly. This class provides additional convenience methods consistent
    with MixedMultiGraph.

    Parameters
    ----------
    *args
        Additional arguments passed to nx.MultiDiGraph.
    **kwargs
        Additional keyword arguments passed to nx.MultiDiGraph.

    Examples
    --------
    >>> G = DirectedMultiGraph()
    >>> key1 = G.add_edge(1, 2, weight=1.0)
    >>> key2 = G.add_edge(1, 2, weight=2.0)  # Parallel edge
    >>> key1 != key2
    True
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        """
        Initialize a directed multi-graph.

        Parameters
        ----------
        *args
            Additional arguments passed to nx.MultiDiGraph.
        **kwargs
            Additional keyword arguments passed to nx.MultiDiGraph.
        """
        super().__init__(*args, **kwargs)

    def _clear_cache(self) -> None:
        """
        Clear any cached computations.

        This method is a placeholder for subclasses to override if they
        maintain cached values that need to be cleared when the graph changes.
        """
        pass

    def remove_node(self, v: T) -> None:
        """
        Remove node v from the graph.

        Also removes all edges incident to v.

        Parameters
        ----------
        v : T
            Node to remove.

        Examples
        --------
        >>> G = DirectedMultiGraph()
        >>> G.add_edge(1, 2)
        0
        >>> G.add_edge(2, 3)
        0
        >>> G.remove_node(2)
        >>> list(G.nodes())
        [1, 3]
        """
        super().remove_node(v)
        self._clear_cache()

    def remove_nodes_from(self, nodes: List[T] | Set[T]) -> None:
        """
        Remove all nodes in 'nodes' from the graph.

        Parameters
        ----------
        nodes : List[T] | Set[T]
            Iterable of nodes to remove.

        Examples
        --------
        >>> G = DirectedMultiGraph()
        >>> G.add_node(1)
        >>> G.add_node(2)
        >>> G.add_node(3)
        >>> G.remove_nodes_from([1, 2])
        >>> list(G.nodes())
        [3]
        """
        for v in nodes:
            self.remove_node(v)

    def remove_edges_from(self, edges: List[Tuple[T, T] | Tuple[T, T, int]]) -> None:
        """
        Remove all edges in 'edges' from the graph.

        Parameters
        ----------
        edges : List[Tuple[T, T] | Tuple[T, T, int]]
            List of edges to remove. Can be (u, v) or (u, v, key) tuples.

        Examples
        --------
        >>> G = DirectedMultiGraph()
        >>> G.add_edge(1, 2)
        0
        >>> G.add_edge(2, 3)
        0
        >>> G.remove_edges_from([(1, 2), (2, 3)])
        """
        for edge in edges:
            if len(edge) == 2:
                u, v = edge
                # Remove all parallel edges between u and v
                if u in self and v in self[u]:
                    for key in list(self[u][v].keys()):
                        super().remove_edge(u, v, key=key)
            elif len(edge) == 3:
                u, v, key = edge
                super().remove_edge(u, v, key=key)
            else:
                raise ValueError(f"Invalid edge format: {edge}")

    def indegree(self, v: T) -> int:
        """
        Return the indegree of vertex v.

        The indegree is the number of directed edges pointing to v.

        Parameters
        ----------
        v : T
            Vertex.

        Returns
        -------
        int
            Indegree of v.

        Examples
        --------
        >>> G = DirectedMultiGraph()
        >>> G.add_edge(1, 2)
        0
        >>> G.add_edge(3, 2)
        0
        >>> G.indegree(2)
        2
        """
        return self.in_degree(v)

    def outdegree(self, v: T) -> int:
        """
        Return the outdegree of vertex v.

        The outdegree is the number of directed edges pointing from v.

        Parameters
        ----------
        v : T
            Vertex.

        Returns
        -------
        int
            Outdegree of v.

        Examples
        --------
        >>> G = DirectedMultiGraph()
        >>> G.add_edge(1, 2)
        0
        >>> G.add_edge(1, 3)
        0
        >>> G.outdegree(1)
        2
        """
        return self.out_degree(v)

    def is_cutedge(self, u: T, v: T, key: Optional[int] = None) -> bool:
        """
        Check if edge (u, v) is a cut-edge.

        A cut-edge is an edge whose removal increases the number of
        weakly connected components.

        Parameters
        ----------
        u : T
            Source node.
        v : T
            Target node.
        key : Optional[int], optional
            Edge key. If None, checks the first edge found. By default None.

        Returns
        -------
        bool
            True if (u, v) is a cut-edge, False otherwise.

        Raises
        ------
        ValueError
            If (u, v) is not an edge in the graph.

        Examples
        --------
        >>> G = DirectedMultiGraph()
        >>> G.add_edge(1, 2)
        0
        >>> G.add_edge(2, 3)
        0
        >>> G.is_cutedge(2, 3)
        True
        """
        if u not in self or v not in self[u]:
            raise ValueError(f"Edge ({u}, {v}) is not in the graph.")

        # Find edge key if not provided
        if key is None:
            if u in self and v in self[u]:
                key = next(iter(self[u][v].keys()))
            else:
                raise ValueError(f"Edge ({u}, {v}) is not in the graph.")

        modified_graph = self.copy()
        modified_graph.remove_edge(u, v, key=key)
        return nx.number_weakly_connected_components(modified_graph) != nx.number_weakly_connected_components(self)

    def is_cutvertex(self, v: T) -> bool:
        """
        Check if v is a cut-vertex.

        A cut-vertex is a vertex whose removal increases the number of
        weakly connected components.

        Parameters
        ----------
        v : T
            Vertex.

        Returns
        -------
        bool
            True if v is a cut-vertex, False otherwise.

        Raises
        ------
        ValueError
            If v is not a vertex in the graph.

        Examples
        --------
        >>> G = DirectedMultiGraph()
        >>> G.add_edge(1, 2)
        0
        >>> G.add_edge(2, 3)
        0
        >>> G.is_cutvertex(2)
        True
        """
        if v not in self:
            raise ValueError(f"Vertex {v} is not in the graph.")

        modified_graph = self.copy()
        modified_graph.remove_node(v)
        return nx.number_weakly_connected_components(modified_graph) != nx.number_weakly_connected_components(self)

    def identify_two_nodes(self, u: T, v: T) -> None:
        """
        Identify two nodes u and v by keeping node u.

        All edges incident to v are moved to u, and v is removed.

        Parameters
        ----------
        u : T
            Node to keep.
        v : T
            Node to identify with u (will be removed).

        Examples
        --------
        >>> G = DirectedMultiGraph()
        >>> G.add_edge(1, 2)
        0
        >>> G.add_edge(2, 3)
        0
        >>> G.identify_two_nodes(1, 2)
        >>> list(G.nodes())
        [1, 3]
        """
        nx.contracted_nodes(self, u, v, self_loops=False, copy=False)
        self._clear_cache()

    def identify_node_set(self, nodes: List[T] | Set[T]) -> None:
        """
        Identify all nodes in the set by keeping the first node.

        Parameters
        ----------
        nodes : List[T] | Set[T]
            Iterable of nodes to identify. The first node will be kept.

        Examples
        --------
        >>> G = DirectedMultiGraph()
        >>> G.add_edge(1, 2)
        0
        >>> G.add_edge(2, 3)
        0
        >>> G.identify_node_set([1, 2, 3])
        >>> len(G.nodes()) <= 3
        True
        """
        nodes_list = list(nodes)
        if len(nodes_list) < 2:
            return

        for i in range(1, len(nodes_list)):
            self.identify_two_nodes(nodes_list[0], nodes_list[i])

    def clear(self) -> None:
        """
        Clear the whole directed multi-graph.

        Removes all nodes and edges.

        Examples
        --------
        >>> G = DirectedMultiGraph()
        >>> G.add_edge(1, 2)
        0
        >>> G.add_edge(2, 3)
        0
        >>> G.clear()
        >>> len(G.nodes())
        0
        """
        super().clear()
        self._clear_cache()
