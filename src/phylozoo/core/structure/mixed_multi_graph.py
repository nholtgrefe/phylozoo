"""
Mixed multi-graph module.

This module provides the MixedMultiGraph class for working with mixed multi-graphs.
"""

from typing import Any, List, Optional, Set, Tuple, TypeVar

import networkx as nx

T = TypeVar('T')


class MixedMultiGraph(nx.MultiGraph):
    """
    Mixed multi-graph with undirected and directed edges.

    This class allows parallel directed edges (multiple directed edges between
    the same nodes) but enforces that undirected edges are simple (no parallel
    undirected edges). Each parallel directed edge can have different parameters
    (weights, attributes, etc.) via edge keys.

    Parameters
    ----------
    directed_edges : Optional[List[Tuple[T, T]] | List[Tuple[T, T, int]]], optional
        List of directed edges. Can be (u, v) tuples or (u, v, key) tuples.
        If keys are not provided, they will be auto-generated. By default None.
    *args
        Additional arguments passed to nx.MultiGraph.
    **kwargs
        Additional keyword arguments passed to nx.MultiGraph.

    Attributes
    ----------
    _directed_edge_keys : Set[Tuple[T, T, int]]
        Set of (u, v, key) tuples identifying all directed edges.

    Examples
    --------
    >>> G = MixedMultiGraph()
    >>> G.add_directed_edge(1, 2)
    0
    >>> G.add_directed_edge(1, 2)  # Parallel directed edge
    1
    >>> G.add_edge(2, 3)  # Undirected edge
    >>> G.add_edge(2, 3)  # Replaces previous undirected edge
    """

    def __init__(
        self, directed_edges: Optional[List[Tuple[T, T] | Tuple[T, T, int]]] = None, *args: Any, **kwargs: Any
    ) -> None:
        """
        Initialize a mixed multi-graph.

        Parameters
        ----------
        directed_edges : Optional[List[Tuple[T, T] | Tuple[T, T, int]]], optional
            List of directed edges. Can be (u, v) or (u, v, key) tuples.
            By default None.
        *args
            Additional arguments passed to nx.MultiGraph.
        **kwargs
            Additional keyword arguments passed to nx.MultiGraph.
        """
        super().__init__(*args, **kwargs)
        self._directed_edge_keys: Set[Tuple[T, T, int]] = set()

        # Load directed edges if given
        if directed_edges:
            for edge in directed_edges:
                if len(edge) == 2:
                    u, v = edge
                    self.add_directed_edge(u, v)
                elif len(edge) == 3:
                    u, v, key = edge
                    self.add_directed_edge(u, v, key=key)
                else:
                    raise ValueError(f"Invalid edge format: {edge}")

    def _clear_cache(self) -> None:
        """
        Clear any cached computations.

        This method is a placeholder for subclasses to override if they
        maintain cached values that need to be cleared when the graph changes.
        """
        pass

    @classmethod
    def from_graph(cls, graph: nx.Graph, directed_edges: Optional[List[Tuple[T, T]]] = None) -> 'MixedMultiGraph':
        """
        Create a MixedMultiGraph from a NetworkX Graph.

        Parameters
        ----------
        graph : nx.Graph
            NetworkX graph to convert.
        directed_edges : Optional[List[Tuple[T, T]]], optional
            List of edges that should be directed. By default None.

        Returns
        -------
        MixedMultiGraph
            New MixedMultiGraph instance.

        Examples
        --------
        >>> G = nx.Graph()
        >>> G.add_edge(1, 2)
        >>> G.add_edge(2, 3)
        >>> M = MixedMultiGraph.from_graph(G, directed_edges=[(1, 2)])
        """
        mg = cls(incoming_graph_data=graph.edges, directed_edges=directed_edges)
        return mg

    def remove_node(self, v: T) -> None:
        """
        Remove node v from the graph.

        Also removes all edges incident to v, including both directed
        and undirected edges.

        Parameters
        ----------
        v : T
            Node to remove.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2)
        0
        >>> G.add_edge(2, 3)
        >>> G.remove_node(2)
        >>> list(G.nodes())
        [1, 3]
        """
        # Remove all directed edges incident to v
        edges_to_remove: List[Tuple[T, T, int]] = []
        for u, w, key in self._directed_edge_keys:
            if v == u or v == w:
                edges_to_remove.append((u, w, key))

        for u, w, key in edges_to_remove:
            self._directed_edge_keys.discard((u, w, key))

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
        >>> G = MixedMultiGraph()
        >>> G.add_node(1)
        >>> G.add_node(2)
        >>> G.add_node(3)
        >>> G.remove_nodes_from([1, 2])
        >>> list(G.nodes())
        [3]
        """
        for v in nodes:
            self.remove_node(v)

    def remove_edge(self, u: T, v: T, key: Optional[int] = None) -> None:
        """
        Remove edge (u, v) from the graph.

        If key is provided, removes the specific edge with that key.
        If key is None and the edge is directed, removes one directed edge.
        If key is None and the edge is undirected, removes the undirected edge.

        Parameters
        ----------
        u : T
            Source node.
        v : T
            Target node.
        key : Optional[int], optional
            Edge key. If None, removes one edge (directed or undirected).
            By default None.

        Raises
        ------
        KeyError
            If the specified edge does not exist.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> key = G.add_directed_edge(1, 2)
        >>> G.remove_edge(1, 2, key=key)
        """
        # Check if this is a directed edge
        if key is not None:
            if (u, v, key) in self._directed_edge_keys:
                self._directed_edge_keys.discard((u, v, key))
        else:
            # Try to find and remove a directed edge first
            if u in self and v in self[u]:
                # Check all keys for this edge
                for k in list(self[u][v].keys()):
                    if (u, v, k) in self._directed_edge_keys:
                        self._directed_edge_keys.discard((u, v, k))
                        key = k
                        break

        super().remove_edge(u, v, key)
        self._clear_cache()

    def remove_edges_from(self, edges: List[Tuple[T, T] | Tuple[T, T, int]]) -> None:
        """
        Remove all edges in 'edges' from the graph.

        Parameters
        ----------
        edges : List[Tuple[T, T] | Tuple[T, T, int]]
            List of edges to remove. Can be (u, v) or (u, v, key) tuples.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2)
        0
        >>> G.add_edge(2, 3)
        >>> G.remove_edges_from([(1, 2), (2, 3)])
        """
        for edge in edges:
            if len(edge) == 2:
                u, v = edge
                self.remove_edge(u, v)
            elif len(edge) == 3:
                u, v, key = edge
                self.remove_edge(u, v, key=key)
            else:
                raise ValueError(f"Invalid edge format: {edge}")

    def remove_directed_edge(self, u: T, v: T, key: Optional[int] = None) -> None:
        """
        Remove directed edge (u, v) from the graph.

        Parameters
        ----------
        u : T
            Source node.
        v : T
            Target node.
        key : Optional[int], optional
            Edge key. If None, removes one directed edge. By default None.

        Raises
        ------
        ValueError
            If no directed edge (u, v) exists with the specified key.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> key = G.add_directed_edge(1, 2)
        >>> G.remove_directed_edge(1, 2, key=key)
        """
        if key is None:
            # Find first directed edge
            found = False
            if u in self and v in self[u]:
                for k in list(self[u][v].keys()):
                    if (u, v, k) in self._directed_edge_keys:
                        key = k
                        found = True
                        break
            if not found:
                raise ValueError(f"No directed edge ({u}, {v}) exists.")
        elif (u, v, key) not in self._directed_edge_keys:
            raise ValueError(f"Directed edge ({u}, {v}, {key}) does not exist or has different direction.")

        self._directed_edge_keys.discard((u, v, key))
        super().remove_edge(u, v, key)
        self._clear_cache()

    def remove_directed_edges_from(self, edges: List[Tuple[T, T] | Tuple[T, T, int]]) -> None:
        """
        Remove all directed edges in 'edges' from the graph.

        Parameters
        ----------
        edges : List[Tuple[T, T] | Tuple[T, T, int]]
            List of directed edges to remove. Can be (u, v) or (u, v, key) tuples.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2)
        0
        >>> G.add_directed_edge(2, 3)
        0
        >>> G.remove_directed_edges_from([(1, 2), (2, 3)])
        """
        for edge in edges:
            if len(edge) == 2:
                u, v = edge
                self.remove_directed_edge(u, v)
            elif len(edge) == 3:
                u, v, key = edge
                self.remove_directed_edge(u, v, key=key)
            else:
                raise ValueError(f"Invalid edge format: {edge}")

    def add_directed_edge(self, u: T, v: T, key: Optional[int] = None, **attr: Any) -> int:
        """
        Add directed edge (u, v) to the graph.

        Parallel directed edges are allowed. Each parallel edge can have
        different attributes (weights, etc.) via the key parameter.

        Parameters
        ----------
        u : T
            Source node.
        v : T
            Target node.
        key : Optional[int], optional
            Edge key. If None, auto-generates a key. By default None.
        **attr
            Edge attributes (e.g., weight, label, etc.).

        Returns
        -------
        int
            The key of the added edge.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> key1 = G.add_directed_edge(1, 2, weight=1.0)
        >>> key2 = G.add_directed_edge(1, 2, weight=2.0)  # Parallel edge
        >>> key1 != key2
        True
        """
        # Ensure nodes exist
        if u not in self:
            self.add_node(u)
        if v not in self:
            self.add_node(v)

        # Auto-generate key if not provided
        if key is None:
            # Find next available key for (u, v)
            if u in self and v in self[u]:
                existing_keys = set(self[u][v].keys())
                key = max(existing_keys) + 1 if existing_keys else 0
            else:
                key = 0

        # Add edge to underlying MultiGraph
        super().add_edge(u, v, key=key, **attr)

        # Track as directed edge
        self._directed_edge_keys.add((u, v, key))

        # If there was an undirected edge, remove it (directed takes precedence)
        # Check if there's an undirected edge (not in _directed_edge_keys)
        if u in self and v in self[u]:
            for k in list(self[u][v].keys()):
                if (u, v, k) not in self._directed_edge_keys and (v, u, k) not in self._directed_edge_keys:
                    # This is an undirected edge, remove it
                    super().remove_edge(u, v, key=k)
                    break

        self._clear_cache()
        return key

    def add_directed_edges_from(self, edges: List[Tuple[T, T] | Tuple[T, T, int]], **attr: Any) -> None:
        """
        Add all directed edges in 'edges' to the graph.

        Parameters
        ----------
        edges : List[Tuple[T, T] | Tuple[T, T, int]]
            List of directed edges. Can be (u, v) or (u, v, key) tuples.
        **attr
            Edge attributes applied to all edges.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edges_from([(1, 2), (2, 3), (3, 1)])
        """
        for edge in edges:
            if len(edge) == 2:
                u, v = edge
                self.add_directed_edge(u, v, **attr)
            elif len(edge) == 3:
                u, v, key = edge
                self.add_directed_edge(u, v, key=key, **attr)
            else:
                raise ValueError(f"Invalid edge format: {edge}")

    def add_edge(self, u: T, v: T, key: Optional[int] = None, **attr: Any) -> int:
        """
        Add undirected edge (u, v) to the graph.

        Undirected edges cannot be parallel. If an undirected edge already
        exists between u and v, it will be replaced.

        Parameters
        ----------
        u : T
            First node.
        v : T
            Second node.
        key : Optional[int], optional
            Edge key. If None, uses key=0. By default None.
        **attr
            Edge attributes.

        Returns
        -------
        int
            The key of the added edge (always 0 for undirected edges).

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_edge(1, 2)
        0
        >>> G.add_edge(1, 2)  # Replaces previous undirected edge
        0
        """
        # Ensure nodes exist
        if u not in self:
            self.add_node(u)
        if v not in self:
            self.add_node(v)

        # Remove any existing undirected edge between u and v
        # (undirected edges cannot be parallel)
        if u in self and v in self[u]:
            for k in list(self[u][v].keys()):
                if (u, v, k) not in self._directed_edge_keys and (v, u, k) not in self._directed_edge_keys:
                    # This is an undirected edge, remove it
                    super().remove_edge(u, v, key=k)
                    break

        # Use key=0 for undirected edges (enforce no parallel)
        if key is None:
            key = 0

        # Add edge to underlying MultiGraph
        super().add_edge(u, v, key=key, **attr)

        # Do NOT add to _directed_edge_keys (it's undirected)

        self._clear_cache()
        return key

    def undirected_edges(self, keys: bool = False) -> List[Tuple[T, T] | Tuple[T, T, int]]:
        """
        Return all undirected edges of the graph.

        Parameters
        ----------
        keys : bool, optional
            If True, return (u, v, key) tuples. If False, return (u, v) tuples.
            By default False.

        Returns
        -------
        List[Tuple[T, T] | Tuple[T, T, int]]
            List of undirected edges.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_edge(1, 2)
        0
        >>> G.add_directed_edge(2, 3)
        0
        >>> G.undirected_edges()
        [(1, 2)]
        """
        und: List[Tuple[T, T] | Tuple[T, T, int]] = []
        for u, v, k in self.edges(keys=True):
            if (u, v, k) not in self._directed_edge_keys and (v, u, k) not in self._directed_edge_keys:
                if keys:
                    und.append((u, v, k))
                else:
                    und.append((u, v))
        return und

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
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2)
        0
        >>> G.add_directed_edge(3, 2)
        0
        >>> G.indegree(2)
        2
        """
        return len([(u, w, k) for u, w, k in self._directed_edge_keys if w == v])

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
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2)
        0
        >>> G.add_directed_edge(1, 3)
        0
        >>> G.outdegree(1)
        2
        """
        return len([(u, w, k) for u, w, k in self._directed_edge_keys if u == v])

    def is_cutedge(self, u: T, v: T, key: Optional[int] = None) -> bool:
        """
        Check if edge (u, v) is a cut-edge.

        A cut-edge is an edge whose removal increases the number of
        connected components.

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
        >>> G = MixedMultiGraph()
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
        return nx.number_connected_components(modified_graph) != nx.number_connected_components(self)

    def is_cutvertex(self, v: T) -> bool:
        """
        Check if v is a cut-vertex.

        A cut-vertex is a vertex whose removal increases the number of
        connected components.

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
        >>> G = MixedMultiGraph()
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
        return nx.number_connected_components(modified_graph) != nx.number_connected_components(self)

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
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2)
        0
        >>> G.add_edge(2, 3)
        0
        >>> G.identify_two_nodes(1, 2)
        >>> list(G.nodes())
        [1, 3]
        """
        nx.contracted_nodes(self, u, v, self_loops=False, copy=False)

        # Remove directed edges that became self-loops or were between u and v
        edges_to_remove: List[Tuple[T, T, int]] = []
        for x, y, k in list(self._directed_edge_keys):
            if x == v or y == v:
                edges_to_remove.append((x, y, k))
            elif (x == u and y == u) or (x == v and y == u) or (x == u and y == v):
                edges_to_remove.append((x, y, k))

        for x, y, k in edges_to_remove:
            self._directed_edge_keys.discard((x, y, k))
            if (x, y, k) in list(self.edges(keys=True)):
                try:
                    super().remove_edge(x, y, key=k)
                except (KeyError, nx.NetworkXError):
                    pass

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
        >>> G = MixedMultiGraph()
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
        Clear the whole mixed graph.

        Removes all nodes and edges.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2)
        0
        >>> G.add_edge(2, 3)
        0
        >>> G.clear()
        >>> len(G.nodes())
        0
        """
        super().clear()
        self._directed_edge_keys.clear()
        self._clear_cache()
