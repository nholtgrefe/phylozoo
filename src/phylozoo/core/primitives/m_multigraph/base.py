"""
Mixed multi-graph module.

This module provides the MixedMultiGraph class for working with mixed multi-graphs.
"""

from typing import Any, Dict, Iterator, List, Set, Tuple, TypeVar, TYPE_CHECKING

import networkx as nx

from phylozoo.utils.identifier_warnings import warn_on_keyword, warn_on_none_value

if TYPE_CHECKING:
    from ..d_multigraph import DirectedMultiGraph

T = TypeVar('T')


class MixedMultiGraph:
    """
    Mixed multi-graph with undirected and directed edges.

    This class uses composition with separate NetworkX graphs:
    - _undirected: nx.MultiGraph for undirected edges
    - _directed: nx.MultiDiGraph for directed edges
    - _combined: nx.MultiGraph combining all edges for connectivity analysis

    This class allows parallel edges for both directed and undirected edges,
    including self-loops. However, edges between the same two nodes must be either
    all directed or all undirected - mixing is not allowed. This mutual exclusivity is enforced
    automatically: adding a directed edge will remove any undirected edges between
    the same nodes, and vice versa. Each parallel edge can have different
    parameters (weights, attributes, etc.) via edge keys.

    Parameters
    ----------
    directed_edges : list[tuple[T, T] | tuple[T, T, int] | dict[str, Any]] | None, optional
        List of directed edges. Can be:
        - (u, v) tuples (key auto-generated)
        - (u, v, key) tuples (explicit key)
        - Dict with 'u', 'v' keys and optional 'key' and edge attributes
        If keys are not provided, they will be auto-generated. By default None.
    undirected_edges : list[tuple[T, T] | tuple[T, T, int] | dict[str, Any]] | None, optional
        List of undirected edges. Can be:
        - (u, v) tuples (key auto-generated)
        - (u, v, key) tuples (explicit key)
        - Dict with 'u', 'v' keys and optional 'key' and edge attributes
        If keys are not provided, they will be auto-generated. By default None.

    Attributes
    ----------
    _undirected : nx.MultiGraph
        NetworkX MultiGraph storing undirected edges.
        **Warning:** Do not modify directly. Use class methods instead.
    _directed : nx.MultiDiGraph
        NetworkX MultiDiGraph storing directed edges.
        **Warning:** Do not modify directly. Use class methods instead.
    _combined : nx.MultiGraph
        Combined undirected view of all edges for connectivity analysis.
        **Warning:** Do not modify directly. Use class methods instead.
    
    Notes
    -----
    The underlying graphs (`_undirected`, `_directed`, `_combined`) are accessible
    but should NOT be modified directly. All modifications must go through the
    class methods (add_edge, remove_edge, etc.) to ensure state synchronization.
    Direct modification will desynchronize the graphs and cause incorrect behavior.

    Examples
    --------
    >>> G = MixedMultiGraph()
    >>> G.add_directed_edge(1, 2)
    0
    >>> G.add_directed_edge(1, 2)  # Parallel directed edge
    1
    >>> G.add_undirected_edge(2, 3)  # Undirected edge
    0
    >>> G.add_undirected_edge(2, 3)  # Parallel undirected edge
    1
    >>> from phylozoo.core.primitives.m_multigraph.features import number_of_connected_components
    >>> number_of_connected_components(G)
    1
    >>> # Initialize with edges (including attributes)
    >>> G2 = MixedMultiGraph(
    ...     undirected_edges=[(1, 2), {'u': 2, 'v': 3, 'weight': 5.0}],
    ...     directed_edges=[(3, 4), {'u': 4, 'v': 5, 'weight': 10.0, 'label': 'test'}]
    ... )
    >>> G2.number_of_edges()
    4
    >>> # Create from NetworkX graphs
    >>> import networkx as nx
    >>> from phylozoo.core.primitives.m_multigraph.conversions import graph_to_mixedmultigraph
    >>> nx_g = nx.Graph()
    >>> nx_g.add_edge(1, 2, weight=1.0)
    >>> G3 = graph_to_mixedmultigraph(nx_g)
    >>> G3.number_of_edges()
    1
    """

    def __init__(
        self,
        directed_edges: list[tuple[T, T] | tuple[T, T, int] | dict[str, Any]] | None = None,
        undirected_edges: list[tuple[T, T] | tuple[T, T, int] | dict[str, Any]] | None = None,
    ) -> None:
        """
        Initialize a mixed multi-graph.

        Parameters
        ----------
        directed_edges : list[tuple[T, T] | tuple[T, T, int] | dict[str, Any]] | None, optional
            List of directed edges. Can be:
            - (u, v) tuples (key auto-generated)
            - (u, v, key) tuples (explicit key)
            - Dict with 'u', 'v' keys and optional 'key' and edge attributes
            By default None.
        undirected_edges : list[tuple[T, T] | tuple[T, T, int] | dict[str, Any]] | None, optional
            List of undirected edges. Can be:
            - (u, v) tuples (key auto-generated)
            - (u, v, key) tuples (explicit key)
            - Dict with 'u', 'v' keys and optional 'key' and edge attributes
            By default None.

        Examples
        --------
        >>> G = MixedMultiGraph(undirected_edges=[(1, 2), (2, 3)])
        >>> G2 = MixedMultiGraph(directed_edges=[(1, 2, 0), {'u': 2, 'v': 3, 'weight': 5.0}])
        """
        self._undirected: nx.MultiGraph = nx.MultiGraph()
        self._directed: nx.MultiDiGraph = nx.MultiDiGraph()
        self._combined: nx.MultiGraph = nx.MultiGraph()

        # Load undirected edges if given (before directed edges to handle mutual exclusivity)
        if undirected_edges:
            for edge in undirected_edges:
                if isinstance(edge, dict):
                    # Dict format: {'u': u, 'v': v, 'key': key, **attr}
                    u = edge.pop('u')
                    v = edge.pop('v')
                    key = edge.pop('key', None)
                    self.add_undirected_edge(u, v, key=key, **edge)
                elif len(edge) == 2:
                    u, v = edge
                    self.add_undirected_edge(u, v)
                elif len(edge) == 3:
                    u, v, key = edge
                    self.add_undirected_edge(u, v, key=key)
                else:
                    raise ValueError(f"Invalid edge format: {edge}")

        # Load directed edges if given
        if directed_edges:
            for edge in directed_edges:
                if isinstance(edge, dict):
                    # Dict format: {'u': u, 'v': v, 'key': key, **attr}
                    u = edge.pop('u')
                    v = edge.pop('v')
                    key = edge.pop('key', None)
                    self.add_directed_edge(u, v, key=key, **edge)
                elif len(edge) == 2:
                    u, v = edge
                    self.add_directed_edge(u, v)
                elif len(edge) == 3:
                    u, v, key = edge
                    self.add_directed_edge(u, v, key=key)
                else:
                    raise ValueError(f"Invalid edge format: {edge}")

    # ========== NetworkX Compatibility Methods ==========

    def nodes_iter(self, data: bool | str = False) -> Iterator[T] | Iterator[tuple[T, Any]]:
        """
        Return an iterator over nodes.

        Parameters
        ----------
        data : bool | str, optional
            If False (default), return iterator over nodes.
            If True, return iterator of (node, data_dict) tuples.
            If string, return iterator of (node, attribute_value) tuples
            for the given attribute name.

        Returns
        -------
        Iterator[T] | Iterator[tuple[T, Any]]
            Iterator over nodes or (node, data) tuples.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_node(1, weight=2.0)
        >>> G.add_node(2, weight=3.0)
        >>> list(G.nodes_iter())
        [1, 2]
        >>> list(G.nodes_iter(data=True))
        [(1, {'weight': 2.0}), (2, {'weight': 3.0})]
        >>> list(G.nodes_iter(data='weight'))
        [(1, 2.0), (2, 3.0)]
        """
        # Use the combined graph for consistency - it has all nodes
        return self._combined.nodes(data=data)

    def edges_iter(self, keys: bool = False, data: bool | str = False) -> Iterator[tuple[T, T] | tuple[T, T, int] | tuple[T, T, Any] | tuple[T, T, dict[str, Any]] | tuple[T, T, int, Any] | tuple[T, T, int, dict[str, Any]]]:
        """
        Return an iterator over edges.

        Parameters
        ----------
        keys : bool, optional
            If True, return edge keys. By default False.
        data : bool | str, optional
            If False (default), no edge data is included.
            If True, return edge data dictionaries.
            If string, return value of that edge attribute.

        Returns
        -------
        Iterator
            Iterator over edges. Format depends on keys and data parameters.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_undirected_edge(1, 2, weight=1.0)
        0
        >>> list(G.edges_iter())
        [(1, 2)]
        >>> list(G.edges_iter(keys=True))
        [(1, 2, 0)]
        >>> list(G.edges_iter(data=True))
        [(1, 2, {'weight': 1.0})]
        >>> list(G.edges_iter(keys=True, data='weight'))
        [(1, 2, 0, 1.0)]
        """
        # Return edges from combined graph (has all edges)
        return self._combined.edges(keys=keys, data=data)
    
    def undirected_edges_iter(self, keys: bool = False, data: bool | str = False) -> Iterator[tuple[T, T] | tuple[T, T, int] | tuple[T, T, Any] | tuple[T, T, dict[str, Any]] | tuple[T, T, int, Any] | tuple[T, T, int, dict[str, Any]]]:
        """
        Return an iterator over undirected edges.

        Parameters
        ----------
        keys : bool, optional
            If True, return edge keys. By default False.
        data : bool | str, optional
            If False (default), no edge data is included.
            If True, return edge data dictionaries.
            If string, return value of that edge attribute.

        Returns
        -------
        Iterator
            Iterator over undirected edges. Format depends on keys and data parameters.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_undirected_edge(1, 2, weight=1.0)
        0
        >>> list(G.undirected_edges_iter())
        [(1, 2)]
        >>> list(G.undirected_edges_iter(keys=True))
        [(1, 2, 0)]
        >>> list(G.undirected_edges_iter(data='weight'))
        [(1, 2, 1.0)]
        """
        return self._undirected.edges(keys=keys, data=data)
    
    def directed_edges_iter(self, keys: bool = False, data: bool | str = False) -> Iterator[tuple[T, T] | tuple[T, T, int] | tuple[T, T, Any] | tuple[T, T, dict[str, Any]] | tuple[T, T, int, Any] | tuple[T, T, int, dict[str, Any]]]:
        """
        Return an iterator over directed edges.

        Parameters
        ----------
        keys : bool, optional
            If True, return edge keys. By default False.
        data : bool | str, optional
            If False (default), no edge data is included.
            If True, return edge data dictionaries.
            If string, return value of that edge attribute.

        Returns
        -------
        Iterator
            Iterator over directed edges. Format depends on keys and data parameters.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2, weight=1.0)
        0
        >>> list(G.directed_edges_iter())
        [(1, 2)]
        >>> list(G.directed_edges_iter(keys=True))
        [(1, 2, 0)]
        >>> list(G.directed_edges_iter(data='weight'))
        [(1, 2, 1.0)]
        """
        return self._directed.edges(keys=keys, data=data)
    

    def neighbors(self, v: T) -> Iterator[T]:
        """
        Return an iterator over neighbors of node v.

        Parameters
        ----------
        v : T
            Node.

        Returns
        -------
        Iterator[T]
            Iterator over neighbors.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_undirected_edge(1, 2)
        0
        >>> G.add_directed_edge(1, 3)
        0
        >>> list(G.neighbors(1))
        [2, 3]
        """
        neighbors_set: set[T] = set()
        if v in self._undirected:
            neighbors_set.update(self._undirected.neighbors(v))
        if v in self._directed:
            neighbors_set.update(self._directed.predecessors(v))
            neighbors_set.update(self._directed.successors(v))
        return iter(neighbors_set)
    
    def incident_parent_edges(self, v: T, keys: bool = False, data: bool | str = False) -> Iterator[tuple[T, T] | tuple[T, T, int] | tuple[T, T, Any] | tuple[T, T, dict[str, Any]] | tuple[T, T, int, Any] | tuple[T, T, int, dict[str, Any]]]:
        """
        Return an iterator over directed edges entering node v (from parent nodes).
        
        Parameters
        ----------
        v : T
            Node.
        keys : bool, optional
            If True, return edge keys. By default False.
        data : bool | str, optional
            If False (default), no edge data is included.
            If True, return edge data dictionaries.
            If string, return value of that edge attribute.
        
        Returns
        -------
        Iterator
            Iterator over incoming edges. Format depends on keys and data parameters.
        
        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2, weight=1.0)
        0
        >>> G.add_directed_edge(3, 2, weight=2.0)
        0
        >>> list(G.incident_parent_edges(2))
        [(1, 2), (3, 2)]
        >>> list(G.incident_parent_edges(2, keys=True, data=True))
        [(1, 2, 0, {'weight': 1.0}), (3, 2, 0, {'weight': 2.0})]
        >>> list(G.incident_parent_edges(2, data='weight'))
        [(1, 2, 1.0), (3, 2, 2.0)]
        """
        return self._directed.in_edges(v, keys=keys, data=data)
    
    def incident_child_edges(self, v: T, keys: bool = False, data: bool | str = False) -> Iterator[tuple[T, T] | tuple[T, T, int] | tuple[T, T, Any] | tuple[T, T, dict[str, Any]] | tuple[T, T, int, Any] | tuple[T, T, int, dict[str, Any]]]:
        """
        Return an iterator over directed edges leaving node v (to child nodes).
        
        Parameters
        ----------
        v : T
            Node.
        keys : bool, optional
            If True, return edge keys. By default False.
        data : bool | str, optional
            If False (default), no edge data is included.
            If True, return edge data dictionaries.
            If string, return value of that edge attribute.
        
        Returns
        -------
        Iterator
            Iterator over outgoing edges. Format depends on keys and data parameters.
        
        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2, weight=1.0)
        0
        >>> G.add_directed_edge(1, 3, weight=2.0)
        0
        >>> list(G.incident_child_edges(1))
        [(1, 2), (1, 3)]
        >>> list(G.incident_child_edges(1, keys=True, data=True))
        [(1, 2, 0, {'weight': 1.0}), (1, 3, 0, {'weight': 2.0})]
        >>> list(G.incident_child_edges(1, data='weight'))
        [(1, 2, 1.0), (1, 3, 2.0)]
        """
        return self._directed.out_edges(v, keys=keys, data=data)
    
    def incident_undirected_edges(self, v: T, keys: bool = False, data: bool | str = False) -> Iterator[tuple[T, T] | tuple[T, T, int] | tuple[T, T, Any] | tuple[T, T, dict[str, Any]] | tuple[T, T, int, Any] | tuple[T, T, int, dict[str, Any]]]:
        """
        Return an iterator over undirected edges incident to node v.
        
        Parameters
        ----------
        v : T
            Node.
        keys : bool, optional
            If True, return edge keys. By default False.
        data : bool | str, optional
            If False (default), no edge data is included.
            If True, return edge data dictionaries.
            If string, return value of that edge attribute.
        
        Returns
        -------
        Iterator
            Iterator over incident undirected edges. Format depends on keys and data parameters.
        
        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_undirected_edge(1, 2, weight=1.0)
        0
        >>> G.add_undirected_edge(2, 3, weight=2.0)
        0
        >>> list(G.incident_undirected_edges(2))
        [(1, 2), (2, 3)]
        >>> list(G.incident_undirected_edges(2, keys=True, data=True))
        [(1, 2, 0, {'weight': 1.0}), (2, 3, 0, {'weight': 2.0})]
        >>> list(G.incident_undirected_edges(2, data='weight'))
        [(1, 2, 1.0), (2, 3, 2.0)]
        """
        return self._undirected.edges(v, keys=keys, data=data)

    def __contains__(self, v: T) -> bool:
        """
        Check if node v is in the graph.

        Parameters
        ----------
        v : T
            Node to check.

        Returns
        -------
        bool
            True if v is in the graph, False otherwise.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_node(1)
        >>> 1 in G
        True
        >>> 2 in G
        False
        """
        return v in self._undirected or v in self._directed

    def __iter__(self) -> Iterator[T]:
        """
        Iterate over nodes.

        Returns
        -------
        Iterator[T]
            Iterator over nodes.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_nodes_from([1, 2, 3])
        >>> list(G)
        [1, 2, 3]
        """
        return iter(self.nodes)

    def __len__(self) -> int:
        """
        Return the number of nodes.

        Returns
        -------
        int
            Number of nodes.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_nodes_from([1, 2, 3])
        >>> len(G)
        3
        """
        return len(set(self._undirected.nodes()) | set(self._directed.nodes()))

    def __getitem__(self, v: T) -> dict[T, dict[int, dict[str, Any]]]:
        """
        Return adjacency dict for node v.

        Parameters
        ----------
        v : T
            Node.

        Returns
        -------
        dict[T, dict[int, dict[str, Any]]]
            Adjacency dict (actually returns NetworkX's AdjacencyView, which is
            dict-like and supports all dict operations).

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_undirected_edge(1, 2, weight=1.0)
        0
        >>> G[1]
        {2: {0: {'weight': 1.0}}}
        """
        # Return combined adjacency (has all edges)
        return self._combined[v]

    # ========== Properties for Attribute Access ==========

    @property
    def combined_graph(self) -> nx.MultiGraph:
        """
        Get combined undirected graph for NetworkX algorithms.

        Returns
        -------
        nx.MultiGraph
            Combined graph treating all edges as undirected.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2)
        0
        >>> G.add_undirected_edge(2, 3)
        0
        >>> nx.is_connected(G.combined_graph)
        True
        """
        return self._combined

    class NodeView:
        """
        Node view that works as both attribute and method, similar to NetworkX's NodeView.
        
        This class provides a set-like interface for nodes while also being callable
        as a method to get iterators or node data.
        """
        def __init__(self, items: set[T], callable_func: callable):
            """
            Initialize a node view.
            
            Parameters
            ----------
            items : set[T]
                Set of nodes.
            callable_func : callable
                Function to call when used as method.
            """
            self._items = items
            self._callable_func = callable_func
        
        def __call__(self, data: bool | str = False):
            """
            Call as method to get iterator or node data.
            
            Parameters
            ----------
            data : bool | str, optional
                If False (default), return iterator over nodes.
                If True, return iterator of (node, data_dict) tuples.
                If string, return iterator of (node, attribute_value) tuples.
            
            Returns
            -------
            Iterator[T] | Iterator[tuple[T, Any]]
                Iterator over nodes or (node, data) tuples.
            """
            return self._callable_func(data)
        
        def __iter__(self):
            """Iterate over nodes."""
            return iter(self._items)
        
        def __contains__(self, item: T) -> bool:
            """Check if node in view."""
            return item in self._items
        
        def __repr__(self) -> str:
            """String representation."""
            return repr(self._items)
        
        def __len__(self) -> int:
            """Number of nodes."""
            return len(self._items)
        
        def __or__(self, other):
            """Union with other set."""
            return self._items | other
        
        def __and__(self, other):
            """Intersection with other set."""
            return self._items & other
        
        def issubset(self, other):
            """Check if this is a subset of other."""
            return self._items.issubset(other)
    
    class EdgeView:
        """
        Edge view that works as both attribute and method, similar to NetworkX's EdgeView.
        
        This class provides a list-like interface for edges while also being callable
        as a method to get iterators with keys or data.
        """
        def __init__(self, items: list[tuple[T, T]], callable_func: callable):
            """
            Initialize an edge view.
            
            Parameters
            ----------
            items : list[tuple[T, T]]
                List of edges.
            callable_func : callable
                Function to call when used as method.
            """
            self._items = items
            self._callable_func = callable_func
        
        def __call__(self, keys: bool = False, data: bool | str = False):
            """
            Call as method to get iterator with keys or data.
            
            Parameters
            ----------
            keys : bool, optional
                If True, return edge keys. By default False.
            data : bool | str, optional
                If False (default), no edge data is included.
                If True, return edge data dictionaries.
                If string, return value of that edge attribute.
            
            Returns
            -------
            Iterator
                Iterator over edges. Format depends on keys and data parameters.
            """
            return self._callable_func(keys, data)
        
        def __iter__(self):
            """Iterate over edges."""
            return iter(self._items)
        
        def __contains__(self, item: tuple[T, T]) -> bool:
            """Check if edge in view."""
            return item in self._items
        
        def __repr__(self) -> str:
            """String representation."""
            return repr(self._items)
        
        def __len__(self) -> int:
            """Number of edges."""
            return len(self._items)
    
    @property
    def nodes(self) -> 'NodeView':
        """
        Get all nodes (works as both attribute and method).

        When accessed as attribute, returns a NodeView object (set-like).
        When called as method, returns an iterator.

        Returns
        -------
        NodeView
            Node view object that's both callable and set-like.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_nodes_from([1, 2, 3])
        >>> G.nodes  # Attribute access
        {1, 2, 3}
        >>> list(G.nodes())  # Method call
        [1, 2, 3]
        >>> list(G.nodes(data=True))  # Method call with data
        [(1, {}), (2, {}), (3, {})]
        """
        nodes_set = set(self._undirected.nodes()) | set(self._directed.nodes())
        return self.NodeView(nodes_set, self.nodes_iter)

    @property
    def edges(self) -> 'EdgeView':
        """
        Get all edges (works as both attribute and method).

        When accessed as attribute, returns an EdgeView object (list-like).
        When called as method, returns an iterator.

        Returns
        -------
        EdgeView
            Edge view object that's both callable and list-like.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_undirected_edge(1, 2)
        0
        >>> G.add_directed_edge(2, 3)
        0
        >>> G.edges  # Attribute access
        [(1, 2), (2, 3)]
        >>> list(G.edges())  # Method call
        [(1, 2), (2, 3)]
        """
        edges_list = list(self._undirected.edges())
        edges_list.extend(self._directed.edges())
        return self.EdgeView(edges_list, self.edges_iter)

    # ========== Node Operations ==========

    def add_node(self, v: T, **attr: Any) -> None:
        """
        Add node v to the graph.

        Parameters
        ----------
        v : T
            Node to add.
        **attr
            Node attributes.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_node(1, weight=2.0)
        >>> 1 in G
        True
        """
        # Warn on Python keyword identifiers
        warn_on_keyword(v, "Node id")
        # Warn on Python keyword attribute names and None values
        for attr_name, attr_value in attr.items():
            warn_on_keyword(attr_name, "Attribute name")
            warn_on_none_value(attr_value, f"Attribute '{attr_name}'")
        
        self._undirected.add_node(v, **attr)
        self._directed.add_node(v, **attr)
        self._combined.add_node(v, **attr)

    def add_nodes_from(self, nodes: list[T] | set[T], **attr: Any) -> None:
        """
        Add nodes from iterable.

        Parameters
        ----------
        nodes : list[T] | set[T]
            Iterable of nodes.
        **attr
            Attributes to add to all nodes.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_nodes_from([1, 2, 3])
        >>> len(G)
        3
        """
        self._undirected.add_nodes_from(nodes, **attr)
        self._directed.add_nodes_from(nodes, **attr)
        self._combined.add_nodes_from(nodes, **attr)

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
        >>> G.add_undirected_edge(2, 3)
        0
        >>> G.remove_node(2)
        >>> list(G.nodes())
        [1, 3]
        """
        self._undirected.remove_node(v)
        self._directed.remove_node(v)
        self._combined.remove_node(v)

    def remove_nodes_from(self, nodes: list[T] | set[T]) -> None:
        """
        Remove all nodes in 'nodes' from the graph.

        Parameters
        ----------
        nodes : list[T] | set[T]
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

    # ========== Directed Edge Operations ==========

    def add_directed_edge(self, u: T, v: T, key: int | None = None, **attr: Any) -> int:
        """
        Add directed edge (u, v) to the graph.

        Parallel directed edges are allowed. Each parallel edge can have
        different attributes (weights, etc.) via the key parameter.
        
        **Mutual Exclusivity:** If there are any undirected edges between u and v,
        they will be automatically removed. Edges between the same two nodes must
        be either all directed or all undirected - mixing is not allowed. This
        ensures consistent graph semantics.

        Parameters
        ----------
        u : T
            Source node.
        v : T
            Target node.
        key : int | None, optional
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
        >>> G.add_undirected_edge(1, 2)  # This removes the directed edges
        >>> G.add_directed_edge(1, 2)  # This removes the undirected edge
        """
        # Warn on Python keyword identifiers
        warn_on_keyword(u, "Node id")
        warn_on_keyword(v, "Node id")
        if key is not None:
            warn_on_keyword(key, "Edge key")
        # Warn on Python keyword attribute names and None values
        for attr_name, attr_value in attr.items():
            warn_on_keyword(attr_name, "Attribute name")
            warn_on_none_value(attr_value, f"Attribute '{attr_name}'")

        # Ensure nodes exist
        if u not in self._undirected:
            self._undirected.add_node(u)
            self._combined.add_node(u)
        if v not in self._undirected:
            self._undirected.add_node(v)
            self._combined.add_node(v)

        # Auto-generate key if not provided
        if key is None:
            if u in self._directed and v in self._directed[u]:
                existing_keys = set(self._directed[u][v].keys())
                key = max(existing_keys) + 1 if existing_keys else 0
            else:
                key = 0
        # If key is provided, use it as-is (no else needed)

        # Remove any undirected edges between u and v FIRST (mutual exclusivity)
        # This must happen BEFORE adding to _combined to avoid key conflicts
        # when both edge types use the same key value (e.g., key=0)
        if self._undirected.has_edge(u, v):
            for k in list(self._undirected[u][v].keys()):
                self._undirected.remove_edge(u, v, key=k)
                self._combined.remove_edge(u, v, key=k)

        # Now add to directed graph and combined graph
        self._directed.add_edge(u, v, key=key, **attr)
        self._combined.add_edge(u, v, key=key, **attr)

        return key

    def add_directed_edges_from(
        self, edges: list[tuple[T, T] | tuple[T, T, int]], **attr: Any
    ) -> None:
        """
        Add all directed edges in 'edges' to the graph.

        Parameters
        ----------
        edges : list[tuple[T, T] | tuple[T, T, int]]
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

    def remove_directed_edge(self, u: T, v: T, key: int | None = None) -> None:
        """
        Remove directed edge (u, v) from the graph.

        Parameters
        ----------
        u : T
            Source node.
        v : T
            Target node.
        key : int | None, optional
            Edge key. If None, removes one directed edge. By default None.

        Raises
        ------
        ValueError
            If no directed edge (u, v) exists with the specified key.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> key = G.add_directed_edge(1, 2)
        0
        >>> G.remove_directed_edge(1, 2, key=key)
        """
        if key is None:
            if not self._directed.has_edge(u, v):
                raise ValueError(f"No directed edge ({u}, {v}) exists.")
            key = next(iter(self._directed[u][v].keys()))

        if not self._directed.has_edge(u, v, key):
            raise ValueError(f"Directed edge ({u}, {v}, {key}) does not exist.")

        self._directed.remove_edge(u, v, key)
        self._combined.remove_edge(u, v, key)

    def remove_directed_edges_from(
        self, edges: list[tuple[T, T] | tuple[T, T, int]]
    ) -> None:
        """
        Remove all directed edges in 'edges' from the graph.

        Parameters
        ----------
        edges : list[tuple[T, T] | tuple[T, T, int]]
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

    @property
    def directed_edges(self) -> 'EdgeView':
        """
        Get all directed edges (works as both attribute and method).

        When accessed as attribute, returns an EdgeView object (list-like).
        When called as method, returns an iterator with optional keys and data.

        Returns
        -------
        EdgeView
            Edge view object that's both callable and list-like.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2, weight=1.0)
        0
        >>> G.add_directed_edge(1, 2, weight=2.0)  # Parallel edge
        1
        >>> G.add_undirected_edge(2, 3)
        0
        >>> G.directed_edges  # Attribute access
        [(1, 2), (1, 2)]
        >>> list(G.directed_edges())  # Method call
        [(1, 2), (1, 2)]
        >>> list(G.directed_edges(keys=True))  # With keys
        [(1, 2, 0), (1, 2, 1)]
        >>> list(G.directed_edges(data=True))  # With data
        [(1, 2, {'weight': 1.0}), (1, 2, {'weight': 2.0})]
        >>> list(G.directed_edges(keys=True, data=True))  # With keys and data
        [(1, 2, 0, {'weight': 1.0}), (1, 2, 1, {'weight': 2.0})]
        """
        edges_list = list(self._directed.edges())
        return self.EdgeView(edges_list, self.directed_edges_iter)

    # ========== Undirected Edge Operations ==========

    def add_undirected_edge(self, u: T, v: T, key: int | None = None, **attr: Any) -> int:
        """
        Add undirected edge (u, v) to the graph.

        Parallel undirected edges are allowed. Each parallel edge can have
        different attributes (weights, etc.) via the key parameter.
        
        **Mutual Exclusivity:** If there are any directed edges between u and v,
        they will be automatically removed. Edges between the same two nodes must
        be either all directed or all undirected - mixing is not allowed. This
        ensures consistent graph semantics.

        Parameters
        ----------
        u : T
            First node.
        v : T
            Second node.
        key : int | None, optional
            Edge key. If None, auto-generates a key. By default None.
        **attr
            Edge attributes.

        Returns
        -------
        int
            The key of the added edge.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> key1 = G.add_undirected_edge(1, 2, weight=1.0)
        >>> key2 = G.add_undirected_edge(1, 2, weight=2.0)  # Parallel edge
        >>> key1 != key2
        True
        >>> G.add_directed_edge(1, 2)  # This removes the undirected edges
        >>> G.add_undirected_edge(1, 2)  # This removes the directed edge
        """
        # Warn on Python keyword identifiers
        warn_on_keyword(u, "Node id")
        warn_on_keyword(v, "Node id")
        if key is not None:
            warn_on_keyword(key, "Edge key")
        # Warn on Python keyword attribute names and None values
        for attr_name, attr_value in attr.items():
            warn_on_keyword(attr_name, "Attribute name")
            warn_on_none_value(attr_value, f"Attribute '{attr_name}'")

        # Ensure nodes exist
        if u not in self._directed:
            self._directed.add_node(u)
            self._combined.add_node(u)
        if v not in self._directed:
            self._directed.add_node(v)
            self._combined.add_node(v)

        # Remove any directed edges between u and v (mutual exclusivity)
        if self._directed.has_edge(u, v):
            for k in list(self._directed[u][v].keys()):
                self._directed.remove_edge(u, v, key=k)
                self._combined.remove_edge(u, v, key=k)

        # Auto-generate key if not provided
        if key is None:
            if u in self._undirected and v in self._undirected[u]:
                existing_keys = set(self._undirected[u][v].keys())
                key = max(existing_keys) + 1 if existing_keys else 0
            else:
                key = 0

        # Add to both undirected and combined graphs
        self._undirected.add_edge(u, v, key=key, **attr)
        self._combined.add_edge(u, v, key=key, **attr)

        return key

    def add_undirected_edges_from(
        self, edges: list[tuple[T, T] | tuple[T, T, int]], **attr: Any
    ) -> None:
        """
        Add all undirected edges in 'edges' to the graph.

        Parameters
        ----------
        edges : list[tuple[T, T] | tuple[T, T, int]]
            List of undirected edges. Can be (u, v) or (u, v, key) tuples.
        **attr
            Edge attributes applied to all edges.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_undirected_edges_from([(1, 2), (2, 3), (3, 1)])
        """
        for edge in edges:
            if len(edge) == 2:
                u, v = edge
                self.add_undirected_edge(u, v, **attr)
            elif len(edge) == 3:
                u, v, key = edge
                self.add_undirected_edge(u, v, key=key, **attr)
            else:
                raise ValueError(f"Invalid edge format: {edge}")

    def remove_edge(self, u: T, v: T, key: int | None = None) -> None:
        """
        Remove edge (u, v) from the graph.

        If key is provided, removes the specific edge with that key.
        If key is None, tries undirected first, then directed.

        Parameters
        ----------
        u : T
            Source node.
        v : T
            Target node.
        key : int | None, optional
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
        0
        >>> G.remove_edge(1, 2, key=key)
        """
        removed = False

        # Try undirected first
        if self._undirected.has_edge(u, v, key):
            self._undirected.remove_edge(u, v, key)
            self._combined.remove_edge(u, v, key)
            removed = True
        # Then try directed
        elif self._directed.has_edge(u, v, key):
            self._directed.remove_edge(u, v, key)
            self._combined.remove_edge(u, v, key)
            removed = True

        if not removed:
            raise KeyError(f"Edge ({u}, {v}, {key}) not found")

    def remove_edges_from(self, edges: list[tuple[T, T] | tuple[T, T, int]]) -> None:
        """
        Remove all edges in 'edges' from the graph.

        Parameters
        ----------
        edges : list[tuple[T, T] | tuple[T, T, int]]
            List of edges to remove. Can be (u, v) or (u, v, key) tuples.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2)
        0
        >>> G.add_undirected_edge(2, 3)
        0
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

    @property
    def undirected_edges(self) -> 'EdgeView':
        """
        Get all undirected edges (works as both attribute and method).

        When accessed as attribute, returns an EdgeView object (list-like).
        When called as method, returns an iterator with optional keys and data.

        Returns
        -------
        EdgeView
            Edge view object that's both callable and list-like.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_undirected_edge(1, 2, weight=1.0)
        0
        >>> G.add_undirected_edge(1, 2, weight=2.0)  # Parallel edge
        1
        >>> G.add_directed_edge(2, 3)
        0
        >>> G.undirected_edges  # Attribute access
        [(1, 2), (1, 2)]
        >>> list(G.undirected_edges())  # Method call
        [(1, 2), (1, 2)]
        >>> list(G.undirected_edges(keys=True))  # With keys
        [(1, 2, 0), (1, 2, 1)]
        >>> list(G.undirected_edges(data=True))  # With data
        [(1, 2, {'weight': 1.0}), (1, 2, {'weight': 2.0})]
        >>> list(G.undirected_edges(keys=True, data=True))  # With keys and data
        [(1, 2, 0, {'weight': 1.0}), (1, 2, 1, {'weight': 2.0})]
        """
        edges_list = list(self._undirected.edges())
        return self.EdgeView(edges_list, self.undirected_edges_iter)

    # ========== Query Operations ==========

    def has_edge(self, u: T, v: T, key: int | None = None) -> bool:
        """
        Check if edge exists (directed or undirected).

        Parameters
        ----------
        u : T
            Source node.
        v : T
            Target node.
        key : int | None, optional
            Edge key. By default None.

        Returns
        -------
        bool
            True if edge exists, False otherwise.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_undirected_edge(1, 2)
        0
        >>> G.has_edge(1, 2)
        True
        >>> G.has_edge(2, 3)
        False
        """
        return self._directed.has_edge(u, v, key) or self._undirected.has_edge(u, v, key)

    def number_of_nodes(self) -> int:
        """
        Return the number of nodes.

        Returns
        -------
        int
            Number of nodes.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_nodes_from([1, 2, 3])
        >>> G.number_of_nodes()
        3
        """
        return len(self)

    def number_of_edges(self) -> int:
        """
        Return the number of edges (directed + undirected).

        Returns
        -------
        int
            Number of edges.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2)
        0
        >>> G.add_undirected_edge(2, 3)
        0
        >>> G.number_of_edges()
        2
        """
        return self._directed.number_of_edges() + self._undirected.number_of_edges()

    def __repr__(self) -> str:
        """
        Return a concise representation.
        
        Returns
        -------
        str
            Representation containing counts of nodes, directed edges, and undirected edges.
        """
        return (
            "MixedMultiGraph("
            f"nodes={self.number_of_nodes()}, "
            f"directed_edges={self._directed.number_of_edges()}, "
            f"undirected_edges={self._undirected.number_of_edges()}"
            ")"
        )

    def degree(self, v: T) -> int:
        """
        Return the total degree of vertex v.

        The degree is the sum of undirected edges and directed edges
        (both incoming and outgoing).

        Parameters
        ----------
        v : T
            Vertex.

        Returns
        -------
        int
            Total degree of v. Returns 0 if v is not in the graph.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_undirected_edge(1, 2)
        0
        >>> G.add_directed_edge(1, 3)
        0
        >>> G.degree(1)
        2
        >>> G.degree(99)  # Node not in graph
        0
        """
        undirected_degree = self.undirected_degree(v)
        directed_degree = self.indegree(v) + self.outdegree(v)
        return undirected_degree + directed_degree

    def undirected_degree(self, v: T) -> int:
        """
        Return the undirected degree of vertex v.

        The undirected degree is the number of undirected edges incident to v.
        Each undirected edge contributes 1 to the degree.

        Parameters
        ----------
        v : T
            Vertex.

        Returns
        -------
        int
            Undirected degree of v. Returns 0 if v is not in the graph.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_undirected_edge(1, 2)
        0
        >>> G.add_undirected_edge(1, 2)  # Parallel edge
        1
        >>> G.add_undirected_edge(1, 3)
        0
        >>> G.undirected_degree(1)
        3
        >>> G.undirected_degree(99)  # Node not in graph
        0
        """
        if v not in self._undirected:
            return 0
        return self._undirected.degree(v)

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
            Indegree of v. Returns 0 if v is not in the graph.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2)
        0
        >>> G.add_directed_edge(3, 2)
        0
        >>> G.indegree(2)
        2
        >>> G.indegree(99)  # Node not in graph
        0
        """
        if v not in self._directed:
            return 0
        return self._directed.in_degree(v)

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
            Outdegree of v. Returns 0 if v is not in the graph.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2)
        0
        >>> G.add_directed_edge(1, 3)
        0
        >>> G.outdegree(1)
        2
        >>> G.outdegree(99)  # Node not in graph
        0
        """
        if v not in self._directed:
            return 0
        return self._directed.out_degree(v)



    # ========== Graph Operations ==========

    def _validate_synchronization(self) -> bool:
        """
        Validate that the internal graphs are synchronized.
        
        Checks that `_combined` contains all edges from both `_undirected`
        and `_directed` graphs. This is useful for debugging if the graphs
        have been modified directly (which should not happen).
        
        Returns
        -------
        bool
            True if graphs are synchronized, False otherwise.
        
        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2)
        0
        >>> G.add_undirected_edge(2, 3)
        0
        >>> G._validate_synchronization()
        True
        >>> # Direct modification (BAD - don't do this!)
        >>> G._undirected.add_edge(99, 100)
        >>> G._validate_synchronization()
        False
        """
        # Check that all undirected edges are in combined
        for u, v, key in self._undirected.edges(keys=True):
            if not self._combined.has_edge(u, v, key):
                return False
        
        # Check that all directed edges are in combined
        for u, v, key in self._directed.edges(keys=True):
            if not self._combined.has_edge(u, v, key):
                return False
        
        # Check that combined doesn't have extra edges (shouldn't happen, but check anyway)
        combined_edges = set(self._combined.edges(keys=True))
        undirected_edges = set(self._undirected.edges(keys=True))
        directed_edges = set(self._directed.edges(keys=True))
        expected_edges = undirected_edges | directed_edges
        
        if combined_edges != expected_edges:
            return False
        
        return True

    def copy(self) -> 'MixedMultiGraph':
        """
        Create a copy of the graph.

        Returns
        -------
        MixedMultiGraph
            A copy of the graph.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2)
        0
        >>> H = G.copy()
        >>> H.add_undirected_edge(2, 3)
        0
        >>> G.number_of_edges()
        1
        >>> H.number_of_edges()
        2
        """
        new_graph = MixedMultiGraph()
        new_graph._undirected = self._undirected.copy()
        new_graph._directed = self._directed.copy()
        new_graph._combined = self._combined.copy()
        return new_graph

    def clear(self) -> None:
        """
        Clear the whole mixed graph.

        Removes all nodes and edges.

        Examples
        --------
        >>> G = MixedMultiGraph()
        >>> G.add_directed_edge(1, 2)
        0
        >>> G.add_undirected_edge(2, 3)
        0
        >>> G.clear()
        >>> len(G.nodes())
        0
        """
        self._undirected.clear()
        self._directed.clear()
        self._combined.clear()

    # Note: identify_two_nodes and identify_node_set are now functions
    # in mm_operations module (NetworkX-style API)
