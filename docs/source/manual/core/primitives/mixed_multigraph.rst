Mixed Multi-Graphs
==================

The :mod:`phylozoo.core.primitives.m_multigraph` module provides the
:class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph` class, a mixed
multigraph that supports both directed and undirected edges with parallel edges of both
types. It serves as the foundation for :class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork`
and enables representation of phylogenetic networks with both tree-like and reticulate
evolutionary relationships.

All classes and functions on this page can be imported from the mixed multigraph submodule:

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph

Working with Mixed Multi-Graphs
--------------------------------

Creating and modifying mixed multi-graphs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mixed multi-graphs support both directed and undirected edges in one structure. You
specify directed and undirected edge lists separately; between any given pair of nodes,
edges must be either all directed or all undirected (mutual exclusivity). Parallel edges
are allowed for both types.

**Creating**

Use ``directed_edges`` and ``undirected_edges`` as lists of ``(u, v)`` tuples or dicts
with ``u`` and ``v``. You can pass only one of them for a graph that is fully directed
or fully undirected.

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph

   empty_graph = MixedMultiGraph()
   tree_graph = MixedMultiGraph(undirected_edges=[(1, 2), (2, 3), (3, 4)])
   network_graph = MixedMultiGraph(directed_edges=[(1, 2), (1, 3), (2, 4), (3, 4)])
   mixed_graph = MixedMultiGraph(
       directed_edges=[(1, 3), (2, 3)],
       undirected_edges=[(3, 4), (4, 5), (4, 6)]
   )

Optional ``nodes`` and ``attributes`` let you attach node or graph-level metadata:

.. code-block:: python

   attributed_graph = MixedMultiGraph(
       directed_edges=[{'u': 1, 'v': 2, 'weight': 1.0, 'type': 'hybrid'}],
       undirected_edges=[{'u': 2, 'v': 3, 'length': 0.5, 'bootstrap': 95}]
   )
   labeled_graph = MixedMultiGraph(
       directed_edges=[(1, 2)],
       undirected_edges=[(2, 3)],
       nodes=[{'id': 1, 'label': 'A'}, {'id': 2, 'label': 'B'}, {'id': 3, 'label': 'C'}]
   )

.. tip::
   For graphs without parallel edges, you can build or manipulate them in NetworkX
   (e.g. :class:`networkx.Graph` or :class:`networkx.DiGraph`) and then convert to
   :class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph`. See the
   **NetworkX conversion** section below.


**Directed edges**

Add or remove directed edges with
:meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.add_directed_edge` and
:meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.remove_directed_edge`;
each edge has a unique key.

.. code-block:: python

   key1 = mixed_graph.add_directed_edge(1, 2, weight=1.0, type='reticulation')
   key2 = mixed_graph.add_directed_edge(1, 2, weight=2.0)
   mixed_graph.remove_directed_edge(1, 2, key=key1)

Batch operations for directed edges:
:meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.add_directed_edges_from` and
:meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.remove_directed_edges_from`.

.. code-block:: python

   mixed_graph.add_directed_edges_from([(1, 4), (2, 4)])
   mixed_graph.remove_directed_edges_from([(1, 4)])

**Undirected edges**

Undirected edges are added with
:meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.add_undirected_edge` and
removed with :meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.remove_edge`
(which applies to undirected edges between the given pair). Parallel undirected edges are supported.

.. code-block:: python

   key3 = mixed_graph.add_undirected_edge(2, 3, length=0.5, bootstrap=95)
   key4 = mixed_graph.add_undirected_edge(2, 3, length=0.7, bootstrap=87)
   mixed_graph.remove_edge(2, 3, key=key3)

Batch operations for undirected edges:
:meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.add_undirected_edges_from` and
:meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.remove_edges_from`.

.. code-block:: python

   mixed_graph.add_undirected_edges_from([(4, 5), (5, 6)])
   mixed_graph.remove_edges_from([(5, 6)])

**Mutual exclusivity**

Between the same pair of nodes, all edges must be either directed or undirected. Adding
a directed edge between a pair that currently has undirected edges (or the reverse)
replaces those edges with the new type.

.. code-block:: python

   graph = MixedMultiGraph()
   graph.add_undirected_edge(1, 2, weight=1.0)
   graph.add_directed_edge(1, 2, weight=2.0)  # Replaces the undirected edge
   assert graph.has_directed_edge(1, 2) and not graph.has_undirected_edge(1, 2)

**Node operations**

Nodes are added with optional attributes; the :attr:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.nodes`
view gives attribute access by node.

.. code-block:: python

   mixed_graph.add_node(5, label='taxon_A', type='leaf', support=100)
   node_attrs = mixed_graph.nodes[5]
   in_degree = mixed_graph.indegree(3)
   out_degree = mixed_graph.outdegree(1)
   undirected_degree = mixed_graph.undirected_degree(3)
   total_degree = mixed_graph.degree(3)

You can add or remove multiple nodes with
:meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.add_nodes_from`,
:meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.remove_node`, and
:meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.remove_nodes_from`, and
:meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.generate_node_ids`
returns an iterator of fresh integer IDs.

.. code-block:: python

   mixed_graph.add_nodes_from([7, 8, 9], type='internal')
   mixed_graph.remove_node(9)

**Accessing graph structure**

Use :meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.number_of_nodes` and
:meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.number_of_edges` for
counts. The :attr:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.directed_edges`
and :attr:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.undirected_edges`
views yield ``(u, v, key)`` tuples; :attr:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.edges`
returns all edges. :meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.neighbors`
returns all nodes adjacent to a vertex (via any edge type).

.. code-block:: python

   num_nodes = mixed_graph.number_of_nodes()
   num_edges = mixed_graph.number_of_edges()
   directed_edges = list(mixed_graph.directed_edges())
   undirected_edges = list(mixed_graph.undirected_edges())
   all_edges = list(mixed_graph.edges())
   neighbors = list(mixed_graph.neighbors(3))

Use :meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.has_edge` for
edge membership (optionally with a key). For edges incident to a vertex, use
:meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.incident_parent_edges`,
:meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.incident_child_edges`, and
:meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.incident_undirected_edges`.
:meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.copy` returns a shallow
copy; :meth:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph.clear` removes all
nodes and edges.

.. code-block:: python

   mixed_graph.has_edge(1, 3, key=0)
   incoming = list(mixed_graph.incident_parent_edges(3, keys=True))
   outgoing = list(mixed_graph.incident_child_edges(1, keys=True))
   undir_at_3 = list(mixed_graph.incident_undirected_edges(3))
   graph_copy = mixed_graph.copy()
   graph_copy.clear()

File Input/Output
^^^^^^^^^^^^^^^^^

Mixed multi-graphs support reading and writing in PhyloZoo DOT format (default), which preserves
the directed/undirected edge distinction:

- **PhyloZoo DOT** (default): Extended DOT format for mixed graphs — see :doc:`DOT format <../../utils/io/formats/dot>`

.. code-block:: python

   # Load from file (auto-detects format by extension)
   mixed_graph = MixedMultiGraph.load("phylogenetic_network.pzdot")

   # Load with explicit format
   mixed_graph = MixedMultiGraph.load("network.txt", format="phylozoo-dot")

   # Save to file
   mixed_graph.save("output.pzdot")

.. seealso::
   The :class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph` class uses the
   :class:`~phylozoo.utils.io.IOMixin` interface, providing consistent file handling across PhyloZoo
   classes. For details on the I/O system, see the :doc:`I/O manual <../../utils/io/overview>`.

Graph Analysis Features
-----------------------

Mixed multi-graphs provide specialized algorithms for phylogenetic network analysis.

**Connectivity**

The :func:`~phylozoo.core.primitives.m_multigraph.features.is_connected` function checks weak connectivity (ignoring edge directions).

The :func:`~phylozoo.core.primitives.m_multigraph.features.number_of_connected_components` function counts the number of connected components.

The :func:`~phylozoo.core.primitives.m_multigraph.features.connected_components` function returns an iterator over all connected components.

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph.features import (
       is_connected, number_of_connected_components, connected_components
   )

   # Check weak connectivity (ignoring edge directions)
   connected = is_connected(mixed_graph)

   # Count connected components
   num_components = number_of_connected_components(mixed_graph)

   # Get all connected components
   for component in connected_components(mixed_graph):
       print(f"Component: {component}")

**Source components**

The :func:`~phylozoo.core.primitives.m_multigraph.features.source_components` function finds source components (nodes with no incoming directed edges), which is important for phylogenetic network analysis.

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph.features import source_components

   # Find source components (no incoming directed edges)
   sources = source_components(mixed_graph)  # List of (nodes, directed_edges, undirected_edges) tuples

**Biconnectivity**

The :func:`~phylozoo.core.primitives.m_multigraph.features.biconnected_components` function returns an iterator over biconnected components.

The :func:`~phylozoo.core.primitives.m_multigraph.features.bi_edge_connected_components` function returns an iterator over bi-edge-connected components.

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph.features import (
       biconnected_components, bi_edge_connected_components
   )

   # Get biconnected components
   for component in biconnected_components(mixed_graph):
       print(f"Biconnected component: {component}")

   # Get bi-edge-connected components
   for component in bi_edge_connected_components(mixed_graph):
       print(f"Bi-edge-connected component: {component}")

**Cut edges and vertices**

The :func:`~phylozoo.core.primitives.m_multigraph.features.cut_edges` function finds all cut edges.

The :func:`~phylozoo.core.primitives.m_multigraph.features.cut_vertices` function finds all cut vertices.

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph.features import cut_edges, cut_vertices

   # Find cut edges
   cuts = cut_edges(mixed_graph)  # Returns list of (u, v, key) tuples

   # Find cut vertices
   cut_nodes = cut_vertices(mixed_graph)  # Returns set of nodes

**Up-down paths**

The :func:`~phylozoo.core.primitives.m_multigraph.features.updown_path_vertices` function finds all vertices that lie on an up-down path between two given vertices (a path that first goes up via directed edges, then down via directed edges).

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph.features import updown_path_vertices

   # Find vertices on up-down paths between x and y
   vertices = updown_path_vertices(mixed_graph, x=1, y=5)

**Parallel edges and self-loops**

The :func:`~phylozoo.core.primitives.m_multigraph.features.has_parallel_edges` function checks if the graph contains any parallel edges.

The :func:`~phylozoo.core.primitives.m_multigraph.features.has_self_loops` function checks for edges connecting a node to itself.

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph.features import (
       has_parallel_edges, has_self_loops
   )

   # Check for parallel edges
   has_parallel = has_parallel_edges(mixed_graph)

   # Check for self-loops
   has_loops = has_self_loops(mixed_graph)

**Isomorphism**

The :func:`~phylozoo.core.primitives.m_multigraph.isomorphism.is_isomorphic` function compares graph structures while respecting the mixed edge types, checking if two graphs have the same topology regardless of node labeling.

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph.isomorphism import is_isomorphic

   graph1 = MixedMultiGraph(
       directed_edges=[(1, 3)],
       undirected_edges=[(3, 4)]
   )
   graph2 = MixedMultiGraph(
       directed_edges=[(2, 4)],
       undirected_edges=[(4, 5)]
   )

   # Basic isomorphism check
   isomorphic = is_isomorphic(graph1, graph2)

   # With node attributes
   graph1.add_node(1, label='root')
   graph2.add_node(2, label='root')
   isomorphic_with_attrs = is_isomorphic(graph1, graph2, node_attrs=['label'])

Graph Transformations
---------------------

The mixed multi-graph module provides functions for transforming graph structures.

**Vertex identification**

The :func:`~phylozoo.core.primitives.m_multigraph.transformations.identify_vertices` function merges multiple vertices into a single vertex, combining their incident edges.

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph.transformations import identify_vertices

   # Merge vertices 1, 2, and 3 into a single vertex
   identify_vertices(mixed_graph, [1, 2, 3])

**Orientation**

The :func:`~phylozoo.core.primitives.m_multigraph.transformations.orient_away_from_vertex` function orients all undirected edges away from a given root vertex, converting the mixed graph to a directed graph.

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph.transformations import orient_away_from_vertex

   # Orient all edges away from root vertex
   directed_graph = orient_away_from_vertex(mixed_graph, root=1)

**Degree-2 node suppression**

The :func:`~phylozoo.core.primitives.m_multigraph.transformations.suppress_degree2_node` function removes a degree-2 node and connects its neighbors directly.

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph.transformations import suppress_degree2_node

   # Suppress a degree-2 node
   suppress_degree2_node(mixed_graph, node=5)

**Parallel edge identification**

The :func:`~phylozoo.core.primitives.m_multigraph.transformations.identify_parallel_edge` function merges all parallel edges between two nodes into a single edge.

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph.transformations import identify_parallel_edge

   # Merge all parallel edges between nodes 1 and 2
   identify_parallel_edge(mixed_graph, u=1, v=2)

**Subgraph extraction**

The :func:`~phylozoo.core.primitives.m_multigraph.transformations.subgraph` function creates a subgraph containing only the specified nodes and their incident edges.

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph.transformations import subgraph

   # Extract subgraph with specific nodes
   sub = subgraph(mixed_graph, nodes=[1, 2, 3, 4])

NetworkX Conversion
-----------------

Convert from various NetworkX graph types. The :func:`~phylozoo.core.primitives.m_multigraph.conversions.graph_to_mixedmultigraph` function converts a NetworkX Graph to a MixedMultiGraph, the :func:`~phylozoo.core.primitives.m_multigraph.conversions.multigraph_to_mixedmultigraph` function converts a NetworkX MultiGraph to a MixedMultiGraph, and the :func:`~phylozoo.core.primitives.m_multigraph.conversions.multidigraph_to_mixedmultigraph` function converts a NetworkX MultiDiGraph to a MixedMultiGraph.
Also, the :func:`~phylozoo.core.primitives.m_multigraph.conversions.directedmultigraph_to_mixedmultigraph` function converts a DirectedMultiGraph to a MixedMultiGraph.

.. code-block:: python

   import networkx as nx
   from phylozoo.core.primitives.m_multigraph.conversions import (
       graph_to_mixedmultigraph,
       multigraph_to_mixedmultigraph,
       multidigraph_to_mixedmultigraph
   )

   # Convert undirected NetworkX Graph
   nx_graph = nx.Graph()
   nx_graph.add_edge(1, 2, weight=1.0)
   mmg_from_graph = graph_to_mixedmultigraph(nx_graph)

   # Convert undirected NetworkX MultiGraph
   nx_multigraph = nx.MultiGraph()
   nx_multigraph.add_edge(1, 2, key=0, weight=1.0)
   mmg_from_multigraph = multigraph_to_mixedmultigraph(nx_multigraph)

   # Convert directed NetworkX MultiDiGraph
   nx_multidigraph = nx.MultiDiGraph()
   nx_multidigraph.add_edge(1, 2, key=0, weight=1.0)
   mmg_from_multidigraph = multidigraph_to_mixedmultigraph(nx_multidigraph)

   # Convert from DirectedMultiGraph
   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   from phylozoo.core.primitives.m_multigraph.conversions import directedmultigraph_to_mixedmultigraph

   dmg = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
   mmg_from_dmg = directedmultigraph_to_mixedmultigraph(dmg)

.. tip::
   The class stores three NetworkX graphs that can be accessed for further use with
   NetworkX algorithms. The undirected edges are stored as a :class:`networkx.MultiGraph`
   in the ``_undirected`` attribute, the directed edges as a :class:`networkx.MultiDiGraph`
   in the ``_directed`` attribute, and the ``_combined`` attribute stores a combined
   undirected view of the graph as a :class:`networkx.MultiGraph` for quick connectivity
   analyses.

See Also
--------

- :doc:`API Reference <../../../../api/core/primitives/index>` - Complete function signatures and detailed examples
- :doc:`Directed Multi-Graph <directed_multigraph>` - Directed-only graphs
- :doc:`Networks <../networks/overview>` - Network classes using mixed multi-graphs
