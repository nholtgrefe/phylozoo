Directed Multi-Graphs
=====================

The :mod:`phylozoo.core.primitives.d_multigraph` module provides the
:class:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph` class, a directed multigraph that supports parallel
edges between nodes. This class serves as the foundation for :class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork`
and enables representation of complex phylogenetic relationships with reticulations,
, and other directed graph structures. It builds on the :class:`networkx.MultiDiGraph` class.

All classes and functions on this page can be imported from the directed multigraph submodule:

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph


Working with Directed Multi-Graphs
----------------------------------

Creating and modifying directed multi-graphs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can build a directed multi-graph from scratch (empty or from an edge list), add or remove
edges and nodes, and query structure such as neighbors and degrees. Each edge has a unique key
so that parallel edges between the same pair of nodes are distinct.

Construct a graph with no arguments for an empty graph, or pass ``edges`` as a list of
``(u, v)`` tuples or dicts with ``u`` and ``v`` keys. Repeated ``(u, v)`` pairs create
parallel edges.

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph

   empty_graph = DirectedMultiGraph()
   graph = DirectedMultiGraph(edges=[(1, 2), (2, 3), (3, 1)])
   parallel_graph = DirectedMultiGraph(edges=[(1, 2), (1, 2), (1, 2)])

You can attach arbitrary attributes to edges (e.g. ``weight``, ``label``) and to the graph
itself via the ``attributes`` argument:

.. code-block:: python

   attributed_graph = DirectedMultiGraph(edges=[
       {'u': 1, 'v': 2, 'weight': 1.0, 'label': 'edge1'},
       {'u': 2, 'v': 3, 'weight': 2.0, 'type': 'reticulation'}
   ])
   annotated_graph = DirectedMultiGraph(
       edges=[(1, 2), (2, 3)],
       attributes={'source': 'simulation', 'version': '1.0'}
   )

.. tip::
   For graphs without parallel edges, you can build or manipulate them in NetworkX
   (e.g. :class:`networkx.DiGraph`) and then convert to
   :class:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph`. See the section
   on "NetworkX Conversion" below.

**Edge operations**

The :meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.add_edge` method adds a new edge and returns the key of the new edge.
Multiple edges between the same nodes are allowed (parallel edges).

.. code-block:: python

   graph = DirectedMultiGraph()
   key1 = graph.add_edge(1, 2, weight=1.0, label='primary')
   key2 = graph.add_edge(1, 2, weight=2.0, label='secondary')
   key3 = graph.add_edge(2, 3, key=10, edge_type='reticulation')

To remove an edge between two nodes, use the key and :meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.remove_edge`.

.. code-block:: python

   graph.remove_edge(1, 2, key=key1)

For adding or removing many edges at once, use
:meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.add_edges_from` and
:meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.remove_edges_from`:

.. code-block:: python

   graph.add_edges_from([(1, 2), (2, 3), (3, 1)])
   graph.remove_edges_from([(2, 3)])  # (u, v) or (u, v, key)

**Node operations**

Add nodes with optional attributes; the :attr:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.nodes`
view lets you read attributes by node.

.. code-block:: python

   graph.add_node(4, label='taxon_A', type='leaf', bootstrap=95)
   node_attrs = graph.nodes[4]

To add or remove multiple nodes at once, use
:meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.add_nodes_from`,
:meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.remove_node`, and
:meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.remove_nodes_from`.
The :meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.generate_node_ids` method
returns an iterator of fresh integer IDs.

.. code-block:: python

   graph.add_nodes_from([5, 6, 7], type='internal')
   graph.remove_node(7)
   ids = list(graph.generate_node_ids(3))

**Accessing graph structure**

Counts and membership are available via :meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.number_of_nodes`
and :meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.number_of_edges`; use ``v in graph`` for node
membership and :meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.has_edge` for edge membership. Iterate
over :attr:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.nodes` and
:attr:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.edges` for all
nodes and ``(u, v, key)`` tuples.

.. code-block:: python

   num_nodes = graph.number_of_nodes()
   num_edges = graph.number_of_edges()
   nodes = list(graph.nodes())
   edges = list(graph.edges())
   edge_data = graph.get_edge_data(1, 2, key=0)

For connectivity and degrees, use
:meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.successors` and
:meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.predecessors` for
neighbors by direction, and :meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.outdegree`,
:meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.indegree`, and
:meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.degree` for in-, out-,
and total degree.

.. code-block:: python

   successors = list(graph.successors(2))
   predecessors = list(graph.predecessors(2))
   out_degree = graph.outdegree(2)
   in_degree = graph.indegree(2)
   total_degree = graph.degree(2)

:meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.neighbors` returns
all nodes adjacent to a vertex (predecessors and successors combined). For edges incident
to a vertex, use
:meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.incident_parent_edges`
(incoming) and
:meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.incident_child_edges`
(outgoing); both support ``keys`` and ``data``. :meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.copy`
returns a shallow copy of the graph, and :meth:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph.clear`
removes all nodes and edges.

.. code-block:: python

   neighbors = list(graph.neighbors(2))
   incoming = list(graph.incident_parent_edges(2, keys=True))
   outgoing = list(graph.incident_child_edges(2, keys=True))
   graph_copy = graph.copy()
   graph_copy.clear()

File Input/Output
^^^^^^^^^^^^^^^^^

Directed multi-graphs support reading and writing in standard graph formats:

- **DOT** (default): GraphViz format for visualization and storage — see :doc:`DOT format <../../utils/io/formats/dot>`
- **Edge list**: Simple text format with one edge per line — see :doc:`Edge list format <../../utils/io/formats/edgelist>`

.. code-block:: python

   # Load from file (auto-detects format by extension)
   graph = DirectedMultiGraph.load("phylogeny.dot")

   # Load with explicit format
   graph = DirectedMultiGraph.load("phylogeny.edgelist", format="edgelist")

   # Save to file
   graph.save("output.dot")
   graph.save("output.edgelist", format="edgelist")

.. seealso::
   The :class:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph` class uses the
   :class:`~phylozoo.utils.io.mixin.IOMixin` interface, providing consistent file handling across PhyloZoo
   classes. For details on the I/O system, see the :doc:`I/O manual <../../utils/io/overview>`.


Graph Analysis Features
-----------------------

The directed multi-graph module provides algorithms for analyzing graph properties.

**Connectivity**

The :func:`~phylozoo.core.primitives.d_multigraph.features.is_connected` function checks if a graph is weakly connected (ignoring edge directions).

The :func:`~phylozoo.core.primitives.d_multigraph.features.number_of_connected_components` function counts the number of connected components.

The :func:`~phylozoo.core.primitives.d_multigraph.features.connected_components` function returns an iterator over all connected components.

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph.features import (
       is_connected, number_of_connected_components, connected_components
   )

   # Check if graph is weakly connected (ignoring edge directions)
   connected = is_connected(graph)

   # Count connected components
   num_components = number_of_connected_components(graph)

   # Get all connected components
   for component in connected_components(graph):
       print(f"Component: {component}")

**Biconnectivity**

The :func:`~phylozoo.core.primitives.d_multigraph.features.biconnected_components` function returns an iterator over biconnected components (maximal sets of nodes where any two nodes are connected by at least two node-disjoint paths).

The :func:`~phylozoo.core.primitives.d_multigraph.features.bi_edge_connected_components` function returns an iterator over bi-edge-connected components (maximal sets of nodes where any two nodes are connected by at least two edge-disjoint paths).

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph.features import (
       biconnected_components, bi_edge_connected_components
   )

   # Get biconnected components
   for component in biconnected_components(graph):
       print(f"Biconnected component: {component}")

   # Get bi-edge-connected components
   for component in bi_edge_connected_components(graph):
       print(f"Bi-edge-connected component: {component}")

**Cut edges and vertices**

The :func:`~phylozoo.core.primitives.d_multigraph.features.cut_edges` function finds all cut edges (edges whose removal increases the number of connected components).

The :func:`~phylozoo.core.primitives.d_multigraph.features.cut_vertices` function finds all cut vertices (nodes whose removal increases the number of connected components).

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph.features import cut_edges, cut_vertices

   # Find cut edges
   cuts = cut_edges(graph)  # Returns list of (u, v, key) tuples

   # Find cut vertices
   cut_nodes = cut_vertices(graph)  # Returns set of nodes

**Parallel edges and self-loops**

The :func:`~phylozoo.core.primitives.d_multigraph.features.has_parallel_edges` function checks if the graph contains any parallel edges (multiple edges between the same pair of nodes).

The :func:`~phylozoo.core.primitives.d_multigraph.features.has_self_loops` function checks for edges that connect a node to itself.

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph.features import (
       has_parallel_edges, has_self_loops
   )

   # Check for parallel edges
   has_parallel = has_parallel_edges(graph)

   # Check for self-loops
   has_loops = has_self_loops(graph)

**Isomorphism**

The :func:`~phylozoo.core.primitives.d_multigraph.isomorphism.is_isomorphic` function checks if two graphs have the same structure regardless of node labeling. This is useful for comparing phylogenetic networks that may have different node labels but identical topologies.

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph.isomorphism import is_isomorphic

   graph1 = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
   graph2 = DirectedMultiGraph(edges=[(4, 5), (5, 6)])

   # Basic isomorphism check
   isomorphic = is_isomorphic(graph1, graph2)  # True

   # With node attributes
   graph1.add_node(1, label='A')
   graph2.add_node(4, label='A')
   isomorphic_with_attrs = is_isomorphic(graph1, graph2, node_attrs=['label'])

   # With edge attributes
   isomorphic_with_edges = is_isomorphic(graph1, graph2, edge_attrs=['weight'])

Graph Transformations
--------------------

The directed multi-graph module provides functions for transforming graph structures.

**Vertex identification**

The :func:`~phylozoo.core.primitives.d_multigraph.transformations.identify_vertices` function merges multiple vertices into a single vertex, combining their incident edges.

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph.transformations import identify_vertices

   # Merge vertices 1, 2, and 3 into a single vertex
   identify_vertices(graph, [1, 2, 3])

**Degree-2 node suppression**

The :func:`~phylozoo.core.primitives.d_multigraph.transformations.suppress_degree2_node` function removes a degree-2 node and connects its neighbors directly.

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph.transformations import suppress_degree2_node

   # Suppress a degree-2 node
   suppress_degree2_node(graph, node=5)

**Parallel edge identification**

The :func:`~phylozoo.core.primitives.d_multigraph.transformations.identify_parallel_edge` function merges all parallel edges between two nodes into a single edge.

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph.transformations import identify_parallel_edge

   # Merge all parallel edges between nodes 1 and 2
   identify_parallel_edge(graph, u=1, v=2)

**Subgraph extraction**

The :func:`~phylozoo.core.primitives.d_multigraph.transformations.subgraph` function creates a subgraph containing only the specified nodes and their incident edges.

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph.transformations import subgraph

   # Extract subgraph with specific nodes
   sub = subgraph(graph, nodes=[1, 2, 3, 4])

NetworkX Conversion
--------------------

Convert from NetworkX graphs to DirectedMultiGraph. The :func:`~phylozoo.core.primitives.d_multigraph.conversions.digraph_to_directedmultigraph` function converts a NetworkX DiGraph to a DirectedMultiGraph, and the :func:`~phylozoo.core.primitives.d_multigraph.conversions.multidigraph_to_directedmultigraph` function converts a NetworkX MultiDiGraph to a DirectedMultiGraph.

.. code-block:: python

   import networkx as nx
   from phylozoo.core.primitives.d_multigraph.conversions import (
       digraph_to_directedmultigraph,
       multidigraph_to_directedmultigraph
   )

   # Convert from NetworkX DiGraph
   nx_digraph = nx.DiGraph()
   nx_digraph.add_edge(1, 2, weight=1.0)
   dmg_from_digraph = digraph_to_directedmultigraph(nx_digraph)

   # Convert from NetworkX MultiDiGraph
   nx_multidigraph = nx.MultiDiGraph()
   nx_multidigraph.add_edge(1, 2, key=0, weight=1.0)
   nx_multidigraph.add_edge(1, 2, key=1, weight=2.0)
   dmg_from_multidigraph = multidigraph_to_directedmultigraph(nx_multidigraph)


.. tip::
   The class stores two NetworkX graphs that can be accessed for further use with NetworkX algorithms. The whole directed multigraph is stored as a :class:`networkx.MultiDiGraph` in the ``_graph`` attribute . The ``_combined`` attribute stores a completely undirected view of the graph as a :class:`networkx.MultiGraph` for quick connectivity analyses.


See Also
--------

- :doc:`API Reference <../../../../api/core/primitives/index>` - Complete function signatures and detailed examples
- :doc:`Mixed Multi-Graph <mixed_multigraph>` - Graphs with both directed and undirected edges
- :doc:`Networks <../networks/overview>` - Network classes using directed multi-graphs
