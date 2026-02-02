Directed Multi-Graphs
=====================

The :mod:`phylozoo.core.primitives.d_multigraph` module provides the
:class:`DirectedMultiGraph` class, a directed multigraph that supports parallel
edges between nodes. This class serves as the foundation for :class:`DirectedPhyNetwork`
and enables representation of complex phylogenetic relationships with reticulations,
multiple inheritance, and other directed graph structures.

All classes and functions on this page can be imported from the directed multigraph submodule:

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph


.. note::
   :class: dropdown

   **Implementation details**

   DirectedMultiGraph is optimized for phylogenetic network operations:

   - Uses NetworkX MultiDiGraph internally for directed edge storage
   - Maintains a combined MultiGraph for efficient connectivity analysis
   - Supports parallel edges with unique keys for identification
   - Edge and node attributes enable rich phylogenetic annotations

Working with Directed Multi-Graphs
-----------------------------------

Creating Directed Multi-Graphs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Directed multi-graphs can be created from various input formats:

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph

   # Empty graph
   empty_graph = DirectedMultiGraph()

   # From edge list
   graph = DirectedMultiGraph(edges=[(1, 2), (2, 3), (3, 1)])

   # With parallel edges (multiple edges between same nodes)
   parallel_graph = DirectedMultiGraph(edges=[(1, 2), (1, 2), (1, 2)])

   # From edge dictionaries with attributes
   attributed_graph = DirectedMultiGraph(edges=[
       {'u': 1, 'v': 2, 'weight': 1.0, 'label': 'edge1'},
       {'u': 2, 'v': 3, 'weight': 2.0, 'type': 'reticulation'}
   ])

   # With graph-level attributes
   annotated_graph = DirectedMultiGraph(
       edges=[(1, 2), (2, 3)],
       attributes={'source': 'simulation', 'version': '1.0'}
   )


Edge Operations
^^^^^^^^^^^^^^^

Directed multi-graphs support comprehensive edge manipulation with parallel edge handling:

.. code-block:: python

   graph = DirectedMultiGraph()

   # Add edges (returns unique key for each edge)
   key1 = graph.add_edge(1, 2, weight=1.0, label='primary')
   key2 = graph.add_edge(1, 2, weight=2.0, label='secondary')  # Parallel edge

   # Add edge with explicit key
   key3 = graph.add_edge(2, 3, key=10, type='reticulation')

   # Remove specific edge by key
   graph.remove_edge(1, 2, key=key1)

   # Update edge attributes
   graph._graph[1][2][key2]['weight'] = 3.0

   # Get all edges between specific nodes
   edges_1_2 = graph.get_edge_data(1, 2)  # Dict of key -> attributes

Node Operations
^^^^^^^^^^^^^^^

Nodes can be added with attributes and their properties accessed:

.. code-block:: python

   # Add node with attributes
   graph.add_node(4, label='taxon_A', type='leaf', bootstrap=95)

   # Access node attributes
   node_attrs = graph.nodes[4]  # Returns attribute dictionary

   # Check node existence
   exists = graph.has_node(4)

   # Get node degree information
   total_degree = graph.degree(4)  # Total edges incident to node


Accessing Graph Structure
^^^^^^^^^^^^^^^^^^^^^^^^^

Directed multi-graphs provide comprehensive access to their structure and properties:

.. code-block:: python

   # Basic properties
   num_nodes = graph.number_of_nodes()
   num_edges = graph.number_of_edges()

   # Node access
   nodes = list(graph.nodes())  # All node identifiers
   node_exists = graph.has_node(1)

   # Edge access
   edges = list(graph.edges())  # All (u, v, key) tuples
   edge_exists = graph.has_edge(1, 2, key=0)

   # Edge data access
   edge_data = graph.get_edge_data(1, 2, key=0)  # Returns attribute dict

   # Connectivity
   successors = list(graph.successors(2))    # Nodes reachable from 2
   predecessors = list(graph.predecessors(2)) # Nodes that reach 2

   # Degrees
   out_degree = graph.out_degree(2)  # Number of outgoing edges
   in_degree = graph.in_degree(2)    # Number of incoming edges

File Input/Output
^^^^^^^^^^^^^^^^^

Directed multi-graphs support serialization in standard graph formats:

- **DOT**: GraphViz format for visualization and storage
- **Edge List**: Simple text format with one edge per line

.. code-block:: python

   # Save to file (DOT format)
   graph.save('phylogeny.dot')

   # Load from file (auto-detects format by extension)
   loaded_graph = DirectedMultiGraph.load('phylogeny.dot')

   # Save to edge list format
   graph.save('phylogeny.edgelist', format='edgelist')

   # Load from edge list
   loaded_from_edgelist = DirectedMultiGraph.load('phylogeny.edgelist', format='edgelist')

   # Convert to string representation
   dot_string = graph.to_string(format='dot')
   edgelist_string = graph.to_string(format='edgelist')

   # Create from string
   graph_from_string = DirectedMultiGraph.from_string(dot_string, format='dot')
   graph_from_edgelist = DirectedMultiGraph.from_string(edgelist_string, format='edgelist')


Graph Analysis Features
-----------------------

The directed multi-graph module provides algorithms for analyzing graph properties.

Connectivity Analysis
^^^^^^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.primitives.d_multigraph.features.is_connected` function checks if a graph is weakly connected (ignoring edge directions).

The :func:`phylozoo.core.primitives.d_multigraph.features.number_of_connected_components` function counts the number of connected components.

The :func:`phylozoo.core.primitives.d_multigraph.features.connected_components` function returns an iterator over all connected components.

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

Biconnectivity Analysis
^^^^^^^^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.primitives.d_multigraph.features.biconnected_components` function returns an iterator over biconnected components (maximal sets of nodes where any two nodes are connected by at least two node-disjoint paths).

The :func:`phylozoo.core.primitives.d_multigraph.features.bi_edge_connected_components` function returns an iterator over bi-edge-connected components (maximal sets of nodes where any two nodes are connected by at least two edge-disjoint paths).

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

Cut Analysis
^^^^^^^^^^^^

The :func:`phylozoo.core.primitives.d_multigraph.features.cut_edges` function finds all cut edges (edges whose removal increases the number of connected components).

The :func:`phylozoo.core.primitives.d_multigraph.features.cut_vertices` function finds all cut vertices (nodes whose removal increases the number of connected components).

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph.features import cut_edges, cut_vertices

   # Find cut edges
   cuts = cut_edges(graph)  # Returns list of (u, v, key) tuples

   # Find cut vertices
   cut_nodes = cut_vertices(graph)  # Returns set of nodes

Parallel Edges and Self-Loops
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.primitives.d_multigraph.features.has_parallel_edges` function checks if the graph contains any parallel edges (multiple edges between the same pair of nodes).

The :func:`phylozoo.core.primitives.d_multigraph.features.has_self_loops` function checks for edges that connect a node to itself.

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph.features import (
       has_parallel_edges, has_self_loops
   )

   # Check for parallel edges
   has_parallel = has_parallel_edges(graph)

   # Check for self-loops
   has_loops = has_self_loops(graph)

Graph Isomorphism
^^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.primitives.d_multigraph.isomorphism.is_isomorphic` function checks if two graphs have the same structure regardless of node labeling. This is useful for comparing phylogenetic networks that may have different node labels but identical topologies.

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

Vertex Identification
^^^^^^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.primitives.d_multigraph.transformations.identify_vertices` function merges multiple vertices into a single vertex, combining their incident edges.

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph.transformations import identify_vertices

   # Merge vertices 1, 2, and 3 into a single vertex
   identify_vertices(graph, [1, 2, 3])

Degree-2 Node Suppression
^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.primitives.d_multigraph.transformations.suppress_degree2_node` function removes a degree-2 node and connects its neighbors directly.

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph.transformations import suppress_degree2_node

   # Suppress a degree-2 node
   suppress_degree2_node(graph, node=5)

Parallel Edges Identification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.primitives.d_multigraph.transformations.identify_parallel_edge` function merges all parallel edges between two nodes into a single edge.

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph.transformations import identify_parallel_edge

   # Merge all parallel edges between nodes 1 and 2
   identify_parallel_edge(graph, u=1, v=2)

Subgraph Extraction
^^^^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.primitives.d_multigraph.transformations.subgraph` function creates a subgraph containing only the specified nodes and their incident edges.

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph.transformations import subgraph

   # Extract subgraph with specific nodes
   sub = subgraph(graph, nodes=[1, 2, 3, 4])

NetworkX Conversion
--------------------

Convert from NetworkX graphs to DirectedMultiGraph. The :func:`phylozoo.core.primitives.d_multigraph.conversions.digraph_to_directedmultigraph` function converts a NetworkX DiGraph to a DirectedMultiGraph, and the :func:`phylozoo.core.primitives.d_multigraph.conversions.multidigraph_to_directedmultigraph` function converts a NetworkX MultiDiGraph to a DirectedMultiGraph.

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

See Also
--------

- :doc:`API Reference <../../../api/core/primitives>` - Complete function signatures and detailed examples
- :doc:`Mixed Multi-Graph <mixed_multigraph>` - Graphs with both directed and undirected edges
- :doc:`Networks (Basic) <../networks/basic>` - Network classes using directed multi-graphs
