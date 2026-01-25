Directed Multi-Graphs
=====================

The :mod:`phylozoo.core.primitives.d_multigraph` module provides the
:class:`DirectedMultiGraph` class, a directed multigraph that supports parallel
edges between nodes. This class serves as the foundation for :class:`DirectedPhyNetwork`
and enables representation of complex phylogenetic relationships with reticulations,
multiple inheritance, and other directed graph structures.

All classes and functions on this page can be imported from the core primitives module:

.. code-block:: python

   from phylozoo.core.primitives import *
   # or directly
   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph

Working with Directed Multi-Graphs
-----------------------------------

Directed multi-graphs extend traditional directed graphs by allowing multiple
parallel edges between the same pair of nodes. This is essential for phylogenetic
networks where a single parent-child relationship may have multiple evolutionary
paths or reticulate events.

.. note::
   :class: dropdown

   **Implementation details**

   DirectedMultiGraph is optimized for phylogenetic network operations:

   - Uses NetworkX MultiDiGraph internally for directed edge storage
   - Maintains a combined MultiGraph for efficient connectivity analysis
   - Supports parallel edges with unique keys for identification
   - Edge and node attributes enable rich phylogenetic annotations
   - Immutable design ensures data integrity in analysis pipelines

   For implementation details, see :mod:`src/phylozoo/core/primitives/d_multigraph/base.py`.

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

Graph Analysis Features
-----------------------

The directed multi-graph module provides algorithms for analyzing graph properties:

**Connectivity Analysis**

The directed multi-graph module provides functions for analyzing graph connectivity.
The :func:`is_connected` function checks if a graph is weakly connected (ignoring
edge directions), while :func:`number_of_connected_components` counts the number
of connected components. The :func:`has_self_loops` function checks for edges
that connect a node to itself.

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph.features import (
       is_connected, number_of_connected_components, has_self_loops
   )

   # Check if graph is weakly connected (ignoring edge directions)
   connected = is_connected(graph)

   # Count connected components
   num_components = number_of_connected_components(graph)

   # Check for self-loops
   has_loops = has_self_loops(graph)

**Graph Isomorphism**

The :func:`is_isomorphic` function checks if two graphs have the same structure
regardless of node labeling. This is useful for comparing phylogenetic networks
that may have different node labels but identical topologies.

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

Graph Conversions
^^^^^^^^^^^^^^^^^

Convert from NetworkX graphs to DirectedMultiGraph:

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

File Input/Output
^^^^^^^^^^^^^^^^^

Directed multi-graphs support serialization in standard graph formats:

- **DOT**: GraphViz format for visualization and storage

.. code-block:: python

   # Save to file
   graph.save('phylogeny.dot')

   # Load from file
   loaded_graph = DirectedMultiGraph.load('phylogeny.dot')

   # Convert to string representation
   dot_string = graph.to_string(format='dot')

   # Create from string
   graph_from_string = DirectedMultiGraph.from_string(dot_string, format='dot')

See Also
--------

- :doc:`API Reference <../../../api/core/primitives>` - Complete function signatures and detailed examples
- :doc:`Mixed Multi-Graph <mixed_multigraph>` - Graphs with both directed and undirected edges
- :doc:`Networks (Basic) <../networks/basic>` - Network classes using directed multi-graphs
- :doc:`Partition <partition>` - Set partitions for hierarchical analysis