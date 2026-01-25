Mixed Multi-Graphs
==================

The :mod:`phylozoo.core.primitives.m_multigraph` module provides the
:class:`MixedMultiGraph` class, a mixed multigraph that supports both directed
and undirected edges with parallel edges of both types. This class serves as the
foundation for :class:`SemiDirectedPhyNetwork` and enables representation of
phylogenetic networks with both tree-like and reticulate evolutionary relationships.

All classes and functions on this page can be imported from the core primitives module:

.. code-block:: python

   from phylozoo.core.primitives import *
   # or directly
   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph

Working with Mixed Multi-Graphs
--------------------------------

Mixed multi-graphs combine the capabilities of directed and undirected graphs,
allowing complex phylogenetic structures where some relationships are directed
(reticulations, hybridization) and others are undirected (tree-like relationships).
This flexibility makes them ideal for semi-directed phylogenetic networks.

.. note::
   :class: dropdown

   **Implementation details**

   MixedMultiGraph is designed for complex phylogenetic network analysis:

   - Separate NetworkX MultiDiGraph and MultiGraph for directed/undirected edges
   - Combined view for efficient connectivity analysis
   - Mutual exclusivity: edges between same nodes must be all directed or all undirected
   - Parallel edges supported for both edge types
   - Rich attribute support for phylogenetic annotations

   For implementation details, see :mod:`src/phylozoo/core/primitives/m_multigraph/base.py`.

Creating Mixed Multi-Graphs
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mixed multi-graphs can be created with separate directed and undirected edge specifications:

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph

   # Empty graph
   empty_graph = MixedMultiGraph()

   # Graph with only undirected edges (tree-like)
   tree_graph = MixedMultiGraph(undirected_edges=[(1, 2), (2, 3), (3, 4)])

   # Graph with only directed edges (reticulate)
   network_graph = MixedMultiGraph(directed_edges=[(1, 2), (1, 3), (2, 4), (3, 4)])

   # Mixed graph (semi-directed network)
   mixed_graph = MixedMultiGraph(
       directed_edges=[(1, 3), (2, 3)],      # Reticulations to node 3
       undirected_edges=[(3, 4), (4, 5), (4, 6)]  # Tree structure below
   )

   # With edge attributes
   attributed_graph = MixedMultiGraph(
       directed_edges=[{'u': 1, 'v': 2, 'weight': 1.0, 'type': 'hybrid'}],
       undirected_edges=[{'u': 2, 'v': 3, 'length': 0.5, 'bootstrap': 95}]
   )

   # With graph attributes
   annotated_graph = MixedMultiGraph(
       directed_edges=[(1, 2)],
       undirected_edges=[(2, 3)],
       attributes={'name': 'phylogenetic_network', 'method': 'inference'}
   )

Mutual Exclusivity Constraint
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A key feature of mixed multi-graphs is the mutual exclusivity constraint: edges
between the same pair of nodes must be either all directed or all undirected:

.. code-block:: python

   graph = MixedMultiGraph()

   # Add undirected edge between nodes 1 and 2
   graph.add_undirected_edge(1, 2, weight=1.0)

   # Adding directed edge between same nodes removes undirected edge
   graph.add_directed_edge(1, 2, weight=2.0)  # Removes the undirected edge

   # Now only directed edge exists
   assert graph.has_directed_edge(1, 2)
   assert not graph.has_undirected_edge(1, 2)

   # Vice versa: undirected edge removes directed edge
   graph.add_undirected_edge(1, 2, weight=3.0)  # Removes the directed edge
   assert graph.has_undirected_edge(1, 2)
   assert not graph.has_directed_edge(1, 2)

This constraint ensures mathematical consistency while allowing flexible network topologies.

Accessing Graph Structure
^^^^^^^^^^^^^^^^^^^^^^^^^

Mixed multi-graphs provide comprehensive access to their structure:

.. code-block:: python

   # Basic properties
   num_nodes = mixed_graph.number_of_nodes()
   num_edges = mixed_graph.number_of_edges()

   # Node access
   nodes = list(mixed_graph.nodes())
   node_exists = mixed_graph.has_node(1)

   # Directed edge access
   directed_edges = list(mixed_graph.directed_edges())  # (u, v, key) tuples
   directed_exists = mixed_graph.has_directed_edge(1, 3, key=0)

   # Undirected edge access
   undirected_edges = list(mixed_graph.undirected_edges())  # (u, v, key) tuples
   undirected_exists = mixed_graph.has_undirected_edge(3, 4, key=0)

   # Combined edge access
   all_edges = list(mixed_graph.edges())  # All edges regardless of direction

   # Connectivity
   neighbors = list(mixed_graph.neighbors(3))     # All adjacent nodes
   successors = list(mixed_graph.successors(1))   # Nodes reachable via directed edges
   predecessors = list(mixed_graph.predecessors(3)) # Nodes that reach via directed edges

Edge Operations
^^^^^^^^^^^^^^^

Mixed multi-graphs support comprehensive edge operations for both edge types:

**Directed Edges**

Directed edges represent relationships with a clear direction, such as parent-child
relationships in phylogenetic trees or reticulation events. You can add multiple
parallel directed edges between the same nodes, each with its own attributes.

.. code-block:: python

   # Add directed edges
   key1 = mixed_graph.add_directed_edge(1, 2, weight=1.0, type='reticulation')
   key2 = mixed_graph.add_directed_edge(1, 2, weight=2.0, type='parallel_reticulation')

   # Remove directed edge
   mixed_graph.remove_directed_edge(1, 2, key=key1)

   # Access directed edge data
   edge_data = mixed_graph.get_directed_edge_data(1, 2, key=key2)

**Undirected Edges**

Undirected edges represent symmetric relationships, such as tree edges in semi-directed
phylogenetic networks. Like directed edges, you can have multiple parallel undirected
edges between the same nodes.

.. code-block:: python

   # Add undirected edges
   key3 = mixed_graph.add_undirected_edge(2, 3, length=0.5, bootstrap=95)
   key4 = mixed_graph.add_undirected_edge(2, 3, length=0.7, bootstrap=87)

   # Remove undirected edge
   mixed_graph.remove_undirected_edge(2, 3, key=key3)

   # Access undirected edge data
   edge_data = mixed_graph.get_undirected_edge_data(2, 3, key=key4)

Node Operations
^^^^^^^^^^^^^^^

Nodes can be added with attributes and their connectivity analyzed:

.. code-block:: python

   # Add node with attributes
   mixed_graph.add_node(5, label='taxon_A', type='leaf', support=100)

   # Access node attributes
   node_attrs = mixed_graph.nodes[5]

   # Degree analysis
   in_degree = mixed_graph.in_degree(3)     # Directed edges pointing to node
   out_degree = mixed_graph.out_degree(1)   # Directed edges from node
   undirected_degree = mixed_graph.degree(3) # All edges incident to node

Graph Analysis Features
-----------------------

Mixed multi-graphs provide specialized algorithms for phylogenetic network analysis:

**Connectivity and Components**

The mixed multi-graph module provides specialized connectivity analysis functions.
The :func:`is_connected` function checks weak connectivity (ignoring edge directions),
while :func:`number_of_connected_components` counts components. The :func:`source_components`
function finds source components (nodes with no incoming directed edges), which is
important for phylogenetic network analysis. The :func:`has_self_loops` function
checks for edges connecting a node to itself.

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph.features import (
       is_connected, number_of_connected_components, has_self_loops, source_components
   )

   # Check weak connectivity (ignoring edge directions)
   connected = is_connected(mixed_graph)

   # Count connected components
   num_components = number_of_connected_components(mixed_graph)

   # Find source components (no incoming directed edges)
   sources = source_components(mixed_graph)  # List of node sets

   # Check for self-loops
   has_loops = has_self_loops(mixed_graph)

**Graph Isomorphism**

The :func:`is_isomorphic` function compares graph structures while respecting the
mixed edge types, checking if two graphs have the same topology regardless of
node labeling.

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph.isomorphism import is_isomorphic

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

Graph Conversions
^^^^^^^^^^^^^^^^^

Convert from various NetworkX graph types:

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

File Input/Output
^^^^^^^^^^^^^^^^^

Mixed multi-graphs support specialized phylogenetic formats:

- **PhyloZoo DOT**: Extended DOT format preserving directed/undirected distinction

.. code-block:: python

   # Save to file
   mixed_graph.save('phylogenetic_network.pzdot')

   # Load from file
   loaded_graph = MixedMultiGraph.load('phylogenetic_network.pzdot')

   # Convert to string
   pzdot_string = mixed_graph.to_string(format='phylozoo-dot')

   # Create from string
   graph_from_string = MixedMultiGraph.from_string(pzdot_string, format='phylozoo-dot')

See Also
--------

- :doc:`API Reference <../../../api/core/primitives>` - Complete function signatures and detailed examples
- :doc:`Directed Multi-Graph <directed_multigraph>` - Directed-only graphs
- :doc:`Networks (Basic) <../networks/basic>` - Network classes using mixed multi-graphs
- :doc:`Circular Ordering <circular_ordering>` - Circular arrangements for quartet analysis