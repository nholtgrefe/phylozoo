Mixed Multi-Graph
=================

The ``MixedMultiGraph`` class provides a mixed multigraph supporting both directed and 
undirected edges. This class is used internally by ``SemiDirectedPhyNetwork`` but can 
also be used directly for graph operations.

Overview
--------

.. figure:: ../../../../_static/images/example_mixed_multigraph.png
   :alt: Example MixedMultiGraph
   :width: 400px
   :align: center
   
   Example of a MixedMultiGraph with both directed and undirected edges.

A ``MixedMultiGraph`` supports:
- Both directed and undirected edges
- Parallel edges for both edge types
- Self-loops
- Edge attributes (weights, labels, etc.)
- Node attributes
- Graph-level attributes

**Important**: Edges between the same two nodes must be either all directed or all 
undirected - mixing is not allowed. This mutual exclusivity is enforced automatically: 
adding a directed edge will remove any undirected edges between the same nodes, and 
vice versa.

The class uses composition with NetworkX:
- ``_undirected``: ``nx.MultiGraph`` for undirected edges
- ``_directed``: ``nx.MultiDiGraph`` for directed edges
- ``_combined``: ``nx.MultiGraph`` combining all edges for connectivity analysis

Creating Graphs
---------------

Create a graph with both edge types:

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
   
   # Empty graph
   G = MixedMultiGraph()
   
   # Graph with undirected edges
   G = MixedMultiGraph(undirected_edges=[(1, 2), (2, 3)])
   
   # Graph with directed edges
   G = MixedMultiGraph(directed_edges=[(1, 2), (2, 3)])
   
   # Graph with both types
   G = MixedMultiGraph(
       directed_edges=[(1, 2)],
       undirected_edges=[(2, 3), (3, 4)]
   )
   
   # Graph with edge attributes
   G = MixedMultiGraph(
       directed_edges=[{'u': 1, 'v': 2, 'weight': 1.0}],
       undirected_edges=[{'u': 2, 'v': 3, 'weight': 2.0, 'label': 'tree_edge'}]
   )
   
   # Graph with graph attributes
   G = MixedMultiGraph(
       directed_edges=[(1, 2)],
       undirected_edges=[(2, 3)],
       attributes={'source': 'file.nex', 'version': '1.0'}
   )

Edge Operations
--------------

Add and remove edges:

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
   
   G = MixedMultiGraph()
   
   # Add directed edge (returns key)
   key1 = G.add_directed_edge(1, 2, weight=1.0)
   key2 = G.add_directed_edge(1, 2, weight=2.0)  # Parallel directed edge
   
   # Add undirected edge
   key3 = G.add_undirected_edge(2, 3, weight=1.5)
   key4 = G.add_undirected_edge(2, 3, weight=2.5)  # Parallel undirected edge
   
   # Remove directed edge
   G.remove_directed_edge(1, 2, key=key1)
   
   # Remove undirected edge
   G.remove_undirected_edge(2, 3, key=key3)
   
   # Check if edge exists
   has_dir = G.has_directed_edge(1, 2, key=key2)  # True
   has_undir = G.has_undirected_edge(2, 3, key=key4)  # True
   
   # Get edge data
   edge_data = G.get_directed_edge_data(1, 2, key=key2)
   # Returns: {'weight': 2.0}

Mutual Exclusivity
------------------

Edges between the same nodes must be either all directed or all undirected:

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
   
   G = MixedMultiGraph()
   
   # Add undirected edge
   G.add_undirected_edge(1, 2)
   
   # Adding a directed edge between same nodes removes undirected edge
   G.add_directed_edge(1, 2)  # Removes the undirected edge
   
   # Now only directed edge exists
   assert G.has_directed_edge(1, 2)
   assert not G.has_undirected_edge(1, 2)
   
   # Adding undirected edge removes directed edge
   G.add_undirected_edge(1, 2)  # Removes the directed edge
   assert G.has_undirected_edge(1, 2)
   assert not G.has_directed_edge(1, 2)

Node Operations
---------------

Work with nodes:

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
   
   G = MixedMultiGraph(
       directed_edges=[(1, 2)],
       undirected_edges=[(2, 3)]
   )
   
   # Add node with attributes
   G.add_node(4, label='node4', type='leaf')
   
   # Get node attributes
   node_data = G.nodes[4]  # Returns dict of attributes
   
   # Check if node exists
   has_node = G.has_node(4)  # True
   
   # Get all nodes
   all_nodes = list(G.nodes())  # [1, 2, 3, 4]
   
   # Get neighbors (both directed and undirected)
   neighbors = list(G.neighbors(2))  # [1, 3]
   
   # Get in-degree and out-degree (directed edges only)
   in_deg = G.in_degree(2)  # 1 (from directed edge 1->2)
   out_deg = G.out_degree(2)  # 0
   
   # Get total degree (both directed and undirected)
   total_deg = G.degree(2)  # 2

Graph Features
--------------

Analyze graph properties:

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
   from phylozoo.core.primitives.m_multigraph.features import (
       is_connected,
       number_of_connected_components,
       has_self_loops,
       source_components,
   )
   
   G = MixedMultiGraph(
       directed_edges=[(1, 2)],
       undirected_edges=[(2, 3), (4, 5)]
   )
   
   # Check connectivity
   connected = is_connected(G)  # False (two components)
   num_components = number_of_connected_components(G)  # 2
   
   # Get source components (components with no incoming directed edges)
   sources = source_components(G)  # List of sets of nodes
   
   # Check for self-loops
   has_loops = has_self_loops(G)  # False
   
   # Get number of edges and nodes
   num_edges = G.number_of_edges()  # 3
   num_nodes = G.number_of_nodes()  # 5

Isomorphism
-----------

Check if two graphs are isomorphic:

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
   from phylozoo.core.primitives.m_multigraph.isomorphism import is_isomorphic
   
   G1 = MixedMultiGraph(
       directed_edges=[(1, 2)],
       undirected_edges=[(2, 3)]
   )
   G2 = MixedMultiGraph(
       directed_edges=[(4, 5)],
       undirected_edges=[(5, 6)]
   )
   
   # Check isomorphism (ignoring node labels)
   isomorphic = is_isomorphic(G1, G2)  # True
   
   # Check with node attributes
   G1.add_node(1, label='A')
   G2.add_node(4, label='A')
   isomorphic = is_isomorphic(G1, G2, node_attrs=['label'])  # True
   
   # Check with edge attributes
   G3 = MixedMultiGraph(undirected_edges=[{'u': 1, 'v': 2, 'weight': 1.0}])
   G4 = MixedMultiGraph(undirected_edges=[{'u': 3, 'v': 4, 'weight': 1.0}])
   isomorphic = is_isomorphic(G3, G4, edge_attrs=['weight'])  # True

Conversions
-----------

Convert from NetworkX graphs:

.. code-block:: python

   import networkx as nx
   from phylozoo.core.primitives.m_multigraph.conversions import (
       graph_to_mixedmultigraph,
       multigraph_to_mixedmultigraph,
       multidigraph_to_mixedmultigraph,
   )
   
   # From NetworkX Graph (becomes undirected edges)
   nx_graph = nx.Graph()
   nx_graph.add_edge(1, 2, weight=1.0)
   mmg = graph_to_mixedmultigraph(nx_graph)
   
   # From NetworkX MultiGraph (becomes undirected edges)
   nx_multigraph = nx.MultiGraph()
   nx_multigraph.add_edge(1, 2, key=0, weight=1.0)
   mmg = multigraph_to_mixedmultigraph(nx_multigraph)
   
   # From NetworkX MultiDiGraph (becomes directed edges)
   nx_multidigraph = nx.MultiDiGraph()
   nx_multidigraph.add_edge(1, 2, key=0, weight=1.0)
   mmg = multidigraph_to_mixedmultigraph(nx_multidigraph)

I/O Operations
--------------

Save and load graphs:

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
   
   G = MixedMultiGraph(
       directed_edges=[(1, 2)],
       undirected_edges=[(2, 3)]
   )
   
   # Save to file (PhyloZoo DOT format)
   G.save('graph.pzdot')
   
   # Load from file
   G2 = MixedMultiGraph.load('graph.pzdot')
   
   # Convert to string
   pzdot_string = G.to_string(format='phylozoo-dot')
   
   # Create from string
   G3 = MixedMultiGraph.from_string(pzdot_string, format='phylozoo-dot')

Example
-------

Here's a complete example:

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
   from phylozoo.core.primitives.m_multigraph.features import (
       is_connected,
       source_components,
   )
   
   # Create a mixed graph
   G = MixedMultiGraph(
       directed_edges=[
           (1, 3),  # Hybrid edge
           (2, 3),  # Hybrid edge (parallel)
       ],
       undirected_edges=[
           (3, 4),  # Tree edge
           (4, 5),  # Tree edge
           (4, 6),  # Tree edge
       ],
       attributes={'name': 'example_mixed_graph'}
   )
   
   # Add node attributes
   G.add_node(5, label='A')
   G.add_node(6, label='B')
   
   # Check properties
   print(f"Nodes: {G.number_of_nodes()}")  # 6
   print(f"Edges: {G.number_of_edges()}")  # 5
   print(f"Connected: {is_connected(G)}")  # True
   
   # Get source components
   sources = source_components(G)
   print(f"Source components: {sources}")
   
   # Save graph
   G.save('example.pzdot')

.. seealso::
   For network classes that use MixedMultiGraph, see :doc:`Networks (Basic) <../networks/basic>`.
   For directed-only graphs, see :doc:`Directed Multi-Graph <directed_multigraph>`.
