Directed Multi-Graph
====================

The ``DirectedMultiGraph`` class provides a directed multigraph where all edges are directed.
This class is used internally by ``DirectedPhyNetwork`` but can also be used directly for 
graph operations.

Overview
--------

.. figure:: ../../../../_static/images/example_directed_multigraph.png
   :alt: Example DirectedMultiGraph
   :width: 400px
   :align: center
   
   Example of a DirectedMultiGraph with parallel edges and cycles.

A ``DirectedMultiGraph`` supports:
- Parallel directed edges (multiple edges between the same nodes)
- Self-loops
- Edge attributes (weights, labels, etc.)
- Node attributes
- Graph-level attributes

The class uses composition with NetworkX:
- ``_graph``: ``nx.MultiDiGraph`` for directed edges
- ``_combined``: ``nx.MultiGraph`` combining all edges as undirected for connectivity analysis

Creating Graphs
---------------

Create a graph from edges:

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   
   # Empty graph
   G = DirectedMultiGraph()
   
   # Graph with edges
   G = DirectedMultiGraph(edges=[(1, 2), (2, 3), (3, 1)])
   
   # Graph with parallel edges
   G = DirectedMultiGraph(edges=[(1, 2), (1, 2)])  # Two edges from 1 to 2
   
   # Graph with edge attributes
   G = DirectedMultiGraph(edges=[
       {'u': 1, 'v': 2, 'weight': 1.0, 'label': 'edge1'},
       {'u': 2, 'v': 3, 'weight': 2.0}
   ])
   
   # Graph with graph attributes
   G = DirectedMultiGraph(
       edges=[(1, 2), (2, 3)],
       attributes={'source': 'file.nex', 'version': '1.0'}
   )

Edge Operations
--------------

Add and remove edges:

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   
   G = DirectedMultiGraph()
   
   # Add edge (returns key)
   key1 = G.add_edge(1, 2, weight=1.0)
   key2 = G.add_edge(1, 2, weight=2.0)  # Parallel edge
   
   # Add edge with explicit key
   key3 = G.add_edge(2, 3, key=5, label='custom')
   
   # Remove edge by key
   G.remove_edge(1, 2, key=key1)
   
   # Check if edge exists
   has_edge = G.has_edge(1, 2, key=key2)  # True
   
   # Get edge data
   edge_data = G.get_edge_data(1, 2, key=key2)
   # Returns: {'weight': 2.0}

Node Operations
---------------

Work with nodes:

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   
   G = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
   
   # Add node with attributes
   G.add_node(4, label='node4', type='leaf')
   
   # Get node attributes
   node_data = G.nodes[4]  # Returns dict of attributes
   
   # Check if node exists
   has_node = G.has_node(4)  # True
   
   # Get all nodes
   all_nodes = list(G.nodes())  # [1, 2, 3, 4]
   
   # Get neighbors
   neighbors = list(G.neighbors(2))  # [3]
   
   # Get in-degree and out-degree
   in_deg = G.in_degree(2)  # 1
   out_deg = G.out_degree(2)  # 1

Graph Features
--------------

Analyze graph properties:

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   from phylozoo.core.primitives.d_multigraph.features import (
       is_connected,
       number_of_connected_components,
       has_self_loops,
   )
   
   G = DirectedMultiGraph(edges=[(1, 2), (2, 3), (4, 5)])
   
   # Check connectivity
   connected = is_connected(G)  # False (two components)
   num_components = number_of_connected_components(G)  # 2
   
   # Check for self-loops
   has_loops = has_self_loops(G)  # False
   
   # Get number of edges and nodes
   num_edges = G.number_of_edges()  # 3
   num_nodes = G.number_of_nodes()  # 5

Isomorphism
-----------

Check if two graphs are isomorphic:

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   from phylozoo.core.primitives.d_multigraph.isomorphism import is_isomorphic
   
   G1 = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
   G2 = DirectedMultiGraph(edges=[(4, 5), (5, 6)])
   
   # Check isomorphism (ignoring node labels)
   isomorphic = is_isomorphic(G1, G2)  # True
   
   # Check with node attributes
   G1.add_node(1, label='A')
   G2.add_node(4, label='A')
   isomorphic = is_isomorphic(G1, G2, node_attrs=['label'])  # True
   
   # Check with edge attributes
   G3 = DirectedMultiGraph(edges=[{'u': 1, 'v': 2, 'weight': 1.0}])
   G4 = DirectedMultiGraph(edges=[{'u': 3, 'v': 4, 'weight': 1.0}])
   isomorphic = is_isomorphic(G3, G4, edge_attrs=['weight'])  # True

Conversions
-----------

Convert from NetworkX graphs:

.. code-block:: python

   import networkx as nx
   from phylozoo.core.primitives.d_multigraph.conversions import (
       digraph_to_directedmultigraph,
       multidigraph_to_directedmultigraph,
   )
   
   # From NetworkX DiGraph
   nx_graph = nx.DiGraph()
   nx_graph.add_edge(1, 2, weight=1.0)
   dmg = digraph_to_directedmultigraph(nx_graph)
   
   # From NetworkX MultiDiGraph
   nx_multigraph = nx.MultiDiGraph()
   nx_multigraph.add_edge(1, 2, key=0, weight=1.0)
   nx_multigraph.add_edge(1, 2, key=1, weight=2.0)
   dmg = multidigraph_to_directedmultigraph(nx_multigraph)

I/O Operations
--------------

Save and load graphs:

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   
   G = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
   
   # Save to file (DOT format)
   G.save('graph.dot')
   
   # Load from file
   G2 = DirectedMultiGraph.load('graph.dot')
   
   # Convert to string
   dot_string = G.to_string(format='dot')
   
   # Create from string
   G3 = DirectedMultiGraph.from_string(dot_string, format='dot')

Example
-------

Here's a complete example:

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   from phylozoo.core.primitives.d_multigraph.features import is_connected
   
   # Create a directed graph with parallel edges
   G = DirectedMultiGraph(
       edges=[
           (1, 2),
           (1, 2),  # Parallel edge
           (2, 3),
           (3, 4),
           (4, 1),  # Creates a cycle
       ],
       attributes={'name': 'example_graph'}
   )
   
   # Add node attributes
   G.add_node(1, label='A')
   G.add_node(2, label='B')
   
   # Check properties
   print(f"Nodes: {G.number_of_nodes()}")  # 4
   print(f"Edges: {G.number_of_edges()}")  # 5
   print(f"Connected: {is_connected(G)}")  # True
   
   # Save graph
   G.save('example.dot')

.. seealso::
   For network classes that use DirectedMultiGraph, see :doc:`Networks (Basic) <../networks/basic>`.
   For mixed graphs (with undirected edges), see :doc:`Mixed Multi-Graph <mixed_multigraph>`.
