Introduction
============

What is PhyloZoo?
-----------------

PhyloZoo is a Python package for working with phylogenetic networks. Phylogenetic networks 
extend phylogenetic trees by allowing for both divergence (splitting) and merging events, 
making them suitable for modeling processes like hybridization, horizontal gene transfer, 
and admixture.

Key Features
------------

* **Network Representation**: Support for both directed and semi-directed phylogenetic networks
* **Network Manipulation**: Tools for creating, modifying, and analyzing networks
* **Sequence Analysis**: Multiple sequence alignment (MSA) support with efficient bootstrapping
* **Network Inference**: Algorithms for inferring networks from data (e.g., SQuaRE)
* **Visualization**: Flexible plotting system for networks
* **Diversity Analysis**: Phylogenetic diversity calculations (Panda module)

Quick Start
-----------

Here's a simple example to get started:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   
   # Create a simple network
   network = DirectedPhyNetwork(
       nodes=[("root", "internal"), ("A", "leaf"), ("B", "leaf")],
       edges=[("root", "A"), ("root", "B")]
   )
   
   print(f"Network has {network.num_nodes} nodes and {network.num_edges} edges")

Core Concepts
-------------

**Directed Networks**: Fully directed phylogenetic networks with a single root, where 
hybrid nodes have in-degree >= 2.

**Semi-Directed Networks**: Networks with directed hybrid edges and undirected tree edges, 
allowing for more flexible representation of evolutionary relationships.

**Quartets**: Four-taxon unrooted trees, fundamental building blocks for network inference.

**Splits**: Bipartitions of taxa, used to represent evolutionary relationships.

For more detailed information, see the :doc:`Library <../library/index>` section.
