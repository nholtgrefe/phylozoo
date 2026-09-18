
[![PyPI](https://img.shields.io/pypi/v/phylozoo)](https://pypi.org/project/phylozoo/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue)](https://github.com/nholtgrefe/phylozoo/blob/master/LICENSE.md)
[![CI](https://github.com/nholtgrefe/phylozoo/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/nholtgrefe/phylozoo/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-stable-blue)](https://nholtgrefe.github.io/phylozoo/)

<img src="https://github.com/nholtgrefe/phylozoo/blob/master/docs/source/_static/phylozoo_compact.png" alt="PhyloZoo logo" width="375" align="left">

PhyloZoo is a Python package for working with phylogenetic networks and related evolutionary
data types. PhyloZoo aims to provide the foundational infrastructure for
phylogenetic network analysis in Python—a common framework that other packages can build on.

<br>

## Key Features

- **Directed & semi-directed networks** — represent phylogenetic networks as fully directed rooted DAGs or as semi-directed/mixed graphs that allow root uncertainty. Both representations are validated on construction to guarantee well-formed phylogenetic objects. Includes a rich library of operations: network classifications, generators, conversions between representations, and much more.
- **Quartets, triplets, splits & distance matrices** — support for quartet and triplet systems, split systems, and pairwise distance matrices: the core building blocks for phylogenetic inference and comparison.
- **Multiple sequence alignments** — store and manipulate sequence data with efficient NumPy-backed arrays, including bootstrapping and site-pattern extraction.
- **Flexible visualization** — plot networks with different layouts and fine-grained control over styling, labels, and coloring via Matplotlib.
- **Standard file formats** — read and write common phylogenetic formats including eNewick, FASTA, NEXUS, and PHYLIP, making it easy to integrate with existing workflows.

## Installation

To install the recommended version that includes vizualization, do:

```bash
pip install phylozoo[viz]
```

## Quickstart

Build a small rooted network with a hybrid node, inspect it, save it and plot it:

```python
import phylozoo as pz

# A rooted network on four taxa with one hybrid node (edges are parent-child tuples)
network = pz.DirectedPhyNetwork(
    edges=[
        ("root", "u1"), ("root", "u2"),
        ("u1", "A"), ("u1", "x"), ("x", "B"),
        ("x", "h"), ("u2", "h"),        # the two edges into the hybrid node h
        ("u2", "D"), ("h", "C"),
    ],
)

print(network.leaves)                            # {'A', 'B', 'C', 'D'}
print(pz.dnetwork.classifications.level(network))  # 1

# Save network to file (eNewick format)
network.save("my_network.enewick")               # ((A,(B,(C)#H1)),(D,#H1));

# Draw it as text (no matplotlib needed)
network.pretty_print()
#         ┌───────────────────────● A
#         │
# ┌───────●       ┌───────────────● B
# │       └───────●
# ○               └┄┄┄┄┄┄>◆───────● C
# │                       ^
# └───────────────────────●───────● D

# Plot the network
from phylozoo.viz import plot

plot(network, show=True)
```

See the [quickstart tutorial](https://nholtgrefe.github.io/phylozoo/tutorials/quickstart.html)
for a longer tour of what PhyloZoo can do.

## Documentation

For detailed documentation, installation instructions, tutorials, and API reference, visit the **[PhyloZoo docs](https://nholtgrefe.github.io/phylozoo/)**.

## Citation

If you use PhyloZoo in your research, please cite:

> Niels Holtgrefe. PhyloZoo: a unified framework for phylogenetic network analysis in Python. bioRxiv:[10.64898/2026.06.09.731120](https://doi.org/10.64898/2026.06.09.731120), 2026.
