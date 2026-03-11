
[![PyPI](https://img.shields.io/pypi/v/phylozoo)](https://pypi.org/project/phylozoo/)
[![License](https://img.shields.io/github/license/nholtgrefe/phylozoo)](https://github.com/nholtgrefe/phylozoo/blob/main/LICENSE.md)
[![Docs](https://img.shields.io/badge/docs-stable-blue)](https://nholtgrefe.github.io/phylozoo/)

# PhyloZoo <img src="https://raw.githubusercontent.com/nholtgrefe/phylozoo/main/docs/source/_static/phylozoo_full.svg" alt="PhyloZoo" width="300" align="right">

PhyloZoo is a Python package for working with phylogenetic networks and related evolutionary
data types. PhyloZoo aims to provide the foundational infrastructure for
phylogenetic network analysis in Python — a common framework that other packages can build on.

## Key Features

- **Directed & semi-directed networks** — represent phylogenetic networks as fully directed rooted DAGs or as semi-directed/mixed graphs that allow root uncertainty. Both representations are validated on construction to guarantee well-formed phylogenetic objects. Includes a rich library of operations: network classifications, generators, conversions between representations, and much more.
- **Quartets, splits & distance matrices** — support for quartet systems, split systems, and pairwise distance matrices: the core building blocks for phylogenetic inference and comparison.
- **Multiple sequence alignments** — store and manipulate sequence data with efficient NumPy-backed arrays, including bootstrapping and site-pattern extraction.
- **Flexible visualization** — plot networks with different layouts and fine-grained control over styling, labels, and coloring via Matplotlib.
- **Standard file formats** — read and write common phylogenetic formats including eNewick, DOT, FASTA, and NEXUS, making it easy to integrate with existing workflows.
- **Performance** — leverages NumPy and optional Numba JIT compilation for computationally intensive algorithms.


## Installation

To install the recommended version that includes vizualization, do:

```bash
pip install phylozoo[viz]
```

## Documentation

For detailed documentation, installation instructions, tutorials, and API reference, visit the **[PhyloZoo docs](https://nholtgrefe.github.io/phylozoo/)**.

## Citation

If you use PhyloZoo in your research, please cite:

> Niels Holtgrefe (2026). *PhyloZoo*. Available at: https://github.com/nholtgrefe/phylozoo
