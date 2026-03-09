"""
Shared layout utilities for PhyloZoo visualization.

This module provides common layout computation (NetworkX/Graphviz dispatch)
and position normalization used by all layout modules. Layouts determine
only node positions (placement); node sizes come from the style.
"""

from __future__ import annotations

from typing import Any, TypeVar

import networkx as nx

from phylozoo.utils.exceptions import (
    PhyloZooImportError,
    PhyloZooLayoutError,
)

T = TypeVar('T')

GRAPHVIZ_LAYOUTS = ('dot', 'neato', 'fdp', 'sfdp', 'twopi', 'circo')

NETWORKX_LAYOUTS = {
    'spring': nx.spring_layout,
    'circular': nx.circular_layout,
    'kamada_kawai': nx.kamada_kawai_layout,
    'planar': nx.planar_layout,
    'random': nx.random_layout,
    'shell': nx.shell_layout,
    'spectral': nx.spectral_layout,
    'spiral': nx.spiral_layout,
    'bipartite': nx.bipartite_layout,
}


def compute_nx_positions(
    G: nx.Graph,
    layout: str = 'spring',
    **kwargs: Any,
) -> dict[Any, tuple[float, float]]:
    """
    Compute node positions using NetworkX or Graphviz.

    Parameters
    ----------
    G : nx.Graph
        NetworkX graph (DiGraph, Graph, MultiGraph, etc.).
    layout : str, optional
        Layout algorithm name. NetworkX: 'spring', 'circular', 'kamada_kawai',
        'planar', 'random', 'shell', 'spectral', 'spiral', 'bipartite'.
        Graphviz: 'dot', 'neato', 'fdp', 'sfdp', 'twopi', 'circo'.
        By default 'spring'.
    **kwargs
        Additional parameters passed to the layout algorithm.

    Returns
    -------
    dict
        Node ID -> (x, y) position mapping.

    Raises
    ------
    PhyloZooLayoutError
        If layout is not supported or computation fails.
    PhyloZooImportError
        If Graphviz layout requested but pygraphviz not installed.
    """
    pos: dict[Any, tuple[float, float]]

    if layout in GRAPHVIZ_LAYOUTS:
        try:
            pos = nx.nx_agraph.graphviz_layout(G, prog=layout, **kwargs)
        except ImportError:
            raise PhyloZooImportError(
                f"Graphviz layout '{layout}' requires pygraphviz. "
                "Install with: pip install pygraphviz"
            )
        except Exception as e:
            raise PhyloZooLayoutError(
                f"Graphviz layout '{layout}' failed: {e}"
            ) from e
    elif layout in NETWORKX_LAYOUTS:
        pos = NETWORKX_LAYOUTS[layout](G, **kwargs)
    else:
        supported = ', '.join(
            sorted(NETWORKX_LAYOUTS.keys()) + list(GRAPHVIZ_LAYOUTS)
        )
        raise PhyloZooLayoutError(
            f"Unsupported layout algorithm: '{layout}'. "
            f"Supported: {supported}"
        )

    return pos


def normalize_positions(
    pos: dict[T, tuple[float, float]],
) -> dict[T, tuple[float, float]]:
    """
    Center positions at origin and scale to unit range.

    Parameters
    ----------
    pos : dict
        Node ID -> (x, y) position mapping.

    Returns
    -------
    dict
        Normalized positions (centered, max dimension = 1).
    """
    if not pos:
        return pos

    xs = [x for x, _ in pos.values()]
    ys = [y for _, y in pos.values()]
    if not xs or not ys:
        return pos

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    width = max_x - min_x if max_x != min_x else 1.0
    height = max_y - min_y if max_y != min_y else 1.0
    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2
    scale = 1.0 / max(width, height) if max(width, height) > 0 else 1.0

    return {
        node: ((x - center_x) * scale, (y - center_y) * scale)
        for node, (x, y) in pos.items()
    }


def compute_layout_center(
    pos: dict[T, tuple[float, float]],
) -> tuple[float, float]:
    """
    Compute the center of a layout from node positions.

    Parameters
    ----------
    pos : dict
        Node ID -> (x, y) position mapping.

    Returns
    -------
    tuple[float, float]
        (center_x, center_y) of the bounding box.
    """
    if not pos:
        return (0.0, 0.0)
    xs = [x for x, _ in pos.values()]
    ys = [y for _, y in pos.values()]
    return ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2)
