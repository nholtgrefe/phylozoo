"""
Centralized matplotlib import with helpful error for optional viz dependency.

Visualization is optional; install with: pip install phylozoo[viz]
"""

from __future__ import annotations

from phylozoo.utils.exceptions import PhyloZooImportError

try:
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    from matplotlib.path import Path as MPath
except ImportError as e:
    raise PhyloZooImportError(
        "Visualization requires matplotlib. Install with: pip install phylozoo[viz]"
    ) from e

__all__ = ["mpatches", "plt", "Circle", "MPath"]
