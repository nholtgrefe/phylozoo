"""
PyQtGraph Example for Network Plotting

This example demonstrates how to use the PyQtGraph backend for interactive
visualization of phylogenetic networks.

Usage:
    python examples/pyqtgraph_example.py

Features:
    - Interactive zooming and panning
    - Real-time interaction
    - High-performance rendering

Requirements:
    pip install pyqtgraph
"""

import sys

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.visualize.network_plot import plot_network_with_layout

# Check for PyQtGraph
try:
    from PyQt5.QtWidgets import QApplication
    HAS_PYQT = True
except ImportError:
    try:
        from PyQt6.QtWidgets import QApplication
        HAS_PYQT = True
    except ImportError:
        HAS_PYQT = False
        print("ERROR: PyQt5 or PyQt6 is required for PyQtGraph.")
        print("Install with: pip install pyqtgraph")
        sys.exit(1)


def main():
    """Create and display a network using PyQtGraph."""
    print("=" * 60)
    print("PyQtGraph Example - Interactive Network Plotting")
    print("=" * 60)
    print()
    
    # Create a network with hybrid node
    print("Creating network...")
    net = DirectedPhyNetwork(
        edges=[
            (8, 5), (8, 6),
            (5, 1), (5, 2),
            (6, 3), (6, 9),
            {'u': 5, 'v': 4, 'gamma': 0.6},
            {'u': 6, 'v': 4, 'gamma': 0.4},
            (4, 10),
        ],
        nodes=[
            (1, {'label': 'A'}),
            (2, {'label': 'B'}),
            (3, {'label': 'C'}),
            (9, {'label': 'D'}),
            (10, {'label': 'E'}),
        ]
    )
    
    # Create QApplication
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    
    print("Opening PyQtGraph window...")
    print("  - Left-click and drag to pan")
    print("  - Right-click and drag to zoom")
    print("  - Mouse wheel to zoom in/out")
    print("  - Close window to exit")
    print()
    
    # Plot with PyQtGraph backend
    win = plot_network_with_layout(
        net,
        backend='pyqtgraph',
        orientation='top-bottom',
        with_labels=True,
        node_size=12,
        leaf_size=15,
        hybrid_color='coral',
        hybrid_edge_color='darkred',
    )
    
    print("Window opened! Interact with the plot.")
    print("Close the window when done.")
    print()
    
    # Run event loop
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()

