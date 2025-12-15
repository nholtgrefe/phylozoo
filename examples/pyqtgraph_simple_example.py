"""
Simple PyQtGraph example for network plotting.

This demonstrates the PyQtGraph backend for interactive visualization.
Run this script to see an interactive plot window.

Note: Requires PyQtGraph to be installed: pip install pyqtgraph
"""

import sys

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.visualize.network_plot import plot_network_with_layout

# Check if PyQtGraph is available
try:
    import pyqtgraph as pg
    from PyQt5.QtWidgets import QApplication
    HAS_PYQTGAPH = True
except ImportError:
    HAS_PYQTGAPH = False
    print("ERROR: PyQtGraph is not installed.")
    print("Install it with: pip install pyqtgraph")
    sys.exit(1)


def main():
    """Create and display a network using PyQtGraph."""
    print("Creating network with PyQtGraph backend...")
    
    # Create a simple network
    net = DirectedPhyNetwork(
        edges=[
            (7, 5), (7, 6),
            (5, 1), (5, 2),
            (6, 3), (6, 4),
        ],
        nodes=[
            (1, {'label': 'A'}),
            (2, {'label': 'B'}),
            (3, {'label': 'C'}),
            (4, {'label': 'D'}),
        ]
    )
    
    # Create QApplication if it doesn't exist
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    
    # Plot with PyQtGraph backend
    print("Opening PyQtGraph window...")
    print("  - Use mouse to zoom and pan")
    print("  - Close window to exit")
    
    win = plot_network_with_layout(
        net,
        backend='pyqtgraph',
        orientation='top-bottom',
        with_labels=True,
        node_size=15,
        leaf_size=18,
    )
    
    # Run event loop
    print("\nWindow opened. Interact with the plot, then close to exit.")
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()

