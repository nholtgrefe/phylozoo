"""
Examples demonstrating PyQtGraph backend for network plotting.

This script shows how to use the PyQtGraph backend for interactive
visualization of phylogenetic networks.

Note: PyQtGraph is an optional dependency. If not installed, the examples
will fall back to matplotlib or raise an ImportError.
"""

import sys

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
from phylozoo.visualize.network_plot import plot_network_with_layout


def example_directed_pyqtgraph():
    """Example: Directed network with PyQtGraph backend."""
    print("Creating directed network with PyQtGraph...")
    
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
    
    try:
        win = plot_network_with_layout(
            net,
            backend='pyqtgraph',
            orientation='top-bottom',
            with_labels=True,
        )
        print("  PyQtGraph window opened. Close the window to continue.")
        print("  Note: The window will stay open until you close it.")
        
        # Start Qt event loop if running as script
        if __name__ == '__main__':
            from PyQt5.QtWidgets import QApplication
            app = QApplication.instance()
            if app is None:
                app = QApplication(sys.argv)
            app.exec_()
        
        return win
    except ImportError as e:
        print(f"  PyQtGraph not available: {e}")
        print("  Install with: pip install pyqtgraph")
        return None


def example_hybrid_network_pyqtgraph():
    """Example: Network with hybrid node using PyQtGraph."""
    print("Creating hybrid network with PyQtGraph...")
    
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
    
    try:
        win = plot_network_with_layout(
            net,
            backend='pyqtgraph',
            orientation='top-bottom',
            hybrid_color='coral',
            hybrid_edge_color='darkred',
            with_labels=True,
        )
        print("  PyQtGraph window opened. Close the window to continue.")
        
        if __name__ == '__main__':
            from PyQt5.QtWidgets import QApplication
            app = QApplication.instance()
            if app is None:
                app = QApplication(sys.argv)
            app.exec_()
        
        return win
    except ImportError as e:
        print(f"  PyQtGraph not available: {e}")
        print("  Install with: pip install pyqtgraph")
        return None


def example_semidirected_pyqtgraph():
    """Example: Semi-directed network with PyQtGraph."""
    print("Creating semi-directed network with PyQtGraph...")
    
    net = SemiDirectedPhyNetwork(
        undirected_edges=[
            (5, 1),
            (6, 2), (6, 3), (6, 4),
            (6, 7), (7, 8),
        ],
        directed_edges=[
            {'u': 6, 'v': 5, 'gamma': 0.6},
            {'u': 7, 'v': 5, 'gamma': 0.4},
        ],
        nodes=[
            (1, {'label': 'A'}),
            (2, {'label': 'B'}),
            (3, {'label': 'C'}),
            (4, {'label': 'D'}),
            (8, {'label': 'E'}),
        ]
    )
    
    try:
        win = plot_network_with_layout(
            net,
            backend='pyqtgraph',
            leaf_color='lightgreen',
            hybrid_color='salmon',
            with_labels=True,
        )
        print("  PyQtGraph window opened. Close the window to continue.")
        
        if __name__ == '__main__':
            from PyQt5.QtWidgets import QApplication
            app = QApplication.instance()
            if app is None:
                app = QApplication(sys.argv)
            app.exec_()
        
        return win
    except ImportError as e:
        print(f"  PyQtGraph not available: {e}")
        print("  Install with: pip install pyqtgraph")
        return None


def example_interactive_comparison():
    """Example: Compare matplotlib and pyqtgraph side by side."""
    print("Creating comparison example...")
    
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
    
    print("  Matplotlib version:")
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 6))
        plot_network_with_layout(
            net,
            ax=ax,
            backend='matplotlib',
            orientation='top-bottom',
        )
        plt.savefig('example_matplotlib_comparison.png', dpi=150, bbox_inches='tight')
        print("    Saved: example_matplotlib_comparison.png")
        plt.close()
    except Exception as e:
        print(f"    Error: {e}")
    
    print("  PyQtGraph version:")
    try:
        win = plot_network_with_layout(
            net,
            backend='pyqtgraph',
            orientation='top-bottom',
        )
        print("    PyQtGraph window opened (interactive)")
        
        if __name__ == '__main__':
            from PyQt5.QtWidgets import QApplication
            app = QApplication.instance()
            if app is None:
                app = QApplication(sys.argv)
            app.exec_()
        
        return win
    except ImportError as e:
        print(f"    PyQtGraph not available: {e}")
        print("    Install with: pip install pyqtgraph")
        return None


def main():
    """Run PyQtGraph examples."""
    print("=" * 60)
    print("PyQtGraph Examples - Interactive Network Plotting")
    print("=" * 60)
    print()
    
    # Check if pyqtgraph is available
    try:
        import pyqtgraph
        print(f"PyQtGraph version: {pyqtgraph.__version__}")
        print()
    except ImportError:
        print("WARNING: PyQtGraph is not installed.")
        print("Install it with: pip install pyqtgraph")
        print("Examples will show ImportError messages.")
        print()
    
    try:
        # Run examples
        print("Example 1: Simple directed network")
        example_directed_pyqtgraph()
        print()
        
        print("Example 2: Network with hybrid node")
        example_hybrid_network_pyqtgraph()
        print()
        
        print("Example 3: Semi-directed network")
        example_semidirected_pyqtgraph()
        print()
        
        print("Example 4: Comparison (matplotlib vs pyqtgraph)")
        example_interactive_comparison()
        print()
        
        print("=" * 60)
        print("All PyQtGraph examples completed!")
        print("=" * 60)
        print()
        print("Note: PyQtGraph windows are interactive - you can zoom, pan,")
        print("      and interact with the plots in real-time.")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()

