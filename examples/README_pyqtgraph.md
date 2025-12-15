# PyQtGraph Backend Examples

This directory contains examples demonstrating the PyQtGraph backend for interactive network visualization.

## Requirements

Install PyQtGraph (and its dependencies):
```bash
pip install pyqtgraph
```

Note: PyQtGraph requires either PyQt5 or PyQt6. If you don't have them installed, PyQtGraph will attempt to install PyQt5 automatically, or you can install it manually:
```bash
pip install pyqt5
# or
pip install pyqt6
```

## Usage

### Simple Example

Run the basic example:
```bash
python examples/pyqtgraph_example.py
```

This will open an interactive window where you can:
- **Pan**: Left-click and drag
- **Zoom**: Right-click and drag, or use mouse wheel
- **Close**: Close the window to exit

### Programmatic Usage

```python
from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.visualize.network_plot import plot_network_with_layout
from PyQt5.QtWidgets import QApplication
import sys

# Create network
net = DirectedPhyNetwork(
    edges=[(7, 5), (7, 6), (5, 1), (5, 2), (6, 3), (6, 4)],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), 
           (3, {'label': 'C'}), (4, {'label': 'D'})]
)

# Create QApplication (required for PyQtGraph)
app = QApplication.instance()
if app is None:
    app = QApplication(sys.argv)

# Plot with PyQtGraph backend
win = plot_network_with_layout(
    net,
    backend='pyqtgraph',
    orientation='top-bottom',
    with_labels=True,
)

# Run event loop
sys.exit(app.exec_())
```

## Features

- **Interactive**: Zoom, pan, and interact with plots in real-time
- **High Performance**: Fast rendering for large networks
- **Auto-detection**: Falls back to matplotlib if PyQtGraph is not available
- **Same API**: Uses the same `plot_network_with_layout()` function with `backend='pyqtgraph'`

## Comparison with Matplotlib

| Feature | Matplotlib | PyQtGraph |
|---------|-----------|-----------|
| Interactive | No | Yes |
| Export to file | Yes | Yes (via screenshot) |
| Performance | Good | Excellent |
| Dependencies | Built-in | Optional |

## Notes

- PyQtGraph windows are interactive and stay open until closed
- The QApplication event loop must be running for the window to be responsive
- If PyQtGraph is not available, the function will automatically fall back to matplotlib

