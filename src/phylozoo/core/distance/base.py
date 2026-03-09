"""
Distance matrix base module.

This module provides the core DistanceMatrix class for working with distance matrices.
"""

from __future__ import annotations

from typing import TypeVar

import numpy as np

from ...utils.exceptions import PhyloZooValueError
from phylozoo.utils.io import IOMixin

T = TypeVar('T')


class DistanceMatrix(IOMixin):
    """
    An immutable distance matrix.
    
    A DistanceMatrix represents pairwise distances between a set of labeled items.
    The matrix is stored as a symmetric numpy array and is immutable after initialization.
    
    Parameters
    ----------
    distance_matrix : numpy.ndarray
        A symmetric square 2D numpy array representing pairwise distances.
        Must be square and symmetric.
    labels : list[T] | None, optional
        List of labels corresponding to the rows/columns of the distance matrix.
        If None, defaults to `[0, 1, 2, ..., n-1]` where n is the matrix size.
        By default None.
    
    Notes
    -----
    The class is immutable after initialization. To create a modified version,
    create a new DistanceMatrix instance with the modified data.

    Supported I/O formats:

    - ``nexus`` (default): ``.nexus``, ``.nex``, ``.nxs``
    - ``phylip``: ``.phy``, ``.phylip``
    - ``csv``: ``.csv``
    
    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> 
    >>> # Create from numpy array
    >>> matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
    >>> dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
    >>> len(dm)
    3
    >>> dm.get_distance('A', 'B')
    1.0
    
    >>> # Default labels (0, 1, 2, ...)
    >>> dm2 = DistanceMatrix(matrix)
    >>> dm2.labels
    (0, 1, 2)
    
    Attributes
    ----------
    np_array : numpy.ndarray
        Read-only access to the underlying numpy array.
    labels : tuple[T, ...]
        Tuple of labels corresponding to the rows/columns (immutable).
    indices : tuple[int, ...]
        Tuple of indices `(0, 1, 2, ..., len(self)-1)` (immutable).
    """
    
    # Minimal class-level defaults so static analyzers recognize attributes.
    # These are immediately replaced in __init__ for each instance.
    _matrix = np.empty((0, 0))
    _labels = ()
    _indices = ()
    _label_to_index = {}
    
    # I/O format configuration
    _default_format = 'nexus'
    _supported_formats = ['nexus', 'phylip', 'csv']

    def __init__(
        self,
        distance_matrix: np.ndarray,
        labels: list[T] | None = None,
    ) -> None:
        """
        Initialize a distance matrix.
        
        Parameters
        ----------
        distance_matrix : numpy.ndarray
            A symmetric square 2D numpy array representing pairwise distances.
        labels : list[T] | None, optional
            List of labels corresponding to the rows/columns. If None, uses
            `[0, 1, 2, ..., n-1]`. By default None.
        
        Raises
        ------
        TypeError
            If distance_matrix is not a numpy ndarray.
        PhyloZooValueError
            If distance_matrix is not square, not symmetric, or dimensions don't
            match the number of labels, or if labels are not unique.
        """
        # Type checking
        if not isinstance(distance_matrix, np.ndarray):
            raise TypeError("Distance matrix must be a numpy ndarray")
        
        # Shape validation
        if len(distance_matrix.shape) != 2:
            raise PhyloZooValueError("Distance matrix must be a 2D array")
        if distance_matrix.shape[0] != distance_matrix.shape[1]:
            raise PhyloZooValueError("Distance matrix must be square")
        
        # Symmetry validation (with tolerance for floating point)
        if not np.allclose(distance_matrix, distance_matrix.T, rtol=1e-10, atol=1e-10):
            raise PhyloZooValueError("Distance matrix must be symmetric")
        
        # Convert to float64 and make contiguous for better performance
        matrix = np.ascontiguousarray(distance_matrix, dtype=np.float64)
        
        # Make matrix read-only for immutability
        matrix.setflags(write=False)
        
        # Store as private attribute (immutability)
        self._matrix = matrix
        n = int(matrix.shape[0])

        # Handle labels (assume labels are hashable)
        if labels is None:
            self._labels = tuple(range(n))
        else:
            labels_list = list(labels)
            if len(labels_list) != n:
                raise PhyloZooValueError(
                    f"Number of labels ({len(labels_list)}) must match "
                    f"matrix size ({n})"
                )
            # Ensure labels are unique to avoid ambiguous lookups
            if len(set(labels_list)) != len(labels_list):
                raise PhyloZooValueError("Labels must be unique")
            # Store as tuple for immutability
            self._labels = tuple(labels_list)

        # Build O(1) lookup mapping for label -> index (labels must be hashable)
        self._label_to_index = {lbl: idx for idx, lbl in enumerate(self._labels)}
        
        # Store indices as tuple
        self._indices = tuple(range(n))
    
    @property
    def np_array(self) -> np.ndarray:
        """
        Get the underlying numpy array (read-only).
        
        Returns
        -------
        numpy.ndarray
            The distance matrix as a read-only numpy array.
        
        Notes
        -----
        The returned array is read-only. To modify, create a new DistanceMatrix.
        
        Examples
        --------
        >>> import numpy as np
        >>> dm = DistanceMatrix(np.array([[0, 1], [1, 0]]), labels=['A', 'B'])
        >>> arr = dm.np_array
        >>> arr[0, 1]
        1.0
        """
        return self._matrix
    
    @property
    def labels(self) -> tuple[T, ...]:
        """
        Get the labels corresponding to rows/columns.
        
        Returns
        -------
        tuple[T, ...]
            Tuple of labels (immutable).
        """
        return self._labels
    
    @property
    def indices(self) -> tuple[int, ...]:
        """
        Get the indices (0, 1, 2, ..., len(self)-1).
        
        Returns
        -------
        tuple[int, ...]
            Tuple of indices (immutable).
        """
        return self._indices
    
    def get_index(self, label: T) -> int:
        """
        Get the index of a label in the distance matrix.
        
        Parameters
        ----------
        label : T
            Label to look up.
        
        Returns
        -------
        int
            Index of the label in the matrix.
        
        Raises
        ------
        PhyloZooValueError
            If label is not found in the distance matrix.
        
        Examples
        --------
        >>> import numpy as np
        >>> dm = DistanceMatrix(np.array([[0, 1], [1, 0]]), labels=['A', 'B'])
        >>> dm.get_index('A')
        0
        >>> dm.get_index('B')
        1
        """
        try:
            return self._label_to_index[label]
        except KeyError as exc:
            raise PhyloZooValueError(f"Label {label} not found in distance matrix") from exc
    
    def get_distance(self, label1: T, label2: T) -> float:
        """
        Get the distance between two labels.
        
        Parameters
        ----------
        label1 : T
            First label.
        label2 : T
            Second label.
        
        Returns
        -------
        float
            Distance between the two labels.
        
        Raises
        ------
        ValueError
            If either label is not found in the distance matrix.
        
        Examples
        --------
        >>> import numpy as np
        >>> dm = DistanceMatrix(np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]]), 
        ...                     labels=['A', 'B', 'C'])
        >>> dm.get_distance('A', 'B')
        1.0
        >>> dm.get_distance('B', 'A')  # Symmetric
        1.0
        >>> dm.get_distance('A', 'C')
        2.0
        """
        idx1 = self.get_index(label1)
        idx2 = self.get_index(label2)
        return float(self._matrix[idx1, idx2])
    
    def copy(self) -> 'DistanceMatrix':
        """
        Create a copy of the distance matrix.
        
        Returns
        -------
        DistanceMatrix
            A new DistanceMatrix instance with copied data (also immutable).
        
        Examples
        --------
        >>> import numpy as np
        >>> dm1 = DistanceMatrix(np.array([[0, 1], [1, 0]]), labels=['A', 'B'])
        >>> dm2 = dm1.copy()
        >>> dm1 is dm2
        False
        >>> dm1.get_distance('A', 'B') == dm2.get_distance('A', 'B')
        True
        """
        # Create new instance with copied matrix
        new_dm = DistanceMatrix.__new__(DistanceMatrix)
        matrix_copy = self._matrix.copy()
        matrix_copy.setflags(write=False)  # Make copy also immutable
        new_dm._matrix = matrix_copy
        new_dm._labels = self._labels  # Tuples are immutable, so safe to share
        new_dm._indices = self._indices  # Tuples are immutable, so safe to share
        new_dm._label_to_index = dict(self._label_to_index)
        return new_dm
    
    def __len__(self) -> int:
        """
        Return the size of the distance matrix.
        
        Returns
        -------
        int
            Number of rows/columns.
        """
        return self._matrix.shape[0]
    
    def __contains__(self, label: T) -> bool:
        """
        Check if a label is in the distance matrix.
        
        Parameters
        ----------
        label : T
            Label to check.
        
        Returns
        -------
        bool
            True if label is in the matrix, False otherwise.
        """
        # Use mapping for consistent/fast membership checks
        return label in self._label_to_index
    
    def __repr__(self) -> str:
        """
        Return string representation of the distance matrix.
        
        Returns
        -------
        str
            String representation.
        """
        return f"DistanceMatrix(size={len(self)}, labels={self._labels})"
    
    def __str__(self) -> str:
        """
        Return human-readable string representation.
        
        For small matrices (up to 10 elements), prints the full upper triangle.
        For larger matrices, truncates the display. Always includes element names.
        
        Returns
        -------
        str
            Human-readable string with matrix contents (upper triangle only).
        """
        n = len(self)
        max_display = 10
        
        # Format label strings
        label_strs = [str(label) for label in self._labels]
        if n > 0:
            max_label_width = max(len(s) for s in label_strs)
        else:
            max_label_width = 0
        
        lines = []
        
        if n == 0:
            lines.append("DistanceMatrix (empty)")
            return "\n".join(lines)
        
        # Header
        lines.append("DistanceMatrix:")
        
        if n <= max_display:
            # Print full upper triangle
            # Calculate column widths for alignment (need to fit both labels and values)
            col_widths = []
            for j in range(n):
                # Width needed for label
                label_width = len(label_strs[j])
                # Width needed for values in this column
                max_val_width = label_width  # Start with label width
                for i in range(j):  # Only upper triangle
                    val_str = f"{self._matrix[i, j]:.6f}".rstrip('0').rstrip('.')
                    max_val_width = max(max_val_width, len(val_str))
                # Ensure minimum width for readability
                col_widths.append(max(max_val_width, 6))
            
            # Print header row with labels
            header = " " * (max_label_width + 2)
            for j in range(n):
                header += f"{label_strs[j]:>{col_widths[j] + 1}}"
            lines.append(header)
            
            # Print each row (upper triangle only)
            for i in range(n):
                row_str = f"{label_strs[i]:>{max_label_width}}  "
                for j in range(n):
                    if j <= i:
                        # Lower triangle or diagonal: print spaces
                        row_str += " " * (col_widths[j] + 1)
                    else:
                        # Upper triangle: print value
                        val_str = f"{self._matrix[i, j]:.6f}".rstrip('0').rstrip('.')
                        row_str += f"{val_str:>{col_widths[j] + 1}}"
                lines.append(row_str)
        else:
            # Truncated display
            num_cols = max_display // 2
            col_width = 10
            
            # Print header with first few and last few columns
            header = " " * (max_label_width + 2)
            # First few columns
            for j in range(num_cols):
                header += f"{label_strs[j]:>{col_width}}"
            header += "  ...  "
            # Last few columns
            for j in range(n - num_cols, n):
                header += f"{label_strs[j]:>{col_width}}"
            lines.append(header)
            
            # Print rows (upper triangle only, truncated)
            for i in range(n):
                if i < num_cols:
                    # First few rows
                    row_str = f"{label_strs[i]:>{max_label_width}}  "
                    # First few columns
                    for j in range(num_cols):
                        if j <= i:
                            row_str += " " * col_width
                        else:
                            val_str = f"{self._matrix[i, j]:.4f}".rstrip('0').rstrip('.')
                            row_str += f"{val_str:>{col_width}}"
                    row_str += "  ...  "
                    # Last few columns (only if in upper triangle)
                    for j in range(n - num_cols, n):
                        if j <= i:
                            row_str += " " * col_width
                        else:
                            val_str = f"{self._matrix[i, j]:.4f}".rstrip('0').rstrip('.')
                            row_str += f"{val_str:>{col_width}}"
                    lines.append(row_str)
                elif i == num_cols:
                    # Ellipsis row
                    row_str = " " * (max_label_width + 2)
                    row_str += " " * (num_cols * col_width) + "  ...  " + " " * (num_cols * col_width)
                    lines.append(row_str)
                elif i >= n - num_cols:
                    # Last few rows
                    row_str = f"{label_strs[i]:>{max_label_width}}  "
                    # First few columns (only if in upper triangle)
                    for j in range(num_cols):
                        if j <= i:
                            row_str += " " * col_width
                        else:
                            val_str = f"{self._matrix[i, j]:.4f}".rstrip('0').rstrip('.')
                            row_str += f"{val_str:>{col_width}}"
                    row_str += "  ...  "
                    # Last few columns
                    for j in range(n - num_cols, n):
                        if j <= i:
                            row_str += " " * col_width
                        else:
                            val_str = f"{self._matrix[i, j]:.4f}".rstrip('0').rstrip('.')
                            row_str += f"{val_str:>{col_width}}"
                    lines.append(row_str)
            
            # Add summary line
            lines.append(f"\n... ({n} elements total, showing {num_cols} rows/columns) ...")
        
        return "\n".join(lines)

