"""
Tests for distance matrix I/O functions.
"""

import os
import tempfile
import pytest
import numpy as np

from phylozoo.core.distance import DistanceMatrix
from phylozoo.core.distance.io import save_nexus


class TestSaveNexus:
    """Test save_nexus function."""
    
    def test_save_nexus_basic(self) -> None:
        """Test basic NEXUS file saving."""
        matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            save_nexus(dm, file_path)
            
            assert os.path.exists(file_path)
            
            # Read and verify content
            with open(file_path, 'r') as f:
                content = f.read()
            
            assert '#NEXUS' in content
            assert 'BEGIN Taxa' in content
            assert 'BEGIN Distances' in content
            assert 'A' in content
            assert 'B' in content
            assert 'C' in content
    
    def test_save_nexus_with_overwrite(self) -> None:
        """Test saving with overwrite=True."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            
            # Save first time
            save_nexus(dm, file_path)
            assert os.path.exists(file_path)
            
            # Save again with overwrite
            save_nexus(dm, file_path, overwrite=True)
            assert os.path.exists(file_path)
    
    def test_save_nexus_without_overwrite_raises_error(self) -> None:
        """Test that saving without overwrite raises error if file exists."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            
            # Save first time
            save_nexus(dm, file_path)
            
            # Try to save again without overwrite
            with pytest.raises(FileExistsError, match="already exists"):
                save_nexus(dm, file_path, overwrite=False)
    
    def test_save_nexus_invalid_extension(self) -> None:
        """Test that invalid file extension raises ValueError."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.txt')
            
            with pytest.raises(ValueError, match="must end with"):
                save_nexus(dm, file_path)
    
    def test_save_nexus_nex_extension(self) -> None:
        """Test that .nex extension is accepted."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nex')
            save_nexus(dm, file_path)
            assert os.path.exists(file_path)
    
    def test_save_nexus_format(self) -> None:
        """Test that NEXUS file has correct format."""
        matrix = np.array([[0, 1.5], [1.5, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            save_nexus(dm, file_path)
            
            with open(file_path, 'r') as f:
                content = f.read()
            
            # Check format labels
            assert 'DIMENSIONS ntax=2' in content
            assert 'FORMAT triangle=LOWER diagonal LABELS' in content
            assert 'MATRIX' in content
    
    def test_save_nexus_lower_triangle(self) -> None:
        """Test that only lower triangle is saved."""
        matrix = np.array([
            [0, 1, 2],
            [1, 0, 1.5],
            [2, 1.5, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            save_nexus(dm, file_path)
            
            with open(file_path, 'r') as f:
                lines = f.readlines()
            
            # Find matrix section
            matrix_start = next(i for i, line in enumerate(lines) if 'MATRIX' in line)
            # Get lines until we hit ';' or 'END'
            matrix_lines = []
            for line in lines[matrix_start + 1:]:
                if line.strip() == ';' or 'END' in line:
                    break
                if line.strip():  # Non-empty line
                    matrix_lines.append(line)
            
            # Check that we have 3 rows (A, B, C)
            assert len(matrix_lines) == 3
    
    def test_save_nexus_creates_directory(self) -> None:
        """Test that save_nexus creates directory if needed."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            subdir = os.path.join(tmpdir, 'subdir')
            file_path = os.path.join(subdir, 'test.nexus')
            
            # Directory doesn't exist yet
            assert not os.path.exists(subdir)
            
            # Save should create directory
            save_nexus(dm, file_path)
            assert os.path.exists(file_path)
            assert os.path.exists(subdir)
    
    def test_save_nexus_precision(self) -> None:
        """Test that distances are saved with appropriate precision."""
        matrix = np.array([[0, 1.123456789], [1.123456789, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            save_nexus(dm, file_path)
            
            with open(file_path, 'r') as f:
                content = f.read()
            
            # Should have 6 decimal places
            assert '1.123457' in content or '1.123456' in content

