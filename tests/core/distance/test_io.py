"""
Tests for distance matrix I/O functions.

This test suite covers the IOMixin integration for DistanceMatrix, including
save, load, convert, and format conversion methods.
"""

import os
import tempfile

import numpy as np
import pytest

from phylozoo.core.distance import DistanceMatrix


class TestDistanceMatrixIO:
    """Test IOMixin methods for DistanceMatrix."""
    
    def test_to_string(self) -> None:
        """Test converting distance matrix to NEXUS string."""
        matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        
        nexus_str = dm.to_string()
        
        assert '#NEXUS' in nexus_str
        assert 'BEGIN Taxa' in nexus_str
        assert 'BEGIN Distances' in nexus_str
        assert 'A' in nexus_str
        assert 'B' in nexus_str
        assert 'C' in nexus_str
    
    def test_save_basic(self) -> None:
        """Test basic NEXUS file saving."""
        matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            dm.save(file_path, overwrite=True)
            
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
    
    def test_save_with_overwrite(self) -> None:
        """Test saving with overwrite=True."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            
            # Save first time
            dm.save(file_path, overwrite=True)
            assert os.path.exists(file_path)
            
            # Save again with overwrite
            dm.save(file_path, overwrite=True)
            assert os.path.exists(file_path)
    
    def test_save_without_overwrite_raises_error(self) -> None:
        """Test that saving without overwrite raises error if file exists."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            
            # Save first time
            dm.save(file_path, overwrite=True)
            
            # Try to save again without overwrite
            with pytest.raises(FileExistsError, match="already exists"):
                dm.save(file_path, overwrite=False)
    
    def test_save_nex_extension(self) -> None:
        """Test that .nex extension is accepted."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nex')
            dm.save(file_path, overwrite=True)
            assert os.path.exists(file_path)
    
    def test_save_format(self) -> None:
        """Test that NEXUS file has correct format."""
        matrix = np.array([[0, 1.5], [1.5, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            dm.save(file_path, overwrite=True)
            
            with open(file_path, 'r') as f:
                content = f.read()
            
            # Check format labels
            assert 'DIMENSIONS ntax=2' in content
            assert 'FORMAT triangle=LOWER diagonal LABELS' in content
            assert 'MATRIX' in content
    
    def test_save_lower_triangle(self) -> None:
        """Test that only lower triangle is saved."""
        matrix = np.array([
            [0, 1, 2],
            [1, 0, 1.5],
            [2, 1.5, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            dm.save(file_path, overwrite=True)
            
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
    
    def test_save_creates_directory(self) -> None:
        """Test that save creates directory if needed."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            subdir = os.path.join(tmpdir, 'subdir')
            file_path = os.path.join(subdir, 'test.nexus')
            
            # Directory doesn't exist yet
            assert not os.path.exists(subdir)
            
            # Save should create directory
            dm.save(file_path, overwrite=True)
            assert os.path.exists(file_path)
            assert os.path.exists(subdir)
    
    def test_save_precision(self) -> None:
        """Test that distances are saved with appropriate precision."""
        matrix = np.array([[0, 1.123456789], [1.123456789, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            dm.save(file_path, overwrite=True)
            
            with open(file_path, 'r') as f:
                content = f.read()
            
            # Should have 6 decimal places
            assert '1.123457' in content or '1.123456' in content
    
    def test_load_basic(self) -> None:
        """Test loading distance matrix from NEXUS file."""
        matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            dm.save(file_path, overwrite=True)
            
            # Load and verify
            dm2 = DistanceMatrix.load(file_path)
            assert dm2.labels == dm.labels
            assert np.allclose(dm2._matrix, dm._matrix)
    
    def test_load_auto_detect_format(self) -> None:
        """Test that load auto-detects format from extension."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nex')
            dm.save(file_path, overwrite=True)
            
            dm2 = DistanceMatrix.load(file_path)
            assert dm2.labels == dm.labels
            assert np.allclose(dm2._matrix, dm._matrix)
    
    def test_from_string(self) -> None:
        """Test creating distance matrix from NEXUS string."""
        nexus_str = '''#NEXUS

BEGIN Taxa;
    DIMENSIONS ntax=3;
    TAXLABELS
        A
        B
        C
    ;
END;

BEGIN Distances;
    DIMENSIONS ntax=3;
    FORMAT triangle=LOWER diagonal LABELS;
    MATRIX
    A 0.000000
    B 1.000000 0.000000
    C 2.000000 1.000000 0.000000
    ;
END;'''
        
        dm = DistanceMatrix.from_string(nexus_str)
        
        assert len(dm) == 3
        assert dm.labels == ('A', 'B', 'C')
        assert dm.get_distance('A', 'B') == 1.0
        assert dm.get_distance('A', 'C') == 2.0
        assert dm.get_distance('B', 'C') == 1.0
    
    def test_round_trip_save_load(self) -> None:
        """Test that save/load round trip preserves data."""
        matrix = np.array([[0, 1.5, 2.3], [1.5, 0, 1.1], [2.3, 1.1, 0]])
        dm = DistanceMatrix(matrix, labels=['X', 'Y', 'Z'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            dm.save(file_path, overwrite=True)
            
            dm2 = DistanceMatrix.load(file_path)
            assert dm2.labels == dm.labels
            assert np.allclose(dm2._matrix, dm._matrix)
    
    def test_round_trip_to_string_from_string(self) -> None:
        """Test that to_string/from_string round trip preserves data."""
        matrix = np.array([[0, 1, 2], [1, 0, 1.5], [2, 1.5, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        
        nexus_str = dm.to_string()
        dm2 = DistanceMatrix.from_string(nexus_str)
        
        assert dm2.labels == dm.labels
        assert np.allclose(dm2._matrix, dm._matrix)
    
    def test_convert(self) -> None:
        """Test converting between file formats."""
        matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path1 = os.path.join(tmpdir, 'test1.nexus')
            file_path2 = os.path.join(tmpdir, 'test2.nex')
            
            dm.save(file_path1, overwrite=True)
            DistanceMatrix.convert(file_path1, file_path2, overwrite=True)
            
            dm2 = DistanceMatrix.load(file_path2)
            assert dm2.labels == dm.labels
            assert np.allclose(dm2._matrix, dm._matrix)
    
    def test_convert_string(self) -> None:
        """Test converting string representations."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        nexus_str1 = dm.to_string()
        nexus_str2 = DistanceMatrix.convert_string(nexus_str1, 'nexus', 'nexus')
        
        dm2 = DistanceMatrix.from_string(nexus_str2)
        assert dm2.labels == dm.labels
        assert np.allclose(dm2._matrix, dm._matrix)
    
    def test_from_string_malformed_raises_error(self) -> None:
        """Test that malformed NEXUS string raises ValueError."""
        with pytest.raises(ValueError):
            DistanceMatrix.from_string("not a nexus file")
    
    def test_from_string_missing_taxa_block(self) -> None:
        """Test that missing Taxa block raises ValueError."""
        nexus_str = '''#NEXUS

BEGIN Distances;
    DIMENSIONS ntax=2;
    FORMAT triangle=LOWER diagonal LABELS;
    MATRIX
    A 0.000000
    B 1.000000 0.000000
    ;
END;'''
        
        with pytest.raises(ValueError, match="Could not find Taxa block"):
            DistanceMatrix.from_string(nexus_str)
    
    def test_from_string_missing_distances_block(self) -> None:
        """Test that missing Distances block raises ValueError."""
        nexus_str = '''#NEXUS

BEGIN Taxa;
    DIMENSIONS ntax=2;
    TAXLABELS
        A
        B
    ;
END;'''
        
        with pytest.raises(ValueError, match="Could not find Distances block"):
            DistanceMatrix.from_string(nexus_str)
    
    def test_from_string_mismatched_dimensions(self) -> None:
        """Test that mismatched dimensions raise ValueError."""
        nexus_str = '''#NEXUS

BEGIN Taxa;
    DIMENSIONS ntax=3;
    TAXLABELS
        A
        B
        C
    ;
END;

BEGIN Distances;
    DIMENSIONS ntax=3;
    FORMAT triangle=LOWER diagonal LABELS;
    MATRIX
    A 0.000000
    B 1.000000 0.000000
    ;
END;'''
        
        with pytest.raises(ValueError, match="Number of matrix rows"):
            DistanceMatrix.from_string(nexus_str)


class TestPhylipFormat:
    """Test cases for PHYLIP format."""
    
    def test_to_phylip(self) -> None:
        """Test converting distance matrix to PHYLIP string."""
        matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        
        phylip_str = dm.to_string(format='phylip')
        
        assert phylip_str.startswith('3\n')
        assert 'A         0.00000 1.00000 2.00000' in phylip_str
        assert 'B         1.00000 0.00000 1.00000' in phylip_str
        assert 'C         2.00000 1.00000 0.00000' in phylip_str
    
    def test_from_phylip(self) -> None:
        """Test parsing PHYLIP string."""
        phylip_str = '''3
A         0.00000 1.00000 2.00000
B         1.00000 0.00000 1.00000
C         2.00000 1.00000 0.00000
'''
        
        dm = DistanceMatrix.from_string(phylip_str, format='phylip')
        
        assert len(dm) == 3
        assert dm.labels == ('A', 'B', 'C')
        assert dm.get_distance('A', 'B') == 1.0
        assert dm.get_distance('A', 'C') == 2.0
        assert dm.get_distance('B', 'C') == 1.0
    
    def test_round_trip_phylip(self) -> None:
        """Test that PHYLIP round trip preserves data."""
        matrix = np.array([[0, 1.5, 2.3], [1.5, 0, 1.1], [2.3, 1.1, 0]])
        dm = DistanceMatrix(matrix, labels=['X', 'Y', 'Z'])
        
        phylip_str = dm.to_string(format='phylip')
        dm2 = DistanceMatrix.from_string(phylip_str, format='phylip')
        
        assert dm2.labels == dm.labels
        assert np.allclose(dm2._matrix, dm._matrix)
    
    def test_save_load_phylip(self) -> None:
        """Test saving and loading PHYLIP files."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.phy')
            dm.save(file_path, format='phylip', overwrite=True)
            
            dm2 = DistanceMatrix.load(file_path)
            assert dm2.labels == dm.labels
            assert np.allclose(dm2._matrix, dm._matrix)
    
    def test_from_phylip_malformed(self) -> None:
        """Test that malformed PHYLIP string raises ValueError."""
        with pytest.raises(ValueError):
            DistanceMatrix.from_string("not phylip", format='phylip')
    
    def test_from_phylip_insufficient_rows(self) -> None:
        """Test that insufficient rows raise ValueError."""
        phylip_str = '''3
A         0.00000 1.00000 2.00000
B         1.00000 0.00000 1.00000
'''
        
        with pytest.raises(ValueError, match="Expected.*lines"):
            DistanceMatrix.from_string(phylip_str, format='phylip')
    
    def test_from_phylip_non_symmetric(self) -> None:
        """Test that non-symmetric matrix raises ValueError."""
        phylip_str = '''2
A         0.00000 1.00000
B         2.00000 0.00000
'''
        
        with pytest.raises(ValueError, match="not symmetric"):
            DistanceMatrix.from_string(phylip_str, format='phylip')


class TestCSVFormat:
    """Test cases for CSV format."""
    
    def test_to_csv(self) -> None:
        """Test converting distance matrix to CSV string."""
        matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        
        csv_str = dm.to_string(format='csv')
        
        assert csv_str.startswith(',A,B,C\n')
        assert 'A,0.000000,1.000000,2.000000' in csv_str
        assert 'B,1.000000,0.000000,1.000000' in csv_str
        assert 'C,2.000000,1.000000,0.000000' in csv_str
    
    def test_to_csv_no_header(self) -> None:
        """Test CSV output without header."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        csv_str = dm.to_string(format='csv', include_header=False)
        
        assert not csv_str.startswith(',')
        assert 'A,0.000000,1.000000' in csv_str
    
    def test_to_csv_tab_delimiter(self) -> None:
        """Test CSV output with tab delimiter."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        csv_str = dm.to_string(format='csv', delimiter='\t')
        
        assert '\t' in csv_str
        assert ',' not in csv_str or csv_str.count(',') < csv_str.count('\t')
    
    def test_from_csv(self) -> None:
        """Test parsing CSV string."""
        csv_str = ''',A,B,C
A,0.0,1.0,2.0
B,1.0,0.0,1.0
C,2.0,1.0,0.0
'''
        
        dm = DistanceMatrix.from_string(csv_str, format='csv')
        
        assert len(dm) == 3
        assert dm.labels == ('A', 'B', 'C')
        assert dm.get_distance('A', 'B') == 1.0
        assert dm.get_distance('A', 'C') == 2.0
        assert dm.get_distance('B', 'C') == 1.0
    
    def test_from_csv_no_header(self) -> None:
        """Test parsing CSV string without header."""
        csv_str = '''A,0.0,1.0
B,1.0,0.0
'''
        
        dm = DistanceMatrix.from_string(csv_str, format='csv', has_header=False)
        
        assert len(dm) == 2
        assert dm.labels == ('A', 'B')
        assert dm.get_distance('A', 'B') == 1.0
    
    def test_from_csv_tab_delimiter(self) -> None:
        """Test parsing CSV string with tab delimiter."""
        csv_str = '''	A	B
A	0.0	1.0
B	1.0	0.0
'''
        
        dm = DistanceMatrix.from_string(csv_str, format='csv', delimiter='\t')
        
        assert len(dm) == 2
        assert dm.labels == ('A', 'B')
        assert dm.get_distance('A', 'B') == 1.0
    
    def test_round_trip_csv(self) -> None:
        """Test that CSV round trip preserves data."""
        matrix = np.array([[0, 1.5, 2.3], [1.5, 0, 1.1], [2.3, 1.1, 0]])
        dm = DistanceMatrix(matrix, labels=['X', 'Y', 'Z'])
        
        csv_str = dm.to_string(format='csv')
        dm2 = DistanceMatrix.from_string(csv_str, format='csv')
        
        assert dm2.labels == dm.labels
        assert np.allclose(dm2._matrix, dm._matrix)
    
    def test_save_load_csv(self) -> None:
        """Test saving and loading CSV files."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.csv')
            dm.save(file_path, format='csv', overwrite=True)
            
            dm2 = DistanceMatrix.load(file_path)
            assert dm2.labels == dm.labels
            assert np.allclose(dm2._matrix, dm._matrix)
    
    def test_from_csv_malformed(self) -> None:
        """Test that malformed CSV string raises ValueError."""
        with pytest.raises(ValueError):
            DistanceMatrix.from_string("not csv", format='csv')
    
    def test_from_csv_non_symmetric(self) -> None:
        """Test that non-symmetric matrix raises ValueError."""
        csv_str = ''',A,B
A,0.0,1.0
B,2.0,0.0
'''
        
        with pytest.raises(ValueError, match="not symmetric"):
            DistanceMatrix.from_string(csv_str, format='csv')
    
    def test_from_csv_mismatched_dimensions(self) -> None:
        """Test that mismatched dimensions raise ValueError."""
        csv_str = ''',A,B,C
A,0.0,1.0
B,1.0,0.0
'''
        
        with pytest.raises(ValueError, match="has.*distances|Header has"):
            DistanceMatrix.from_string(csv_str, format='csv')


class TestFormatConversion:
    """Test cases for format conversion."""
    
    def test_convert_nexus_to_phylip(self) -> None:
        """Test converting from NEXUS to PHYLIP."""
        matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            nexus_file = os.path.join(tmpdir, 'test.nexus')
            phylip_file = os.path.join(tmpdir, 'test.phy')
            
            dm.save(nexus_file, overwrite=True)
            DistanceMatrix.convert(nexus_file, phylip_file, overwrite=True)
            
            dm2 = DistanceMatrix.load(phylip_file)
            assert dm2.labels == dm.labels
            assert np.allclose(dm2._matrix, dm._matrix)
    
    def test_convert_phylip_to_csv(self) -> None:
        """Test converting from PHYLIP to CSV."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            phylip_file = os.path.join(tmpdir, 'test.phy')
            csv_file = os.path.join(tmpdir, 'test.csv')
            
            dm.save(phylip_file, format='phylip', overwrite=True)
            DistanceMatrix.convert(phylip_file, csv_file, overwrite=True)
            
            dm2 = DistanceMatrix.load(csv_file)
            assert dm2.labels == dm.labels
            assert np.allclose(dm2._matrix, dm._matrix)
    
    def test_convert_string_nexus_to_phylip(self) -> None:
        """Test converting string from NEXUS to PHYLIP."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        nexus_str = dm.to_string(format='nexus')
        phylip_str = DistanceMatrix.convert_string(nexus_str, 'nexus', 'phylip')
        
        dm2 = DistanceMatrix.from_string(phylip_str, format='phylip')
        assert dm2.labels == dm.labels
        assert np.allclose(dm2._matrix, dm._matrix)
    
    def test_convert_string_phylip_to_csv(self) -> None:
        """Test converting string from PHYLIP to CSV."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        
        phylip_str = dm.to_string(format='phylip')
        csv_str = DistanceMatrix.convert_string(phylip_str, 'phylip', 'csv')
        
        dm2 = DistanceMatrix.from_string(csv_str, format='csv')
        assert dm2.labels == dm.labels
        assert np.allclose(dm2._matrix, dm._matrix)

