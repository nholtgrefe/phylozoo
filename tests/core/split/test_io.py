"""
Tests for split system I/O functions.

This test suite covers the IOMixin integration for SplitSystem and WeightedSplitSystem,
including save, load, convert, and format conversion methods.
"""

import os
import tempfile

import pytest

from phylozoo.core.split import Split, SplitSystem, WeightedSplitSystem


class TestSplitSystemIO:
    """Test IOMixin methods for SplitSystem."""
    
    def test_to_string(self) -> None:
        """Test converting split system to NEXUS string."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        system = SplitSystem([split1, split2])
        
        nexus_str = system.to_string()
        
        assert '#NEXUS' in nexus_str
        assert 'BEGIN TAXA' in nexus_str
        assert 'BEGIN SPLITS' in nexus_str
        assert 'FORMAT LABELS=YES' in nexus_str
        assert '1' in nexus_str
        assert '2' in nexus_str
        assert '3' in nexus_str
        assert '4' in nexus_str
        assert 'WEIGHTS=YES' not in nexus_str  # Should not have weights
    
    def test_save_basic(self) -> None:
        """Test basic NEXUS file saving."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        system = SplitSystem([split1, split2])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            system.save(file_path, overwrite=True)
            
            assert os.path.exists(file_path)
            
            # Read and verify content
            with open(file_path, 'r') as f:
                content = f.read()
            
            assert '#NEXUS' in content
            assert 'BEGIN TAXA' in content
            assert 'BEGIN SPLITS' in content
            assert '1' in content
            assert '2' in content
    
    def test_save_with_overwrite(self) -> None:
        """Test saving with overwrite=True."""
        split1 = Split({1, 2}, {3, 4})
        system = SplitSystem([split1])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            
            # Save first time
            system.save(file_path, overwrite=True)
            assert os.path.exists(file_path)
            
            # Save again with overwrite
            system.save(file_path, overwrite=True)
            assert os.path.exists(file_path)
    
    def test_save_without_overwrite_raises_error(self) -> None:
        """Test that saving without overwrite raises error if file exists."""
        split1 = Split({1, 2}, {3, 4})
        system = SplitSystem([split1])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            
            # Save first time
            system.save(file_path, overwrite=True)
            
            # Try to save again without overwrite
            with pytest.raises(FileExistsError, match="already exists"):
                system.save(file_path, overwrite=False)
    
    @pytest.mark.parametrize("ext", [".nex", ".nxs"])
    def test_save_nex_extensions(self, ext: str) -> None:
        """Test that .nex and .nxs extensions are accepted."""
        split1 = Split({1, 2}, {3, 4})
        system = SplitSystem([split1])
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, f'test{ext}')
            system.save(file_path, overwrite=True)
            assert os.path.exists(file_path)
    
    def test_from_string(self) -> None:
        """Test parsing NEXUS string to SplitSystem."""
        nexus_str = """#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=4;
    TAXLABELS
        1
        2
        3
        4
    ;
END;

BEGIN SPLITS;
    DIMENSIONS NSPLITS=2;
    FORMAT LABELS=YES;
    MATRIX
        [1] (1 2) (3 4)
        [2] (1 3) (2 4)
    ;
END;"""
        
        system = SplitSystem.from_string(nexus_str, format='nexus')
        
        assert len(system) == 2
        assert system.elements == {1, 2, 3, 4}
    
    def test_load(self) -> None:
        """Test loading SplitSystem from file."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        system = SplitSystem([split1, split2])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            system.save(file_path, overwrite=True)
            
            # Load it back
            loaded_system = SplitSystem.load(file_path)
            
            assert len(loaded_system) == len(system)
            assert loaded_system.elements == system.elements

    def test_load_from_multi_block_nexus(self) -> None:
        """SplitSystem.load from NEXUS with DISTANCES and SPLITS uses only SPLITS."""
        multi_block_nexus = """#NEXUS

BEGIN TAXA;
    DIMENSIONS ntax=4;
    TAXLABELS
        1
        2
        3
        4
    ;
END;

BEGIN DISTANCES;
    DIMENSIONS ntax=4;
    FORMAT triangle=LOWER;
    MATRIX
        1 0
        2 1 0
        3 2 1 0
        4 3 2 1 0
    ;
END;

BEGIN SPLITS;
    DIMENSIONS NSPLITS=2;
    FORMAT LABELS=YES;
    MATRIX
        [1] (1 2) (3 4)
        [2] (1 3) (2 4)
    ;
END;
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'multi.nexus')
            with open(file_path, 'w') as f:
                f.write(multi_block_nexus)
            loaded = SplitSystem.load(file_path)
        assert len(loaded) == 2
        assert loaded.elements == {1, 2, 3, 4}
    
    def test_round_trip(self) -> None:
        """Test round-trip conversion: save and load."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        split3 = Split({1, 4}, {2, 3})
        system = SplitSystem([split1, split2, split3])
        
        nexus_str = system.to_string()
        system2 = SplitSystem.from_string(nexus_str)
        
        assert len(system2) == len(system)
        assert system2.elements == system.elements
    
    def test_empty_split_system(self) -> None:
        """Test I/O with empty split system."""
        system = SplitSystem()
        
        nexus_str = system.to_string()
        assert '#NEXUS' in nexus_str
        assert 'ntax=0' in nexus_str
        assert 'DIMENSIONS NSPLITS=0' in nexus_str
        
        system2 = SplitSystem.from_string(nexus_str)
        assert len(system2) == 0
    
    def test_convert(self) -> None:
        """Test format conversion between files."""
        split1 = Split({1, 2}, {3, 4})
        system = SplitSystem([split1])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, 'input.nexus')
            output_file = os.path.join(tmpdir, 'output.nexus')
            
            system.save(input_file, overwrite=True)
            SplitSystem.convert(input_file, output_file, overwrite=True)
            
            assert os.path.exists(output_file)
            loaded = SplitSystem.load(output_file)
            assert len(loaded) == len(system)
    
    def test_convert_string(self) -> None:
        """Test string format conversion."""
        nexus_str = """#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=4;
    TAXLABELS
        1
        2
        3
        4
    ;
END;

BEGIN SPLITS;
    DIMENSIONS NSPLITS=1;
    FORMAT LABELS=YES;
    MATRIX
        [1] (1 2) (3 4)
    ;
END;"""
        
        # Convert string to string (should be same for same format)
        converted = SplitSystem.convert_string(nexus_str, 'nexus', 'nexus')
        assert '#NEXUS' in converted
        assert 'BEGIN SPLITS' in converted
    
    @pytest.mark.parametrize("method", ["from_string", "to_string"])
    def test_invalid_format_raises_error(self, method: str) -> None:
        """Test that invalid format raises PhyloZooFormatError for from_string and to_string."""
        from phylozoo.utils.exceptions import PhyloZooFormatError
        with pytest.raises(PhyloZooFormatError, match="not supported"):
            if method == "from_string":
                SplitSystem.from_string("invalid", format='invalid_format')
            else:
                SplitSystem([Split({1, 2}, {3, 4})]).to_string(
                    format='invalid_format'
                )
    
    def test_from_string_malformed_nexus_raises_error(self) -> None:
        """Test that malformed NEXUS string raises PhyloZooParseError."""
        from phylozoo.utils.exceptions import PhyloZooParseError
        with pytest.raises(PhyloZooParseError, match="Could not find"):
            SplitSystem.from_string("invalid nexus", format='nexus')
    
    def test_string_labels(self) -> None:
        """Test I/O with string labels."""
        split1 = Split({'A', 'B'}, {'C', 'D'})
        system = SplitSystem([split1])
        
        nexus_str = system.to_string()
        system2 = SplitSystem.from_string(nexus_str)
        
        assert len(system2) == len(system)
        assert system2.elements == system.elements


class TestWeightedSplitSystemIO:
    """Test IOMixin methods for WeightedSplitSystem."""
    
    def test_to_string(self) -> None:
        """Test converting weighted split system to NEXUS string."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 0.8, split2: 0.6}
        system = WeightedSplitSystem(weights)
        
        nexus_str = system.to_string()
        
        assert '#NEXUS' in nexus_str
        assert 'BEGIN TAXA' in nexus_str
        assert 'BEGIN SPLITS' in nexus_str
        assert 'FORMAT LABELS=YES WEIGHTS=YES' in nexus_str
        assert '0.800000' in nexus_str
        assert '0.600000' in nexus_str
    
    def test_save_basic(self) -> None:
        """Test basic NEXUS file saving."""
        split1 = Split({1, 2}, {3, 4})
        weights = {split1: 2.5}
        system = WeightedSplitSystem(weights)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            system.save(file_path, overwrite=True)
            
            assert os.path.exists(file_path)
            
            # Read and verify content
            with open(file_path, 'r') as f:
                content = f.read()
            
            assert '#NEXUS' in content
            assert 'BEGIN TAXA' in content
            assert 'BEGIN SPLITS' in content
            assert 'WEIGHTS=YES' in content
            assert '2.500000' in content
    
    def test_from_string(self) -> None:
        """Test parsing NEXUS string to WeightedSplitSystem."""
        nexus_str = """#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=4;
    TAXLABELS
        1
        2
        3
        4
    ;
END;

BEGIN SPLITS;
    DIMENSIONS NSPLITS=2;
    FORMAT LABELS=YES WEIGHTS=YES;
    MATRIX
        [1] (1 2) (3 4) 0.8
        [2] (1 3) (2 4) 0.6
    ;
END;"""
        
        system = WeightedSplitSystem.from_string(nexus_str, format='nexus')
        
        assert len(system) == 2
        assert system.elements == {1, 2, 3, 4}
        
        # Check weights
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        assert system.get_weight(split1) == 0.8
        assert system.get_weight(split2) == 0.6
    
    def test_from_string_without_weights_uses_default(self) -> None:
        """Test parsing NEXUS without WEIGHTS=YES uses default weight."""
        nexus_str = """#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=4;
    TAXLABELS
        1
        2
        3
        4
    ;
END;

BEGIN SPLITS;
    DIMENSIONS NSPLITS=1;
    FORMAT LABELS=YES;
    MATRIX
        [1] (1 2) (3 4) 1.0
    ;
END;"""
        
        system = WeightedSplitSystem.from_string(nexus_str, format='nexus')
        
        assert len(system) == 1
        split1 = Split({1, 2}, {3, 4})
        assert system.get_weight(split1) == 1.0
    
    def test_load(self) -> None:
        """Test loading WeightedSplitSystem from file."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 0.8, split2: 0.6}
        system = WeightedSplitSystem(weights)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.nexus')
            system.save(file_path, overwrite=True)
            
            # Load it back
            loaded_system = WeightedSplitSystem.load(file_path)
            
            assert len(loaded_system) == len(system)
            assert loaded_system.elements == system.elements
            assert loaded_system.get_weight(split1) == 0.8
            assert loaded_system.get_weight(split2) == 0.6
    
    def test_round_trip(self) -> None:
        """Test round-trip conversion: save and load."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 0.8, split2: 0.6}
        system = WeightedSplitSystem(weights)
        
        nexus_str = system.to_string()
        system2 = WeightedSplitSystem.from_string(nexus_str)
        
        assert len(system2) == len(system)
        assert system2.elements == system.elements
        assert system2.get_weight(split1) == 0.8
        assert system2.get_weight(split2) == 0.6
    
    def test_empty_weighted_split_system(self) -> None:
        """Test I/O with empty weighted split system."""
        system = WeightedSplitSystem()
        
        nexus_str = system.to_string()
        assert '#NEXUS' in nexus_str
        assert 'ntax=0' in nexus_str
        assert 'DIMENSIONS NSPLITS=0' in nexus_str
        
        system2 = WeightedSplitSystem.from_string(nexus_str)
        assert len(system2) == 0
    
    def test_convert(self) -> None:
        """Test format conversion between files."""
        split1 = Split({1, 2}, {3, 4})
        weights = {split1: 2.5}
        system = WeightedSplitSystem(weights)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, 'input.nexus')
            output_file = os.path.join(tmpdir, 'output.nexus')
            
            system.save(input_file, overwrite=True)
            WeightedSplitSystem.convert(input_file, output_file, overwrite=True)
            
            assert os.path.exists(output_file)
            loaded = WeightedSplitSystem.load(output_file)
            assert len(loaded) == len(system)
            assert loaded.get_weight(split1) == 2.5
    
    def test_from_string_missing_weight_raises_error(self) -> None:
        """Test that missing weight when WEIGHTS=YES raises PhyloZooParseError."""
        from phylozoo.utils.exceptions import PhyloZooParseError
        nexus_str = """#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=4;
    TAXLABELS
        1
        2
        3
        4
    ;
END;

BEGIN SPLITS;
    DIMENSIONS NSPLITS=1;
    FORMAT LABELS=YES WEIGHTS=YES;
    MATRIX
        [1] (1 2) (3 4)
    ;
END;"""
        
        with pytest.raises(PhyloZooParseError, match="Weight required but missing"):
            WeightedSplitSystem.from_string(nexus_str, format='nexus')
    
    def test_from_string_invalid_weight_raises_error(self) -> None:
        """Test that invalid weight raises PhyloZooParseError."""
        from phylozoo.utils.exceptions import PhyloZooParseError
        nexus_str = """#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=4;
    TAXLABELS
        1
        2
        3
        4
    ;
END;

BEGIN SPLITS;
    DIMENSIONS NSPLITS=1;
    FORMAT LABELS=YES WEIGHTS=YES;
    MATRIX
        [1] (1 2) (3 4) invalid
    ;
END;"""
        
        with pytest.raises(PhyloZooParseError, match="Could not parse weight"):
            WeightedSplitSystem.from_string(nexus_str, format='nexus')
    
    def test_from_string_negative_weight_raises_error(self) -> None:
        """Test that negative weight raises ValueError."""
        nexus_str = """#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=4;
    TAXLABELS
        1
        2
        3
        4
    ;
END;

BEGIN SPLITS;
    DIMENSIONS NSPLITS=1;
    FORMAT LABELS=YES WEIGHTS=YES;
    MATRIX
        [1] (1 2) (3 4) -0.5
    ;
END;"""
        
        from phylozoo.utils.exceptions import PhyloZooValueError
        with pytest.raises(PhyloZooValueError, match="Weight must be positive"):
            WeightedSplitSystem.from_string(nexus_str, format='nexus')
    
    def test_string_labels(self) -> None:
        """Test I/O with string labels."""
        split1 = Split({'A', 'B'}, {'C', 'D'})
        weights = {split1: 1.5}
        system = WeightedSplitSystem(weights)
        
        nexus_str = system.to_string()
        system2 = WeightedSplitSystem.from_string(nexus_str)
        
        assert len(system2) == len(system)
        assert system2.elements == system.elements
        assert system2.get_weight(split1) == 1.5


class TestFormatConversion:
    """Test format conversion between SplitSystem and WeightedSplitSystem."""
    
    def test_convert_split_system_to_weighted(self) -> None:
        """Test converting SplitSystem file to WeightedSplitSystem."""
        split1 = Split({1, 2}, {3, 4})
        system = SplitSystem([split1])
        
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, 'input.nexus')
            output_file = os.path.join(tmpdir, 'output.nexus')
            
            system.save(input_file, overwrite=True)
            
            # Load as WeightedSplitSystem (should work with default weights)
            loaded = WeightedSplitSystem.load(input_file)
            assert len(loaded) == 1
            # Note: weights might be 1.0 if parser adds default weight
    
    def test_convert_weighted_to_split_system(self) -> None:
        """Test converting WeightedSplitSystem file to SplitSystem."""
        split1 = Split({1, 2}, {3, 4})
        weights = {split1: 0.8}
        system = WeightedSplitSystem(weights)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, 'input.nexus')
            output_file = os.path.join(tmpdir, 'output.nexus')
            
            system.save(input_file, overwrite=True)
            
            # Load as SplitSystem (should ignore weights)
            loaded = SplitSystem.load(input_file)
            assert len(loaded) == 1

