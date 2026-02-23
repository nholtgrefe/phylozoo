"""
Tests for phylozoo.utils.io module.

This test suite covers:
- Utility functions (read_file_safely, write_file_safely, ensure_directory_exists)
- FormatRegistry class
- IOMixin class and all its methods
"""

import os
import tempfile
from pathlib import Path

import pytest

from phylozoo.utils.io import (
    FormatRegistry,
    IOMixin,
    ensure_directory_exists,
    read_file_safely,
    write_file_safely,
)
from phylozoo.utils.exceptions import (
    PhyloZooFormatError,
    PhyloZooIOError,
    PhyloZooValueError,
)


class TestUtilityFunctions:
    """Test cases for utility functions."""

    def test_read_file_safely(self) -> None:
        """Test reading a file safely."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write("test content")
            temp_path = f.name
        
        try:
            content = read_file_safely(temp_path)
            assert content == "test content"
        finally:
            os.unlink(temp_path)

    def test_read_file_safely_file_not_found(self) -> None:
        """Test that read_file_safely raises FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError, match="File not found"):
            read_file_safely("nonexistent_file.txt")

    def test_read_file_safely_not_a_file(self) -> None:
        """Test that read_file_safely raises PhyloZooValueError for directories."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with pytest.raises(PhyloZooValueError, match="Path is not a file"):
                read_file_safely(tmpdir)

    def test_write_file_safely(self) -> None:
        """Test writing to a file safely."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            temp_path = f.name
        
        try:
            write_file_safely(temp_path, "test content")
            assert read_file_safely(temp_path) == "test content"
        finally:
            os.unlink(temp_path)

    def test_write_file_safely_creates_file(self) -> None:
        """Test that write_file_safely creates file if it doesn't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "new_file.txt"
            write_file_safely(filepath, "new content")
            assert filepath.exists()
            assert read_file_safely(filepath) == "new content"

    def test_ensure_directory_exists(self) -> None:
        """Test that ensure_directory_exists creates directories."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "subdir" / "file.txt"
            directory = ensure_directory_exists(filepath)
            assert directory.exists()
            assert directory == filepath.parent

    def test_ensure_directory_exists_nested(self) -> None:
        """Test that ensure_directory_exists creates nested directories."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "level1" / "level2" / "level3" / "file.txt"
            directory = ensure_directory_exists(filepath)
            assert directory.exists()
            assert (directory / "level3").exists() is False  # filepath.parent is level3
            assert filepath.parent.exists()


class TestFormatRegistry:
    """Test cases for FormatRegistry class."""

    def test_register_format(self) -> None:
        """Test registering a format."""
        class TestClass:
            pass
        
        def reader(string: str, **kwargs):
            return TestClass()
        
        def writer(obj: TestClass, **kwargs):
            return "test"
        
        FormatRegistry.register(
            TestClass, 'test_format',
            reader=reader,
            writer=writer,
            extensions=['.test'],
            default=True
        )
        
        # Verify registration
        assert FormatRegistry.get_reader(TestClass, 'test_format') == reader
        assert FormatRegistry.get_writer(TestClass, 'test_format') == writer

    def test_register_multiple_formats(self) -> None:
        """Test registering multiple formats for the same class."""
        class TestClass:
            pass
        
        def reader1(string: str, **kwargs):
            return TestClass()
        
        def writer1(obj: TestClass, **kwargs):
            return "format1"
        
        def reader2(string: str, **kwargs):
            return TestClass()
        
        def writer2(obj: TestClass, **kwargs):
            return "format2"
        
        FormatRegistry.register(
            TestClass, 'format1',
            reader=reader1,
            writer=writer1,
            extensions=['.fmt1']
        )
        
        FormatRegistry.register(
            TestClass, 'format2',
            reader=reader2,
            writer=writer2,
            extensions=['.fmt2']
        )
        
        assert FormatRegistry.get_writer(TestClass, 'format1') == writer1
        assert FormatRegistry.get_writer(TestClass, 'format2') == writer2

    def test_get_writer_not_registered(self) -> None:
        """Test that get_writer raises PhyloZooFormatError for unregistered format."""
        class TestClass:
            pass
        
        with pytest.raises(PhyloZooFormatError, match="No writer registered"):
            FormatRegistry.get_writer(TestClass, 'nonexistent')

    def test_get_reader_not_registered(self) -> None:
        """Test that get_reader raises PhyloZooFormatError for unregistered format."""
        class TestClass:
            pass
        
        with pytest.raises(PhyloZooFormatError, match="No reader registered"):
            FormatRegistry.get_reader(TestClass, 'nonexistent')

    def test_detect_format_from_extension(self) -> None:
        """Test format detection from file extension."""
        class TestClass:
            pass
        
        def reader(string: str, **kwargs):
            return TestClass()
        
        def writer(obj: TestClass, **kwargs):
            return "test"
        
        FormatRegistry.register(
            TestClass, 'test_format',
            reader=reader,
            writer=writer,
            extensions=['.test', '.TEST'],
            default=True
        )
        
        assert FormatRegistry.detect_format('file.test', TestClass) == 'test_format'
        assert FormatRegistry.detect_format('file.TEST', TestClass) == 'test_format'

    def test_detect_format_falls_back_to_default(self) -> None:
        """Test that detect_format falls back to default if extension not found."""
        class TestClass:
            pass
        
        def reader(string: str, **kwargs):
            return TestClass()
        
        def writer(obj: TestClass, **kwargs):
            return "test"
        
        FormatRegistry.register(
            TestClass, 'test_format',
            reader=reader,
            writer=writer,
            extensions=['.test'],
            default=True
        )
        
        # Unknown extension should use default
        assert FormatRegistry.detect_format('file.unknown', TestClass) == 'test_format'

    def test_detect_format_no_default(self) -> None:
        """Test that detect_format raises PhyloZooFormatError if no default and unknown extension."""
        class TestClass:
            pass
        
        def reader(string: str, **kwargs):
            return TestClass()
        
        def writer(obj: TestClass, **kwargs):
            return "test"
        
        FormatRegistry.register(
            TestClass, 'test_format',
            reader=reader,
            writer=writer,
            extensions=['.test']
            # No default=True
        )
        
        with pytest.raises(PhyloZooFormatError, match="Could not detect format"):
            FormatRegistry.detect_format('file.unknown', TestClass)


class TestIOMixin:
    """Test cases for IOMixin class."""

    def setup_method(self) -> None:
        """Set up test fixtures."""
        # Create a test class
        class TestNetwork(IOMixin):
            _default_format = 'format1'
            _supported_formats = ['format1', 'format2']
            
            def __init__(self, data: str):
                self.data = data
            
            def is_tree(self) -> bool:
                return self.data == "tree"
        
        self.TestNetwork = TestNetwork
        
        # Register format handlers
        def to_format1(network: TestNetwork, **kwargs) -> str:
            return f"FORMAT1:{network.data}"
        
        def from_format1(string: str, **kwargs) -> TestNetwork:
            if string.startswith("FORMAT1:"):
                return TestNetwork(string[8:])
            return TestNetwork(string)
        
        def to_format2(network: TestNetwork, **kwargs) -> str:
            # Tree-only format
            if not network.is_tree():
                raise ValueError("Format2 only supports trees")
            return f"FORMAT2:{network.data}"
        
        def from_format2(string: str, **kwargs) -> TestNetwork:
            if string.startswith("FORMAT2:"):
                return TestNetwork(string[8:])
            return TestNetwork(string)
        
        FormatRegistry.register(
            TestNetwork, 'format1',
            reader=from_format1,
            writer=to_format1,
            extensions=['.fmt1'],
            default=True
        )
        
        FormatRegistry.register(
            TestNetwork, 'format2',
            reader=from_format2,
            writer=to_format2,
            extensions=['.fmt2']
        )

    def test_to_string_default_format(self) -> None:
        """Test converting to string with default format."""
        net = self.TestNetwork("test")
        result = net.to_string()
        assert result == "FORMAT1:test"

    def test_to_string_specified_format(self) -> None:
        """Test converting to string with specified format."""
        net = self.TestNetwork("test")
        result = net.to_string(format='format1')
        assert result == "FORMAT1:test"

    def test_to_string_unsupported_format(self) -> None:
        """Test that to_string raises PhyloZooFormatError for unsupported format."""
        net = self.TestNetwork("test")
        with pytest.raises(PhyloZooFormatError, match="not supported"):
            net.to_string(format='nonexistent')

    def test_to_string_with_kwargs(self) -> None:
        """Test that to_string passes kwargs to writer."""
        class TestClass(IOMixin):
            _default_format = 'test'
            _supported_formats = ['test']
            
            def __init__(self, data: str):
                self.data = data
        
        def writer(obj: TestClass, prefix: str = "", **kwargs) -> str:
            return f"{prefix}{obj.data}"
        
        FormatRegistry.register(
            TestClass, 'test',
            reader=lambda s, **kwargs: TestClass(s),
            writer=writer,
            extensions=['.test'],
            default=True
        )
        
        obj = TestClass("data")
        result = obj.to_string(prefix="PREFIX:")
        assert result == "PREFIX:data"

    def test_save_auto_detect_format(self) -> None:
        """Test saving with auto-detected format."""
        net = self.TestNetwork("test")
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fmt1') as f:
            temp_path = f.name
        
        try:
            net.save(temp_path, overwrite=True)
            content = read_file_safely(temp_path)
            assert content == "FORMAT1:test"
        finally:
            os.unlink(temp_path)

    def test_save_explicit_format(self) -> None:
        """Test saving with explicit format."""
        net = self.TestNetwork("test")
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            temp_path = f.name
        
        try:
            net.save(temp_path, format='format1', overwrite=True)
            content = read_file_safely(temp_path)
            assert content == "FORMAT1:test"
        finally:
            os.unlink(temp_path)

    def test_save_file_exists_error(self) -> None:
        """Test that save raises FileExistsError if file exists and overwrite=False."""
        net = self.TestNetwork("test")
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fmt1') as f:
            temp_path = f.name
        
        try:
            net.save(temp_path, overwrite=True)
            with pytest.raises(FileExistsError, match="already exists"):
                net.save(temp_path, overwrite=False)
        finally:
            os.unlink(temp_path)

    def test_save_overwrite(self) -> None:
        """Test that save overwrites existing file when overwrite=True."""
        net1 = self.TestNetwork("test1")
        net2 = self.TestNetwork("test2")
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fmt1') as f:
            temp_path = f.name
        
        try:
            net1.save(temp_path, overwrite=True)
            net2.save(temp_path, overwrite=True)
            content = read_file_safely(temp_path)
            assert content == "FORMAT1:test2"
        finally:
            os.unlink(temp_path)

    def test_save_creates_directory(self) -> None:
        """Test that save creates directory if it doesn't exist."""
        net = self.TestNetwork("test")
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = Path(tmpdir) / "subdir" / "file.fmt1"
            net.save(filepath, overwrite=True)
            assert filepath.exists()
            content = read_file_safely(filepath)
            assert content == "FORMAT1:test"

    def test_save_format_requirement_error(self) -> None:
        """Test that save raises error if format requirements not met."""
        net = self.TestNetwork("network")  # Not a tree
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fmt2') as f:
            temp_path = f.name
        
        try:
            with pytest.raises(ValueError, match="Format2 only supports trees"):
                net.save(temp_path, overwrite=True)
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)

    def test_load_auto_detect_format(self) -> None:
        """Test loading with auto-detected format."""
        net = self.TestNetwork("test")
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fmt1') as f:
            f.write("FORMAT1:loaded")
            temp_path = f.name
        
        try:
            loaded = self.TestNetwork.load(temp_path)
            assert loaded.data == "loaded"
        finally:
            os.unlink(temp_path)

    def test_load_explicit_format(self) -> None:
        """Test loading with explicit format."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write("FORMAT1:loaded")
            temp_path = f.name
        
        try:
            loaded = self.TestNetwork.load(temp_path, format='format1')
            assert loaded.data == "loaded"
        finally:
            os.unlink(temp_path)

    def test_load_file_not_found(self) -> None:
        """Test that load raises FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError, match="File not found"):
            self.TestNetwork.load("nonexistent.fmt1")

    def test_from_string_default_format(self) -> None:
        """Test creating instance from string with default format."""
        string = "FORMAT1:parsed"
        net = self.TestNetwork.from_string(string)
        assert net.data == "parsed"

    def test_from_string_specified_format(self) -> None:
        """Test creating instance from string with specified format."""
        string = "FORMAT1:parsed"
        net = self.TestNetwork.from_string(string, format='format1')
        assert net.data == "parsed"

    def test_from_string_unsupported_format(self) -> None:
        """Test that from_string raises PhyloZooFormatError for unsupported format."""
        with pytest.raises(PhyloZooFormatError, match="not supported"):
            self.TestNetwork.from_string("test", format='nonexistent')

    def test_from_string_with_kwargs(self) -> None:
        """Test that from_string passes kwargs to reader."""
        class TestClass(IOMixin):
            _default_format = 'test'
            _supported_formats = ['test']
            
            def __init__(self, data: str, extra: str = ""):
                self.data = data
                self.extra = extra
        
        def reader(string: str, extra: str = "", **kwargs) -> TestClass:
            return TestClass(string, extra=extra)
        
        FormatRegistry.register(
            TestClass, 'test',
            reader=reader,
            writer=lambda obj, **kwargs: obj.data,
            extensions=['.test'],
            default=True
        )
        
        obj = TestClass.from_string("data", extra="extra_value")
        assert obj.data == "data"
        assert obj.extra == "extra_value"

    def test_convert_file_to_file(self) -> None:
        """Test converting between file formats."""
        net = self.TestNetwork("tree")
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fmt1') as f1:
            temp_path1 = f1.name
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fmt2') as f2:
            temp_path2 = f2.name
        
        try:
            # Save in format1
            net.save(temp_path1, overwrite=True)
            
            # Convert to format2
            self.TestNetwork.convert(temp_path1, temp_path2, overwrite=True)
            
            # Load and verify
            loaded = self.TestNetwork.load(temp_path2)
            assert loaded.data == "tree"
            
            # Verify format2 content
            content = read_file_safely(temp_path2)
            assert content == "FORMAT2:tree"
        finally:
            if os.path.exists(temp_path1):
                os.unlink(temp_path1)
            if os.path.exists(temp_path2):
                os.unlink(temp_path2)

    def test_convert_explicit_formats(self) -> None:
        """Test converting with explicit format specifications."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f1:
            f1.write("FORMAT1:test")
            temp_path1 = f1.name
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f2:
            temp_path2 = f2.name
        
        try:
            self.TestNetwork.convert(
                temp_path1, temp_path2,
                input_format='format1',
                output_format='format1',
                overwrite=True
            )
            
            content = read_file_safely(temp_path2)
            assert content == "FORMAT1:test"
        finally:
            if os.path.exists(temp_path1):
                os.unlink(temp_path1)
            if os.path.exists(temp_path2):
                os.unlink(temp_path2)

    def test_convert_format_requirement_error(self) -> None:
        """Test that convert raises error if format requirements not met."""
        net = self.TestNetwork("network")  # Not a tree
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fmt1') as f1:
            temp_path1 = f1.name
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fmt2') as f2:
            temp_path2 = f2.name
        
        try:
            net.save(temp_path1, overwrite=True)
            
            with pytest.raises(ValueError, match="Format2 only supports trees"):
                self.TestNetwork.convert(temp_path1, temp_path2, overwrite=True)
        finally:
            if os.path.exists(temp_path1):
                os.unlink(temp_path1)
            if os.path.exists(temp_path2):
                os.unlink(temp_path2)

    def test_convert_file_not_found(self) -> None:
        """Test that convert raises FileNotFoundError for missing input file."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fmt2') as f:
            temp_path = f.name
        
        try:
            with pytest.raises(FileNotFoundError):
                self.TestNetwork.convert("nonexistent.fmt1", temp_path, overwrite=True)
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)

    def test_convert_file_exists_error(self) -> None:
        """Test that convert raises FileExistsError if output exists and overwrite=False."""
        net = self.TestNetwork("test")
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fmt1') as f1:
            temp_path1 = f1.name
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fmt1') as f2:
            temp_path2 = f2.name
        
        try:
            net.save(temp_path1, overwrite=True)
            net.save(temp_path2, overwrite=True)
            
            with pytest.raises(FileExistsError, match="already exists"):
                self.TestNetwork.convert(temp_path1, temp_path2, overwrite=False)
        finally:
            if os.path.exists(temp_path1):
                os.unlink(temp_path1)
            if os.path.exists(temp_path2):
                os.unlink(temp_path2)

    def test_convert_string(self) -> None:
        """Test converting string representations between formats."""
        string1 = "FORMAT1:test"
        string2 = self.TestNetwork.convert_string(string1, 'format1', 'format1')
        assert string2 == "FORMAT1:test"

    def test_convert_string_different_formats(self) -> None:
        """Test converting string between different formats."""
        # Convert tree from format1 to format2
        string1 = "FORMAT1:tree"
        string2 = self.TestNetwork.convert_string(string1, 'format1', 'format2')
        assert string2 == "FORMAT2:tree"

    def test_convert_string_format_requirement_error(self) -> None:
        """Test that convert_string raises error if format requirements not met."""
        string1 = "FORMAT1:network"  # Not a tree
        with pytest.raises(ValueError, match="Format2 only supports trees"):
            self.TestNetwork.convert_string(string1, 'format1', 'format2')

    def test_convert_string_with_kwargs(self) -> None:
        """Test that convert_string passes kwargs to reader and writer."""
        class TestClass(IOMixin):
            _default_format = 'test'
            _supported_formats = ['test']
            
            def __init__(self, data: str, prefix: str = ""):
                self.data = data
                self.prefix = prefix
        
        def reader(string: str, prefix: str = "", **kwargs) -> TestClass:
            return TestClass(string, prefix=prefix)
        
        def writer(obj: TestClass, suffix: str = "", **kwargs) -> str:
            return f"{obj.prefix}{obj.data}{suffix}"
        
        FormatRegistry.register(
            TestClass, 'test',
            reader=reader,
            writer=writer,
            extensions=['.test'],
            default=True
        )
        
        result = TestClass.convert_string(
            "data", 'test', 'test',
            prefix="PREFIX:",
            suffix=":SUFFIX"
        )
        assert result == "PREFIX:data:SUFFIX"

    def test_round_trip_save_load(self) -> None:
        """Test that save/load round trip preserves data."""
        net = self.TestNetwork("test_data")
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fmt1') as f:
            temp_path = f.name
        
        try:
            net.save(temp_path, overwrite=True)
            loaded = self.TestNetwork.load(temp_path)
            assert loaded.data == net.data
        finally:
            os.unlink(temp_path)

    def test_round_trip_to_string_from_string(self) -> None:
        """Test that to_string/from_string round trip preserves data."""
        net = self.TestNetwork("test_data")
        string = net.to_string()
        loaded = self.TestNetwork.from_string(string)
        assert loaded.data == net.data

