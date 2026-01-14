"""
IO utilities module.

This module provides shared utilities for file I/O operations that can be
reused across different classes and modules. It includes:

1. **Utility Functions**: Low-level file operations (read, write, validation)
2. **FormatRegistry**: Central registry for I/O format handlers
3. **IOMixin**: Mixin class providing unified I/O methods for phylozoo classes

The IOMixin provides a consistent API for saving/loading objects to/from files
with automatic format detection and extensible format support.

Format-Specific Requirements
----------------------------
Format writers can enforce class-specific requirements. For example, a writer
for a tree-only format can check if the object is a tree and raise an error
if not. This allows different formats to have different constraints:

    >>> def to_newick(network):
    ...     # Newick format only supports trees
    ...     if not network.is_tree():
    ...         raise ValueError("Newick format only supports trees, not networks")
    ...     return tree_string
    ... 
    >>> # Register tree-only format
    >>> FormatRegistry.register(
    ...     Network, 'newick',
    ...     reader=from_newick,
    ...     writer=to_newick,  # Will validate tree requirement
    ...     extensions=['.newick', '.nwk']
    ... )

Writers can also accept optional parameters via **kwargs for format-specific
options (e.g., precision, include_labels, etc.).

Examples
--------
Basic usage with IOMixin:

    >>> from phylozoo.utils.io import IOMixin, FormatRegistry
    >>> 
    >>> class MyNetwork(IOMixin):
    ...     _default_format = 'enewick'
    ...     _supported_formats = ['enewick', 'newick']
    ...     
    ...     def __init__(self, data):
    ...         self.data = data
    ...     
    ...     def is_tree(self):
    ...         return False  # Example method
    ... 
    >>> # Register format handlers
    >>> def to_enewick(network, **kwargs):
    ...     return f"Network({network.data});"
    ... 
    >>> def from_enewick(string, **kwargs):
    ...     return MyNetwork("parsed_data")
    ... 
    >>> def to_newick(network, **kwargs):
    ...     # Example: tree-only format with validation
    ...     if not network.is_tree():
    ...         raise ValueError("Newick only supports trees")
    ...     return f"Tree({network.data});"
    ... 
    >>> FormatRegistry.register(
    ...     MyNetwork, 'enewick',
    ...     reader=from_enewick,
    ...     writer=to_enewick,
    ...     extensions=['.enewick', '.eNewick'],
    ...     default=True
    ... )
    >>> 
    >>> FormatRegistry.register(
    ...     MyNetwork, 'newick',
    ...     reader=from_enewick,  # Same reader for simplicity
    ...     writer=to_newick,
    ...     extensions=['.newick', '.nwk']
    ... )
    >>> 
    >>> # Save and load with auto-detection
    >>> net = MyNetwork("test")
    >>> net.save('network.enewick')  # Auto-detects format
    >>> net2 = MyNetwork.load('network.enewick')
    >>> 
    >>> # Format conversion
    >>> MyNetwork.convert('network.enewick', 'network.newick')
    >>> 
    >>> # String conversion
    >>> string = net.to_string(format='enewick')
    >>> net3 = MyNetwork.from_string(string, format='enewick')
    >>> 
    >>> # Format-specific parameters
    >>> net.save('network.enewick', precision=8, include_labels=True)

Format conversion:

    >>> # Convert between formats
    >>> MyNetwork.convert('input.enewick', 'output.newick')
    >>> 
    >>> # Convert string representations
    >>> enewick_str = "(A,B);"
    >>> newick_str = MyNetwork.convert_string(
    ...     enewick_str, 'enewick', 'newick'
    ... )
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, TypeVar

from phylozoo.utils.exceptions import (
    PhyloZooFormatError,
    PhyloZooIOError,
    PhyloZooValueError,
)

T = TypeVar('T')


def read_file_safely(filepath: str | Path, encoding: str = 'utf-8') -> str:
    """
    Read a file safely with error handling.
    
    Parameters
    ----------
    filepath : str | Path
        Path to the file to read.
    encoding : str, optional
        File encoding, by default 'utf-8'.
    
    Returns
    -------
    str
        File contents as string.
    
    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    PhyloZooValueError
        If the path is not a file.
    PhyloZooIOError
        If the file cannot be read.
    
    Examples
    --------
    >>> from phylozoo.utils.io import read_file_safely
    >>> content = read_file_safely("data.txt")
    >>> len(content) > 0
    True
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    if not path.is_file():
        raise PhyloZooValueError(f"Path is not a file: {filepath}")
    
    try:
        return path.read_text(encoding=encoding)
    except Exception as e:
        raise PhyloZooIOError(f"Error reading file {filepath}: {e}") from e


def write_file_safely(filepath: str | Path, content: str, encoding: str = 'utf-8') -> None:
    """
    Write content to a file safely with error handling.
    
    Parameters
    ----------
    filepath : str | Path
        Path to the file to write.
    content : str
        Content to write to the file.
    encoding : str, optional
        File encoding, by default 'utf-8'.
    
    Raises
    ------
    PhyloZooIOError
        If the file cannot be written.
    
    Examples
    --------
    >>> from phylozoo.utils.io import write_file_safely
    >>> write_file_safely("output.txt", "Hello, world!")
    """
    path = Path(filepath)
    try:
        path.write_text(content, encoding=encoding)
    except Exception as e:
        raise PhyloZooIOError(f"Error writing file {filepath}: {e}") from e


def ensure_directory_exists(filepath: str | Path) -> Path:
    """
    Ensure the directory containing the filepath exists, creating it if necessary.
    
    Parameters
    ----------
    filepath : str | Path
        Path to a file (directory will be created if needed).
    
    Returns
    -------
    Path
        Path object to the directory.
    
    Examples
    --------
    >>> from phylozoo.utils.io import ensure_directory_exists
    >>> dir_path = ensure_directory_exists("output/data.txt")
    >>> dir_path.exists()
    True
    """
    path = Path(filepath)
    directory = path.parent
    directory.mkdir(parents=True, exist_ok=True)
    return directory


class FormatRegistry:
    """
    Central registry for I/O format handlers.
    
    This class provides a registry system for registering format readers and writers
    for different object types. Formats can be automatically detected from file
    extensions, and multiple formats can be registered for the same class.
    
    Examples
    --------
    >>> from phylozoo.utils.io import FormatRegistry
    >>> 
    >>> # Register a format
    >>> FormatRegistry.register(
    ...     MyClass, 'json',
    ...     reader=from_json,
    ...     writer=to_json,
    ...     extensions=['.json']
    ... )
    >>> 
    >>> # Get handlers
    >>> writer = FormatRegistry.get_writer(MyClass, 'json')
    >>> reader = FormatRegistry.get_reader(MyClass, 'json')
    >>> 
    >>> # Detect format from file extension
    >>> format = FormatRegistry.detect_format('data.json', MyClass)
    >>> format
    'json'
    """
    
    _readers: dict[tuple[type, str], Callable] = {}
    _writers: dict[tuple[type, str], Callable] = {}
    _default_formats: dict[type, str] = {}
    _extensions: dict[str, str] = {}  # .nex -> nexus, .newick -> newick
    
    @classmethod
    def register(
        cls,
        obj_type: type,
        format: str,
        reader: Callable,
        writer: Callable,
        extensions: list[str] | None = None,
        default: bool = False,
    ) -> None:
        """
        Register format handlers for an object type.
        
        Parameters
        ----------
        obj_type : type
            The class type to register handlers for.
        format : str
            Format name (e.g., 'enewick', 'nexus', 'json').
        reader : Callable
            Function that takes a string and returns an instance of obj_type.
            Signature: `reader(string: str, **kwargs) -> obj_type`
        writer : Callable
            Function that takes an instance of obj_type and returns a string.
            Signature: `writer(obj: obj_type, **kwargs) -> str`
        extensions : list[str] | None, optional
            List of file extensions that map to this format (e.g., ['.enewick', '.eNewick']).
            If None, no automatic detection is available. By default None.
        default : bool, optional
            If True, set this format as the default for obj_type. By default False.
        
        Examples
        --------
        >>> FormatRegistry.register(
        ...     MyNetwork, 'enewick',
        ...     reader=from_enewick,
        ...     writer=to_enewick,
        ...     extensions=['.enewick', '.eNewick'],
        ...     default=True
        ... )
        """
        key = (obj_type, format)
        cls._readers[key] = reader
        cls._writers[key] = writer
        
        if extensions:
            for ext in extensions:
                cls._extensions[ext.lower()] = format
        
        if default:
            cls._default_formats[obj_type] = format
    
    @classmethod
    def get_writer(cls, obj_type: type, format: str) -> Callable:
        """
        Get writer function for a type and format.
        
        Parameters
        ----------
        obj_type : type
            The class type.
        format : str
            Format name.
        
        Returns
        -------
        Callable
            Writer function.
        
        Raises
        ------
        PhyloZooFormatError
            If no writer is registered for the type and format.
        
        Examples
        --------
        >>> writer = FormatRegistry.get_writer(MyNetwork, 'enewick')
        >>> string = writer(my_network_instance)
        """
        key = (obj_type, format)
        if key not in cls._writers:
            raise PhyloZooFormatError(
                f"No writer registered for {obj_type.__name__} format '{format}'. "
                f"Available formats: {cls._get_available_formats(obj_type)}"
            )
        return cls._writers[key]
    
    @classmethod
    def get_reader(cls, obj_type: type, format: str) -> Callable:
        """
        Get reader function for a type and format.
        
        Parameters
        ----------
        obj_type : type
            The class type.
        format : str
            Format name.
        
        Returns
        -------
        Callable
            Reader function.
        
        Raises
        ------
        PhyloZooFormatError
            If no reader is registered for the type and format.
        
        Examples
        --------
        >>> reader = FormatRegistry.get_reader(MyNetwork, 'enewick')
        >>> network = reader(string_content)
        """
        key = (obj_type, format)
        if key not in cls._readers:
            raise PhyloZooFormatError(
                f"No reader registered for {obj_type.__name__} format '{format}'. "
                f"Available formats: {cls._get_available_formats(obj_type)}"
            )
        return cls._readers[key]
    
    @classmethod
    def detect_format(cls, filepath: str | Path, obj_type: type) -> str:
        """
        Detect format from file extension.
        
        Parameters
        ----------
        filepath : str | Path
            Path to the file.
        obj_type : type
            The class type to detect format for.
        
        Returns
        -------
        str
            Detected format name.
        
        Raises
        ------
        PhyloZooFormatError
            If format cannot be detected and no default is set.
        
        Examples
        --------
        >>> format = FormatRegistry.detect_format('network.enewick', MyNetwork)
        >>> format
        'enewick'
        """
        path = Path(filepath)
        ext = path.suffix.lower()
        
        if ext in cls._extensions:
            format = cls._extensions[ext]
            # Verify this format is registered for this type
            if (obj_type, format) in cls._readers:
                return format
        
        # Fall back to default format
        if obj_type in cls._default_formats:
            return cls._default_formats[obj_type]
        
        # If no default, try to find any registered format
        available = cls._get_available_formats(obj_type)
        if available:
            raise PhyloZooFormatError(
                f"Could not detect format from extension '{ext}' for {obj_type.__name__}. "
                f"Available formats: {available}. "
                f"Please specify format explicitly."
            )
        else:
            raise PhyloZooFormatError(
                f"No formats registered for {obj_type.__name__}. "
                f"Please register a format handler first."
            )
    
    @classmethod
    def _get_available_formats(cls, obj_type: type) -> list[str]:
        """Get list of available formats for a type."""
        formats = []
        for (registered_type, format) in cls._readers.keys():
            if registered_type == obj_type:
                formats.append(format)
        return sorted(set(formats))


class IOMixin:
    """
    Mixin providing unified I/O methods for phylozoo classes.
    
    Classes that inherit from this mixin get consistent methods for saving/loading
    to/from files and converting to/from strings. Format handlers must be registered
    using FormatRegistry before use.
    
    Classes should define:
    - `_default_format`: str - Default format name
    - `_supported_formats`: list[str] - List of supported format names
    
    Examples
    --------
    >>> from phylozoo.utils.io import IOMixin, FormatRegistry
    >>> 
    >>> class MyNetwork(IOMixin):
    ...     _default_format = 'enewick'
    ...     _supported_formats = ['enewick', 'newick']
    ...     
    ...     def __init__(self, data):
    ...         self.data = data
    ... 
    >>> # Register format (typically done in module __init__ or io.py)
    >>> FormatRegistry.register(
    ...     MyNetwork, 'enewick',
    ...     reader=lambda s: MyNetwork(s.strip()),
    ...     writer=lambda n: str(n.data),
    ...     extensions=['.enewick'],
    ...     default=True
    ... )
    >>> 
    >>> # Save to file (auto-detects format from extension)
    >>> net = MyNetwork("test_data")
    >>> net.save('network.enewick')
    >>> 
    >>> # Load from file
    >>> net2 = MyNetwork.load('network.enewick')
    >>> 
    >>> # Convert to/from string
    >>> string = net.to_string()
    >>> net3 = MyNetwork.from_string(string)
    """
    
    _default_format: str = 'default'
    _supported_formats: list[str] = []
    
    def to_string(self, format: str | None = None, **kwargs) -> str:
        """
        Convert instance to string representation.
        
        Parameters
        ----------
        format : str | None, optional
            Format name. If None, uses default format.
        **kwargs
            Additional arguments passed to format writer.
        
        Returns
        -------
        str
            String representation in specified format.
        
        Raises
        ------
        PhyloZooFormatError
            If format is not registered or not supported.
        
        Examples
        --------
        >>> net = MyNetwork("data")
        >>> string = net.to_string(format='enewick')
        >>> string
        'data'
        """
        format = format or self._default_format
        
        if format not in self._supported_formats:
            raise PhyloZooFormatError(
                f"Format '{format}' not supported. "
                f"Supported formats: {self._supported_formats}"
            )
        
        writer = FormatRegistry.get_writer(type(self), format)
        return writer(self, **kwargs)
    
    def save(
        self,
        filepath: str | Path,
        format: str | None = None,
        overwrite: bool = False,
        **kwargs
    ) -> None:
        """
        Save instance to file.
        
        Parameters
        ----------
        filepath : str | Path
            Path to output file.
        format : str | None, optional
            Format name. If None, detected from file extension.
        overwrite : bool, optional
            If True, overwrite existing file. If False, raise error if file exists.
            By default False.
        **kwargs
            Additional arguments passed to format writer.
        
        Raises
        ------
        FileExistsError
            If file exists and overwrite is False.
        PhyloZooFormatError
            If format cannot be detected or is not supported.
        PhyloZooIOError
            If file cannot be written.
        
        Examples
        --------
        >>> net = MyNetwork("data")
        >>> net.save('network.enewick')  # Auto-detects format
        >>> net.save('network.txt', format='enewick')  # Explicit format
        >>> net.save('network.enewick', overwrite=True)  # Overwrite existing
        """
        # Detect format if not specified
        if format is None:
            format = FormatRegistry.detect_format(filepath, type(self))
        
        # Check overwrite
        path = Path(filepath)
        if not overwrite and path.exists():
            raise FileExistsError(
                f"File {filepath} already exists. Use overwrite=True to overwrite."
            )
        
        # Ensure directory exists
        ensure_directory_exists(filepath)
        
        # Convert to string and write
        content = self.to_string(format=format, **kwargs)
        write_file_safely(filepath, content)
    
    @classmethod
    def load(cls, filepath: str | Path, format: str | None = None, **kwargs):
        """
        Load instance from file.
        
        Parameters
        ----------
        filepath : str | Path
            Path to input file.
        format : str | None, optional
            Format name. If None, detected from file extension.
        **kwargs
            Additional arguments passed to format reader.
        
        Returns
        -------
        Instance of cls
            Loaded instance.
        
        Raises
        ------
        FileNotFoundError
            If file does not exist.
        PhyloZooFormatError
            If format cannot be detected or is not supported.
        PhyloZooIOError
            If file cannot be read.
        
        Examples
        --------
        >>> net = MyNetwork.load('network.enewick')  # Auto-detects format
        >>> net = MyNetwork.load('network.txt', format='enewick')  # Explicit format
        """
        # Detect format if not specified
        if format is None:
            format = FormatRegistry.detect_format(filepath, cls)
        
        # Read file
        content = read_file_safely(filepath)
        
        # Use format reader
        return cls.from_string(content, format=format, **kwargs)
    
    @classmethod
    def from_string(cls, string: str, format: str | None = None, **kwargs):
        """
        Create instance from string.
        
        Parameters
        ----------
        string : str
            String representation.
        format : str | None, optional
            Format name. If None, uses default format.
        **kwargs
            Additional arguments passed to format reader.
        
        Returns
        -------
        Instance of cls
            Parsed instance.
        
        Raises
        ------
        PhyloZooFormatError
            If format is not registered or not supported.
        
        Examples
        --------
        >>> string = "Network(data);"
        >>> net = MyNetwork.from_string(string, format='enewick')
        >>> net.data
        'data'
        """
        format = format or cls._default_format
        
        if format not in cls._supported_formats:
            raise PhyloZooFormatError(
                f"Format '{format}' not supported. "
                f"Supported formats: {cls._supported_formats}"
            )
        
        reader = FormatRegistry.get_reader(cls, format)
        return reader(string, **kwargs)
    
    @classmethod
    def convert(
        cls,
        input_file: str | Path,
        output_file: str | Path,
        input_format: str | None = None,
        output_format: str | None = None,
        overwrite: bool = False,
        **kwargs
    ) -> None:
        """
        Convert a file from one format to another.
        
        This is a convenience method that loads from one format and saves
        to another. Useful for format migration or interoperability.
        
        Parameters
        ----------
        input_file : str | Path
            Path to input file.
        output_file : str | Path
            Path to output file.
        input_format : str | None, optional
            Input format. If None, auto-detected from input_file extension.
        output_format : str | None, optional
            Output format. If None, auto-detected from output_file extension.
        overwrite : bool, optional
            If True, overwrite existing output file. Default False.
        **kwargs
            Additional arguments passed to both reader and writer.
            Can include format-specific options.
        
        Returns
        -------
        None
        
        Raises
        ------
        FileNotFoundError
            If input file does not exist.
        FileExistsError
            If output file exists and overwrite is False.
        PhyloZooFormatError
            If format cannot be detected or object doesn't meet format requirements
            (e.g., trying to save a network to a tree-only format).
        PhyloZooIOError
            If file cannot be read or written.
        
        Examples
        --------
        >>> # Convert eNewick to Newick (auto-detects formats)
        >>> Network.convert('network.enewick', 'network.newick')
        >>> 
        >>> # Convert with explicit formats
        >>> Network.convert('data.txt', 'data.nex', 
        ...                 input_format='enewick', 
        ...                 output_format='nexus')
        >>> 
        >>> # Convert with format-specific options
        >>> Network.convert('network.enewick', 'network.newick',
        ...                 include_internal_labels=True)
        >>> 
        >>> # Will raise error if network is not a tree and format requires trees
        >>> Network.convert('network.enewick', 'tree.newick')
        Traceback (most recent call last):
            ...
        ValueError: Newick format only supports trees, not networks
        """
        # Load from input format
        obj = cls.load(input_file, format=input_format, **kwargs)
        
        # Save to output format
        obj.save(output_file, format=output_format, overwrite=overwrite, **kwargs)
    
    @classmethod
    def convert_string(
        cls,
        content: str,
        input_format: str,
        output_format: str,
        **kwargs
    ) -> str:
        """
        Convert a string representation from one format to another.
        
        Useful for format transformations without file I/O.
        
        Parameters
        ----------
        content : str
            Input string in input_format.
        input_format : str
            Format of input string.
        output_format : str
            Desired output format.
        **kwargs
            Additional arguments passed to reader/writer.
            Can include format-specific options.
        
        Returns
        -------
        str
            Converted string in output_format.
        
        Raises
        ------
        ValueError
            If format is not registered or object doesn't meet format requirements
            (e.g., trying to convert a network to a tree-only format).
        
        Examples
        --------
        >>> # Convert eNewick string to Newick string
        >>> enewick_str = "(A,B);"
        >>> newick_str = Network.convert_string(
        ...     enewick_str, 'enewick', 'newick'
        ... )
        >>> 
        >>> # Convert with format-specific options
        >>> newick_str = Network.convert_string(
        ...     enewick_str, 'enewick', 'newick',
        ...     include_internal_labels=True
        ... )
        """
        # Parse from input format
        obj = cls.from_string(content, format=input_format, **kwargs)
        
        # Convert to output format
        return obj.to_string(format=output_format, **kwargs)

