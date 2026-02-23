"""
Central registry for I/O format handlers.

Provides FormatRegistry for registering format readers and writers
for different object types. Formats can be detected from file extensions.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

from phylozoo.utils.exceptions import PhyloZooFormatError


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

    _readers: dict[tuple[type, str], Callable[..., Any]] = {}
    _writers: dict[tuple[type, str], Callable[..., Any]] = {}
    _default_formats: dict[type, str] = {}
    _extensions: dict[str, str] = {}

    @classmethod
    def register(
        cls,
        obj_type: type,
        format: str,
        reader: Callable[..., Any],
        writer: Callable[..., Any],
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
            List of file extensions that map to this format.
            If None, no automatic detection is available. By default None.
        default : bool, optional
            If True, set this format as the default for obj_type. By default False.
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
    def get_writer(cls, obj_type: type, format: str) -> Callable[..., Any]:
        """Get writer function for a type and format."""
        key = (obj_type, format)
        if key not in cls._writers:
            raise PhyloZooFormatError(
                f"No writer registered for {obj_type.__name__} format '{format}'. "
                f"Available formats: {cls._get_available_formats(obj_type)}"
            )
        return cls._writers[key]

    @classmethod
    def get_reader(cls, obj_type: type, format: str) -> Callable[..., Any]:
        """Get reader function for a type and format."""
        key = (obj_type, format)
        if key not in cls._readers:
            raise PhyloZooFormatError(
                f"No reader registered for {obj_type.__name__} format '{format}'. "
                f"Available formats: {cls._get_available_formats(obj_type)}"
            )
        return cls._readers[key]

    @classmethod
    def detect_format(cls, filepath: str | Path, obj_type: type) -> str:
        """Detect format from file extension."""
        path = Path(filepath)
        ext = path.suffix.lower()

        if ext in cls._extensions:
            format = cls._extensions[ext]
            if (obj_type, format) in cls._readers:
                return format

        if obj_type in cls._default_formats:
            return cls._default_formats[obj_type]

        available = cls._get_available_formats(obj_type)
        if available:
            raise PhyloZooFormatError(
                f"Could not detect format from extension '{ext}' for {obj_type.__name__}. "
                f"Available formats: {available}. "
                f"Please specify format explicitly."
            )
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
