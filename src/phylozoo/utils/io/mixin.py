"""
Mixin providing unified I/O methods for phylozoo classes.

Classes that inherit from IOMixin get save/load, to_string/from_string,
and convert/convert_string. Format handlers must be registered via FormatRegistry.
"""

from __future__ import annotations

from pathlib import Path

from phylozoo.utils.exceptions import PhyloZooFormatError

from .file_ops import ensure_directory_exists, read_file_safely, write_file_safely
from .registry import FormatRegistry


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
    ...     _supported_formats = ['enewick']
    ...     def __init__(self, data):
    ...         self.data = data
    >>>
    >>> FormatRegistry.register(
    ...     MyNetwork, 'enewick',
    ...     reader=lambda s: MyNetwork(s.strip()),
    ...     writer=lambda n: str(n.data),
    ...     extensions=['.enewick'],
    ...     default=True
    ... )
    >>> net = MyNetwork("test_data")
    >>> net.save('network.enewick')
    >>> net2 = MyNetwork.load('network.enewick')
    """

    _default_format: str = 'default'
    _supported_formats: list[str] = []

    def to_string(self, format: str | None = None, **kwargs: object) -> str:
        """
        Convert instance to string representation in the given format.

        Parameters
        ----------
        format : str or None, optional
            Output format name. If None, uses the class's default format.
        **kwargs
            Additional format-specific options passed to the writer.

        Returns
        -------
        str
            String representation of the instance in the requested format.

        Raises
        ------
        PhyloZooFormatError
            If the requested format is not in the class's supported formats.
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
        **kwargs: object,
    ) -> None:
        """
        Save instance to a file.

        Parameters
        ----------
        filepath : str or Path
            Path to the output file.
        format : str or None, optional
            Output format name. If None, inferred from the file extension.
        overwrite : bool, optional
            If True, overwrite existing file. Default is False.
        **kwargs
            Additional format-specific options passed to the writer.

        Raises
        ------
        FileExistsError
            If the file already exists and overwrite is False.
        PhyloZooFormatError
            If the format is not supported.
        """
        if format is None:
            format = FormatRegistry.detect_format(filepath, type(self))
        path = Path(filepath)
        if not overwrite and path.exists():
            raise FileExistsError(
                f"File {filepath} already exists. Use overwrite=True to overwrite."
            )
        ensure_directory_exists(filepath)
        content = self.to_string(format=format, **kwargs)
        write_file_safely(filepath, content)

    @classmethod
    def load(cls, filepath: str | Path, format: str | None = None, **kwargs: object):  # noqa: ANN206
        """
        Load instance from a file.

        Parameters
        ----------
        filepath : str or Path
            Path to the input file.
        format : str or None, optional
            Input format name. If None, inferred from the file extension.
        **kwargs
            Additional format-specific options passed to the reader.

        Returns
        -------
        instance of cls
            A new instance loaded from the file.

        Raises
        ------
        PhyloZooFormatError
            If the format is not supported.
        """
        if format is None:
            format = FormatRegistry.detect_format(filepath, cls)
        content = read_file_safely(filepath)
        return cls.from_string(content, format=format, **kwargs)

    @classmethod
    def from_string(cls, string: str, format: str | None = None, **kwargs: object):  # noqa: ANN206
        """
        Create instance from a string representation.

        Parameters
        ----------
        string : str
            String content in the specified format.
        format : str or None, optional
            Input format name. If None, uses the class's default format.
        **kwargs
            Additional format-specific options passed to the reader.

        Returns
        -------
        instance of cls
            A new instance parsed from the string.

        Raises
        ------
        PhyloZooFormatError
            If the format is not supported.
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
        **kwargs: object,
    ) -> None:
        """
        Convert a file from one format to another.

        Parameters
        ----------
        input_file : str or Path
            Path to the input file.
        output_file : str or Path
            Path to the output file.
        input_format : str or None, optional
            Input format name. If None, inferred from the input file extension.
        output_format : str or None, optional
            Output format name. If None, inferred from the output file extension.
        overwrite : bool, optional
            If True, overwrite existing output file. Default is False.
        **kwargs
            Additional format-specific options passed to reader and writer.

        Raises
        ------
        FileExistsError
            If the output file already exists and overwrite is False.
        PhyloZooFormatError
            If either format is not supported.
        """
        obj = cls.load(input_file, format=input_format, **kwargs)
        obj.save(output_file, format=output_format, overwrite=overwrite, **kwargs)

    @classmethod
    def convert_string(
        cls,
        content: str,
        input_format: str,
        output_format: str,
        **kwargs: object,
    ) -> str:
        """
        Convert a string representation from one format to another.

        Parameters
        ----------
        content : str
            String content in the input format.
        input_format : str
            Input format name.
        output_format : str
            Output format name.
        **kwargs
            Additional format-specific options passed to reader and writer.

        Returns
        -------
        str
            String representation in the output format.

        Raises
        ------
        PhyloZooFormatError
            If either format is not supported.
        """
        obj = cls.from_string(content, format=input_format, **kwargs)
        return obj.to_string(format=output_format, **kwargs)
