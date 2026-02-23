"""
Low-level file operations for phylozoo I/O.

Provides read_file_safely, write_file_safely, and ensure_directory_exists
used by the IOMixin and format handlers.
"""

from __future__ import annotations

from pathlib import Path

from phylozoo.utils.exceptions import PhyloZooIOError, PhyloZooValueError


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
