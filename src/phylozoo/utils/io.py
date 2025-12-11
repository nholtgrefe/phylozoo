"""
IO utilities module.

This module provides shared utilities for file I/O operations that can be
reused across different classes and modules.
"""

from pathlib import Path
from typing import List


def check_file_exists(filepath: str | Path) -> Path:
    """
    Check if a file exists and return Path object.
    
    Parameters
    ----------
    filepath : str | Path
        Path to the file to check.
    
    Returns
    -------
    Path
        Path object if file exists.
    
    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the path exists but is not a file.
    
    Examples
    --------
    >>> from phylozoo.utils.io import check_file_exists
    >>> path = check_file_exists("data.fasta")
    >>> path
    PosixPath('data.fasta')
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    if not path.is_file():
        raise ValueError(f"Path is not a file: {filepath}")
    return path


def validate_file_extension(filepath: str | Path, allowed_extensions: list[str]) -> None:
    """
    Validate that a file has an allowed extension.
    
    Parameters
    ----------
    filepath : str | Path
        Path to the file to validate.
    allowed_extensions : list[str]
        List of allowed extensions (e.g., ['.fasta', '.fa', '.fas']).
    
    Raises
    ------
    ValueError
        If the file extension is not in the allowed list.
    
    Examples
    --------
    >>> from phylozoo.utils.io import validate_file_extension
    >>> validate_file_extension("data.fasta", ['.fasta', '.fa'])
    >>> validate_file_extension("data.txt", ['.fasta', '.fa'])
    Traceback (most recent call last):
        ...
    ValueError: File extension '.txt' not allowed...
    """
    path = Path(filepath)
    if path.suffix.lower() not in [ext.lower() for ext in allowed_extensions]:
        raise ValueError(
            f"File extension '{path.suffix}' not allowed. "
            f"Allowed extensions: {allowed_extensions}"
        )


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
    IOError
        If the file cannot be read.
    
    Examples
    --------
    >>> from phylozoo.utils.io import read_file_safely
    >>> content = read_file_safely("data.txt")
    >>> len(content) > 0
    True
    """
    path = check_file_exists(filepath)
    try:
        return path.read_text(encoding=encoding)
    except Exception as e:
        raise IOError(f"Error reading file {filepath}: {e}") from e


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
    IOError
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
        raise IOError(f"Error writing file {filepath}: {e}") from e


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

