"""
PhyloZoo I/O utilities (phylozoo.utils.io).

Provides FormatRegistry, IOMixin, and file operations. Shared format structure
(NEXUS, PHYLIP) is in phylozoo.utils.io.formats; class-specific readers/writers
are in core/distance/io, core/sequence/io, core/split/io, etc.
"""

from .file_ops import ensure_directory_exists, read_file_safely, write_file_safely
from .mixin import IOMixin
from .registry import FormatRegistry

__all__ = [
    'FormatRegistry',
    'IOMixin',
    'read_file_safely',
    'write_file_safely',
    'ensure_directory_exists',
]
