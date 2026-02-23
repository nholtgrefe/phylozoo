"""Shared format structure (NEXUS, PHYLIP). Class-specific (de)serialization lives in core/*/io."""

from . import nexus
from . import phylip

__all__ = ['nexus', 'phylip']
