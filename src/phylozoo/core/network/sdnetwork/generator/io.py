"""
Semi-directed generator I/O module.

A :class:`~phylozoo.core.network.sdnetwork.generator.base.SemiDirectedGenerator`
is serialised through its underlying
:class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph`, so it
supports exactly the mixed-multigraph DOT formats:

- **phylozoo-dot** (extensions: .pzdot, *default*)
- **dot** (extensions: .dot, .gv) — standard Graphviz, with undirected edges
  marked ``dir=none``.

Writing delegates to the graph writer; reading parses the graph and rebuilds the
generator (so a round-tripped generator is validated on construction). The handlers
are registered with FormatRegistry when this module is imported (via the package
``__init__``), giving the generator the usual ``save``/``load``/``to_string``/
``from_string``/``convert`` methods from :class:`~phylozoo.utils.io.IOMixin`.
"""

from __future__ import annotations

from typing import Any

from phylozoo.utils.io import FormatRegistry

from ....primitives.m_multigraph.io import (
    from_dot as _mm_from_dot,
    from_phylozoo_dot as _mm_from_phylozoo_dot,
    to_dot as _mm_to_dot,
    to_phylozoo_dot as _mm_to_phylozoo_dot,
)
from .base import SemiDirectedGenerator


def to_phylozoo_dot(generator: SemiDirectedGenerator, **kwargs: Any) -> str:
    """Serialise a generator as a phylozoo-dot string (via its MixedMultiGraph)."""
    return _mm_to_phylozoo_dot(generator.graph, **kwargs)


def from_phylozoo_dot(string: str, **kwargs: Any) -> SemiDirectedGenerator:
    """Parse a phylozoo-dot string into a (validated) SemiDirectedGenerator."""
    return SemiDirectedGenerator(_mm_from_phylozoo_dot(string, **kwargs))


def to_dot(generator: SemiDirectedGenerator, **kwargs: Any) -> str:
    """Serialise a generator as a standard Graphviz DOT string (via its MixedMultiGraph)."""
    return _mm_to_dot(generator.graph, **kwargs)


def from_dot(string: str, **kwargs: Any) -> SemiDirectedGenerator:
    """Parse a standard Graphviz DOT string into a (validated) SemiDirectedGenerator."""
    return SemiDirectedGenerator(_mm_from_dot(string, **kwargs))


# Register format handlers with FormatRegistry
FormatRegistry.register(
    SemiDirectedGenerator,
    "phylozoo-dot",
    reader=from_phylozoo_dot,
    writer=to_phylozoo_dot,
    extensions=[".pzdot"],
    default=True,
)

FormatRegistry.register(
    SemiDirectedGenerator,
    "dot",
    reader=from_dot,
    writer=to_dot,
    extensions=[".dot", ".gv"],
)
