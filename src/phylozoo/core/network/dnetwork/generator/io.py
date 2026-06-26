"""
Directed generator I/O module.

A :class:`~phylozoo.core.network.dnetwork.generator.base.DirectedGenerator` is
serialised through its underlying
:class:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph`, so it
supports exactly the directed-multigraph formats:

- **dot** (extensions: .dot, .gv, *default*) — standard Graphviz ``digraph``.
- **edgelist** (extensions: .el)

Writing delegates to the graph writer; reading parses the graph and rebuilds the
generator (so a round-tripped generator is validated on construction). The handlers
are registered with FormatRegistry when this module is imported (via the package
``__init__``), giving the generator the usual ``save``/``load``/``to_string``/
``from_string``/``convert`` methods from :class:`~phylozoo.utils.io.IOMixin`.
"""

from __future__ import annotations

from typing import Any

from phylozoo.utils.io import FormatRegistry

from ....primitives.d_multigraph.io import (
    from_dot as _dm_from_dot,
    from_edgelist as _dm_from_edgelist,
    to_dot as _dm_to_dot,
    to_edgelist as _dm_to_edgelist,
)
from .base import DirectedGenerator


def to_dot(generator: DirectedGenerator, **kwargs: Any) -> str:
    """Serialise a generator as a Graphviz DOT string (via its DirectedMultiGraph)."""
    return _dm_to_dot(generator.graph, **kwargs)


def from_dot(string: str, **kwargs: Any) -> DirectedGenerator:
    """Parse a Graphviz DOT string into a (validated) DirectedGenerator."""
    return DirectedGenerator(_dm_from_dot(string, **kwargs))


def to_edgelist(generator: DirectedGenerator, **kwargs: Any) -> str:
    """Serialise a generator as an edge-list string (via its DirectedMultiGraph)."""
    return _dm_to_edgelist(generator.graph, **kwargs)


def from_edgelist(string: str, **kwargs: Any) -> DirectedGenerator:
    """Parse an edge-list string into a (validated) DirectedGenerator."""
    return DirectedGenerator(_dm_from_edgelist(string, **kwargs))


# Register format handlers with FormatRegistry
FormatRegistry.register(
    DirectedGenerator,
    "dot",
    reader=from_dot,
    writer=to_dot,
    extensions=[".dot", ".gv"],
    default=True,
)

FormatRegistry.register(
    DirectedGenerator,
    "edgelist",
    reader=from_edgelist,
    writer=to_edgelist,
    extensions=[".el"],
)
