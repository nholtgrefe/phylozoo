# AGENT_OVERVIEW.md — High-level architecture and rationale

This overview helps an AI agent understand the "big picture" of `phylozoo` so changes are consistent with architectural intent.

Major components

- `src/phylozoo/core/` — Core algorithms and data structures (splits, networks, primitives).
- `src/phylozoo/inference/` — Inference and estimation code (network inference algorithms).
- `src/phylozoo/plotting/` — Plotting utilities for networks and examples.
- `src/phylozoo/utils/` — Misc utilities used across modules.
- `examples/` — Runnable scripts showing how the package is used end-to-end.
- `tests/` — Pytest-based test suite; tests exercise core algorithms and examples.

Data flows and service boundaries

- Data typically flows from input parsers/loaders (e.g. `load_data.py` or example scripts) into core data structures (Split, Network classes), then into inference algorithms and finally into plotting or serialization layers.
- Core data types are immutable when possible (e.g., `DistanceMatrix`). Prefer creating new instances for modified state rather than mutating in-place.

Design rationale

- Performance: NumPy and Numba are used for computational hotspots. This motivates design choices like contiguous arrays, explicit dtypes, and avoiding Python-level loops in hot paths.
- Immutability: Several data containers are intentionally immutable to prevent subtle bugs and to make reasoning about state easier.
- src-layout: The package follows a src/ layout so installation is straightforward with editable installs and tests can also add `src/` to `sys.path` for ad-hoc runs.

Key files to inspect for architecture intents

- `pyproject.toml` — dependency and tool configuration, Python version requirement.
- `tests/conftest.py` — inserts `src/` into `sys.path` for test runs; tests are configured in `pyproject.toml`.
- `examples/` — run these after editable install to confirm runtime behavior.
- `docs/` — Sphinx config and docs building scripts.

When modifying architecture

- Small, localized changes are preferred. For larger design changes, create a design note in `sandbox/findings.txt` and request maintainer review before major refactors.
- Always include tests and a clear migration/compatibility note if public APIs are modified.
