# AGENT_PROJECT_SPECIFICS.md — Project-specific rules, PR checklist and gotchas

This file contains the maintainer's exact preferences and gotchas for `phylozoo`.

Core rules (copy/paste from .instructions)

- ALWAYS use type hints for all parameters, returns and class attributes.
- ALWAYS use NumPy-style docstrings and mirror examples in tests.
- ALWAYS write tests for new behavior and include docstrings in tests.
- Prefer SemiDirectedPhyNetwork in tests unless specifically testing MixedPhyNetwork.
- NEVER commit changes without maintainer approval (use branches/PRs only).

Tooling and config notes

- `pyproject.toml` configures pytest, mypy, black, ruff and project dependencies.
- Target Python versions: >=3.10 (project uses pyproject for exact targets).
- Black line-length: 100; Ruff uses the same line length.

Warnings and CI

- Project declares filtered warnings in `pyproject.toml`; do not reintroduce noisy warnings.
- Use pytest markers `slow` and `integration` when appropriate.

Test path hack

- `tests/conftest.py` inserts `src/` into `sys.path`. For robust development prefer `pip install -e .`.

Immutable containers and numpy

- Several core containers (e.g., `DistanceMatrix`) make their underlying arrays read-only and store labels as tuples. Avoid in-place changes to `.np_array`.
 - Several core containers (e.g., `DistanceMatrix`) make their underlying arrays read-only and store labels as tuples. Avoid in-place changes to `.np_array`.
 - `DistanceMatrix` uses a dict for O(1) label lookups. Labels MUST be hashable (int, str, tuple, frozenset). If you need to use composite/mutable labels, convert them to a hashable canonical form first.

Numba guidance

- Numba-compiled functions must receive proper dtypes and contiguous arrays.
- Provide pure-Python wrappers for readability and testing.
- If you must change numba-compiled code, include a micro-benchmark and tests to show no performance regression.

PR checklist for AI agents (must follow before creating a PR)

1. Add/modify tests covering the change.
2. Run `pytest -q` and ensure tests pass locally.
3. Run `mypy src/phylozoo` and fix type issues.
4. Run `black src/phylozoo tests` and `ruff check src/phylozoo tests`.
5. Add a short, focused commit message.
6. Open a PR with a clear description and links to any benchmarks or failing/passing examples.

Agent VCS hard rule

- Agents MUST NOT perform any VCS operations: do not create branches, make commits, push, merge, rebase, or otherwise interact with the repository's version control. Produce changes as patches/diffs or suggested edits for a human maintainer to apply.

Legacy Python compatibility note

- The original `.instructions` states broader Python compatibility (>=3.7 supporting 3.7-3.11). The authoritative requirement for tooling and CI is `pyproject.toml` (>=3.10). When in doubt, follow `pyproject.toml`.

Performance enforcement

- The `.instructions` emphasises: "Always perform optimizations of code efficiency." Prefer correctness and tests first. For performance-sensitive work include micro-benchmarks, document results in `sandbox/findings.txt`, and attach those findings to the PR description.

If anything in this file conflicts with the maintainer's intent, ask before changing behavior.