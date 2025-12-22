## Quick guide for AI coding agents working on PhyloZoo

This repo is a Python package (src layout) for phylogenetic network analysis. The guidance below is intentionally prescriptive and contains project-specific rules copied from the project's `.instructions` file so agents behave consistently with the maintainer's preferences.

-----------------------------------------------------------------
Core project layout & quick commands
-----------------------------------------------------------------

- Repo layout to know:
  - `src/phylozoo/` — main package. Many modules live under `core/`, `plotting/`, `inference/`.
  - `examples/` — runnable scripts demonstrating common usage (useful for smoke tests).
  - `tests/` — pytest suite. `tests/conftest.py` inserts `src/` on `sys.path` for tests (preferred: install editable to avoid path hacks).
  - `docs/` — Sphinx docs; has a Makefile and `autobuild.sh`.

- Important files/configs:
  - `pyproject.toml` — build, dependencies, and test/dev extras. Python >= 3.10.
  - `pyproject.toml[tool.pytest.ini_options]` — pytest discovery, addopts, markers (`slow`, `integration`).
  - `pyproject.toml[tool.mypy|black|ruff]` — mypy/black/ruff settings (line length 100, py310/311 target).

- Dev/run commands (copyable):
  - Install editable: `pip install -e .` or `pip install -e "[dev]"` for dev deps.
  - Run tests: `pytest` (uses config in `pyproject.toml`).
  - Type check: `mypy src/phylozoo`.
  - Format: `black src/phylozoo tests`.
  - Lint: `ruff check src/phylozoo tests`.
  - Build docs: `make -C docs html` or `docs/autobuild.sh`.

-----------------------------------------------------------------
Project-specific rules (from `.instructions`) — follow exactly unless you get explicit approval
-----------------------------------------------------------------

Commits
- Make commit messages very short but still descriptive.

Type hinting
- ALWAYS use type hints for all function parameters, return types, and class attributes.
- Prefer modern built-in generics and PEP 604 unions (e.g., `list[int]`, `dict[str, float]`, `int | None`) — Python >= 3.10 is required.
- Use `TypeVar` for generic types when appropriate.
- Prefer explicit types over `Any`—only use `Any` when absolutely necessary.
- Use `-> None` for functions that don't return values.
- Type hint all class methods, including `__init__` and `__repr__`.

Documentation
- ALWAYS use NumPy-style docstrings for classes, functions, and methods. Include relevant sections (Parameters, Returns, Raises, Examples, Notes, See Also).
- Example docstring template (copy into new functions):

```python
def function_name(param1: Type, param2: Type) -> ReturnType:
    """
    Brief description of the function.
    
    Parameters
    ----------
    param1 : Type
        Description of param1.
    param2 : Type
        Description of param2.
    
    Returns
    -------
    ReturnType
        Description of what is returned.
    
    Raises
    ------
    ValueError
        When and why this exception is raised.
    
    Examples
    --------
    >>> example_usage()
    expected_output
    
    Notes
    -----
    Additional notes or implementation details.
    """
```

- When adding or updating examples in docstrings, mirror those examples in tests.

Testing
- ALWAYS write tests for new methods and classes. Use pytest and place tests in `tests/`.
- Test naming: `test_*.py` and `test_*` functions. Group related tests in `Test*` classes when appropriate.
- Include docstrings in test functions describing intent.
- Use type hints in test functions too.
- Aim to cover success cases and edge cases. If a failing test shows desired behavior ambiguity, ask the maintainer rather than silently changing expected behavior.
- Prefer using `SemiDirectedPhyNetwork` in tests unless testing `MixedPhyNetwork` specifically.

Workflow preferences
- Prefer making direct code changes over running command-line operations. Run commands only when necessary (installing deps, running tests).

Code style & structure
- Follow existing code patterns in the repo. Keep functions small and focused; prefer functions for workflows and methods for basic functionality.
- Use meaningful names, add error handling where appropriate.

Project-specific notes
- This is a phylogenetic networks analysis package. Many modules use `numpy` and `numba` for speed.
- Python version: project requires >=3.10 (pyproject.toml).
- Testing framework: pytest. Type checking: mypy. Formatting: black (line-length 100). Linting: ruff.
- Use Numba where it provides a measurable speedup, but respect numba patterns (avoid rewriting JIT-compatible code into forms numba can't compile).

Performance
- Optimize code efficiency when appropriate, but prefer correctness and test coverage first. For performance-sensitive changes add a micro-benchmark or test.

Findings & experimentation
- Keep a `findings.txt` in `sandbox/` documenting experiments and failed tries so human reviewers can reuse insights.

Important rules
- NEVER commit changes. (Follow the maintainer's instruction: run experiments and propose edits, but do not push commits on your own.)
- When in doubt about behavior, ask the maintainer rather than changing tests to match current behavior.

Agent VCS restriction
- Agents MUST NOT create branches, commits, push, merge, or perform any VCS operations. All suggested changes should be produced as patches, diffs, or proposed edits for a human maintainer to review and commit.

-----------------------------------------------------------------
Project patterns & gotchas (practical examples)
-----------------------------------------------------------------

- Test path hack: `tests/conftest.py` inserts `src/` into `sys.path` so tests can import `phylozoo` without installation. Prefer `pip install -e .` for development to avoid relying on this hack.
- Immutable distance containers: `src/phylozoo/core/distance/base.py::DistanceMatrix` sets the underlying numpy array to read-only and stores labels as tuples. Example usage:

```python
from phylozoo.core.distance.base import DistanceMatrix
import numpy as np
dm = DistanceMatrix(np.array([[0,1],[1,0]]), labels=["A","B"])
assert dm.np_array.flags.writeable is False
```

Note: `DistanceMatrix` builds an internal dictionary for O(1) label lookups. Labels must be hashable (e.g., `int`, `str`, `tuple`). If you pass mutable labels convert them to hashable equivalents first.

- Numba caution: many modules use numba and expect specific array layouts and types. If adding vectorized numpy code, ensure it does not break numba-compiled paths.

Legacy Python compatibility note
- The project's historical `.instructions` mentions supporting Python >=3.7 (3.7–3.11); however `pyproject.toml` currently requires Python >=3.10. Treat `pyproject.toml` as authoritative for tooling and CI, but be aware older code or docs may reference wider compatibility.

Performance note
- The original instructions also state: "Always perform optimizations of code efficiency." Prefer correctness and tests first, but include micro-benchmarks for performance changes and document findings in `sandbox/findings.txt`.

-----------------------------------------------------------------
How to propose changes (PR checklist for AI agents)
-----------------------------------------------------------------

1. Make minimal edits in code; add tests under `tests/` for new behavior.
2. Run `pytest -q` locally and ensure tests pass.
3. Run `mypy src/phylozoo` and fix type issues introduced by the change.
4. Run `black src/phylozoo tests` and `ruff check src/phylozoo tests`.
5. Add a short, focused commit message (the maintainer prefers concise messages).
6. Push a branch and create a PR with description and links to failing/passing examples and any micro-benchmarks.

-----------------------------------------------------------------
If you'd like a longer `AGENT.md` with templates (unit-test template, numba-safe checklist, micro-benchmark harness), tell me which templates to include and I'll add them.

-----------------------------------------------------------------
Notes & contact
-----------------------------------------------------------------

If parts of this file are unclear or you want a different tone/level of detail, say which section to expand or remove and I will iterate.
