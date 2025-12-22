# AGENT.md — General coding style & workflow (for AI agents)

This document contains the general, project-agnostic rules an AI coding agent should follow when making edits in `phylozoo`.
Keep changes small and easy for a human reviewer to validate.

- Commits
  - Make commit messages very short but descriptive.
  - Do not push changes directly; propose via a branch and PR.

- Type hints
  - Always provide complete type hints for function parameters, return types and class attributes.
  - Prefer built-in generics and PEP 604 unions (e.g. `list[int]`, `int | None`).

- Documentation
  - Use NumPy-style docstrings for modules, classes and functions.
  - Include Parameters, Returns, Raises, Examples, and Notes where relevant.

- Tests
  - Add tests for any new public behavior. Use pytest and place tests in `tests/`.
  - Test names: `test_*.py`, test functions `test_*`.
  - Include docstrings in tests explaining intent.

- Formatting & linting
  - Run `black` and `ruff` locally before proposing changes.
  - Respect project settings: line-length 100, target Python 3.10/3.11.

- Type checking
  - Run `mypy src/phylozoo` for API-affecting changes.

- Workflow
  - Prefer changing code over shell operations. Run commands only to validate changes (install/edit/run tests).
  - Install editable for development: `pip install -e .` or `pip install -e "[dev]"`.

- Safety and findings
  - Keep a `findings.txt` in `sandbox/` for experiments and failed attempts.
  - Never exfiltrate secrets or make network calls.

- VCS restriction for agents
  - Agents MUST NOT create branches, make commits, push, merge, or perform any git/VCS operations. Produce suggested edits as patches/diffs for a human maintainer to review and apply.

- Chat responses
    - Be brief in explanations and to the point. For example, when explaining code changes, focus on the key modifications without excessive detail.

Note: The project's `.instructions` historically mentions Python >=3.7 compatibility; the current `pyproject.toml` requires Python >=3.10 — treat `pyproject.toml` as authoritative.

If you follow these guidelines, a human reviewer will be able to understand and validate automated edits quickly.