# How `pip install -e .` Works

## What Does `-e` Mean?

The `-e` flag stands for **"editable"** mode (also called "development mode"). When you run:

```bash
pip install -e .
```

Instead of copying your package files to the Python site-packages directory, pip creates a **link** (symlink or `.pth` file) that points to your source code location.

## What Gets Installed?

When you run `pip install -e .`, pip:

1. **Reads your `pyproject.toml`** to understand the package structure
2. **Creates a link** in your Python environment pointing to your source code
3. **Installs dependencies** listed in `dependencies` or `[project.optional-dependencies]`
4. **Makes the package importable** as `import phylozoo`

## Does Everything Update Automatically?

### ✅ **YES - Code Changes Update Automatically**

**Python code changes are immediately available** without reinstalling:

```bash
# 1. Install in editable mode
pip install -e .

# 2. Edit src/phylozoo/splits.py (add a new function)
# 3. In Python, the new function is immediately available:
python
>>> import phylozoo
>>> from phylozoo.splits import Split
>>> # Your changes are already there!
```

**Why?** Because Python imports directly from your source directory via the link.

### ⚠️ **SOMETIMES - You Need to Restart**

Some things require restarting:

1. **Python interpreter** - If you have a Python session open, restart it:
   ```python
   # Old session - won't see changes
   >>> import phylozoo
   # Close and restart Python
   ```

2. **Jupyter notebooks** - Restart the kernel:
   ```python
   # In Jupyter: Kernel → Restart Kernel
   ```

3. **Long-running scripts** - Restart the script

4. **Cached imports** - Python caches imports, but usually picks up changes

### ❌ **NO - These Don't Update Automatically**

1. **New files/directories** - If you add a new module:
   ```bash
   # After adding src/phylozoo/new_module.py
   # You might need to reinstall:
   pip install -e . --force-reinstall --no-deps
   ```

2. **Package metadata** - Changes to `pyproject.toml` (name, version, dependencies):
   ```bash
   # After changing dependencies in pyproject.toml:
   pip install -e .
   ```

3. **New dependencies** - If you add dependencies to `pyproject.toml`:
   ```bash
   # After adding new dependencies:
   pip install -e .
   ```

## How It Works Technically

### Normal Install (without `-e`):
```
Your code → Copied to → site-packages/phylozoo/
                              ↓
                        Python imports from here
```

### Editable Install (with `-e`):
```
Your code → Link created → site-packages/phylozoo.egg-link
     ↓                          ↓
Python imports directly from your source directory
```

The `.egg-link` file contains the path to your source code, so Python knows where to find it.

## What Gets Created?

After `pip install -e .`, you'll see:

1. **In your environment:**
   - `phylozoo.egg-link` file pointing to your source
   - Package metadata in `*.egg-info/`

2. **In your project:**
   - `src/phylozoo.egg-info/` directory (package metadata)

## When to Reinstall

You need to reinstall (`pip install -e .`) when:

- ✅ Adding new dependencies to `pyproject.toml`
- ✅ Changing package name or version
- ✅ Adding new subpackages/modules (sometimes)
- ✅ Changing `package-dir` or package structure in `pyproject.toml`

You **don't** need to reinstall when:

- ✅ Editing existing Python files
- ✅ Adding new functions/classes to existing files
- ✅ Fixing bugs in existing code
- ✅ Changing docstrings

## Example Workflow

```bash
# 1. Initial setup (once)
pip install -e ".[dev]"

# 2. Make changes to code
vim src/phylozoo/splits.py

# 3. Test immediately (no reinstall needed!)
pytest tests/test_splits.py

# 4. If you add a new dependency to pyproject.toml:
vim pyproject.toml  # Add new dependency
pip install -e .     # Reinstall to get new dependency

# 5. Continue coding - changes are live!
```

## Checking Your Installation

See where your package is installed:

```bash
python -c "import phylozoo; print(phylozoo.__file__)"
```

This will show the path to your source code (not site-packages), confirming editable mode.

## Troubleshooting

### Changes not showing up?

1. **Restart Python/Jupyter** - Most common issue
2. **Check you're in the right environment:**
   ```bash
   which python  # Should show your venv
   ```
3. **Reinstall:**
   ```bash
   pip install -e . --force-reinstall
   ```

### Import errors after adding new files?

```bash
# Reinstall to register new modules
pip install -e . --force-reinstall --no-deps
```

### Package not found?

```bash
# Make sure you're in the project root
cd /home/nholtgreve/Documents/Code/phylozoo
pip install -e .
```

## Summary

- **Code changes**: ✅ Update automatically (restart Python if needed)
- **New dependencies**: ❌ Need `pip install -e .`
- **Package structure changes**: ❌ Need `pip install -e .`
- **New files**: ⚠️ Usually work, sometimes need reinstall

The editable install is perfect for development because you can code and test immediately without reinstalling constantly!

