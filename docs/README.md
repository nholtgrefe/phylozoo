# Documentation

This directory contains the Sphinx documentation for phylozoo.

## Building Documentation Locally

### Quick Build

To build the documentation:

```bash
cd docs
make html
```

Then open `build/html/index.html` in your browser.

### Automatic Rebuild (Recommended for Development)

To automatically rebuild documentation when files change:

```bash
cd docs
sphinx-autobuild source build/html --open-browser
```

This will:
- Start a local server (usually at `http://127.0.0.1:8000`)
- Watch for changes and rebuild automatically
- Auto-reload the browser when docs update

### Manual Build

To build documentation manually:

```bash
cd docs
sphinx-build -b html source build/html
```

Then open `build/html/index.html` in your browser.

## Structure

- `source/` - Source RST files
  - `manual/` - User manual (installation, introduction, etc.)
  - `tutorials/` - Extended tutorials
  - `api/` - API reference
  - `develop/` - Developer guides
  - `reference/` - Reference material (bibliography, changelog)
- `build/` - Built HTML output
- `source/conf.py` - Sphinx configuration
- `requirements.txt` - Documentation dependencies

