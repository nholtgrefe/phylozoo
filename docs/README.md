# Documentation

This directory contains the Sphinx documentation for phylozoo.

## Building Documentation

### Automatic Rebuild (Recommended)

To automatically rebuild documentation when files change:

```bash
sphinx-autobuild source build/html --open-browser
```

This will:
- Start a local server (usually at `http://127.0.0.1:8000`)
- Watch for changes and rebuild automatically
- Auto-reload the browser when docs update

### Manual Build

To build documentation manually:

```bash
sphinx-build -b html source build/html
```

Then open `build/html/index.html` in your browser.

## Structure

- `source/` - Source RST files
- `build/` - Built HTML output
- `source/conf.py` - Sphinx configuration
- `source/api/` - API documentation

