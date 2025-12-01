# Scripts

This directory contains utility scripts for working with phylozoo.

## Available Scripts

- `example_script.py` - Example script demonstrating script structure

## Running Scripts

To run a script:

```bash
python scripts/example_script.py -i input.txt -o output.txt
```

Or from the project root:

```bash
python -m scripts.example_script -i input.txt -o output.txt
```

## Creating New Scripts

When creating new scripts:

1. Use descriptive filenames
2. Include a docstring at the top explaining what the script does
3. Use argparse for command-line arguments
4. Use type hints for all functions
5. Include a `main()` function that can be called directly
6. Add proper error handling
7. Use logging instead of print statements for production scripts

Example structure:

```python
"""
Script description.

This script does [something].
"""
import argparse
import logging
from typing import Optional

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main(input_file: Optional[str] = None) -> None:
    """
    Main function for the script.

    Parameters
    ----------
    input_file : Optional[str], optional
        Path to input file, by default None
    """
    logger.info("Starting script")
    # Script implementation here
    logger.info("Script complete")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script description")
    parser.add_argument("-i", "--input", type=str, help="Input file path")
    args = parser.parse_args()

    main(input_file=args.input)
```

