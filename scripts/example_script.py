"""
Example script for phylozoo.

This is a placeholder script demonstrating the structure for utility scripts.
"""
import argparse
from typing import Optional


def main(input_file: Optional[str] = None, output_file: Optional[str] = None) -> None:
    """
    Main function for the example script.

    Parameters
    ----------
    input_file : Optional[str], optional
        Path to input file, by default None
    output_file : Optional[str], optional
        Path to output file, by default None
    """
    print("Example script for phylozoo")
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")
    # Placeholder implementation
    print("Script functionality not yet implemented")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Example phylozoo script")
    parser.add_argument(
        "-i", "--input", type=str, help="Input file path", default=None
    )
    parser.add_argument(
        "-o", "--output", type=str, help="Output file path", default=None
    )
    args = parser.parse_args()

    main(input_file=args.input, output_file=args.output)

