"""
Script to generate and plot all level-4 semi-directed generators.

This script:
1. Generates all level-4 semi-directed generators
2. Prints the number of generators found
3. Plots all generators in a single figure using subplots
"""

import math

import matplotlib.pyplot as plt

from phylozoo.core.network.sdnetwork.generator import all_level_k_generators
from phylozoo.visualize import plot_mixed_multigraph


def main() -> None:
    """Generate and plot all level-4 semi-directed generators."""
    print("Generating all level-4 semi-directed generators...")
    
    # Generate all level-4 generators
    generators = all_level_k_generators(3)
    
    num_generators = len(generators)
    print(f"\nFound {num_generators} level-4 semi-directed generators")
    
    if num_generators == 0:
        print("No generators to plot.")
        return
    
    # Calculate grid dimensions for subplots
    # Try to make it roughly square
    cols = math.ceil(math.sqrt(num_generators))
    rows = math.ceil(num_generators / cols)
    
    # Create figure with subplots
    # Use squeeze=False to always get a 2D array, then flatten
    fig, axes_2d = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows), squeeze=False)
    
    # Flatten to 1D list
    axes = [axes_2d[i, j] for i in range(rows) for j in range(cols)]
    
    # Plot each generator
    for idx, generator in enumerate(generators):
        ax = axes[idx]
        
        # Plot the generator's graph
        plot_mixed_multigraph(
            generator.graph,
            layout='spring',
            backend='matplotlib',
            ax=ax,
        )
        
        # Set title with generator index
        ax.set_title(f'Generator {idx + 1}', fontsize=10, fontweight='bold')
        ax.axis('off')  # Remove axes for cleaner look
    
    # Hide unused subplots
    for idx in range(num_generators, len(axes)):
        axes[idx].axis('off')
    
    # Adjust layout
    plt.suptitle(
        f'All Level-4 Semi-Directed Generators ({num_generators} total)',
        fontsize=14,
        fontweight='bold',
        y=0.995
    )
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    
    # Show the plot
    print("\nDisplaying plot...")
    plt.show()


if __name__ == '__main__':
    main()

