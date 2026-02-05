#!/usr/bin/env python3
"""
2D Visualization of Critical Cellular Automaton
Reads output from the Fortran simulation and creates an animated MP4 video
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation


def parse_output_file(filename, grid_size=(20, 20)):
    nx, ny = grid_size
    iterations = {}

    try:
        with open(filename, "r") as f:
            current_iteration = None

            for line in f:
                line = line.strip()
                if not line:
                    continue
                if "Iteration:" in line:
                    parts = line.split(":")
                    if len(parts) >= 2:
                        try:
                            current_iteration = int(parts[1].strip())
                            if current_iteration not in iterations:
                                iterations[current_iteration] = {}
                        except ValueError:
                            continue
                else:
                    if current_iteration is not None:
                        try:
                            parts = line.split()
                            if len(parts) >= 3:
                                i = int(parts[0]) - 1  # Convert to 0-indexed
                                j = int(parts[1]) - 1
                                value = int(parts[2])

                                iterations[current_iteration][(i, j)] = value
                        except (ValueError, IndexError):
                            continue
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)

    grid_states = []
    sorted_iterations = sorted(iterations.keys())

    for iteration in sorted_iterations:
        grid = np.zeros((nx, ny), dtype=int)
        for (i, j), value in iterations[iteration].items():
            if 0 <= i < nx and 0 <= j < ny:
                grid[i, j] = value
        grid_states.append(grid)

    return grid_states, sorted_iterations


def create_visualization(
    grid_states, output_file="animation.mp4", fps=30, grid_size=(20, 20)
):
    nx, ny = grid_size
    fig, ax = plt.subplots(figsize=(8, 8))

    # Create initial plot with colorbar outside the update loop
    grid = grid_states[0]
    display_grid = np.clip(grid, 0, 4).astype(float)
    im = ax.imshow(display_grid, cmap="hot", vmin=0, vmax=4, origin="lower")
    cbar = plt.colorbar(im, ax=ax, label="Grain Count")

    def update(frame):
        """Update function for animation."""
        grid = grid_states[frame]

        # Update image data instead of clearing
        display_grid = np.clip(grid, 0, 4).astype(float)
        im.set_data(display_grid)

        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        total_grains = np.sum(grid)
        ax.set_title(
            f"Iteration {frame + 1} | Total Grains: {total_grains}",
            fontsize=12,
            fontweight="bold",
        )

        return (im,)

    print(f"Creating animation with {len(grid_states)} frames...")
    anim = FuncAnimation(
        fig,
        update,
        frames=len(grid_states),
        interval=1000 // fps,
        blit=False,
        repeat=True,
    )

    # Save as MP4
    print(f"Saving animation to '{output_file}'...")
    try:
        from matplotlib.animation import FFMpegWriter

        writer = FFMpegWriter(fps=fps, bitrate=1800)
        anim.save(output_file, writer=writer, dpi=100)
        print(f"✓ Animation saved successfully to '{output_file}'")
    except Exception as e:
        print(f"Error saving animation: {e}")
        print("Make sure you have ffmpeg installed:")
        print("  macOS: brew install ffmpeg")
        print("  Ubuntu: sudo apt-get install ffmpeg")
        sys.exit(1)

    plt.close(fig)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Visualize 2D critical cellular automaton simulation as heatmap video"
    )
    parser.add_argument(
        "--input",
        "-i",
        default="output.txt",
        help="Input txt file from Fortran simulation (default: output.txt)",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="animation.mp4",
        help="Output MP4 file (default: animation.mp4)",
    )
    parser.add_argument(
        "--fps", type=int, default=30, help="Frames per second (default: 30)"
    )
    parser.add_argument(
        "--grid-size",
        type=int,
        nargs=2,
        default=[40, 40],
        metavar=("NX", "NY"),
        help="Grid dimensions (default: 40 40)",
    )

    args = parser.parse_args()
    grid_size = tuple(args.grid_size)

    if not Path(args.input).exists():
        print(f"Error: Input file '{args.input}' not found.")
        sys.exit(1)

    print(f"Reading simulation data from '{args.input}'...")
    grid_states, iterations = parse_output_file(args.input, grid_size=grid_size)

    if not grid_states:
        print("Error: No valid data found in input file.")
        sys.exit(1)

    print(f"Parsed {len(grid_states)} simulation frames")
    print(f"Grid size: {grid_size[0]} x {grid_size[1]}")

    create_visualization(
        grid_states, output_file=args.output, fps=args.fps, grid_size=grid_size
    )


if __name__ == "__main__":
    main()
