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
from matplotlib.colors import ListedColormap


def parse_output_file(filename, grid_size=(20, 20)):
    nx, ny = grid_size
    grid_states = []
    current_grid = np.zeros((nx, ny), dtype=int)
    sorted_iterations = []

    try:
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith("#INIT"):
                    current_grid = np.zeros((nx, ny), dtype=int)
                elif line.startswith("#D"):
                    grid_states.append(current_grid.copy())
                    try:
                        iteration = int(line[2:])
                        sorted_iterations.append(iteration)
                        if iteration % 1000 == 0:
                            print(f"Found iteration: {iteration}")
                    except ValueError:
                        pass
                else:
                    try:
                        parts = line.split(",")
                        if len(parts) >= 3:
                            i = int(parts[0]) - 1
                            j = int(parts[1]) - 1
                            value = int(parts[2])

                            if 0 <= i < nx and 0 <= j < ny:
                                current_grid[i, j] = value
                    except (ValueError, IndexError):
                        continue

            if len(grid_states) == 0 or not np.array_equal(
                grid_states[-1], current_grid
            ):
                grid_states.append(current_grid.copy())

    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)

    return grid_states, sorted_iterations


def create_visualization(
    grid_states, output_file="animation.mp4", fps=30, grid_size=(40, 40)
):
    nx, ny = grid_size
    fig, ax = plt.subplots(figsize=(8, 8))

    max_value = max(np.max(grid) for grid in grid_states)

    grid = grid_states[0]
    im = ax.imshow(grid, cmap="Greys", vmin=0, vmax=max_value, origin="lower")
    cbar = plt.colorbar(im, ax=ax, label="Grain Count")

    def update(frame):
        """Update function for animation."""
        grid = grid_states[frame]

        im.set_data(grid)

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

    print(f"Saving animation to '{output_file}'...")
    try:
        from matplotlib.animation import FFMpegWriter

        writer = FFMpegWriter(
            fps=fps, bitrate=2000, extra_args=["-preset", "ultrafast"]
        )
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
        default="output.csv",
        help="Input txt file from Fortran simulation (default: output.csv)",
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
