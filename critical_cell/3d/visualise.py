#!/usr/bin/env python3
"""
3D Visualization of Critical Cellular Automaton
Reads output from the Fortran simulation and creates an animated MP4 video
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D


def parse_output_file(filename, grid_size=(20, 20, 20)):
    nx, ny, nz = grid_size
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
                            if len(parts) >= 4:
                                i = int(parts[0]) - 1  # Convert to 0-indexed
                                j = int(parts[1]) - 1
                                k = int(parts[2]) - 1
                                value = int(parts[3])

                                iterations[current_iteration][(i, j, k)] = value
                        except (ValueError, IndexError):
                            continue
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)

    grid_states = []
    sorted_iterations = sorted(iterations.keys())

    for iteration in sorted_iterations:
        grid = np.zeros((nx, ny, nz), dtype=int)
        for (i, j, k), value in iterations[iteration].items():
            if 0 <= i < nx and 0 <= j < ny and 0 <= k < nz:
                grid[i, j, k] = value
        grid_states.append(grid)

    return grid_states, sorted_iterations


def create_visualization(
    grid_states, output_file="animation.mp4", fps=30, grid_size=(20, 20, 20)
):
    nx, ny, nz = grid_size
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection="3d")

    all_values = np.concatenate([grid.flatten() for grid in grid_states])
    vmin, vmax = all_values.min(), all_values.max()

    def update(frame):
        """Update function for animation."""
        ax.clear()

        grid = grid_states[frame]
        i_coords, j_coords, k_coords = np.where(grid > 0)

        if len(i_coords) > 0:
            values = grid[i_coords, j_coords, k_coords]
            normalized = np.log1p(values) / np.log1p(vmax)
            # normalized = (normalized - np.min(normalized)) / (
            #     np.max(normalized) - np.min(normalized)
            # )
            normalized = normalized**2
            sizes = 50 + normalized * 400
            colors_rgb = plt.cm.hot(normalized)
            colors_rgba = colors_rgb.copy()
            colors_rgba[:, 3] = normalized

            scatter = ax.scatter(
                i_coords,
                j_coords,
                k_coords,
                c=colors_rgba,
                s=sizes,
                edgecolors="none",
            )

        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_xlim(0, nx)
        ax.set_ylim(0, ny)
        ax.set_zlim(0, nz)

        rotation_angle = (frame / len(grid_states)) * 360
        ax.view_init(elev=20, azim=rotation_angle)

        total_grains = np.sum(grid)
        ax.set_title(
            f"Iteration {frame + 1} | Total Grains: {total_grains}",
            fontsize=12,
            fontweight="bold",
        )

        return (ax,)

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
        description="Visualize critical cellular automaton simulation as 3D rotating video"
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
        nargs=3,
        default=[20, 20, 20],
        metavar=("NX", "NY", "NZ"),
        help="Grid dimensions (default: 20 20 20)",
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
    print(f"Grid size: {grid_size[0]} x {grid_size[1]} x {grid_size[2]}")

    create_visualization(
        grid_states, output_file=args.output, fps=args.fps, grid_size=grid_size
    )


if __name__ == "__main__":
    main()
