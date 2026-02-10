#!/usr/bin/env python3
"""
2D Visualization of Critical Cellular Automaton
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation


def parse_output_file(filename, grid_size=(20, 20)):
    """Parse output file efficiently using NumPy."""
    nx, ny = grid_size

    # First pass: count frames
    with open(filename, "r") as f:
        num_frames = sum(1 for line in f if line.startswith("#D"))

    if num_frames == 0:
        return [], []

    # Pre-allocate storage
    grid_states = np.zeros((num_frames + 1, nx, ny), dtype=np.int32)
    iterations = []
    current_frame = 0

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("#INIT"):
                current_frame = 0
                grid_states[0, :, :] = 0
            elif line.startswith("#D"):
                current_frame += 1
                # Copy previous state
                grid_states[current_frame] = grid_states[current_frame - 1]
                try:
                    iteration = int(line[2:])
                    iterations.append(iteration)
                    if iteration % 1000 == 0:
                        print(f"Found iteration: {iteration}")
                except ValueError:
                    iterations.append(current_frame)
            else:
                try:
                    i, j, value = map(int, line.split(",")[:3])
                    i, j = i - 1, j - 1  # Convert to 0-indexed

                    if 0 <= i < nx and 0 <= j < ny:
                        grid_states[current_frame, i, j] = value
                except (ValueError, IndexError):
                    continue

    return grid_states[: current_frame + 1], iterations


def create_visualization(
    grid_states, output_file="animation.mp4", fps=30, grid_size=(40, 40)
):
    """Create optimized animation with proper blitting."""
    fig, ax = plt.subplots(figsize=(8, 8))

    max_value = np.max(grid_states)

    # Initial setup
    im = ax.imshow(
        grid_states[0],
        cmap="Greys",
        vmin=0,
        vmax=max_value,
        origin="lower",
        interpolation="nearest",  # Faster rendering
    )
    cbar = plt.colorbar(im, ax=ax, label="Grain Count")

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    title = ax.set_title("", fontsize=12, fontweight="bold")

    # Pre-calculate total grains for all frames
    total_grains = grid_states.sum(axis=(1, 2))

    def update(frame):
        """Update function with minimal operations."""
        im.set_array(grid_states[frame])
        title.set_text(f"Iteration {frame + 1} | Total Grains: {total_grains[frame]}")
        return im, title

    print(f"Creating animation with {len(grid_states)} frames...")
    anim = FuncAnimation(
        fig,
        update,
        frames=len(grid_states),
        interval=1000 // fps,
        blit=True,  # Now properly returns artists
        repeat=True,
    )

    print(f"Saving animation to '{output_file}'...")
    writer = FFMpegWriter(fps=fps, bitrate=2000, extra_args=["-preset", "ultrafast"])

    try:
        anim.save(output_file, writer=writer, dpi=100)
        print(f"✓ Animation saved successfully to '{output_file}'")
    except Exception as e:
        print(f"Error saving animation: {e}")
        print("Make sure you have ffmpeg installed")
        sys.exit(1)

    plt.close(fig)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Visualize 2D critical cellular automaton"
    )
    parser.add_argument("--input", "-i", default="output.csv")
    parser.add_argument("--output", "-o", default="animation.mp4")
    parser.add_argument("--fps", type=int, default=30)
    parser.add_argument(
        "--grid-size", type=int, nargs=2, default=[40, 40], metavar=("NX", "NY")
    )

    args = parser.parse_args()
    grid_size = tuple(args.grid_size)

    if not Path(args.input).exists():
        print(f"Error: Input file '{args.input}' not found.")
        sys.exit(1)

    print(f"Reading simulation data from '{args.input}'...")
    grid_states, iterations = parse_output_file(args.input, grid_size=grid_size)

    if len(grid_states) == 0:
        print("Error: No valid data found in input file.")
        sys.exit(1)

    print(f"Parsed {len(grid_states)} simulation frames")
    print(f"Grid size: {grid_size[0]} x {grid_size[1]}")

    create_visualization(
        grid_states, output_file=args.output, fps=args.fps, grid_size=grid_size
    )


if __name__ == "__main__":
    main()
