#!/usr/bin/env python3
"""Visualization for HDF5 simulation output"""

from curses import meta

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation


class SimulationData:
    """Class to handle HDF5 simulation data"""

    def __init__(self, filename):
        self.filename = filename
        self.file = h5py.File(filename, "r")

        # Load metadata
        meta = self.file["metadata"]
        self.nx = meta["nx"][()]
        self.nz = meta["nz"][()]
        self.dx = meta["dx"][()]
        self.dz = meta["dz"][()]
        self.dt = meta["dt"][()]
        self.nt = meta["nt"][()]
        self.save_interval = meta["save_interval"][()]

        # Load source info
        self.source_x = meta["source_x"][:]  # Changed from source_i
        self.source_z = meta["source_z"][:]  # Changed from source_k
        self.source_freq = meta["source_frequency"][
            :
        ]  # Changed from source_f0
        self.source_amp = meta["source_amplitude"][:]  # New field

        self.source_f0 = meta["source_f0"][:]
        self.source_t0 = meta["source_t0"][:]

        # Get available fields
        self.fields = list(self.file["fields"].keys())

        # Get number of frames
        self.n_frames = len(
            self.file[f"fields/{self.fields[0]}"].keys(),
        )

        print(f"Loaded simulation data from {filename}")
        print(f"  Grid: {self.nx} x {self.nz}")
        print(f"  Spacing: dx={self.dx}m, dz={self.dz}m")
        print(f"  Time step: dt={self.dt}s")
        print(f"  Total steps: {self.nt}")
        print(f"  Frames: {self.n_frames}")
        print(f"  Available fields: {', '.join(self.fields)}")
        print(f"  Sources: {len(self.source_i)}")

    def get_frame(self, field, frame_number):
        """Get a specific frame for a field"""
        dataset_name = f"frame_{frame_number:06d}"
        data = self.file[f"fields/{field}/{dataset_name}"][:]

        # Reshape from flat to 2D
        return data.reshape(self.nz, self.nx)

    def get_frame_time(self, field, frame_number):
        """Get the time of a specific frame"""
        dataset_name = f"frame_{frame_number:06d}"
        return self.file[f"fields/{field}/{dataset_name}"].attrs[
            "time"
        ]

    def get_frame_timestep(self, field, frame_number):
        """Get the timestep of a specific frame"""
        dataset_name = f"frame_{frame_number:06d}"
        return self.file[f"fields/{field}/{dataset_name}"].attrs[
            "timestep"
        ]

    def get_materials(self):
        """Get material properties"""
        materials = self.file["materials"]
        return {
            "vp": materials["vp"][:].reshape(self.nz, self.nx),
            "vs": materials["vs"][:].reshape(self.nz, self.nx),
            "rho": materials["rho"][:].reshape(self.nz, self.nx),
            "lambda": materials["lambda"][:].reshape(
                self.nz,
                self.nx,
            ),
            "mu": materials["mu"][:].reshape(self.nz, self.nx),
        }

    def close(self):
        """Close the HDF5 file"""
        self.file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


def visualize_field(
    data,
    field="vmag",
    save_video=False,
    output="animation.mp4",
):
    """Create animation of a field"""
    print(f"Creating animation for field: {field}")

    # Determine color scale from all frames
    print("Determining color scale...")
    max_vals = []
    for frame in range(min(20, data.n_frames)):
        field_data = data.get_frame(field, frame)
        max_vals.append(np.abs(field_data).max())
    vmax = np.percentile(max_vals, 95)
    print(f"Color scale: ±{vmax:.2e}")

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 10))

    # Initial plot
    first_frame = data.get_frame(field, 0)
    im = ax.imshow(
        first_frame,
        cmap="seismic",
        aspect="auto",
        vmin=-vmax,
        vmax=vmax,
        interpolation="bilinear",
        extent=[0, data.nx * data.dx, data.nz * data.dz, 0],
    )

    cbar = plt.colorbar(im, ax=ax, label=f"{field} (m/s)")
    ax.set_xlabel("X (m)", fontsize=12)
    ax.set_ylabel("Z (m)", fontsize=12)

    # Plot sources
    for x, z in zip(data.source_x, data.source_z):
        ax.plot(
            x * data.dx,
            z * data.dz,
            "r*",
            markersize=15,
            markeredgecolor="white",
            markeredgewidth=1,
        )

    title = ax.set_title(f"{field} - t=0.000s (step 0)", fontsize=14)

    def update(frame):
        """Update function for animation"""
        field_data = data.get_frame(field, frame)
        im.set_data(field_data)

        time = data.get_frame_time(field, frame)
        timestep = data.get_frame_timestep(field, frame)
        title.set_text(f"{field} - t={time:.3f}s (step {timestep})")

        if frame % 10 == 0:
            print(f"Processing frame {frame}/{data.n_frames}")

        return [im, title]

    # Create animation
    print("Creating animation...")
    anim = FuncAnimation(
        fig,
        update,
        frames=data.n_frames,
        interval=50,
        blit=True,
        repeat=True,
    )

    if save_video:
        print(f"Saving video to {output}...")
        writer = FFMpegWriter(fps=30, bitrate=5000)
        anim.save(output, writer=writer, dpi=150)
        print("Video saved!")
    else:
        plt.tight_layout()
        plt.show()

    return anim


def plot_comparison(
    data,
    fields=["vmag", "divergence", "curl"],
    frame=50,
):
    """Plot multiple fields side by side"""
    n_fields = len(fields)
    fig, axes = plt.subplots(1, n_fields, figsize=(6 * n_fields, 5))

    if n_fields == 1:
        axes = [axes]

    time = data.get_frame_time(fields[0], frame)
    timestep = data.get_frame_timestep(fields[0], frame)

    for ax, field in zip(axes, fields):
        field_data = data.get_frame(field, frame)
        vmax = np.abs(field_data).max()

        im = ax.imshow(
            field_data,
            cmap="seismic",
            aspect="auto",
            vmin=-vmax,
            vmax=vmax,
            extent=[0, data.nx * data.dx, data.nz * data.dz, 0],
        )
        ax.set_title(f"{field}")
        ax.set_xlabel("X (m)")
        ax.set_ylabel("Z (m)")
        plt.colorbar(im, ax=ax)

        # Plot sources
        for i, k in zip(data.source_i, data.source_k):
            ax.plot(i * data.dx, k * data.dz, "r*", markersize=10)

    fig.suptitle(f"t={time:.3f}s (step {timestep})", fontsize=14)
    plt.tight_layout()
    plt.savefig(f"comparison_frame_{frame:06d}.png", dpi=200)
    print(f"Saved comparison_frame_{frame:06d}.png")
    plt.show()


def plot_materials(data):
    """Plot material properties"""
    materials = data.get_materials()

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Vp
    im1 = axes[0].imshow(
        materials["vp"],
        aspect="auto",
        cmap="viridis",
        extent=[0, data.nx * data.dx, data.nz * data.dz, 0],
    )
    axes[0].set_title("P-wave velocity (Vp)")
    axes[0].set_xlabel("X (m)")
    axes[0].set_ylabel("Z (m)")
    plt.colorbar(im1, ax=axes[0], label="m/s")

    # Vs
    im2 = axes[1].imshow(
        materials["vs"],
        aspect="auto",
        cmap="viridis",
        extent=[0, data.nx * data.dx, data.nz * data.dz, 0],
    )
    axes[1].set_title("S-wave velocity (Vs)")
    axes[1].set_xlabel("X (m)")
    plt.colorbar(im2, ax=axes[1], label="m/s")

    # Density
    im3 = axes[2].imshow(
        materials["rho"],
        aspect="auto",
        cmap="viridis",
        extent=[0, data.nx * data.dx, data.nz * data.dz, 0],
    )
    axes[2].set_title("Density (ρ)")
    axes[2].set_xlabel("X (m)")
    plt.colorbar(im3, ax=axes[2], label="kg/m³")

    plt.tight_layout()
    plt.savefig("materials.png", dpi=200)
    print("Saved materials.png")
    plt.show()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Visualize HDF5 simulation data",
    )
    parser.add_argument(
        "--file",
        default="simulation_output.h5",
        help="HDF5 file",
    )
    parser.add_argument(
        "--field",
        default="vmag",
        help="Field to visualize",
    )
    parser.add_argument(
        "--video",
        action="store_true",
        help="Save video",
    )
    parser.add_argument(
        "--output",
        default="animation.mp4",
        help="Output video file",
    )
    parser.add_argument(
        "--compare",
        type=int,
        help="Create comparison plot at frame",
    )
    parser.add_argument(
        "--materials",
        action="store_true",
        help="Plot material properties",
    )

    args = parser.parse_args()

    with SimulationData(args.file) as data:
        if args.materials:
            plot_materials(data)
        elif args.compare is not None:
            plot_comparison(data, frame=args.compare)
        else:
            visualize_field(
                data,
                args.field,
                save_video=args.video,
                output=args.output,
            )
