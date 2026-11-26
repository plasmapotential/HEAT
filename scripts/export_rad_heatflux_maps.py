#!/usr/bin/env python3
"""Export HEAT radiated heat-flux point clouds and reference-style grids.

Geometry and the $MW/m^2$ cell values come straight from each run's
0.000000000/paraview/HF_allSources_all_mesh.vtp. That mesh pairs every
triangle with its own flux value, so it is the source of truth; matching
triangle order between the cached STLs and HF_rad_all.csv is not reliable
(it previously undercounted total radiated power by roughly 15%).

Each run gets:
- a text file with columns x, y, z, heat flux (one row per mesh triangle)
- a phi_deg/Z_m/R_m/q_MW_m2/area_m2 grid text file, matching the format and
  1-degree x 1-cm binning of the reference maps in
  DIAG_HEATmap/sparc_000001_..._CPU_v10_more_cores/0.000000000/paraview/
- a PNG of that grid, in the same two-panel layout as the reference maps
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

VTP_MESH_NAME = "HF_allSources_all_mesh.vtp"
# File-name labels matching the reference grids in DIAG_HEATmap/.../paraview
# (Uniform_HF_rad_photonic_load_map.txt, TRANSP_..., PreTQ_...).
GRID_FILE_LABELS = {
    "preTQ": "PreTQ",
    "uniform": "Uniform",
    "transp": "TRANSP",
}
# Bin widths used by the reference phi/Z grid text files.
GRID_PHI_WIDTH_DEG = 1.0
GRID_Z_WIDTH_M = 0.01
NO_SURFACE = "#bdbdbd"
# Uniform and TRANSP stay on the reference figure's 0–0.012 scale. PreTQ
# peaks near 0.041 MW/m², so that cap paints most of its power as one red.
COLOR_MAX = {
    "PreTQ": 0.05,
    "Uniform": 0.012,
    "TRANSP": 0.012,
}
CASE_TITLES = {
    "PreTQ": "Pre-TQ photonic load",
    "Uniform": "Uniform photonic load",
    "TRANSP": "TRANSP photonic load",
}


def load_vtp_mesh(vtp_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Triangle vertices (mm, n×3×3) and per-triangle $MW/m^2$ cell values."""
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(vtp_path))
    reader.Update()
    polydata = reader.GetOutput()

    points_mm = vtk_to_numpy(polydata.GetPoints().GetData()).astype(np.float64)
    connectivity = vtk_to_numpy(polydata.GetPolys().GetData()).reshape(-1, 4)
    if not np.all(connectivity[:, 0] == 3):
        raise ValueError(f"Expected an all-triangle mesh in {vtp_path}")
    triangle_indices = connectivity[:, 1:4]
    vertices_mm = points_mm[triangle_indices]

    flux = vtk_to_numpy(polydata.GetCellData().GetArray(0))
    if flux.shape[0] != vertices_mm.shape[0]:
        raise ValueError(
            f"Cell count {vertices_mm.shape[0]} does not match flux values "
            f"{flux.shape[0]} in {vtp_path}"
        )
    return vertices_mm, flux.astype(np.float64)


def triangle_areas_m2(vertices_mm: np.ndarray) -> np.ndarray:
    edge_a = vertices_mm[:, 1] - vertices_mm[:, 0]
    edge_b = vertices_mm[:, 2] - vertices_mm[:, 0]
    return 0.5 * np.linalg.norm(np.cross(edge_a, edge_b), axis=1) * 1e-6


def surface_samples(vertices_mm: np.ndarray, area: np.ndarray, flux: np.ndarray):
    """Vertices and edge midpoints, so a face paints every bin it covers.

    Centroid-only binning leaves empty rows between the ~10 mm faces on the
    top and bottom of the vessel. Those rows are wall in the VTP mesh.
    """
    edge_midpoints = np.stack([
        0.5 * (vertices_mm[:, 0] + vertices_mm[:, 1]),
        0.5 * (vertices_mm[:, 1] + vertices_mm[:, 2]),
        0.5 * (vertices_mm[:, 2] + vertices_mm[:, 0]),
    ], axis=1)
    points = np.concatenate([vertices_mm, edge_midpoints], axis=1).reshape(-1, 3)
    samples_per_face = points.shape[0] // flux.shape[0]
    share = np.repeat(area / samples_per_face, samples_per_face)
    values = np.repeat(flux, samples_per_face)
    return points, share, values


def grid_file_label(run_name: str) -> str:
    for key, label in GRID_FILE_LABELS.items():
        if f"_{key}_" in f"_{run_name}_":
            return label
    return run_name


def export_reference_grid(
    run_dir: Path,
    points_mm: np.ndarray,
    sample_area: np.ndarray,
    sample_flux: np.ndarray,
    name: str,
) -> None:
    """Write a phi_deg/Z_m/R_m/q_MW_m2/area_m2 grid like the reference maps in
    DIAG_HEATmap/sparc_000001_..._CPU_v10_more_cores/0.000000000/paraview/.

    That reference used 1-degree x 1-cm bins and dropped bins with no wall
    surface, rather than interpolating. This mirrors both choices.
    """
    radius_m = np.hypot(points_mm[:, 0], points_mm[:, 1]) / 1000.0
    height_m = points_mm[:, 2] / 1000.0
    phi_deg = np.degrees(np.arctan2(points_mm[:, 1], points_mm[:, 0]))

    phi_edges = np.arange(-180.0, 180.0 + GRID_PHI_WIDTH_DEG, GRID_PHI_WIDTH_DEG)
    z_min = np.floor(height_m.min() / GRID_Z_WIDTH_M) * GRID_Z_WIDTH_M
    z_max = np.ceil(height_m.max() / GRID_Z_WIDTH_M) * GRID_Z_WIDTH_M
    z_edges = np.arange(z_min, z_max + GRID_Z_WIDTH_M, GRID_Z_WIDTH_M)

    weight = np.histogram2d(phi_deg, height_m, bins=[phi_edges, z_edges], weights=sample_area)[0]
    weighted_flux = np.histogram2d(
        phi_deg, height_m, bins=[phi_edges, z_edges], weights=sample_area * sample_flux
    )[0]
    weighted_radius = np.histogram2d(
        phi_deg, height_m, bins=[phi_edges, z_edges], weights=sample_area * radius_m
    )[0]

    phi_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    phi_grid, z_grid = np.meshgrid(phi_centers, z_centers, indexing="ij")

    filled = weight > 0
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_flux = weighted_flux / weight
        mean_radius = weighted_radius / weight

    rows = np.column_stack([
        phi_grid[filled],
        z_grid[filled],
        mean_radius[filled],
        mean_flux[filled],
        weight[filled],
    ])

    output_dir = run_dir / "0.000000000" / "paraview"
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / f"{grid_file_label(name)}_HF_rad_photonic_load_map.txt"
    with output_path.open("w") as grid_file:
        grid_file.write("phi_deg Z_m R_m q_MW_m2 area_m2\n")
        for phi, z, r, q, a in rows:
            grid_file.write(f"{phi:.6f} {z:.6f} {r:.6f} {q:.8e} {a:.8e}\n")
    print(f"Wrote {output_path}  rows {len(rows)} of {weight.size} grid cells")
    plot_reference_grid(output_path)


def plot_reference_grid(grid_path: Path) -> None:
    """Draw the two-panel map used for the reference photonic-load figure."""
    grid = np.loadtxt(grid_path, skiprows=1)
    phi, height, radius, flux, area = grid.T
    label = grid_path.name.split("_HF_rad")[0]
    title = CASE_TITLES.get(label, label)
    total_power = float(np.sum(flux * area))
    peak_flux = float(np.max(flux))

    phi_edges = np.arange(-180.0, 180.0 + GRID_PHI_WIDTH_DEG, GRID_PHI_WIDTH_DEG)
    z_min = np.floor(height.min() / GRID_Z_WIDTH_M) * GRID_Z_WIDTH_M
    z_max = np.ceil(height.max() / GRID_Z_WIDTH_M) * GRID_Z_WIDTH_M
    z_edges = np.arange(z_min, z_max + GRID_Z_WIDTH_M, GRID_Z_WIDTH_M)
    phi_index = np.clip(np.searchsorted(phi_edges, phi, side="right") - 1, 0, len(phi_edges) - 2)
    z_index = np.clip(np.searchsorted(z_edges, height, side="right") - 1, 0, len(z_edges) - 2)
    unfolded = np.full((len(phi_edges) - 1, len(z_edges) - 1), np.nan)
    unfolded[phi_index, z_index] = flux

    radius_edges = np.linspace(1.75, 2.70, 381)
    radius_weight = np.histogram2d(radius, height, bins=[radius_edges, z_edges], weights=area)[0]
    radius_flux = np.histogram2d(
        radius, height, bins=[radius_edges, z_edges], weights=area * flux
    )[0]
    with np.errstate(invalid="ignore", divide="ignore"):
        poloidal = radius_flux / radius_weight
    poloidal[radius_weight <= 0] = np.nan
    radius_centers = 0.5 * (radius_edges[:-1] + radius_edges[1:])
    height_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    radius_grid, height_grid = np.meshgrid(radius_centers, height_centers, indexing="ij")
    on_wall = np.isfinite(poloidal)

    colormap = plt.get_cmap("turbo").copy()
    colormap.set_bad(NO_SURFACE)
    color_max = COLOR_MAX.get(label, 0.012)
    norm = plt.Normalize(vmin=0.0, vmax=color_max)
    # Equal aspect would make the poloidal section a thin strip. Draw it wider
    # than that, still close to the unfolded map.
    figure = plt.figure(figsize=(16.8, 7.4))
    z_span = 3.30
    radius_span = 0.95
    axes_bottom = 0.12
    axes_height = 0.70
    equal_width = axes_height * (figure.get_figheight() / figure.get_figwidth()) * (
        radius_span / z_span
    )
    right_width = equal_width * 1.85
    left_margin = 0.055
    panel_gap = 0.012
    left_width = 0.66
    axes = [
        figure.add_axes([left_margin, axes_bottom, left_width, axes_height]),
        figure.add_axes(
            [left_margin + left_width + panel_gap, axes_bottom, right_width, axes_height],
            sharey=None,
        ),
    ]
    figure.suptitle(
        f"{title}   ·   total {total_power:.3f} MW   ·   peak {peak_flux:.4f} MW/m²\n"
        "vessel shells 2453117_objects / 2453117_objects001",
        fontsize=13,
    )
    image = axes[0].pcolormesh(
        phi_edges,
        z_edges,
        np.ma.masked_invalid(unfolded.T),
        cmap=colormap,
        norm=norm,
        shading="flat",
        rasterized=True,
    )
    axes[0].set(
        title="Unfolded vessel wall",
        xlabel="Toroidal angle φ (degrees)",
        ylabel="Height Z (m)",
        xlim=(-180, 180),
        ylim=(-1.65, 1.65),
        facecolor=NO_SURFACE,
    )
    axes[0].legend(
        handles=[Patch(facecolor=NO_SURFACE, edgecolor="#4d4d4d", label="no wall surface")],
        loc="upper right",
        framealpha=0.92,
        fontsize=9,
    )
    axes[1].scatter(
        radius_grid[on_wall],
        height_grid[on_wall],
        c=poloidal[on_wall],
        s=36,
        marker="s",
        linewidths=0,
        cmap=colormap,
        norm=norm,
        rasterized=True,
    )
    axes[1].set(
        title="Poloidal section, mean over φ",
        xlabel="Major radius R (m)",
        xlim=(1.75, 2.70),
        ylim=(-1.65, 1.65),
        facecolor=NO_SURFACE,
    )
    axes[1].tick_params(labelleft=False)
    colorbar = figure.colorbar(
        image,
        cax=figure.add_axes([
            left_margin + left_width + panel_gap + right_width + 0.012,
            axes_bottom,
            0.012,
            axes_height,
        ]),
    )
    colorbar.set_ticks(np.linspace(0.0, color_max, 6 if color_max >= 0.05 else 7))
    colorbar.set_label("Area-weighted radiative load q (MW/m²)")

    map_output = grid_path.with_suffix(".png")
    figure.savefig(map_output, dpi=160)
    plt.close(figure)
    print(f"Wrote {map_output}")


def export_run(run_dir: Path) -> None:
    vtp_path = run_dir / "0.000000000" / "paraview" / VTP_MESH_NAME
    if not vtp_path.is_file():
        raise FileNotFoundError(f"Missing mesh: {vtp_path}")

    output_dir = run_dir / "heatflux_maps"
    output_dir.mkdir(exist_ok=True)
    name = run_dir.name.removeprefix("sparc_000001_")

    print(f"Reading {vtp_path}")
    vertices_mm, flux = load_vtp_mesh(vtp_path)
    area = triangle_areas_m2(vertices_mm)

    text_output = output_dir / f"{name}_HF_rad_xyz_heatflux.txt"
    if not text_output.is_file():
        centroids_mm = vertices_mm.mean(axis=1)
        np.savetxt(
            text_output,
            np.column_stack([centroids_mm, flux]),
            header="x_mm y_mm z_mm heat_flux_MW_m2",
            comments="# ",
            fmt="%.10e",
        )
        print(f"Wrote {text_output}")

    points_mm, sample_area, sample_flux = surface_samples(vertices_mm, area, flux)
    export_reference_grid(run_dir, points_mm, sample_area, sample_flux, name)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Export x/y/z/heat-flux text data and reference-style phi/Z grids "
            "from HEAT runs."
        )
    )
    parser.add_argument("run_dirs", type=Path, nargs="+", help="HEAT run directories")
    args = parser.parse_args()
    for run_dir in args.run_dirs:
        export_run(run_dir)


if __name__ == "__main__":
    main()
