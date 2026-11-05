# /// script
# requires-python = ">=3.12"
# dependencies = ["build123d", "matplotlib", "meshio", "numpy"]
# ///
"""Plots of fea_hook.py results: section stress maps, per-layer peaks, mass trade-off.

Usage (repo root):
  uv run scripts/fea_plots.py CASE [CASE ...] [--out DIR]
      [--label NAME=TEXT ...] [--maps N] [--failed NAME] [--proposed NAME]
      [--check CASE ...] [--workdir DIR]

CASE is a fea_hook.py case: a name under --workdir (the same default as
fea_hook.py) or a path to the case directory. Writes section_szz.png and
surface_szz.png (the first --maps cases on one colour scale), layer_szz.png
and tradeoff.png, and solver_check.png from the --check cases (which need only
report.json), to --out (default: plots under --workdir). scripts/fea_study.py
calls it with --out docs/img/fea for the figures in docs/fea.md.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
import meshio
import numpy as np
from matplotlib.collections import LineCollection

_spec = importlib.util.spec_from_file_location("fea_hook", Path(__file__).with_name("fea_hook.py"))
fea_hook = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fea_hook)

# VTK 10-node tet (corners 0-3, mid-edge nodes on 01 12 20 03 13 23) split into 8 linear tets
SUBTETS = np.array([[0, 4, 6, 7], [4, 1, 5, 8], [6, 5, 2, 9], [7, 8, 9, 3],
                    [4, 9, 5, 8], [4, 9, 8, 7], [4, 9, 7, 6], [4, 9, 6, 5]])


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cases", nargs="+", type=Path,
                    help="fea_hook.py case names (under --workdir) or case directories")
    ap.add_argument("--out", type=Path, help="output directory for the PNGs (default: plots under --workdir)")
    ap.add_argument("--workdir", type=Path, default=fea_hook.DEFAULT_WORKDIR,
                    help=f"where case names are looked up (default {fea_hook.DEFAULT_WORKDIR})")
    ap.add_argument("--label", action="append", default=[], metavar="NAME=TEXT", help="legend text for a case")
    ap.add_argument("--maps", type=int, default=2, help="cases shown in the section map (default 2)")
    ap.add_argument("--failed", metavar="NAME", help="case to mark with the fracture seen on the failed part")
    ap.add_argument("--proposed", metavar="NAME", help="case to highlight as the proposed design")
    ap.add_argument("--check", action="append", type=Path, default=[], metavar="CASE",
                    help="case whose per-layer profiles from every solver go into solver_check.png (repeatable)")
    args = ap.parse_args()

    def locate(case: Path) -> Path:
        return case if case.is_absolute() or len(case.parts) > 1 else args.workdir / case

    args.cases = [locate(c) for c in args.cases]
    args.check = [locate(c) for c in args.check]
    args.out = args.out or args.workdir / "plots"
    for case in args.cases:
        if not (case / "report.json").is_file():
            ap.error(f"{case}: no report.json (run fea_hook.py first)")
        if not sorted((case / "mesh").glob("results*.vtu")):
            ap.error(f"{case}: no Elmer results VTU under mesh/")
    for case in args.check:
        if not (case / "report.json").is_file():
            ap.error(f"{case}: no report.json (run fea_hook.py first)")
    labels = {}
    for item in args.label:
        name, sep, text = item.partition("=")
        if not sep or not text:
            ap.error(f"--label expects NAME=TEXT, got {item!r}")
        labels[name] = text
    names = [c.name for c in args.cases] + [c.name for c in args.check]
    flagged = {n for n in (args.failed, args.proposed) if n}
    unknown = (set(labels) | flagged) - set(names)
    if unknown:
        ap.error(f"--label/--failed/--proposed name(s) not among the cases: {sorted(unknown)}")
    if not 1 <= args.maps <= len(args.cases):
        ap.error(f"--maps must be between 1 and {len(args.cases)}")
    args.labels = labels
    return args


def load_case(case: Path, label: str) -> dict:
    report = json.loads((case / "report.json").read_text(encoding="utf-8"))
    res = meshio.read(max((case / "mesh").glob("results*.vtu")))
    return {
        "name": case.name,
        "label": label,
        "report": report,
        "ns": fea_hook.load_design(report["overrides"]),
        "points": res.points,
        "tets": res.cells_dict["tetra10"],
        "szz": np.asarray(res.point_data["stress_zz"], dtype=float).ravel(),
        "mass_g": report["mass_g"],
        "peak": max(s for _, s in report["elmer"]["layer_szz_peak"]),
    }


def section(points: np.ndarray, tets10: np.ndarray, values: np.ndarray, x0: float):
    """Cut the mesh with the plane x = x0: (y,z) vertices, triangles, values and outline segments."""
    tets = tets10[:, SUBTETS].reshape(-1, 4)
    d = points[:, 0] - x0
    dt = d[tets]
    cut = tets[(dt.min(axis=1) < 0) & (dt.max(axis=1) > 0)]
    vid: dict[tuple[int, int], int] = {}
    pts, vals, tris = [], [], []

    def vertex(a: int, b: int) -> int:
        key = (a, b) if a < b else (b, a)
        if key not in vid:
            t = d[a] / (d[a] - d[b])
            vid[key] = len(pts)
            pts.append(points[a, 1:] + t * (points[b, 1:] - points[a, 1:]))
            vals.append(values[a] + t * (values[b] - values[a]))
        return vid[key]

    for tet in cut:
        pos = [int(n) for n in tet if d[n] > 0]
        neg = [int(n) for n in tet if d[n] <= 0]
        if len(pos) == 2:
            q = [vertex(pos[0], neg[0]), vertex(pos[0], neg[1]), vertex(pos[1], neg[1]), vertex(pos[1], neg[0])]
            tris += [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
        else:
            lone, others = (pos[0], neg) if len(pos) == 1 else (neg[0], pos)
            tris.append([vertex(lone, o) for o in others])
    pts, tris = np.array(pts), np.array(tris)
    edges = np.sort(np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]]), axis=1)
    uniq, count = np.unique(edges, axis=0, return_counts=True)
    return pts, tris, np.array(vals), pts[uniq[count == 1]]


def plot_sections(cases: list[dict], failed: str | None, out: Path) -> None:
    cuts = [section(c["points"], c["tets"], c["szz"], x0=1e-3) for c in cases]
    vmax = max(np.abs(v).max() for _, _, v, _ in cuts)
    fig, axes = plt.subplots(1, len(cases), figsize=(4.6 * len(cases) + 1.2, 5.0), squeeze=False,
                             layout="constrained")
    for ax, c, (pts, tris, vals, outline) in zip(axes[0], cases, cuts):
        ns = c["ns"]
        im = ax.tripcolor(pts[:, 0], pts[:, 1], tris, vals, shading="gouraud", cmap="RdBu_r",
                          vmin=-vmax, vmax=vmax, rasterized=True)
        ax.add_collection(LineCollection(outline, colors="k", linewidths=0.6))
        rod = plt.Circle((0, ns["rod_cz"]), ns["rod_r"], fill=False, ls="--", lw=0.8, color="0.45")
        ax.add_patch(rod)
        ax.text(0, ns["rod_cz"] + 5, "rod", ha="center", va="center", color="0.45", fontsize=9)
        ax.axhline(0, color="0.3", lw=2.5)
        ax.text(-44, 1.0, "ceiling", fontsize=8, color="0.3", va="bottom")
        ax.axhline(ns["arc_cz"], color="0.35", lw=0.7, ls=":")
        ax.text(44, ns["arc_cz"] + 0.6, "stem / cradle\ntangent", fontsize=7.5, ha="right", va="bottom", color="0.3")
        i = int(np.argmax(vals))
        ax.plot(*pts[i], marker="x", color="k", ms=7, mew=1.5)
        ax.annotate(f"section max {vals[i]:.2f} MPa", pts[i], xytext=(pts[i][0] + 7, pts[i][1] - 9), fontsize=8.5,
                    arrowprops={"arrowstyle": "-", "lw": 0.6})
        if c["name"] == failed:
            y0, y1 = -ns["stem_out"], -ns["r_in"]
            zz = ns["arc_cz"] + 0.9 * np.array([0, 1, -1, 1, -1, 1, -1, 0])
            ax.plot(np.linspace(y0 - 1.5, y1 + 1.5, len(zz)), zz, color="#d62728", lw=2.2)
            ax.annotate("fracture on the\nprinted part", (y0 - 1.5, ns["arc_cz"]), xytext=(-46, ns["arc_cz"] - 14),
                        color="#d62728", fontsize=8.5, arrowprops={"arrowstyle": "->", "color": "#d62728"})
        ax.set_title(f"{c['label']}\n{c['mass_g']:.0f} g PLA, peak layer stress {c['peak']:.2f} MPa", fontsize=10)
        ax.set_xlim(-47, 47)
        ax.set_ylim(min(-60, pts[:, 1].min() - 3), 6)
        ax.set_aspect("equal")
        ax.set_xlabel("y (mm)   wall side  \u2190  \u2192  room side")
        ax.set_ylabel("z (mm)")
    cb = fig.colorbar(im, ax=axes[0].tolist(), shrink=0.8, pad=0.02)
    cb.set_label("\u03c3zz (MPa), + pulls the print layers apart")
    load = cases[0]["report"]["load_N"]
    fig.suptitle(f"Mid-plane section (x = 0), {load:.1f} N ({load / 9.81:.1f} kg) on the cradle, Elmer", fontsize=11)
    fig.savefig(out / "section_szz.png", dpi=150)
    plt.close(fig)


def styles(cases: list[dict], failed: str | None, proposed: str | None) -> dict[str, dict]:
    muted = iter(plt.cm.tab10(np.linspace(0, 1, 10))[[0, 1, 4, 5, 6, 7, 8, 9]])
    out = {}
    for c in cases:
        if c["name"] == failed:
            out[c["name"]] = {"color": "#d62728", "lw": 2.4, "zorder": 3}
        elif c["name"] == proposed:
            out[c["name"]] = {"color": "#2ca02c", "lw": 2.4, "zorder": 3}
        else:
            out[c["name"]] = {"color": next(muted), "lw": 1.1, "zorder": 2}
    return out


def plot_layers(cases: list[dict], failed: str | None, proposed: str | None, out: Path) -> None:
    ns = cases[0]["ns"]
    st = styles(cases, failed, proposed)
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    for c in cases:
        z, s = np.array(c["report"]["elmer"]["layer_szz_peak"]).T
        ax.plot(s, z, **st[c["name"]], label=f"{c['label']}: {c['peak']:.2f} MPa, {c['mass_g']:.0f} g")
    ax.axhline(-ns["plate_t"], color="0.5", lw=0.8, ls="--")
    ax.text(0.02, -ns["plate_t"] - 0.4, "plate underside", fontsize=8, color="0.4", va="top",
            transform=ax.get_yaxis_transform())
    ax.axhline(ns["arc_cz"], color="0.5", lw=0.8, ls=":")
    ax.text(0.02, ns["arc_cz"] - 0.4, "stem / cradle tangent", fontsize=8, color="0.4", va="top",
            transform=ax.get_yaxis_transform())
    ax.axhline(ns["cradle_z"], color="0.5", lw=0.8, ls="--")
    ax.text(0.02, ns["cradle_z"] - 0.4, "cradle bottom", fontsize=8, color="0.4", va="top",
            transform=ax.get_yaxis_transform())
    if failed:
        ax.axhspan(ns["arc_cz"] - 1.5, ns["arc_cz"] + 1.5, color="#d62728", alpha=0.12, lw=0)
        ax.text(0.98, ns["arc_cz"] + 1.7, "fracture on the printed part", fontsize=8, color="#d62728",
                ha="right", transform=ax.get_yaxis_transform())
    ax.set_xlim(left=0)
    ax.set_xlabel("peak \u03c3zz in each 1 mm layer band (MPa)")
    ax.set_ylabel("z (mm), ceiling at 0")
    load = cases[0]["report"]["load_N"]
    ax.set_title(f"Layer-separating stress down the hook, {load / 9.81:.1f} kg on the cradle", fontsize=11)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)
    fig.savefig(out / "layer_szz.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def outer_faces(tets10: np.ndarray) -> np.ndarray:
    """Corner triangles of the faces that belong to exactly one tet (the part's surface)."""
    faces = np.concatenate([tets10[:, [0, 2, 1]], tets10[:, [0, 1, 3]], tets10[:, [1, 2, 3]], tets10[:, [0, 3, 2]]])
    _, inverse, count = np.unique(np.sort(faces, axis=1), axis=0, return_inverse=True, return_counts=True)
    return faces[count[inverse.ravel()] == 1]


def plot_surfaces(cases: list[dict], out: Path) -> None:
    from matplotlib.colors import Normalize
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    elev, azim = -22.0, 58.0
    e, a = np.radians(elev + 25), np.radians(azim - 20)
    light = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
    vmax = max(np.abs(c["szz"]).max() for c in cases)
    norm, cmap = Normalize(-vmax, vmax), plt.get_cmap("RdBu_r")
    fig = plt.figure(figsize=(5.0 * len(cases) + 1.0, 5.2), layout="constrained")
    for k, c in enumerate(cases):
        ax = fig.add_subplot(1, len(cases), k + 1, projection="3d")
        faces = outer_faces(c["tets"])
        tri = c["points"][faces]
        normal = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        normal /= np.linalg.norm(normal, axis=1, keepdims=True)
        rgba = cmap(norm(c["szz"][faces].mean(axis=1)))
        rgba[:, :3] *= (0.55 + 0.45 * np.abs(normal @ light))[:, None]
        ax.add_collection3d(Poly3DCollection(tri, facecolors=rgba, edgecolors=rgba, linewidths=0.25))
        lo, hi = c["points"].min(axis=0), c["points"].max(axis=0)
        ax.set(xlim=(lo[0], hi[0]), ylim=(lo[1], hi[1]), zlim=(lo[2], hi[2]))
        ax.set_box_aspect(hi - lo)
        ax.view_init(elev=elev, azim=azim)
        ax.set_axis_off()
        ax.set_title(f"{c['label']}\npeak layer stress {c['peak']:.2f} MPa", fontsize=10)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=fig.axes, shrink=0.7, pad=0.02)
    cb.set_label("\u03c3zz (MPa), + pulls the print layers apart")
    load = cases[0]["report"]["load_N"]
    fig.suptitle(f"Seen from below on the room side, {load / 9.81:.1f} kg on the cradle, Elmer", fontsize=11)
    fig.savefig(out / "surface_szz.png", dpi=150)
    plt.close(fig)


def plot_check(checks: list[tuple[str, dict]], out: Path) -> None:
    """Per-layer peak sigma_zz from every solver result in the given reports: solver and mesh agreement."""
    solver_name = {"elmer": "Elmer", "ccx": "CalculiX"}
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    colors = iter(plt.cm.tab10(np.linspace(0, 1, 10)))
    for label, report in checks:
        color = next(colors)
        for solver, ls in (("elmer", "-"), ("ccx", "--")):
            if solver not in report:
                continue
            z, s = np.array(report[solver]["layer_szz_peak"]).T
            ax.plot(s, z, ls=ls, lw=1.4, color=color,
                    label=f"{label}, {solver_name[solver]}, h {report['mesh']['h_hook']:g} mm: {s.max():.3f} MPa")
    ax.set_xlim(left=0)
    ax.set_xlabel("peak \u03c3zz in each 1 mm layer band (MPa)")
    ax.set_ylabel("z (mm), ceiling at 0")
    load = checks[0][1]["load_N"]
    ax.set_title(f"Solver and mesh check, {load / 9.81:.1f} kg on the cradle\n"
                 "solid: Elmer, dashed: CalculiX; h = element size in the hook", fontsize=11)
    ax.legend(fontsize=7.5, loc="lower right")
    ax.grid(alpha=0.3)
    fig.savefig(out / "solver_check.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_tradeoff(cases: list[dict], failed: str | None, proposed: str | None, out: Path) -> None:
    st = styles(cases, failed, proposed)
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for c in sorted(cases, key=lambda c: c["mass_g"]):
        ax.plot(c["mass_g"], c["peak"], "o", color=st[c["name"]]["color"], ms=7,
                label=f"{c['label']}: {c['peak']:.2f} MPa, {c['mass_g']:.0f} g")
    ax.legend(fontsize=7.5, loc="upper right")
    ax.set_xlabel("mass in PLA, solid (g)")
    ax.set_ylabel("peak \u03c3zz across layers (MPa)")
    load = cases[0]["report"]["load_N"]
    ax.set_title(f"Design variants: peak layer stress vs mass, {load / 9.81:.1f} kg", fontsize=11)
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.3)
    fig.savefig(out / "tradeoff.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    plt.switch_backend("Agg")
    args.out.mkdir(parents=True, exist_ok=True)
    cases = [load_case(c, args.labels.get(c.name, c.name)) for c in args.cases]
    plot_sections(cases[: args.maps], args.failed, args.out)
    plot_surfaces(cases[: args.maps], args.out)
    plot_layers(cases, args.failed, args.proposed, args.out)
    plot_tradeoff(cases, args.failed, args.proposed, args.out)
    written = ["section_szz.png", "surface_szz.png", "layer_szz.png", "tradeoff.png"]
    if args.check:
        checks = [(args.labels.get(c.name, c.name), json.loads((c / "report.json").read_text(encoding="utf-8")))
                  for c in args.check]
        plot_check(checks, args.out)
        written.append("solver_check.png")
    for c in cases:
        print(f"{c['name']:12s} {c['mass_g']:6.1f} g  peak layer stress {c['peak']:.3f} MPa")
    print(f"wrote {', '.join(written)} to {args.out}")


if __name__ == "__main__":
    main()
