# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy"]
# ///
"""Run the FEA study in docs/fea.md and redraw docs/img/fea.

Runs scripts/fea_hook.py for the design variants (Elmer, 1.2 mm elements in
the hook) and the solver, support, load and mesh checks, then
scripts/fea_plots.py, and prints the numbers docs/fea.md quotes. Cases go to
fea_hook.py's work directory. The 1.0 mm CalculiX check needs about 9 GB of
memory and runs alone at the end; the whole study takes about 40 minutes.

Usage (repo root):
  uv run scripts/fea_study.py --ccx "C:/Program Files/FreeCAD 1.1/bin/ccx.exe" [--jobs 3]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("fea_hook", SCRIPTS / "fea_hook.py")
fea_hook = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fea_hook)

AS_PRINTED = {"stem_t": 8.0, "hook_w": 24.0, "stem_fillet_r": 5.0, "plate_y_stem_side": 28.0, "plate_y_open_side": 28.0}


def plate(half_width: float) -> dict[str, float]:
    return {"plate_y_stem_side": half_width, "plate_y_open_side": half_width}


# (case name, legend label, design parameter overrides, fea_hook.py arguments)
VARIANTS = [
    ("as_printed", "as printed: stem 8, hook 24, R5", AS_PRINTED),
    ("current", "current: stem 18, hook 28, R8", {}),
    ("s12_w28_r5", "stem 12, hook 28, R5", {"stem_t": 12.0, "stem_fillet_r": 5.0, **plate(33.0)}),
    ("s14_w28_r5", "stem 14, hook 28, R5", {"stem_t": 14.0, "stem_fillet_r": 5.0, **plate(35.0)}),
    ("s16_w28_r5", "stem 16, hook 28, R5", {"stem_t": 16.0, "stem_fillet_r": 5.0, **plate(37.0)}),
    ("s14_w28_r8", "stem 14, hook 28, R8", {"stem_t": 14.0, **plate(37.0)}),
    ("s16_w28_r8", "stem 16, hook 28, R8", {"stem_t": 16.0, **plate(39.0)}),
    ("s16_w30_r8", "stem 16, hook 30, R8", {"stem_t": 16.0, "hook_w": 30.0, **plate(39.0)}),
    ("s18_w30_r8", "stem 18, hook 30, R8", {"hook_w": 30.0}),
]
CASES = [(name, label, ov, ["--h-hook", "1.2"]) for name, label, ov in VARIANTS] + [
    ("as_printed_h2", "as printed", AS_PRINTED, ["--h-hook", "2.0", "--solver", "both"]),
    ("as_printed_h2_bonded", "as printed, ceiling held", AS_PRINTED, ["--h-hook", "2.0", "--ceiling", "bonded"]),
    ("current_h2", "current", {}, ["--h-hook", "2.0", "--solver", "both"]),
    ("current_h2_bonded", "current, ceiling held", {}, ["--h-hook", "2.0", "--ceiling", "bonded"]),
    ("current_h2_c15", "current, contact +-15 deg", {}, ["--h-hook", "2.0", "--contact-deg", "15"]),
    ("current_h1.4", "current", {}, ["--h-hook", "1.4", "--solver", "ccx"]),
]
HEAVY = [("current_h1.0", "current", {}, ["--h-hook", "1.0", "--solver", "ccx"])]
SOLVER_CHECK = ["as_printed_h2", "current", "current_h2", "current_h1.4", "current_h1.0"]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ccx", type=Path, required=True, help="CalculiX ccx executable")
    ap.add_argument("--jobs", type=int, default=3, help="cases run at once (default 3)")
    ap.add_argument("--workdir", type=Path, default=fea_hook.DEFAULT_WORKDIR,
                    help=f"parent of the case directories (default {fea_hook.DEFAULT_WORKDIR})")
    ap.add_argument("--out", type=Path, default=SCRIPTS.parent / "docs" / "img" / "fea",
                    help="where the figures go (default docs/img/fea)")
    args = ap.parse_args()
    if not args.ccx.is_file():
        ap.error(f"--ccx: {args.ccx} is not a file")
    if args.jobs < 1:
        ap.error("--jobs must be at least 1")
    return args


def run_case(case: tuple, args: argparse.Namespace) -> str | None:
    """Run one fea_hook.py case; return the tail of its output if it failed."""
    name, _, overrides, extra = case
    cmd = ["uv", "run", str(SCRIPTS / "fea_hook.py"), "--name", name, "--workdir", str(args.workdir),
           "--ccx", str(args.ccx), *extra]
    for key, value in overrides.items():
        cmd += ["--set", f"{key}={value!r}"]
    p = subprocess.run(cmd, cwd=SCRIPTS.parent, capture_output=True, text=True, encoding="utf-8", errors="replace",
                       check=False)
    print(f"{name}: {'done' if p.returncode == 0 else 'FAILED'}", flush=True)
    return None if p.returncode == 0 else (p.stdout + p.stderr)[-2000:]


def summarize(args: argparse.Namespace) -> None:
    reports = {name: json.loads((args.workdir / name / "report.json").read_text(encoding="utf-8"))
               for name, *_ in CASES + HEAVY}
    print("\nDesign variants, Elmer, 1.2 mm (peak layer stress at the stem root and where it turns into the cradle)")
    print("| case | mass (g) | root (MPa) | shoulder (MPa) | deflection (mm) | plate s1 (MPa) | elements |")
    print("| --- | --- | --- | --- | --- | --- | --- |")
    for name, label, _ in VARIANTS:
        r = reports[name]
        e = r["elmer"]
        print(f"| {label} | {r['mass_g']:.1f} | {e['layer_szz_root']['MPa']:.2f} | {e['layer_szz_shoulder']['MPa']:.2f}"
              f" | {e['max_deflection_mm']:.2f} | {e['plate_s1']['MPa']:.2f} | {r['mesh']['tets']} |")
    print("\nChecks (peak layer stress, root / shoulder, and the support reaction)")
    for name, label, _, _ in CASES[len(VARIANTS):] + HEAVY + [CASES[1]]:
        r = reports[name]
        for solver in ("elmer", "ccx"):
            if solver in r:
                s = r[solver]
                print(f"  {name:22s} {solver:5s} h {r['mesh']['h_hook']:<3} {s['layer_szz_root']['MPa']:.3f} /"
                      f" {s['layer_szz_shoulder']['MPa']:.3f} MPa, reaction {s['reaction_N'][2]:.2f} of {r['load_N']:.2f} N")
        if "elmer_vs_ccx_szz_hook" in r:
            print(f"  {'':22s} Elmer - CalculiX sigma_zz in the hook: RMS {100 * r['elmer_vs_ccx_szz_hook']['rms_diff_rel_peak']:.1f} %"
                  f" of the peak, max {r['elmer_vs_ccx_szz_hook']['max_abs_diff_MPa']:.2f} MPa")


def main() -> None:
    args = parse_args()
    with ThreadPoolExecutor(args.jobs) as pool:
        errors = dict(zip([c[0] for c in CASES], pool.map(lambda c: run_case(c, args), CASES), strict=True))
    errors.update({c[0]: run_case(c, args) for c in HEAVY})
    failed = {name: err for name, err in errors.items() if err}
    for name, err in failed.items():
        print(f"\n--- {name} ---\n{err}")
    if failed:
        raise SystemExit(f"{len(failed)} case(s) failed: {', '.join(failed)}")

    labels = {name: label for name, label, *_ in CASES + HEAVY}
    plot = ["uv", "run", str(SCRIPTS / "fea_plots.py"), "--out", str(args.out), "--workdir", str(args.workdir),
            "--maps", "2", "--failed", "as_printed", "--proposed", "current"]
    for name in dict.fromkeys([n for n, *_ in VARIANTS] + SOLVER_CHECK):
        plot += ["--label", f"{name}={labels[name]}"]
    plot += [n for n, *_ in VARIANTS]
    for name in SOLVER_CHECK:
        plot += ["--check", name]
    subprocess.run(plot, cwd=SCRIPTS.parent, check=True)
    summarize(args)


if __name__ == "__main__":
    main()
