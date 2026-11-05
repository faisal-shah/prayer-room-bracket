# /// script
# requires-python = ">=3.12"
# dependencies = ["build123d", "gmsh", "meshio", "numpy", "scipy"]
# ///
"""Linear-static FEA of the ceiling rod hook with Elmer, cross-checked by CalculiX.

The geometry is built by executing scripts/ceiling_rod_hook.py, optionally with
design parameters overridden (--set NAME=VALUE), so the analysis always runs on
the model that is exported for printing.

Model
- Material: isotropic PLA, E = 2.3 GPa (printed, Prusament TDS), nu = 0.35.
  Stresses do not depend on E (one material, force loads, fixed supports);
  deflections scale with 1/E.
- Supports: the counterbore floors under the pan heads are fixed. With
  --ceiling bonded the plate's ceiling face is also held in z, which brackets
  the real, compression-only ceiling contact from the stiff side.
- Load: the hook's share of the rod and curtain weight as a radial bearing
  pressure p0 * cos(90 deg * phi / phi0) on the cradle bottom, |phi| <= phi0.
  Its resultant acts through the rod axis, so stresses scale linearly with it.
- Layer peel stress: printed with the ceiling face on the bed, the layers are
  z-planes and sigma_zz pulls them apart. The report lists the peak sigma_zz
  per 1 mm layer band down the hook (layer_szz_peak) and its two maxima,
  above and below mid-stem: at the root fillet (layer_szz_root) and where the
  stem turns into the cradle (layer_szz_shoulder).

Mesh: Gmsh, 10-node tets, --h-hook in the hook, --h-plate elsewhere. Elmer
(StressSolver) and CalculiX (C3D10) solve the same mesh and loads; with
--solver both the report includes their difference. scripts/fea_study.py runs
the cases that docs/fea.md cites.

Usage (repo root):
  uv run scripts/fea_hook.py --name current
  uv run scripts/fea_hook.py --name check --h-hook 2.0 --solver both --ccx PATH/TO/ccx.exe
  uv run scripts/fea_hook.py --name wide --set hook_w=30.0
  uv run scripts/fea_hook.py --name as_printed --set stem_t=8.0 --set hook_w=24.0 \\
      --set plate_y_stem_side=28.0 --set plate_y_open_side=28.0 --set stem_fillet_r=5.0
"""
from __future__ import annotations

import argparse
import ast
import itertools
import json
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
DESIGN = REPO / "scripts" / "ceiling_rod_hook.py"
E_PLA = 2300.0      # MPa (mm-N-MPa units throughout); printed PLA, Prusament TDS (ISO 527)
NU_PLA = 0.35
PLA_G_PER_MM3 = 1.24e-3
G = 9.81
DEFAULT_WORKDIR = Path(tempfile.gettempdir()) / "prayer-room-bracket-fea"
SEAT, CRADLE, CEILING = 1, 2, 3  # physical surface tags = Elmer boundary ids

# Dunavant degree-5 rule on the reference triangle: area coordinates, weights sum to 1
_A, _B, _C, _D = 0.059715871789770, 0.470142064105115, 0.797426985353087, 0.101286507323456
TRI_L = np.array([[1 / 3, 1 / 3, 1 / 3], [_A, _B, _B], [_B, _A, _B], [_B, _B, _A],
                  [_C, _D, _D], [_D, _C, _D], [_D, _D, _C]])
TRI_W = np.array([0.225] + [0.132394152788506] * 3 + [0.125939180544827] * 3)
GMSH_TO_C3D10 = [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--name", default="case", help="case name (sub-directory of --workdir)")
    ap.add_argument("--set", dest="overrides", action="append", default=[], metavar="NAME=VALUE",
                    help="override a design parameter of ceiling_rod_hook.py (Python literal)")
    ap.add_argument("--load-kg", type=float, default=4.0, help="mass carried by the hook, kg (default 4)")
    ap.add_argument("--contact-deg", type=float, default=30.0,
                    help="half-angle of the bearing pressure on the cradle, deg (default 30)")
    ap.add_argument("--h-hook", type=float, default=1.2, help="element size in the hook, mm (default 1.2)")
    ap.add_argument("--h-plate", type=float, default=3.0, help="element size elsewhere, mm (default 3.0)")
    ap.add_argument("--solver", choices=["elmer", "ccx", "both"], default="elmer")
    ap.add_argument("--elmer-linsys", choices=sorted(ELMER_LINSYS), default="cg",
                    help="Elmer linear solver: cg (iterative, default) or umfpack (direct, slower)")
    ap.add_argument("--ceiling", choices=["free", "bonded"], default="free",
                    help="free: only the screw seats are fixed; bonded: the ceiling face is also held in z "
                         "(bilateral, an upper bound on ceiling support)")
    ap.add_argument("--ccx", type=Path, default=shutil.which("ccx"), help="CalculiX ccx executable")
    ap.add_argument("--workdir", type=Path, default=DEFAULT_WORKDIR,
                    help=f"parent of the case directories (default {DEFAULT_WORKDIR})")
    ap.add_argument("--vtu", type=Path, help="also write the Elmer result VTU to this path")
    args = ap.parse_args()
    if args.load_kg <= 0:
        ap.error("--load-kg must be > 0")
    if not 0 < args.contact_deg <= 90:
        ap.error("--contact-deg must be in (0, 90]")
    if not 0 < args.h_hook <= args.h_plate:
        ap.error("need 0 < --h-hook <= --h-plate")
    if args.solver in ("ccx", "both") and (args.ccx is None or not Path(args.ccx).is_file()):
        ap.error("CalculiX not found: pass --ccx PATH")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.name):
        ap.error("--name may only contain letters, digits, '_', '.', '-'")
    args.overrides = parse_overrides(args.overrides, ap)
    return args


def parse_overrides(items: list[str], ap: argparse.ArgumentParser) -> dict[str, object]:
    out: dict[str, object] = {}
    for item in items:
        name, sep, raw = item.partition("=")
        if not sep or not name.isidentifier():
            ap.error(f"--set expects NAME=VALUE, got {item!r}")
        try:
            out[name] = ast.literal_eval(raw)
        except (ValueError, SyntaxError):
            ap.error(f"--set {name}: {raw!r} is not a Python literal")
    return out


def load_design(overrides: dict[str, object]) -> dict:
    """Execute the design script with parameters replaced; return its namespace."""
    src = DESIGN.read_text(encoding="utf-8")
    params, marker, rest = src.partition("# ---- derived ----")
    if not marker:
        raise SystemExit(f"{DESIGN.name}: '# ---- derived ----' marker not found")
    for name, value in overrides.items():
        pat = re.compile(rf"^{re.escape(name)}(\s*=\s*)[^#\n]*?(\s*#.*)?$", re.MULTILINE)
        params, n = pat.subn(lambda m, name=name, value=value: f"{name}{m.group(1)}{value!r}{m.group(2) or ''}",
                             params, count=1)
        if n != 1:
            raise SystemExit(f"--set {name}: not a design parameter of {DESIGN.name}")
    ns: dict = {"__name__": "fea_hook"}
    exec(compile(params + marker + rest, str(DESIGN), "exec"), ns)  # noqa: S102 - the repo's own design script
    return ns


def classify_faces(gmsh, ns: dict) -> tuple[list[int], list[int], list[int]]:
    """Counterbore floors (screw seats), the cradle's inner cylinder and the plate's ceiling face."""
    z_seat = -ns["plate_t"] + ns["head_bore_depth"]
    seat, cradle, ceiling = [], [], []
    for _, tag in gmsh.model.getEntities(2):
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, tag)
        if max(abs(z0), abs(z1)) < 0.05:
            ceiling.append(tag)
            continue
        if max(abs(z0 - z_seat), abs(z1 - z_seat)) < 0.05 and any(
            abs((x0 + x1) / 2 - sx) < 0.05 and abs((y0 + y1) / 2 - sy) < 0.05 for sx, sy in ns["screws"]
        ):
            seat.append(tag)
            continue
        lo, hi = gmsh.model.getParametrizationBounds(2, tag)
        uv = [c for fu in (0.2, 0.5, 0.8) for fv in (0.2, 0.5, 0.8)
              for c in (lo[0] + fu * (hi[0] - lo[0]), lo[1] + fv * (hi[1] - lo[1]))]
        p = np.asarray(gmsh.model.getValue(2, tag, uv)).reshape(-1, 3)
        if np.all(np.abs(np.hypot(p[:, 1], p[:, 2] - ns["arc_cz"]) - ns["r_in"]) < 1e-3):
            cradle.append(tag)
    if len(seat) != len(ns["screws"]) or len(cradle) != 1 or len(ceiling) != 1:
        raise SystemExit(f"face classification failed: {len(seat)} seats, {len(cradle)} cradle, "
                         f"{len(ceiling)} ceiling faces")
    return seat, cradle, ceiling


def build_mesh(step: Path, msh: Path, ns: dict, h_hook: float, h_plate: float) -> dict:
    import gmsh

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("bracket")
        gmsh.model.occ.importShapes(str(step))
        gmsh.model.occ.synchronize()
        vols = gmsh.model.getEntities(3)
        if len(vols) != 1:
            raise SystemExit(f"expected one solid in {step.name}, got {len(vols)}")
        seat, cradle, ceiling = classify_faces(gmsh, ns)

        bb = gmsh.model.getBoundingBox(3, vols[0][1])
        fld = gmsh.model.mesh.field
        box = fld.add("Box")
        half_x = ns["hook_w"] / 2 + ns["stem_fillet_r"] + 2
        for key, val in {"VIn": h_hook, "VOut": h_plate, "XMin": -half_x, "XMax": half_x,
                         "YMin": bb[1] - 1, "YMax": bb[4] + 1, "ZMin": bb[2] - 1,
                         "ZMax": -ns["plate_t"] + 2, "Thickness": 6}.items():
            fld.setNumber(box, key, val)
        fields = [box]
        for sx, sy in ns["screws"]:
            ball = fld.add("Ball")
            for key, val in {"VIn": min(1.2, h_plate), "VOut": h_plate, "Radius": ns["head_bore_d"] / 2 + 3,
                             "Thickness": 5, "XCenter": sx, "YCenter": sy, "ZCenter": -ns["plate_t"] / 2}.items():
                fld.setNumber(ball, key, val)
            fields.append(ball)
        fmin = fld.add("Min")
        fld.setNumbers(fmin, "FieldsList", fields)
        fld.setAsBackgroundMesh(fmin)
        for key, val in {"Mesh.MeshSizeExtendFromBoundary": 0, "Mesh.MeshSizeFromPoints": 0,
                         "Mesh.MeshSizeFromCurvature": 8, "Mesh.MeshSizeMax": h_plate,
                         # the floor stops curvature sizing from collapsing on degenerate fillet edges
                         "Mesh.MeshSizeMin": h_hook / 4,
                         "Mesh.Algorithm3D": 10, "Mesh.ElementOrder": 2, "Mesh.Optimize": 1,
                         "Mesh.MshFileVersion": 2.2}.items():
            gmsh.option.setNumber(key, val)

        gmsh.model.addPhysicalGroup(3, [vols[0][1]], 1, name="bracket")
        gmsh.model.addPhysicalGroup(2, seat, SEAT, name="screw_seats")
        gmsh.model.addPhysicalGroup(2, cradle, CRADLE, name="cradle")
        gmsh.model.addPhysicalGroup(2, ceiling, CEILING, name="ceiling")
        gmsh.model.mesh.generate(3)
        gmsh.write(str(msh))

        tags, coords, _ = gmsh.model.mesh.getNodes()
        index = np.full(int(tags.max()) + 1, -1, dtype=np.int64)
        index[tags.astype(np.int64)] = np.arange(len(tags))
        etypes, _, enodes = gmsh.model.mesh.getElements(3)
        if list(etypes) != [11]:
            raise SystemExit(f"expected only 10-node tets, got element types {list(etypes)}")
        tris = []
        for tag in cradle:
            _, _, tnodes = gmsh.model.mesh.getElements(2, tag)
            tris.append(index[np.asarray(tnodes[0], dtype=np.int64)].reshape(-1, 6))

        def face_nodes(faces):
            return np.unique(np.concatenate(
                [index[np.asarray(gmsh.model.mesh.getNodes(2, t, includeBoundary=True)[0], dtype=np.int64)]
                 for t in faces]))

        return {
            "X": coords.reshape(-1, 3),
            "tets": index[np.asarray(enodes[0], dtype=np.int64)].reshape(-1, 10),
            "cradle_tris": np.concatenate(tris),
            "seat_nodes": face_nodes(seat),
            "ceiling_nodes": face_nodes(ceiling),
        }
    finally:
        gmsh.finalize()


def bearing_traction(P: np.ndarray, ns: dict, phi0: float) -> np.ndarray:
    """Traction for p0 = 1: radial outward from the cradle axis, cosine in phi from the bottom."""
    dy, dz = P[:, 1], P[:, 2] - ns["arc_cz"]
    r = np.hypot(dy, dz)
    phi = np.arctan2(dy, -dz)
    p = np.where(np.abs(phi) < phi0, np.cos(0.5 * np.pi * phi / phi0), 0.0)
    return np.stack([np.zeros_like(p), p * dy / r, p * dz / r], axis=1)


def consistent_loads(X: np.ndarray, tris: np.ndarray, ns: dict, phi0: float) -> np.ndarray:
    """Nodal forces (p0 = 1) from integrating the traction over the curved 6-node triangles."""
    F = np.zeros_like(X)
    xe = X[tris]                                                   # (T, 6, 3)
    for (l0, l1, l2), w in zip(TRI_L, TRI_W):
        N = np.array([l0 * (2 * l0 - 1), l1 * (2 * l1 - 1), l2 * (2 * l2 - 1), 4 * l0 * l1, 4 * l1 * l2, 4 * l2 * l0])
        dxi = np.array([1 - 4 * l0, 4 * l1 - 1, 0, 4 * (l0 - l1), 4 * l2, -4 * l2])
        deta = np.array([1 - 4 * l0, 0, 4 * l2 - 1, -4 * l1, 4 * l1, 4 * (l0 - l2)])
        P = np.einsum("k,tkd->td", N, xe)
        dA = np.linalg.norm(np.cross(np.einsum("k,tkd->td", dxi, xe), np.einsum("k,tkd->td", deta, xe)), axis=1)
        t = bearing_traction(P, ns, phi0) * (0.5 * w * dA)[:, None]
        np.add.at(F, tris, N[None, :, None] * t[:, None, :])
    return F


def elmer_tables(ns: dict, phi0: float, p0: float, n: int = 121) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Traction on the cradle cylinder tabulated against y (unique on the lower half of the cradle)."""
    r = ns["r_in"]
    y = np.linspace(-r * np.sin(phi0), r * np.sin(phi0), n)
    phi = np.arcsin(y / r)
    p = p0 * np.cos(0.5 * np.pi * phi / phi0)
    y = np.concatenate([[-2 * r], y, [2 * r]])
    f2 = np.concatenate([[0.0], p * np.sin(phi), [0.0]])
    f3 = np.concatenate([[0.0], -p * np.cos(phi), [0.0]])
    f2[[1, -2]] = f3[[1, -2]] = 0.0
    return y, f2, f3


ELMER_LINSYS = {
    "cg": """  Linear System Solver = Iterative
  Linear System Iterative Method = CG
  Linear System Preconditioning = ILU0
  Linear System Max Iterations = 20000
  Linear System Convergence Tolerance = 1.0e-9
  Linear System Residual Output = 500""",
    "umfpack": """  Linear System Solver = Direct
  Linear System Direct Method = umfpack""",
}


def write_elmer_case(case: Path, y: np.ndarray, f2: np.ndarray, f3: np.ndarray, linsys: str, ceiling: str) -> None:
    t2 = "\n".join(f"      {a:.6f} {b:.9e}" for a, b in zip(y, f2))
    t3 = "\n".join(f"      {a:.6f} {b:.9e}" for a, b in zip(y, f3))
    ceiling_bc = f"""
Boundary Condition 3
  Name = "ceiling"
  Target Boundaries(1) = {CEILING}
  Displacement 3 = 0.0
End
""" if ceiling == "bonded" else ""
    sif = f"""Header
  CHECK KEYWORDS Warn
  Mesh DB "." "mesh"
End

Simulation
  Max Output Level = 4
  Coordinate System = Cartesian 3D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals(1) = 1
End

Solver 1
  Equation = Linear elasticity
  Procedure = "StressSolve" "StressSolver"
  Variable = -dofs 3 Displacement
  Calculate Stresses = True
  Calculate Loads = True
  Displace Mesh = False
{ELMER_LINSYS[linsys]}
End

Solver 2
  Exec Solver = After Simulation
  Equation = ResultOutput
  Procedure = "ResultOutputSolve" "ResultOutputSolver"
  Output File Name = "results"
  Vtu Format = True
End

Equation 1
  Active Solvers(1) = 1
End

Material 1
  Name = "PLA"
  Youngs Modulus = {E_PLA}
  Poisson Ratio = {NU_PLA}
End

Body 1
  Target Bodies(1) = 1
  Equation = 1
  Material = 1
End

Boundary Condition 1
  Name = "screw_seats"
  Target Boundaries(1) = {SEAT}
  Displacement 1 = 0.0
  Displacement 2 = 0.0
  Displacement 3 = 0.0
End

Boundary Condition 2
  Name = "cradle"
  Target Boundaries(1) = {CRADLE}
  Force 2 = Variable Coordinate 2
    Real
{t2}
    End
  Force 3 = Variable Coordinate 2
    Real
{t3}
    End
End
{ceiling_bc}"""
    (case / "case.sif").write_text(sif, encoding="ascii")
    (case / "ELMERSOLVER_STARTINFO").write_text("case.sif\n1\n", encoding="ascii")


def run(cmd: list[str], cwd: Path, log: str) -> str:
    exe = shutil.which(cmd[0])
    if exe is None:
        raise SystemExit(f"{cmd[0]} not found on PATH")
    with open(cwd / log, "w", encoding="utf-8", errors="replace") as fh:
        rc = subprocess.run([exe, *cmd[1:]], cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, check=False).returncode
    text = (cwd / log).read_text(encoding="utf-8", errors="replace")
    if rc:
        raise SystemExit(f"{Path(cmd[0]).name} exited {rc}; see {cwd / log}")
    return text


def solve_elmer(case: Path, mesh: dict, support: np.ndarray) -> dict:
    run(["ElmerGrid", "14", "2", "mesh.msh", "-autoclean"], case, "elmergrid.log")
    log = run(["ElmerSolver"], case, "elmersolver.log")
    if "ERROR" in log.upper().replace("ERROR TOLERANCE", ""):
        raise SystemExit(f"ElmerSolver reported an error; see {case / 'elmersolver.log'}")
    import meshio

    vtus = sorted((case / "mesh").glob("results*.vtu")) + sorted(case.glob("results*.vtu"))
    if not vtus:
        raise SystemExit(f"no Elmer results VTU in {case}")
    res = meshio.read(vtus[-1])
    pd = {k.lower(): np.asarray(v, dtype=float) for k, v in res.point_data.items()}
    # same component order as CalculiX: xx yy zz xy yz zx
    S = np.column_stack([pd[f"stress_{c}"].ravel() for c in ("xx", "yy", "zz", "xy", "yz", "xz")])
    order = remap(res.points, mesh["X"])
    out = {"S": S[order], "U": pd["displacement"][order], "vtu": vtus[-1]}
    if "displacement loads" in pd:
        out["reaction"] = pd["displacement loads"][order][support].sum(axis=0)
    return out


def remap(points: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Index into `points` for each node of X (same coordinates, possibly reordered)."""
    from scipy.spatial import cKDTree

    dist, idx = cKDTree(points).query(X)
    if len(points) != len(X) or dist.max() > 1e-6 or len(np.unique(idx)) != len(X):
        raise SystemExit(f"solver mesh does not match the Gmsh mesh (max node offset {dist.max():.3g} mm)")
    return idx


def solve_ccx(case: Path, mesh: dict, loads: np.ndarray, ccx: Path, ceiling: str) -> dict:
    X, tets = mesh["X"], mesh["tets"]
    lines = ["*NODE, NSET=NALL"]
    lines += [f"{i + 1}, {x:.9f}, {y:.9f}, {z:.9f}" for i, (x, y, z) in enumerate(X)]
    lines.append("*ELEMENT, TYPE=C3D10, ELSET=EALL")
    lines += [f"{e + 1}, " + ", ".join(str(n + 1) for n in t[GMSH_TO_C3D10]) for e, t in enumerate(tets)]
    sets = {"SEAT": mesh["seat_nodes"] + 1, "CEIL": mesh["ceiling_nodes"] + 1}
    sets["SUPPORT"] = np.union1d(sets["SEAT"], sets["CEIL"]) if ceiling == "bonded" else sets["SEAT"]
    for name, ids in sets.items():
        lines.append(f"*NSET, NSET={name}")
        lines += [", ".join(map(str, ids[i:i + 16])) for i in range(0, len(ids), 16)]
    lines += ["*MATERIAL, NAME=PLA", "*ELASTIC", f"{E_PLA}, {NU_PLA}",
              "*SOLID SECTION, ELSET=EALL, MATERIAL=PLA", "*STEP", "*STATIC", "*BOUNDARY", "SEAT, 1, 3, 0."]
    if ceiling == "bonded":
        lines.append("CEIL, 3, 3, 0.")
    lines.append("*CLOAD")
    for n in np.nonzero(np.abs(loads).sum(axis=1) > 0)[0]:
        lines += [f"{n + 1}, {d + 1}, {loads[n, d]:.9e}" for d in (1, 2) if loads[n, d] != 0.0]
    lines += ["*NODE FILE", "U", "*EL FILE", "S", "*NODE PRINT, NSET=SUPPORT, TOTALS=ONLY", "RF", "*END STEP"]
    (case / "ccx.inp").write_text("\n".join(lines) + "\n", encoding="ascii")
    run([str(ccx), "-i", "ccx"], case, "ccx.log")
    S = np.full((len(X), 6), np.nan)
    U = np.full((len(X), 3), np.nan)
    block = None
    num = re.compile(r"[-+]?\d\.\d+E[-+]\d+")
    with open(case / "ccx.frd", encoding="ascii", errors="replace") as fh:
        for line in fh:
            if line.startswith(" -4"):
                block = line.split()[1]
            elif line.startswith(" -3"):
                block = None
            elif line.startswith(" -1") and block in ("DISP", "STRESS"):
                m = re.match(r" -1\s*(\d+)(.*)", line)
                vals = [float(v) for v in num.findall(m.group(2))]
                target = S if block == "STRESS" else U
                target[int(m.group(1)) - 1] = vals[: target.shape[1]]
    dat = (case / "ccx.dat").read_text(encoding="ascii", errors="replace")
    tot = re.search(r"total force.*?\n\s*\n?\s*([-+\dEe. ]+)\n", dat, re.DOTALL)
    out = {"S": S, "U": U}
    if tot:
        out["reaction"] = np.array([float(v) for v in tot.group(1).split()])
    return out


def principal(S: np.ndarray) -> np.ndarray:
    T = np.empty((len(S), 3, 3))
    T[:, 0, 0], T[:, 1, 1], T[:, 2, 2] = S[:, 0], S[:, 1], S[:, 2]
    T[:, 0, 1] = T[:, 1, 0] = S[:, 3]
    T[:, 1, 2] = T[:, 2, 1] = S[:, 4]
    T[:, 0, 2] = T[:, 2, 0] = S[:, 5]
    return np.linalg.eigvalsh(T)[:, ::-1]


def summarize(X: np.ndarray, S: np.ndarray, U: np.ndarray, ns: dict) -> dict:
    """Peak stresses in the hook (below the plate) and in the plate away from the fixed seats."""
    s1 = principal(S)[:, 0]
    vm = np.sqrt(0.5 * ((S[:, 0] - S[:, 1]) ** 2 + (S[:, 1] - S[:, 2]) ** 2 + (S[:, 2] - S[:, 0]) ** 2)
                 + 3 * (S[:, 3] ** 2 + S[:, 4] ** 2 + S[:, 5] ** 2))
    hook = X[:, 2] < -ns["plate_t"] - 0.25
    near_seat = np.zeros(len(X), dtype=bool)
    for sx, sy in ns["screws"]:
        near_seat |= np.hypot(X[:, 0] - sx, X[:, 1] - sy) < ns["head_bore_d"] / 2 + 2
    regions = {"hook": hook, "plate": ~hook & ~near_seat}
    out: dict = {"max_deflection_mm": float(np.linalg.norm(U, axis=1).max())}
    for rname, mask in regions.items():
        idx = np.nonzero(mask)[0]
        for qname, q in {"szz": S[:, 2], "sxx": S[:, 0], "syy": S[:, 1], "s1": s1, "von_mises": vm}.items():
            i = idx[np.argmax(q[idx])]
            out[f"{rname}_{qname}"] = {"MPa": round(float(q[i]), 3), "at": [round(float(c), 2) for c in X[i]]}
    out["hook_szz"]["dz_above_tangent_mm"] = round(out["hook_szz"]["at"][2] - ns["arc_cz"], 2)
    # peak layer-normal stress per 1 mm print layer band (layers are z-planes when printed plate-down)
    z_edges = np.arange(-ns["plate_t"] - 1.0, X[:, 2].min() - 1.0, -1.0)
    profile = []
    for z_hi, z_lo in itertools.pairwise(z_edges):
        band = hook & (X[:, 2] <= z_hi) & (X[:, 2] > z_lo)
        if band.any():
            profile.append([round(float((z_hi + z_lo) / 2), 1), round(float(S[band, 2].max()), 3)])
    out["layer_szz_peak"] = profile
    # the stem's two peaks: at the root fillet, and where it turns into the cradle; split at mid-stem
    z_mid = (-ns["plate_t"] + ns["arc_cz"]) / 2
    for key, rows in {"root": [p for p in profile if p[0] > z_mid], "shoulder": [p for p in profile if p[0] <= z_mid]}.items():
        z, s = max(rows, key=lambda p: p[1])
        out[f"layer_szz_{key}"] = {"MPa": s, "z": z}
    return out


def main() -> None:
    args = parse_args()
    case = args.workdir / args.name
    case.mkdir(parents=True, exist_ok=True)
    ns = load_design(args.overrides)
    step = case / "bracket.step"
    ns["export_step"](ns["result"], str(step))
    mesh = build_mesh(step, case / "mesh.msh", ns, args.h_hook, args.h_plate)

    phi0 = np.radians(args.contact_deg)
    force = args.load_kg * G
    unit = consistent_loads(mesh["X"], mesh["cradle_tris"], ns, phi0)
    p0 = force / -unit[:, 2].sum()
    loads = unit * p0

    report = {
        "name": args.name,
        "overrides": args.overrides,
        "load_N": round(force, 3),
        "contact_deg": args.contact_deg,
        "ceiling": args.ceiling,
        "p0_MPa": round(p0, 4),
        "mesh": {"nodes": len(mesh["X"]), "tets": len(mesh["tets"]), "h_hook": args.h_hook, "h_plate": args.h_plate},
        "volume_mm3": round(ns["result"].volume, 1),
        "mass_g": round(ns["result"].volume * PLA_G_PER_MM3, 1),
    }
    print(f"[{args.name}] {report['mesh']['tets']} tets, {report['mesh']['nodes']} nodes, "
          f"load {force:.2f} N, p0 {p0:.4f} MPa")

    if args.solver in ("elmer", "both"):
        write_elmer_case(case, *elmer_tables(ns, phi0, p0), args.elmer_linsys, args.ceiling)
        support = mesh["seat_nodes"]
        if args.ceiling == "bonded":
            support = np.union1d(support, mesh["ceiling_nodes"])
        t0 = time.perf_counter()
        res = solve_elmer(case, mesh, support)
        report["elmer"] = summarize(mesh["X"], res["S"], res["U"], ns)
        report["elmer"]["solve_s"] = round(time.perf_counter() - t0, 1)
        if "reaction" in res:
            report["elmer"]["reaction_N"] = [round(float(v), 3) for v in res["reaction"]]
        if args.vtu:
            shutil.copyfile(res["vtu"], args.vtu)
        elmer_S = res["S"]
    if args.solver in ("ccx", "both"):
        t0 = time.perf_counter()
        res = solve_ccx(case, mesh, loads, args.ccx, args.ceiling)
        report["ccx"] = summarize(mesh["X"], res["S"], res["U"], ns)
        report["ccx"]["solve_s"] = round(time.perf_counter() - t0, 1)
        if "reaction" in res:
            report["ccx"]["reaction_N"] = [round(float(v), 3) for v in res["reaction"]]
        if args.solver == "both":
            hook = mesh["X"][:, 2] < -ns["plate_t"] - 0.25
            d = elmer_S[hook, 2] - res["S"][hook, 2]
            ref = np.abs(res["S"][hook, 2]).max()
            report["elmer_vs_ccx_szz_hook"] = {"max_abs_diff_MPa": round(float(np.abs(d).max()), 4),
                                               "rms_diff_rel_peak": round(float(np.sqrt(np.mean(d ** 2)) / ref), 5)}

    (case / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
