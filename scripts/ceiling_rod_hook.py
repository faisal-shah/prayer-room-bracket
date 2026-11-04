# /// script
# requires-python = ">=3.12"
# dependencies = ["build123d"]
# ///
"""Ceiling-mount mid-span J-hook for a 1 in curtain rod.

Coordinate system: X runs along the rod, +Y points into the room (-Y toward
the window wall), Z is up with the ceiling plane at z = 0. The whole part
hangs below z = 0 and the rod's seated centreline is the X axis at z = rod_cz.

Layout: a crossbar plate lies along the rod on the ceiling with a pan-head
screw at each end, both directly above the rod's centreline, so the hook load
sits between the screws and nothing is cantilevered. The open J hangs from the
middle of the plate: with stem_on = "wall" the stem drops on the wall side of
the rod (hidden behind it), curves under the rod and rises on the room side
to a lip, so the rod enters from the room side. The cradle is set
`cradle_bias` higher than the straight line between the two end brackets so
the rod is lifted onto the hook and bears on it.

Print upright, ceiling face on the bed, with supports under the crook (the
part is sized for PLA: 8 mm walls keep the stem near 3 MPa under load).

Run standalone from the repo root (`uv run scripts/ceiling_rod_hook.py`) to
write output/ceiling_rod_hook.{step,stl,3mf}; the build123d MCP session runs
it with execute_file() and exports through its own tools. The STEP is in
design coordinates; the STL and 3MF are in print orientation.
"""
from math import cos, radians, sin

from build123d import *

# ---- measured on site (mm; converted from the inch figures) ----
rod_d = 25.67                 # rod diameter at the support spot (1.0105 in)
rod_seat_at_brackets = 49.8   # ceiling -> rod underside at the end brackets (2.093 - 0.1315 in, right bracket; left is 50.15)
rod_sag_at_spot = 56.55       # ceiling -> rod underside at the support spot today, unsupported (2.2265 in); reference only
wall_to_rod_far_side = 100.3  # window wall -> room-side face of the rod (3.95 in); reference only

# ---- design (mm) ----
stem_on = "wall"              # "wall" or "room": which side of the rod the stem drops on; the hook opens toward the other side
cradle_bias = 2.0             # cradle sits this much HIGHER than the straight rod line, so the rod pre-loads onto the hook
cradle_clearance = 0.3        # radial clearance between rod and cradle
wall_t = 8.0                  # hook wall thickness (stem and cradle)
hook_w = 24.0                 # hook width along the rod
lip_angle = 15.0              # lip ends this many degrees below the rod centreline (0 = half cup)
plate_t = 8.0                 # ceiling plate thickness
plate_x = 84.0                # plate length along the rod (15 mm from each screw centre to the end)
plate_y_stem_side = 28.0      # plate reach from the rod centreline toward the stem side (covers stem + fillet)
plate_y_open_side = 28.0      # plate reach toward the open side; equal to the stem side centres the hook on the plate
plate_corner_r = 8.0          # plate corner radius in plan view
plate_edge_r = 1.5            # round on the plate's visible (underside) perimeter
stem_fillet_r = 5.0           # fillet where the stem meets the plate
hook_round_r = 2.0            # round on every exposed hook edge (stem sides, cradle sides, lip tip)
screw_hole_d = 6.2            # clearance for the #10-12 screws in Everbilt ribbed anchor packs (#12 major dia 5.5 mm)
head_bore_d = 12.0            # counterbore for the pan head (#12 pan head is 10.8 mm across, #10 is 9.5 mm)
head_bore_depth = 2.5         # pan head 2.8-3.3 mm tall sits 0.3-0.8 mm proud; 5.5 mm of plate stays under it
screw_x = 27.0                # screws at (+/-screw_x, 0): above the rod centreline, either side of the hook,
                              # 15 mm of material to the plate end and 28 mm to each side

# ---- derived ----
if stem_on not in ("room", "wall"):
    raise ValueError(f'stem_on must be "room" or "wall", got {stem_on!r}')
stem_side = 1 if stem_on == "room" else -1  # +Y is the room side
rod_r = rod_d / 2
r_in = rod_r + cradle_clearance             # cradle inner radius
r_out = r_in + wall_t                       # cradle outer radius
cradle_z = -(rod_seat_at_brackets - cradle_bias)   # lowest point of the cradle surface (rod underside when seated)
arc_cz = cradle_z + r_in                    # cradle arc centre height
rod_cz = cradle_z + rod_r                   # rod centre height when seated
stem_top_z = -plate_t / 2                   # stem is buried half-way into the plate so the fuse is clean
arc_start = 0.0 if stem_side > 0 else 180.0          # arc leaves the stem here (deg from +Y, CCW) ...
arc_size = -stem_side * (180.0 - lip_angle)          # ... and sweeps under the rod to the lip
lip_end = arc_start + arc_size
plate_yc = stem_side * (plate_y_stem_side - plate_y_open_side) / 2
plate_y = plate_y_stem_side + plate_y_open_side
screws = [(screw_x, 0.0), (-screw_x, 0.0)]  # (x, y) of each screw; heads on the underside
eps = 1e-4


def _on_plane_z(e, z):
    bb = e.bounding_box()
    return abs(bb.max.Z - z) < eps and abs(bb.min.Z - z) < eps


def _tangent_seam(e):
    """The two edges along X where the straight stem faces meet the cradle arcs."""
    bb = e.bounding_box()
    return _on_plane_z(e, arc_cz) and bb.size.X > hook_w - eps


# ---- base 1: ceiling plate, top face on the ceiling at z = 0, underside perimeter rounded ----
with BuildPart() as plate_bp:
    with BuildSketch(Plane.XY.offset(-plate_t)):
        with Locations((0, plate_yc)):
            RectangleRounded(plate_x, plate_y, plate_corner_r)
    extrude(amount=plate_t)
plate = plate_bp.part
plate = plate.fillet(plate_edge_r, ShapeList([e for e in plate.edges() if _on_plane_z(e, -plate_t)]))

# ---- base 2: J hook, profile in the YZ plane (u = Y, v = Z) extruded along the rod ----
arc_c = (0.0, arc_cz)
lip_outer = (r_out * cos(radians(lip_end)), arc_cz + r_out * sin(radians(lip_end)))
lip_inner = (r_in * cos(radians(lip_end)), arc_cz + r_in * sin(radians(lip_end)))
with BuildPart() as hook_bp:
    with BuildSketch(Plane.YZ):
        with BuildLine():
            Line((stem_side * r_out, stem_top_z), (stem_side * r_out, arc_cz))   # stem, outer face
            CenterArc(arc_c, r_out, start_angle=arc_start, arc_size=arc_size)     # cradle, outer arc
            Line(lip_outer, lip_inner)                                            # lip tip (radial)
            CenterArc(arc_c, r_in, start_angle=arc_start, arc_size=arc_size)      # cradle, inner arc
            Line((stem_side * r_in, arc_cz), (stem_side * r_in, stem_top_z))     # stem, rod-side face
            Line((stem_side * r_in, stem_top_z), (stem_side * r_out, stem_top_z))  # buried top
        make_face()
    extrude(amount=hook_w / 2, both=True)
hook = hook_bp.part

# round every exposed hook edge; skip the buried top and the tangent seams (no edge to round there)
hook_edges = ShapeList([e for e in hook.edges() if not _on_plane_z(e, stem_top_z) and not _tangent_seam(e)])
hook = hook.fillet(hook_round_r, hook_edges)

# ---- fuse, then blend the stem into the plate ----
body = plate + hook
junction_edges = ShapeList(
    [
        e
        for e in body.edges()
        if _on_plane_z(e, -plate_t)
        and abs(e.bounding_box().center().X) < hook_w / 2 + eps
        and stem_side * e.bounding_box().center().Y > r_in - eps
        and stem_side * e.bounding_box().center().Y < r_out + eps
    ]
)
body = body.fillet(stem_fillet_r, junction_edges)

# ---- secondary features: counterbored screw holes for pan heads, on the underside ----
# Built as a cutter solid and subtracted, which keeps the step the same on
# build123d 0.11 (MCP session) and 0.12 (uv run) without add()/insert().
# Pan heads seat flat on the counterbore floor, so nothing wedges the plate open.
with BuildPart() as cutter_bp:
    with Locations(*[Location((x, y, -plate_t), (180, 0, 0)) for x, y in screws]):
        CounterBoreHole(
            radius=screw_hole_d / 2,
            counter_bore_radius=head_bore_d / 2,
            counter_bore_depth=head_bore_depth,
            depth=plate_t,
            mode=Mode.ADD,
        )
result = body - cutter_bp.part

# reference rod in its seated position (not part of the result)
rod = Cylinder(rod_r, 140.0, rotation=(0, 90, 0)).move(Location((0, 0, rod_cz)))

# print orientation: ceiling face on the bed, hook rising as a crook (needs supports under the crook)
print_part = result.rotate(Axis.X, 180)
print_part = print_part.move(Location((0, 0, -print_part.bounding_box().min.Z)))

if __name__ == "__main__":
    export_step(result, "output/ceiling_rod_hook.step")          # design coordinates: ceiling at z = 0, hook hanging below
    export_stl(print_part, "output/ceiling_rod_hook.stl")        # print orientation, ready to slice
    mesher = Mesher()
    mesher.add_shape(print_part)
    mesher.write("output/ceiling_rod_hook.3mf")
    print(f"hook drop to cradle {-cradle_z:.1f} mm, rod centre {-rod_cz:.2f} mm below ceiling")
