# Seam alignment on the loop — High Level Design

## Purpose and scope

The seam placer joins the seams of neighboring layers into strings and smooths each string with
a fitted curve. This document covers the parts of `SeamPlacer` that decide whether the seam line
runs evenly on a smooth surface:

- where the painted (Enforced) candidates of a loop begin and end, and which one is central
  (`process_perimeter_polygon`);
- how a string is collected (`find_seam_string`);
- how the final seam position of a string member is derived from the curve (last step of
  `align_seam_points`), including how a painted (Enforced) seam keeps to painted corners.

Candidate selection, curve fitting and seam placement on the printed loop (`place_seam`) are
referred to only where they shape the result. Alignment runs for the Aligned, Aligned Back and Rear
seam positions; Nearest and Random do not align. A perimeter outside every aligned string keeps its
candidate as the seam.

## Concepts

Aligned and Aligned Back also apply a bounded length optimization after every original string has
been fitted. A scalar cubic B-spline, with the original string's knot spacing, represents signed
displacements along consistently oriented layer contours. Zero coefficients reproduce the original
positions. Each objective evaluation uses actual points on those contours, not a spatial curve
followed by projection.

The objective combines finite penalties derived from the original angle, visibility, overhang and
paint-center preferences with 3D polyline length. The experimental length multiplier defaults to
0.3 after conservative sphere trials; `ORCA_SEAM_CHAIN_LENGTH_WEIGHT=0` bypasses the entire stage. This scalar cost is not an exact
reformulation of the original lexicographic candidate selection.

Exact painted intervals impose linear bounds on the displacement coefficients. Full angle snapping
at 55 degrees fixes the corresponding displacement to zero; the solver eliminates these equalities.
Legacy Enforced candidates without exact intervals retain their original positions when their
allowed painted range cannot be reconstructed reliably. Outside all painted facets' heights, with
no weak modifiers, an Enforced mark is provably only a proximity preference: it stays movable with
the ordinary angle-based anchor cost, without the painted-center preference or its extra weight.
No green interval is invented where paint does not cross the layer.
Fixed points remain part of the length calculation.
Strong Precise Seam, Blocked-selected strings and embedded joins keep the original string.

The solver uses a damped quadratic model, feasible active-set QP steps and backtracking on the actual
objective. QP directions are solved in the active constraints' null space to avoid ill-conditioning
from nearly coincident bounds. At mesh vertices, feasible local and coherent spline probes test the
actual cost at several scales down to 0.001 mm. A complete poll without significant improvement
certifies a stop at this search resolution, not exact KKT stationarity or a global minimum.
All loops have finite budgets. Only a converged, feasible improvement with non-increasing
length replaces a string; errors, exhausted budgets and unsuccessful searches restore the complete
original string. Support metadata is sampled at moved positions; fixed positions keep their metadata.
Smooth scalar displacement couples layers but does not guarantee a continuously differentiable
spatial path on a faceted surface. Final placement on the printed loops remains part of validation.

- **Candidate (C).** The seam candidate chosen for a perimeter in a string. Candidates are the
  polygon vertices plus points inserted near paint (see Painted ranges), so the chosen candidate
  can move in discrete steps from layer to layer.
- **Fitted point (F).** The value of the string's fitted curve at the layer height. It need not
  lie on the loop: where the loop size changes quickly between layers, as near the bottom of a
  sphere, it lies up to a few millimeters off it.
- **Measured angle.** `local_ccw_angle` of a candidate. It is measured at polygon vertices on arms
  at least one nozzle diameter long; inserted points have angle 0.
- **Blend.** `target = t·C + (1 − t)·F`. The angle weight `t = min(1, (|angle C| / 55°)^3)` pulls
  the seam onto sharp corners; at 55° and above the seam is the candidate itself.
- **Enforced run.** The maximal sequence of consecutive Enforced candidates of the perimeter that
  contains C, or the whole loop when every candidate is Enforced. With exact painted ranges it
  covers a range exactly; otherwise it approximates the painted area with the candidate spacing.

## Painted ranges

The paint is a set of small triangles on the model surface. Its **contour** is the set of their
edges that belong to one painted triangle only, taken from the painted triangles with their T-joints
split (`get_facets_strict`), so that edges inside the painted area pair up. Each contour edge keeps
the direction from it into the paint. The contour is built once per object and stored in an AABB
tree.

On a layer, the contour edges crossing the slicing plane give the exact ends of the painted band.
Each crossing within the line width of the loop is moved to its closest loop point; the direction
into the paint, against the loop direction there (the mean of both edges at a vertex), tells whether
paint lies ahead or behind. Along the loop, ends where paint begins and ends must alternate; each
pair bounds a **painted range**. Such a loop has **exact ends**:

- a vertex is Enforced if it lies in a range (within 1 µm);
- points are inserted at the range ends, at the range middle and every 0.2 mm from the range begin,
  and one in the middle of each gap between ranges, so that ranges with no vertex between them stay
  separate patches;
- the central candidate of the longest patch is the middle of its range, unless the patch has a
  sharp corner, which keeps priority as before.

Otherwise, which includes paint not crossing the layer (a dot between layers, the top and bottom of
a band, a fully painted closed surface), ends that do not alternate and a contour running along the
loop, a point is Enforced if paint lies within the line width of it, and edges near paint get points
every 0.2 mm from their start vertex. Such a run reaches past the paint by up to the line width,
depending on the vertices, and its central candidate is the middle one by count; on a smooth band
both shift the seam off the middle of the paint.

## Seam strings

A string starts at a seam and is extended upward layer by layer while the next layer has a seam
within reach. After the upward search ends, the downward search begins at the layer immediately
below the start. This keeps a matching neighboring seam in the same aligned string.

## Final position

For each member of a string:

1. **Sharp candidate** (`t = 1`): the seam is C. This is the cheapest check and comes first.
2. **Painted corner pull** (Enforced C only): `t = max(t, 0.4·k)` with the corner influence
   `k ∈ [0, 1]` below. Retaining the pull near measured corners protects corners; fading it on
   smooth contours reduces transmission of discrete candidate steps into the seam. If C's own angle
   reaches `A1`, `k = 1` without a search.
3. **Place on the loop.** The target is mapped onto the loop, within the allowed arc around C.

### Corner influence

`k` is the largest contribution `a_j · w_j` among candidates `j` of the same perimeter within `R1`
of C along the loop, `d_j` being the distance along the loop:

- `a_j` rises smoothly from 0 at a measured angle of `A0` to 1 at `A1`, so weak bends count partly
  and smooth contours not at all;
- `w_j` is 1 up to `R0` and falls smoothly to 0 at `R1`, so a corner entering or leaving the
  neighborhood does not switch `t` abruptly.

`R0` is the larger of the perimeter's line width and the largest configured nozzle diameter: the
measured angle of a vertex reflects geometry within at least one nozzle diameter, so a candidate
inserted next to a corner usually still sees that corner. Within `R0` of a corner of at least
`A1`, `t` and the target are the same as with a fixed minimum of 0.4.

### Placing the target on the loop

The target is generally off the loop. Its closest loop point is a polygon vertex over a whole range
of target positions, so it snaps to vertices and steps between layers when the vertices differ from
layer to layer.

Instead the loop gets a continuous normal field. Each vertex has the normalized sum of the normals
of its two edges, taken to the nearest points at a different position, so that duplicated points
(shared ends of extrusion paths, a repeated closing point) get one normal. Along an edge the normal
blends linearly between its vertex normals. A **solution** is a loop point whose normal line passes
through the target; on each edge this is a quadratic equation in the edge parameter.

The seam is the solution closest to the target over the whole loop. If it lies outside the allowed
arc, the seam is the end of that arc closer to the solution, so a solution leaving the arc stays at
the end it leaves through. Without any solution, which takes a degenerate loop, the seam is C.

The search walks outward from C on both sides, the shorter walked side first. An edge farther from
the target than the best solution so far is skipped without solving.

Properties:

- On a regular polygon every normal line passes through the center. For a target other than the
  center, the closest solution is the nearer of the two intersections of the loop with the line
  through the center and the target, wherever the vertices are; the seam equals it when it lies
  inside the allowed arc. With an even number of vertices the two intersections are symmetric about
  the center, so the nearer one is on the ray from the center through the target.
- When the closest solution is unique and lies inside the allowed arc, the result depends on the
  target and that arc only, not on which candidate inside the arc was chosen. A candidate step
  therefore does not pass into the seam directly. At an exact distance tie, the first solution found
  from C is kept.
- A selected solution varies continuously while its root branch remains nondegenerate, and so does
  the seam when the solution leaves or enters the allowed arc. The seam may jump when the closest
  solution changes (two solutions equally close, or the closest one disappearing together with
  another where the target crosses the curve of centers of the normal field), or when a solution
  outside the arc becomes equally close to both arc ends. Changes of the candidate types between
  layers change the allowed arc.

### Allowed arc

- For an Enforced C on a loop with exact ends: the Enforced run itself, which ends at the paint.
- For another Enforced C: the Enforced run, extended past each open end by the oversampling distance
  (0.2 mm) but never into a Blocked candidate. The run ends sit on candidates and move in sampling
  steps between layers, and Rear tends to choose the candidate right at an end; without the
  extension a target just past an end would copy those steps into the seam.
- For another C: the part of the loop between the Blocked candidates around it.
- For a Blocked C, or a loop with no limiting candidate: the whole loop.

## Numerical constants

| Constant | Value | Meaning |
| --- | --- | --- |
| Blend minimum near corners | 0.4 | Pull toward C on a painted corner |
| `A0`, `A1` | 15°, 30° | Measured angle where corner contribution starts and becomes full |
| `R0` | max(line width, largest nozzle diameter); 0.5 mm if both are zero | Full corner protection along the loop |
| `R1` | 2·`R0` | Corner influence ends |
| Run extension | 0.2 mm (`enforcer_oversampling_distance`); 0 with exact ends | How far the seam may go past an end of the Enforced run |
| Painted range tolerance | 1 µm (`enforced_end_tolerance`) | Loop positions this close are one point |
| Range end reach | line width | How far from the loop a contour crossing may lie |
| Solver threshold | 1e-12 | Absolute threshold on the equation coefficients that tells degenerate equations apart; not a distance tolerance |

## Known limitations

A painted seam can end up farther from a physical corner than with a fixed minimum of 0.4 when:

- the corner lies beyond `R0` of C along the loop; in the fade band up to `R1` a candidate step
  can also move the target more than before, because `t` and C change together;
- the bend is measured below `A1`, or its measured angle depends on which vertices the polygon
  has (a measured angle is not invariant to inserting or removing collinear vertices);
- the seam is held at an end of the allowed arc.

The pull `t ≥ 0.4` bounds how far the target lies from C only relative to `|C − F|`; with a
distant F the remainder can be millimeters.

The closest solution is searched over the whole loop, by distance to the target, not to C. Where
another part of the same loop is close in space but far along the loop, as across a thin wall or a
U-shaped recess, a target between them can land on that other part.

Near a corner the vertex normal is an average of two directions, so for a target well off the loop
the seam lands where the averaged normals point; the farther the target, the more the result
depends on the loop shape around C rather than on the tangent at C.

The Enforced run is built from candidate types, so a gap in the paint or a narrow blocked area
between two Enforced candidates is not detected; without exact ends the seam can land there or up
to 0.2 mm past a run end.

With exact ends, paint near the loop that does not cross the layer is ignored on that loop. A range
end is the closest loop point to the crossing; where the surface is far from vertical, as near the
bottom of a sphere, the crossing can lie farther than the line width from the loop, and the loop
falls back to the line-width rule. Blocked candidates still use the line-width rule.

A perimeter left out of the aligned strings, for example in a string shorter than six layers,
keeps its candidate as the seam.

Seam metadata stays with C: placement takes the overhang and the angle for inner walls from C,
not from the placed position. `place_seam` still prefers non-bridge paths, places inner walls and
staggers them on its own; the placed point is its input, not its output.

## Implementation and verification

- Code: `src/libslic3r/GCode/SeamPlacer.cpp` — `append_enforcer_contour`, `exact_enforced_intervals`,
  `process_perimeter_polygon`, `find_seam_string`, `align_seam_points` and the helpers
  `enforced_corner_influence` and `place_on_loop` in `SeamPlacerImpl` (declared in
  `SeamPlacer.hpp` for the tests).
- Tests: `tests/fff_print/test_seam_placer.cpp`:
  - a painted face crossed by a loop 0.2 mm in front of it: the Enforced candidates span the face
    exactly (1 µm) and the central one is its middle; over 10 Aligned layers whose vertex before the
    face moves, the seam stays at that middle;
  - a painted circle whose vertices shift irregularly between layers: the seam stays within
    ±0.3 rad of the point where the circle meets the painted face, its step along the loop between
    layers stays below a quarter of the candidate spacing near paint (0.05 mm) and the change of that
    step below 0.15 of it (0.03 mm); with a fixed 0.4 pull and an unchanged fitted point, a candidate
    step of one spacing alone moves the target by 0.4 of it;
  - a 41-layer square whose string starts at layer 39: layers 1 through 39 join the string, including
    layer 38 right below the start;
  - a painted square in Rear and Aligned modes: the seam lies within 0.45 mm of a corner;
  - `place_on_loop` directly: the ray property on a regular 24-gon; the same result from
    neighboring candidates on an asymmetric contour; continuity along target sweeps on a short loop,
    across the point where the closest solution would leave a bounded search, and at a duplicated
    corner point; the ends of an Enforced run with and without a Blocked neighbor and with exact
  ends, and Blocked
    candidates around a Neutral one; a run whose extension ends on the closing edge, also with the
    loop after another one and the run wrapping around the array end; the same arc end held while
    the target passes halfway around the loop; a solution leaving the arc through the end farther
    from the target stays at that end;
  - `enforced_corner_influence` directly: full protection within `R0`, the middle of the fade band,
    no influence from `R1` on, and a partial measured angle.
- Not covered by tests: a painted corner with differing C and F through the whole pipeline, weak
  and lost corner measurements, mixed bridge loops, Arachne and inner-wall staggering.
