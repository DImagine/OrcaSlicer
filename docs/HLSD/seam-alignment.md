# Seam alignment on the loop — High Level Design

## Purpose and scope

The seam placer joins the seams of neighboring layers into strings and smooths each string with
a fitted curve. This document covers two parts of `SeamPlacer` that decide whether the seam line
runs evenly on a smooth surface:

- how a string is collected (`find_seam_string`);
- how the final seam position of a string member is derived from the curve (last step of
  `align_seam_points`), including how a painted (Enforced) seam keeps to painted corners.

Candidate selection, curve fitting and seam placement on the printed loop (`place_seam`) are
unchanged and are referred to only where they shape the result. Alignment runs for the Aligned,
Aligned Back and Rear seam positions; Nearest and Random do not align. A perimeter outside every
aligned string keeps its candidate as the seam.

## Concepts

- **Candidate (C).** The seam candidate chosen for a perimeter in a string. Candidates are the
  polygon vertices plus points inserted every 0.2 mm along edges near paint, so the chosen
  candidate moves in discrete steps from layer to layer.
- **Fitted point (F).** The value of the string's fitted curve at the layer height. It need not
  lie on the loop: where the loop size changes quickly between layers, as near the bottom of a
  sphere, it lies up to a few millimeters off it.
- **Measured angle.** `local_ccw_angle` of a candidate. It is measured at polygon vertices on arms
  at least one nozzle diameter long; inserted points have angle 0.
- **Blend.** `target = t·C + (1 − t)·F`. The angle weight `t = min(1, (|angle C| / 55°)^3)` pulls
  the seam onto sharp corners; at 55° and above the seam is the candidate itself.
- **Enforced run.** The maximal sequence of consecutive Enforced candidates of the perimeter that
  contains C, or the whole loop when every candidate is Enforced. It approximates the painted
  area along the loop with the candidate spacing.

## Seam strings

A string starts at a seam and is extended upward layer by layer while the next layer has a seam
within reach; then it is extended downward from the start. When the upward search ends, the
downward search must begin at the layer right below the start. It used to step once more before
searching, so that layer was left out of the string. Being outside every string, it kept its
candidate while both neighbors were aligned, which showed as a one-layer step of the seam line.

## Final position

For each member of a string:

1. **Sharp candidate** (`t = 1`): the seam is C. This is the cheapest check and comes first.
2. **Painted corner pull** (Enforced C only). A fixed minimum `t ≥ 0.4` toward the candidate used
   to keep a painted seam on painted corners, but on a smooth surface it also passed 40% of every
   candidate step into the seam, which showed as sideways jumps. The minimum is now kept only near
   measured corners: `t = max(t, 0.4·k)` with the corner influence `k ∈ [0, 1]` below. If C's own
   angle reaches `A1`, `k = 1` without a search.
3. **Place on the loop.** The target is mapped onto the loop within the allowed arc around C.

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
`A1`, `t` and the target are the same as with the fixed minimum.

### Placing the target on the loop

The target is generally off the loop. Its closest loop point is a polygon vertex over a whole range
of target positions, so the closest point snaps to vertices and steps between layers when the
vertices differ from layer to layer.

Instead the loop gets a continuous normal field: each vertex has the normalized sum of its two
edge normals, and along an edge the normal blends linearly between its vertex normals. The seam
is the loop point whose normal line passes through the target; of several such points, the one
closest to the target. On each edge this is a quadratic equation in the edge parameter.

Properties:

- On a regular polygon every normal line passes through the center, so the seam is where the ray
  from the center through the target meets the loop, wherever the vertices are.
- The result depends on the target and the allowed arc only, not on which candidate inside that
  arc was chosen. A candidate step therefore no longer passes into the seam directly.
- As the target moves, the seam moves continuously, except where two solutions are equally close
  to the target (the same kind of switch as for the closest loop point).

### Allowed arc

The search walks away from C in both directions:

- for an Enforced C, along the Enforced run and up to the oversampling distance (0.2 mm) past an
  open end of it, never into a Blocked candidate. The run ends sit on candidates and move in
  sampling steps between layers, and Rear tends to choose the candidate right at an end; without
  the extension a target just past an end would copy those steps into the seam;
- for another C, up to the first Blocked candidate (no limit for a Blocked C);
- in all cases at most `2·|target − C| + 0.2 mm` along the loop, and at most one turn.

Without a solution in the allowed arc the seam is the end of the arc closer to the target.

## Numerical constants

| Constant | Value | Meaning |
| --- | --- | --- |
| Blend minimum near corners | 0.4 | Pull toward C on a painted corner |
| `A0`, `A1` | 15°, 30° | Measured angle where corner contribution starts and becomes full |
| `R0` | max(line width, largest nozzle diameter); 0.5 mm if both are zero | Full corner protection along the loop |
| `R1` | 2·`R0` | Corner influence ends |
| Run extension | 0.2 mm (`enforcer_oversampling_distance`) | How far the seam may go past an end of the Enforced run |
| Search window | twice the distance from C to the target, plus 0.2 mm | Longest arc searched on each side of C |

## Known limitations

A painted seam can end up farther from a physical corner than with the fixed minimum when:

- the corner lies beyond `R0` of C along the loop; in the fade band up to `R1` a candidate step
  can also move the target more than before, because `t` and C change together;
- the bend is measured below `A1`, or its measured angle depends on which vertices the polygon
  has (a measured angle is not invariant to inserting or removing collinear vertices);
- the target is placed at an end of the allowed arc.

The pull `t ≥ 0.4` bounds how far the target lies from C only relative to `|C − F|`; with a
distant F the remainder can be millimeters.

Near a corner the vertex normal is an average of two directions, so for a target well off the loop
the seam lands on the side the averaged normals point to; the farther the target, the more the
result depends on the loop shape around C rather than on the tangent at C.

The Enforced run is built from candidate types, so a gap in the paint or a narrow blocked area
between two Enforced candidates is not detected; the seam can land there or up to 0.2 mm past a
run end.

A perimeter left out of the aligned strings, for example in a string shorter than six layers,
keeps its candidate as the seam.

Seam metadata stays with C: placement takes the overhang and the angle for inner walls from C,
not from the placed position. `place_seam` still prefers non-bridge paths, places inner walls and
staggers them on its own; the placed point is its input, not its output.

## Implementation and verification

- Code: `src/libslic3r/GCode/SeamPlacer.cpp` — `find_seam_string`, `align_seam_points` and the
  helpers `enforced_corner_influence` and `place_on_loop` in `SeamPlacerImpl` (declared in
  `SeamPlacer.hpp` for the tests).
- Tests: `tests/fff_print/test_seam_placer.cpp`:
  - a painted circle whose vertices shift irregularly between layers: the seam stays within the
    painted angle range, and its position along the loop changes little between layers (first and
    second difference);
  - a 41-layer square whose string starts near the top: every layer below the start joins it;
  - a painted square in Rear and Aligned modes: the seam lies within 0.45 mm of a corner;
  - `place_on_loop` directly: the ray property on a regular polygon, the same result from
    neighboring candidates on an asymmetric contour, continuity along a sweep of the target on a
    short loop, the end of an Enforced run with and without a Blocked neighbor;
  - `enforced_corner_influence` directly: full protection within `R0`, the fade band, no influence
    from `R1` on, and a partial measured angle.
- Not covered by tests: weak and lost corner measurements, mixed bridge loops, Arachne and
  inner-wall staggering.
