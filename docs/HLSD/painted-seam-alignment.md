# Painted seam alignment — High Level Design

## Purpose and scope

The seam placer aligns the seams of neighboring layers into strings and smooths each string with a fitted curve. This document describes how the final seam position of a **painted** (Enforced) seam point is derived from that curve: a painted seam should run as an even line on a smooth surface and should stay on a painted corner. It covers only the last step of alignment in `SeamPlacer::align_seam_points`. Candidate selection, string search, curve fitting and seam placement on the printed loop are unchanged and are referred to only where they shape the result.

The step applies to the Aligned, Aligned Back and Rear seam positions, the modes that run alignment. Nearest and Random do not align. Seam points that are not Enforced keep the plain blend described below.

## Concepts

- **Candidate (C).** The seam candidate chosen for a perimeter in a string. Candidates are the polygon vertices plus points inserted every 0.2 mm along edges near paint, so the chosen candidate moves in discrete steps from layer to layer.
- **Fitted point (F).** The value of the string's fitted curve at the layer height.
- **Measured angle.** `local_ccw_angle` of a candidate. It is measured at polygon vertices on arms at least one nozzle diameter long; inserted points have angle 0.
- **Blend.** `target = t·C + (1 − t)·F`. The angle weight `t = min(1, (|angle C| / 55°)^3)` pulls the seam onto sharp corners; at 55° and above the seam is the candidate itself.
- **Enforced run.** The maximal sequence of consecutive Enforced candidates of the perimeter that contains C, or the whole loop when every candidate is Enforced. It approximates the painted area along the loop with the candidate spacing.

## Policy

A painted seam point needs two things: it should stay on the paint, and it should snap to a painted corner. A fixed minimum `t ≥ 0.4` toward the candidate gives some of both, but on a smooth surface it also passes 40% of every candidate step into the seam, which shows as sideways jumps of the seam line. The step therefore keeps that minimum only near measured corners and keeps the seam on the paint by clamping it to the Enforced run instead.

For an Enforced point, in order of cost:

1. **Sharp candidate** (`t = 1`): the seam is C.
2. **Single painted point** (both neighbors of C are not Enforced): the seam is C.
3. **Corner influence.** `k ∈ [0, 1]` is the largest corner contribution among candidates of the same perimeter within `R1` of C along the loop; if C's own angle reaches `A1`, `k = 1` without a search. The weight becomes `t = max(t, 0.4·k)`.
4. **Clamp.** The blended target is kept within 0.2 mm of the Enforced run.

### Corner influence

Each candidate `j` at distance `d_j` along the loop contributes `a_j · w_j`:

- `a_j` rises smoothly from 0 at a measured angle of `A0` to 1 at `A1`, so weak bends count partly and smooth contours not at all;
- `w_j` is 1 up to `R0` and falls smoothly to 0 at `R1`, so a corner entering or leaving the neighborhood does not switch `t` abruptly.

`R0` is the larger of the perimeter's line width and the largest configured nozzle diameter: the measured angle of a vertex reflects geometry within at least one nozzle diameter, so a candidate inserted next to a corner still sees that corner. Within `R0` of a corner of at least `A1`, `t` and the blended target are the same as with the fixed minimum.

### Clamping to the Enforced run

The seam position is the point closest to the target within the oversampling distance (0.2 mm) of the run's edges. A target within that distance is kept as it is; a target farther away moves toward the closest run point until it is 0.2 mm from it. The slack matters at the run ends: they sit on candidates, so they move in sampling steps from layer to layer, and Rear tends to place the candidate right at an end. Clamping exactly to the end would copy those steps into the seam; the slack keeps a seam whose target stays near the end on the smooth target. The clamp is continuous in the target.

The closest run point is chosen among distances within 1 µm of the minimum by its distance to C along the run (the shorter arc on a closed loop), then by edge order from the perimeter start. Without this rule numerical noise could move the seam along the run.

## Numerical constants

| Constant | Value | Meaning |
| --- | --- | --- |
| Blend minimum near corners | 0.4 | Pull toward C on a painted corner |
| `A0`, `A1` | 15°, 30° | Measured angle where corner contribution starts and becomes full |
| `R0` | max(line width, largest nozzle diameter); 0.5 mm if both are zero | Full corner protection along the loop |
| `R1` | 2·`R0` | Corner influence ends |
| Run slack | 0.2 mm (`enforcer_oversampling_distance`) | How far the seam may stay off the Enforced run |
| Equal distance | 1 µm | Tie tolerance when choosing the closest run point |

## Known limitations

The seam can end up farther from a physical corner than with a fixed minimum in these cases:

- the corner lies between `R0` and `R1` or beyond `R1` of C along the loop; in the fade band a candidate step can also move the seam more than before, because `t` and C change together;
- the bend is measured below `A1`, or its measured angle depends on which vertices the polygon has (a measured angle is not invariant to inserting or removing collinear vertices);
- the clamp holds the target near an end of the run, the seam is a single painted point, or the closest run point switches to another branch of a U-shaped run;
- seam placement on the printed loop prefers non-bridge paths, places inner walls and staggers them on its own; the clamped point is its input, not its output.

The change of the blended target before clamping is at most `max(0.4 − t_angle, 0)·|C − F|`. No such bound exists for the clamped or printed seam.

The Enforced run is built from candidate types, so a gap in the paint or a narrow blocked area between two Enforced candidates is not detected, and the seam may stay up to 0.2 mm off the run; the seam can land in such places.

The candidate keeps its role for seam metadata: placement takes the overhang and the angle for inner walls from C, not from the clamped position.

## Implementation and verification

- Code: `src/libslic3r/GCode/SeamPlacer.cpp` — `align_seam_points`, with the helpers `enforced_corner_influence`, `is_single_enforced_point` and `clamp_to_enforced_run` in `SeamPlacerImpl`.
- Tests: `tests/fff_print/test_seam_placer.cpp` — a painted circle whose vertices shift irregularly between layers keeps a seam line without kinks on the paint; a painted square keeps the seam at a corner in Rear and Aligned modes.
- Not covered by tests: the fade band, weak and lost corner measurements, U-shaped runs, mixed bridge loops, Arachne and inner-wall staggering.
