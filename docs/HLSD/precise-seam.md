# Precise Seam — High Level Design

## Purpose and scope

Precise Seam lets a helper volume decide where the seam of an object goes. The
user attaches a mesh to an object as a Precise Seam modifier. On every layer,
the part of the external perimeter that lies inside the modifier's slice
determines where the seam must, may or must not be placed. The helper is a
persistent model object rather than paint on the surface, so it keeps working
when the design changes. A body swept along a path on the surface can guide the
seam along any trajectory.

The modifier is non-printing geometry. It takes no part in object slicing,
region assignment, filament selection or brim adhesion, and it affects only seam
placement during G-code export. Objects without Precise Seam volumes follow the
regular seam placement unchanged. In objects that have them, every external
perimeter is normalized first (see [Perimeter preparation](#perimeter-preparation));
perimeters that no modifier reaches then use the regular seam placement.

Precise Seam does not replace the seam placer. It feeds it: a modifier inserts
the points it needs into the perimeter and changes the enforced/blocked type of
seam candidates, the same typing mechanism as seam painting, and the configured
seam position then chooses among them.

This document covers the volume types and their priority, the data flow
through `SeamPlacer`, the extraction of perimeter segments from modifier
slices, both kinds of modifiers, numeric tolerances, diagnostics, known
limitations, project storage, invalidation and the user interface.

## Modifier types

Precise Seam adds six `ModelVolumeType` values after `SUPPORT_ENFORCER`, strong
types first and weak types after them. `is_precise_seam()`,
`is_precise_seam_strong()` and `is_precise_seam_weak()` are range checks that
depend on this order.

| Type | Group | Effect on an intersected perimeter |
| --- | --- | --- |
| `PRECISE_SEAM_CENTER` | strong | seam at the midpoint, by arc length, of the intersection |
| `PRECISE_SEAM_LEFT` | strong | seam at the first point of the intersection |
| `PRECISE_SEAM_RIGHT` | strong | seam at the last point of the intersection |
| `PRECISE_SEAM_ENFORCED` | weak | intersection marked as enforced |
| `PRECISE_SEAM_BLOCKED` | weak | intersection marked as blocked |
| `PRECISE_SEAM_NEUTRAL` | weak | intersection reset to neutral |

A strong modifier fixes a single point. The perimeter gets exactly one enforced
seam candidate there, and every other candidate is blocked. A weak modifier
retypes, and where needed adds, the candidates inside its intersection, like
painting does.

An **intersection** is a continuous part of the external perimeter's centerline
that lies inside the modifier's slice on that layer. It is a portion of the
perimeter, never a chord through the object. The centerline lies half an
extrusion width inside the model surface and depends on print settings, so a
modifier must reach clearly past the surface to cross it unambiguously.

### Terms

- **Segment:** an intersection as the code represents it (`PerimeterSegment`).
  User-facing texts call it an intersection.
- **Fragment:** a piece of the perimeter returned by clipping, before it is tied
  to the source contour.
- **Interval:** the bound part of one source edge, given by the edge index and a
  parameter range on that edge.
- **Zone:** a weak segment with its type (Enforced, Blocked or Neutral).
- **Boundary:** an end of a zone, inserted into the perimeter polygon.
- **Candidate:** a seam candidate of the seam placer, built from the points of
  the processed perimeter polygon (painted enforcers may add more).

First and last are taken along the perimeter oriented counter-clockwise as seen
from above. On an outer wall seen from outside, Left is therefore the left end
of the intersection. On the wall of a hole seen from inside the hole, the two
ends are swapped. Mirroring an object does not mirror the mode: perimeters stay
counter-clockwise, so Left remains the left end seen from outside, and the seam
moves to the other end of the modifier instead of following the mirrored model.

## Priority

The order of volumes in the object is the priority order, highest first.
`ModelObject::sort_volumes()` keeps every strong modifier before every weak one
and preserves the user's order within each group. The object list lets the user
drag a modifier only within its own group. A type change that crosses the group
boundary moves the volume to the end of its new group, with the lowest priority
there.

- **Strong:** modifiers are tried in priority order on each perimeter. The first
  one that yields a usable segment decides the seam. Within that modifier the
  longest segment wins; lengths are never compared across modifiers. Once a
  strong point is placed, no later strong modifier and no weak modifier is
  processed for that perimeter.
- **Weak:** every weak modifier applies. They are applied from the lowest
  priority to the highest, so the highest one overwrites overlapping zones.

A modifier without a usable segment on a perimeter passes the turn to the next
one. This includes a modifier that crosses the perimeter but whose fragments
were all discarded (see [Binding fragments](#binding-fragments-to-source-edges)).
The seam on that layer then comes from a lower-priority modifier, or from weak
or regular placement, and the user is told through the "unable to process
intersection" warning.

## Data flow

1. **Invalidation.** `Print::apply()` treats a change of Precise Seam volumes as
   a change of seam placement and invalidates G-code export; the object is not
   resliced (see [Print invalidation](#print-invalidation)).
2. **Modifier slices.** `SeamPlacer::init()` collects each object's Precise Seam
   volumes once, slices every volume separately and caches its regions with
   their bounding boxes.
3. **Perimeters.** Seam candidates are gathered in parallel over the layers.
   For objects with Precise Seam volumes, each external perimeter polygon is
   normalized and prepared once for all modifiers.
4. **Extraction.** For each modifier, the perimeter is clipped against the
   modifier's regions on that layer. The clipped fragments are bound back to the
   source edges of the perimeter and assembled into segments.
5. **Strong, then weak.** Strong modifiers try to insert one seam point into the
   perimeter polygon. If none succeeds, weak modifiers insert their zone
   boundaries and subdivide enforced edges.
6. **Candidates.** The seam placer builds candidates from the modified polygon.
   Painting assigns types first, weak zones overwrite them, and a strong point
   makes its candidate the only enforced one.
7. **Selection and restoration.** The configured seam position chooses the
   seams and aligns them. Afterwards the exact strong points are restored.
8. **Warnings.** After all objects are processed, one combined user warning is
   issued if any problem was found.

## Modifier slices

`init_precise_seam_data()` collects the Precise Seam volumes of each object:
strong ones in priority order and weak ones in reverse, so that weak zones can
be applied with last-write-wins. Each volume is sliced separately with
`PrintObject::slice_single_volume_regions()`, at the object's layer heights and
with the same centered transformation as the object. The slices keep every
region's outer contour together with its holes as an `ExPolygon`. Volumes are
not merged, so each keeps its own priority, and a modifier may have several
regions on one layer.

`prepare_modifier_slices()` moves the slices into `ModifierRegionsCache`,
pairing each region with the bounding box of its exterior. Empty layers keep
their slots, so the cache is indexed by object layer; `Layer::id()` includes raft
layers, which are subtracted. The cache is filled before candidates are gathered
and is only read afterwards, shared by both modifier kinds and all worker
threads without locking.

## Perimeter preparation

The seam placer works on external perimeter loops, including the walls of
holes. For objects with Precise Seam volumes, consecutive duplicate points and
the repeated closing point of each extrusion loop are removed: adjacent
extrusion paths share endpoints, and the resulting zero-length edges would
prevent point insertion at their junctions. Distinct visits to one point of a
self-touching contour are kept. Objects without Precise Seam volumes keep their
original points, so ordinary seam candidates are unaffected.

Each polygon is made counter-clockwise. A single `PreparedPerimeter` is then
built for all modifiers of that perimeter. It holds a validity check (at least
three points, no consecutive or closing duplicates), the bounding box, and the
clipping line: the polygon as an open polyline with its first point repeated at
the end. The preparation borrows the polygon and is used only while the polygon
is unchanged: strong processing returns immediately after inserting its point,
and weak processing collects all segments before it inserts anything. An
invalid perimeter receives no Precise Seam processing.

## Segment extraction

`extract_perimeter_segments()` turns one modifier's regions on one layer into
segments of the perimeter, each with its geometry and its position on the
source contour. Both modifier kinds consume these segments; the extractor is
told the modifier type so that it prepares only the data that type needs.

### Clipping

Regions whose bounding box does not overlap the perimeter's are skipped. The
clipping line is intersected with each remaining region by `intersection_pl()`,
which clips an open path against an `ExPolygon` with its holes attached, using
the nonzero rule. Clipping an open line yields only pieces of the perimeter, so
a modifier crossing the whole object produces two separate pieces rather than a
chord through the body. Holes in a modifier and several regions of one modifier
simply produce more pieces. The line is cut at vertex zero, so a piece crossing
that vertex arrives as two fragments.

### Binding fragments to source edges

Clipper returns coordinates only. Insertion needs the source edge of every
point, and coordinates alone are ambiguous where a contour visits the same
point twice. Each fragment is therefore bound to the source edges it covers,
producing intervals: an edge index with a parameter range on that edge.

- **Exact path.** For fragments with interior points, the second point is used
  as an anchor that must equal a source vertex exactly. Clipping keeps the
  vertices of an open path unchanged, including collinear ones. The following
  points must match successive source vertices in either direction; later
  occurrences of the anchor are tried if a sequence does not match. Only the two
  end cuts are projected onto their edges.
- **Projection path.** Two-point fragments, and fragments the exact path cannot
  match, are bound by projection. The first source edge that holds both points
  of the first pair, with distinct parameters, establishes the edge and
  direction. Every following pair must continue on the same edge or cross to the
  neighboring edge at their actual shared vertex, in the same direction. A pair
  continuing on the same edge reuses the previous pair's parameter for their
  shared point, so the two projections of one point cannot differ.
- **Failure.** A fragment that cannot be bound continuously is rolled back and
  discarded. Earlier fragments and other fragments are unaffected. The failure
  is counted, logged and reported to the user (see
  [Diagnostics](#diagnostics-and-warnings)).

Two rare rounding cases are handled only after both paths have failed, so the
normal path never pays for them:

- **Cut beside a vertex.** When a modifier boundary crosses within about one
  coordinate unit of a source vertex, Clipper can place the cut at the vertex's
  height but a few units beside it. The end pair then collapses to the vertex's
  parameter or misses both neighboring edges. An end cut closer than the
  snapping radius to a vertex of the fragment's own chain is snapped to that
  vertex: either its neighbor in the fragment (the cut is a rounded copy of it
  and is dropped) or a vertex that shares a source edge with that neighbor. The
  neighbor wins whenever it is within the radius. There is no search for the
  nearest vertex elsewhere on the perimeter. Ends that are themselves source
  vertices and ambiguous choices are left unchanged. Binding is then retried
  once with the same strict rules, so a wrong candidate can only fail again.
- **Contact.** A fragment that still fails but is shorter than the snapping
  radius is accepted as a contact and binds nothing. Insertion would collapse it
  onto one point anyway.

Both outcomes are recoveries, not failures: they show no user warning but leave
a log marker.

### Assembling segments

The intervals are sorted by edge and parameter. Intervals on the same occurrence
of an edge are united when they overlap or meet, by parameter or at the same
integer point; equal coordinates on different edges are never united. A
parameter of 1 is stored as parameter 0 of the next edge, so intervals on
adjacent edges meet exactly at their shared vertex. Consecutive intervals that
meet form one `PerimeterSegment`, and the last segment is joined with the first
when they meet at vertex zero, undoing the artificial cut of the clipping line.

Each segment keeps its polyline, the source edge of every polyline edge, and its
begin and end positions on the source contour. These data are retained for all
modifier types, including for debugging.

### Full containment

A modifier that covers the whole perimeter does not intersect it in any useful
sense. Both consumers skip that modifier for the perimeter and report full
containment to the user. Geometry and bindings stay in the result, but no strong
data are prepared.

The perimeter is fully contained when the united intervals cover every source
edge from parameter 0 to 1. A modifier boundary that merely touches the
perimeter counts as well:

- At a vertex or on an axis-aligned edge, clipping splits the line exactly at the
  touch, the pieces meet at one point, and the coverage is complete.
- On an inclined edge the touching point is usually not representable on the
  integer grid. The boundary pokes a few units across and leaves a real gap, so
  a single segment covers everything except that gap.

Weak insertion would collapse such a segment's boundaries onto one vertex and
turn the intended zone into a single candidate, and strong would put the seam at
the touch. A single segment is therefore also full containment in the cases
where insertion collapses it, exactly up to edges shorter than 2 µm:

- the uncovered length from its end to its begin is below 1 µm, or
- the gap spans one vertex, or starts at a vertex and ends on the next edge, and
  both ends lie within 1 µm of the vertex that ends the first gap edge, since
  each end then snaps onto it from its own edge.

A cheap filter runs first: both cases bring the segment's ends within 2 µm of
each other. A narrow band, an outside contact, a sharp spike and nearly touching
parts of a contour may pass the filter, but their uncovered part is long and
not around one vertex, so they remain ordinary segments. A band narrower than
1 µm then becomes a single candidate, which is its expected result.

## Strong modifiers

For a strong modifier, the extractor prepares each segment's target point
before anything is inserted, together with the source edge it lies on:

- **Left:** the segment's first point.
- **Right:** the segment's last point.
- **Center:** the point at half the segment's arc length.

Arc length is the sum of Euclidean edge lengths, not the chord or a vertex count.
Lengths are measured only when needed: for Center and when several segments of
one modifier compete.

`insert_strong_seam_point()` selects the longest segment of the first modifier
that has one. Exactly equal lengths are resolved by the prepared target points:
greater bed Y first, then smaller X; a complete tie keeps the first segment.
Slice coordinates already include instance rotation and have the bed axes;
centering and XY translation do not change this order. Nearly equal lengths are
not treated as equal, so exact ties occur mainly on axis-aligned geometry.
Geometrically equal segments, such as a symmetric modifier crossing both faces
of a thin wall, differ only by rounding noise that varies between layers, so
the chosen face may alternate. This is accepted deliberately: such a modifier is
ambiguous by itself, and the user is warned and expected to make it cross the
perimeter once. More than one segment raises the "multiple intersections"
warning, counted after joining across vertex zero, and the longest segment is
still used.

The selected point is inserted on its source edge. A point within 1 µm of an
existing vertex is snapped to that vertex. Helper points are added 1 µm on both
sides of it, except on an adjacent edge shorter than 2 µm, which already bounds
the distance.

When the candidates are built, the candidate at the inserted point is the only
enforced one and becomes the central enforcer; every other candidate is blocked.
Every seam position mode then picks it: Aligned and Aligned Back prefer the
central enforcer, while Back, Random and Nearest rank enforced candidates above
blocked ones. Alignment and random placement can still move the final position
along an edge, so after alignment `restore_precise_seam_positions()` writes the
exact point and its index back into every perimeter that has a strong seam.
Inner walls take their seam from the external one as usual, including
staggering.

## Weak modifiers

`collect_weak_modifier_segments()` extracts the segments of every weak modifier
before the polygon is modified, so all positions refer to the same contour. Each
segment becomes a zone with a type and two boundaries, kept in application
order, lowest priority first. The boundaries carry their positions on the source
contour; these remain as provenance after insertion and are not indices into the
modified polygon. The consumer receives only these ready boundaries; it does not
inspect modifier geometry.

`prepare_weak_modifier_segments()` then changes the polygon:

1. **Boundary insertion.** Insertion events are sorted by decreasing source edge
   and parameter, and the polygon is modified from its end towards its start. A
   pending boundary's source index therefore stays valid, and no arc lengths are
   measured for sorting. Vertex zero has the canonical position `(0, 0)` and is
   processed last, and a point on the closing edge is appended rather than
   inserted at index zero. A boundary within 1 µm of either endpoint of its
   current edge, an original vertex or a boundary inserted earlier, is snapped to
   that point, so coincident boundaries share a vertex. A zone narrower than
   1 µm collapses into a single vertex.
2. **Helper points.** A helper point is added 1 µm outside every boundary,
   unless the edge there is shorter than 2 µm, which already bounds it. Random
   placement picks a position along the edge that follows a candidate; the
   helpers keep that edge 1 µm long, so a zone can neither extend nor intrude
   further than that. Coincident boundaries share their helpers.
3. **Enforced subdivision.** Zone types are resolved for the polygon's edges in
   priority order. The edges of a zone are those from its left boundary up to,
   but not including, its right boundary. Enforced edges longer than
   `SeamPlacer::enforcer_oversampling_distance` (0.2 mm) are subdivided into
   steps of at most that length; shorter edges and existing vertices are kept.
   The seam placer then takes the middle candidate, by count, of the longest
   enforced patch, measured in candidates and across the closing edge,
   independently of where the contour starts; the same rule applies to painted
   seams. Subdivision brings that candidate near the geometric middle of a zone
   whose source vertices are sparse or evenly spaced. Densely spaced source
   vertices in one part of the zone pull it toward that part, so a weak Enforced
   zone does not guarantee an arc-length midpoint the way a strong Center
   modifier does.

When candidates are built, painting assigns their types first.
`apply_weak_modifiers_to_perimeter()` then overwrites the types of the
candidates between the boundaries of each zone, both boundaries included,
lowest priority first. Blocked and Enforced zones therefore take precedence over
painting, and Neutral clears painting inside its zone. Several weak segments per
modifier, holes and through-body intersections need no special handling.

## Numeric tolerances

Coordinates are integers in scaled units: 1 nm by default, and 10 nm when a bed
larger than 2147 mm switches `SCALING_FACTOR`. Both Precise Seam tolerances are
deliberately defined in units rather than physical distances. Clipper rounds to
whole units at any scale, so the on-edge tolerance must follow the unit. The
snapping radius following the unit keeps a similar margin over single-precision
candidate coordinates. They are centred on the object, so their step grows with
the object's size, which only large beds allow (about 0.25 µm per step 3 m from
the object's centre). Distances quoted in this document in nanometers and
micrometers assume the default unit; on large printers they are ten times
larger. The enforced subdivision step is a physical distance and stays 0.2 mm.
The design separates two scales: the rounding error of clipping, and the
distance below which points cannot be told apart after insertion.

| Value | Role |
| --- | --- |
| `MACHINE_PRECISION_SQUARED` (2.5 units², about 1.6 nm) | A point lies on an edge if it is this close. It absorbs Clipper's integer rounding of cuts and never bridges a real gap: a one-unit uncovered gap stays a gap. |
| `TOLERANCE_LINEAR` (1000 units, 1 µm) | Insertion snaps points this close to an existing vertex, and helper points are placed this far from boundaries. The same radius bounds the rounding fallback, contacts and the sub-micron full-containment rule, so those decisions match what insertion would produce anyway. |
| `enforcer_oversampling_distance` (0.2 mm) | Maximum step of enforced subdivision. |

Raising the on-edge tolerance would not help with cuts beside a vertex: more
points past a vertex would be clamped to its parameter and collapse. Lowering it
would reject ordinary rounded cuts. The snapping radius is kept far above
clipping precision for robustness: seam candidates hold single-precision
coordinates, whose step is about 8 to 15 nm at typical object coordinates, and
weak boundaries and the strong point are located among the candidates by those
coordinates, so distinct points must stay clearly distinct. 1 µm is also far
below printing precision.

Continuity never depends on two independently computed projections of the same
point: a continuing pair reuses the previous parameter, and touching fragments
are also joined at an equal integer point. Other exact comparisons of computed
parameters, such as rejecting a pair whose two points project to the same
parameter, only decide whether a pair is degenerate. The strong tie rule is the
one exact comparison whose outcome changes a policy decision: only identical
lengths count as a tie.

Boundary contacts are accepted as clipping returns them, without offsets or
special tangency rules. Users should cross the perimeter unambiguously; the
rules above only keep sub-micron results consistent and prevent a zone or a
seam from changing drastically because of rounding.

## Diagnostics and warnings

One `PreciseSeamWarnings` instance is shared by all objects and layers of a
`SeamPlacer::init()` call. After all objects are processed, `SeamPlacer::init()`
issues at most one non-critical warning with the ID `SlicingPreciseSeamWarning`.
It is a single line, "Precise Seam: <causes>. Seam placement may differ from
expected.", because the export warnings dialog shows only the first line of each
warning. Repeated warning events replace the notification instead of appending
to it. The causes are:

- **unable to process intersection:** at least one fragment was discarded by
  binding. Other segments remain usable.
- **multiple intersections with a perimeter detected:** a strong modifier had
  more than one segment on a perimeter. The longest was used.
- **perimeter is fully contained inside modifier and was ignored.**

Multiple weak segments, modifier holes and through-body intersections produce
no warning.

The log carries compact markers for investigating a saved project, not geometry
dumps:

- `[PreciseSeamIntersectionFailed]` for a discarded fragment, with object,
  modifier, layer, height, fragment and failing pair, the failure reason and
  point counts.
- `[PreciseSeamFragmentRecovered]` for a recovery, with `outcome=bound` or
  `outcome=contact`, the same location fields and the original failure reason.

Failures and recoveries have separate atomic counters. The first 10 of each per
`init()` call are logged in detail; later ones are only counted, without
formatting a message. If a limit is exceeded, one summary marker reports the
total and the number omitted. The limit matters for recoveries too: clipping is
deterministic, so a prismatic model can repeat one recovery on every layer.
Parallel processing determines which records come first, and cancellation may
omit the summaries. Callers without shared warning state, such as direct calls
in tests, log every event. The limits never affect discarding, recovery or the
user warning.

## Known limitations

The policies above already define what happens with several strong
intersections and with full containment. The remaining limitations are:

- **The modifier must reach the perimeter centerline.** Boundaries that only
  graze it, within print-setting-dependent distances, are the user's
  responsibility; sub-micron results behave like points (see
  [Numeric tolerances](#numeric-tolerances)). Several near-touches on inclined
  edges can leave several segments separated by gaps of a few units; their zones
  then cover nearly the whole perimeter instead of being skipped.
- **Self-touching perimeters.** Extraction keeps distinct visits of one
  coordinate apart through its source-edge bindings, but the consumers locate
  inserted points by coordinates. A weak zone is typed and subdivided from the
  first vertex with its boundary coordinate, while boundary helpers are added at
  every such vertex. A strong point marks every candidate at its coordinate as
  enforced, and the last one is restored after alignment. If a boundary or a
  strong point falls exactly on a repeated coordinate, a zone may therefore start
  from another visit, or the seam may start at another visit of the same point.
  Carrying visit identity through insertion, refinement, candidates and
  restoration would touch the whole pipeline, so it is not done for this rare
  geometry. Overlapping source visits are likewise outside the binding contract.

## Integration with the application

The remaining sections describe how Precise Seam volumes relate to other seam
settings, how they are stored, how their changes reach the print, and how the
user works with them.

### Other seam settings

- Precise Seam takes part only in outer and hole perimeter seam placement. In
  spiral vase mode the seam placer is not used for perimeters, so the modifiers
  have no effect.
- Scarf seams, the seam gap and wiping start from the chosen point exactly as
  they would from an ordinary seam.
- Seam painting left on a volume from before it became a Precise Seam modifier
  is ignored: the seam gizmo neither shows nor edits it.

### Model storage and 3MF compatibility

Projects must stay readable by earlier releases, and a Precise Seam volume must
not change a print there. Both 3MF writers therefore store it as an ordinary
parameter modifier: `modifier_part` in the Bambu-format part subtype, and
`ParameterModifier` together with the legacy `modifier` flag in the
Prusa-format volume metadata. The seam mode is written separately under
`precise_seam_type`, using the names from `ModelVolume::type_to_string()`
(`precise_seam_center` and so on).

On load, the mode is applied after all other volume metadata, regardless of XML
key order, and only when the base type is a modifier. A missing or unknown mode
leaves an ordinary modifier, and seam metadata on any other base type is
ignored. Files that stored the seam mode directly as the volume type still load.
A project saved again by an earlier release loses the seam mode for good: the
volumes stay ordinary modifiers without settings.

A Precise Seam volume keeps any per-volume settings it had as a part or
modifier, but they are inactive and the object list shows no settings item for
it. The writers prefix these keys with `precise_seam_config:`, so an earlier
reader drops them as unknown options and loads a modifier without settings,
which has no effect on the print. The current reader restores the keys only when
the volume ends up as a Precise Seam type, so the settings return when the user
changes the type back. Configuration values are XML-escaped in both writers, for
every volume type.

### Print invalidation

`Print::apply()` compares the Precise Seam volumes of each object by type, ID
and transformation. Adding, removing, moving, reordering or retyping one cancels
background processing and invalidates only `psGCodeExport`; the sliced layers
are kept. `model_volume_list_update_supports_and_seams()` then brings the
support and Precise Seam volumes of the print's model copy in line with the new
model in one pass. A volume may switch between these two families, since neither
affects object slicing; such a switch also changes the support volumes, so the
support step is invalidated as well. A conversion to or from a part or an
ordinary modifier changes the solid and modifier volume lists and reslices the
object as before.

### User interface

- *Add Precise Seam* in the object menu creates a Center modifier from a
  primitive or a loaded mesh. Text and SVG volumes cannot become Precise Seam
  modifiers: the menu does not offer it, and `ObjectList::set_volume_type()`
  refuses the change.
- *Change Type* has a single *Precise Seam* entry. It converts other volumes to
  Center and keeps the mode of volumes that are already Precise Seam. The
  *Precise Seam Type* submenu appears only when every selected item is a Precise
  Seam volume, including settings rows that resolve to one, and sets the chosen
  mode on all of them.
- Each mode has its own icon in the object list and its own color in the 3D
  view, at 60% opacity: warm orange, gold and dark orange for Center, Left and
  Right; green, red and gray for Enforced, Blocked and Neutral.
- Drag and drop in the object list maps visible rows to volume indices while
  skipping hidden cut connectors, and refreshes the object's row-to-volume map.
- Precise Seam volumes have no filament and cannot be pasted into SLA objects.
  Python plugins see them as `ModelVolumeType` values and through the
  `is_precise_seam*()` methods.

## Implementation and verification

- [PreciseSeam.cpp](../../src/libslic3r/GCode/PreciseSeam.cpp) implements the
  modifier cache, perimeter preparation, segment extraction and binding, strong
  selection and insertion, weak-zone preparation and application, and position
  restoration. [PreciseSeam.hpp](../../src/libslic3r/GCode/PreciseSeam.hpp)
  declares the contracts; [PreciseSeamInternal.hpp](../../src/libslic3r/GCode/PreciseSeamInternal.hpp)
  exposes the binding internals to tests.
- [SeamPlacer.cpp](../../src/libslic3r/GCode/SeamPlacer.cpp) fills the cache,
  normalizes perimeters, calls both consumers while gathering candidates,
  restores strong positions after alignment and issues the warning.
- [Model.hpp](../../src/libslic3r/Model.hpp) defines the types and their order,
  [PrintApply.cpp](../../src/libslic3r/PrintApply.cpp) handles invalidation, and
  [PrintObjectSlice.cpp](../../src/libslic3r/PrintObjectSlice.cpp) slices single
  volumes into structured regions. [bbs_3mf.cpp](../../src/libslic3r/Format/bbs_3mf.cpp)
  and [3mf.cpp](../../src/libslic3r/Format/3mf.cpp) store them.
- [GUI_Factories.cpp](../../src/slic3r/GUI/GUI_Factories.cpp) and
  [GUI_ObjectList.cpp](../../src/slic3r/GUI/GUI_ObjectList.cpp) provide the menus,
  type changes and ordering; [3DScene.cpp](../../src/slic3r/GUI/3DScene.cpp)
  defines the colors.
- [Segment extraction tests](../../tests/libslic3r/test_precise_seam.cpp) cover
  clipping and binding: holes and components, contour origin and reversal,
  repeated coordinates, collinear vertices and rounding, rollback and the
  diagnostic limits, the rounding fallback on synthetic and real Clipper
  fragments, contacts, and full containment including touches and sub-micron
  gaps on inclined edges and around vertices.
- [Precise Seam tests](../../tests/fff_print/test_precise_seam.cpp) cover the
  consumers: strong targets in every mode, including a midpoint on an existing
  vertex or the closing edge, longest-arc selection and tie order in bed axes,
  priorities, weak boundaries that coincide or share an edge, enforced
  subdivision, every warning, raft layer indexing and structured slices.
- [Seam placer tests](../../tests/fff_print/test_seam_placer.cpp) cover
  enforced-patch selection independent of the contour start, fully painted
  contours, duplicate removal, and `Print::apply()` synchronization through type
  changes and restored model snapshots.
- [3MF tests](../../tests/libslic3r/test_precise_seam_3mf.cpp) cover the round
  trip of every mode and of inactive settings, attribute escaping, and the
  metadata combinations that restore a seam mode.
  [Plugin tests](../../tests/slic3rutils/test_precise_seam_plugin.cpp) cover the
  Python bindings.
- There is no automated end-to-end test through `SeamPlacer::init()` with a
  modifier that intersects the object.
