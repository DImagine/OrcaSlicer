# Precise Seam — High Level Design

## Purpose and scope

Precise Seam places the seam where a helper volume intersects the external
wall. The user attaches a mesh to an object as a Precise Seam modifier, and on
every layer the seam placer reads the modifier's slice to decide where the seam
of each external perimeter may, must or must not go. The same mesh keeps
working after the model changes, so the seam does not have to be repainted
after every design revision, and a swept helper body can guide the seam along
any path.

The modifier is non-printing geometry. It does not take part in slicing, region
assignment, filament selection or brim adhesion. It affects only seam
placement, which runs during G-code export.

## Volume types and priority

Precise Seam adds six `ModelVolumeType` values after `SUPPORT_ENFORCER`. The
strong types come first and the weak types follow. `is_precise_seam()`,
`is_precise_seam_strong()` and `is_precise_seam_weak()` are range checks that
depend on this order.

| Type | Group | Effect on the perimeter |
| --- | --- | --- |
| `PRECISE_SEAM_CENTER` | strong | seam at the arc-length midpoint of the intersection |
| `PRECISE_SEAM_LEFT` | strong | seam at the first point of the intersection |
| `PRECISE_SEAM_RIGHT` | strong | seam at the last point of the intersection |
| `PRECISE_SEAM_ENFORCED` | weak | intersection marked as enforced |
| `PRECISE_SEAM_BLOCKED` | weak | intersection marked as blocked |
| `PRECISE_SEAM_NEUTRAL` | weak | intersection reset to neutral |

A strong modifier fixes one point. A weak modifier only changes the
enforced/blocked type of seam candidates, and the configured seam position then
chooses among them. First and last are taken along the perimeter made
counter-clockwise seen from above. On an outer wall seen from outside, Left is
the left end of the intersection. On the wall of a hole seen from inside the
hole, the two ends are swapped.

The order of volumes in the object is the priority order, highest first.
`ModelObject::sort_volumes()` keeps every strong modifier before every weak one
and preserves the user's order within each group. The object list lets the user
drag a modifier only within its own group. A type change that crosses a group
boundary moves the volume to the end of its new group, where it has the lowest
priority. Strong modifiers are tried in this order, and the first one that
yields a seam on a perimeter wins. Weak modifiers are applied from the lowest
priority to the highest, so the highest one overwrites any overlapping zone.

## Model storage and 3MF compatibility

Projects must stay readable by earlier releases, and the modifier must not
change a print there. Both 3MF writers therefore store a Precise Seam volume as
an ordinary parameter modifier: `modifier_part` in the Bambu-format part
subtype, and `ParameterModifier` together with the legacy `modifier` flag in
the Prusa-format volume metadata. The seam mode is written separately under
`precise_seam_type`, using the names from `ModelVolume::type_to_string()`
(`precise_seam_center` and so on).

On load, the mode applies after all other volume metadata, regardless of XML
key order, and only when the base type is a modifier. Missing or unknown modes
leave an ordinary modifier. Seam metadata on any other base type is ignored.
Files that stored the seam mode directly as the volume type still load.

A Precise Seam volume keeps any per-volume settings it had as a part or
modifier, but they are inactive and the object list shows no settings item for
it. The writers prefix these keys with `precise_seam_config:`, so an earlier
reader drops them as unknown options. The volume therefore loads there as a
modifier without settings and has no effect on the print. The current reader
restores the keys only when the volume ends up as a Precise Seam type, so the
settings return when the user changes the type back. Configuration values are
XML-escaped in both writers, for every volume type.

## Print invalidation

`Print::apply()` compares the Precise Seam volumes of each object by type, ID
and transformation. Adding, removing, moving, reordering or retyping one
cancels background processing and invalidates only `psGCodeExport`; the sliced
layers are kept. `model_volume_list_update_supports_and_seams()` then brings
the support and Precise Seam volumes of the print's model copy in line with the
new model in one pass. A volume may switch between the two families, since
neither affects slicing. A conversion to or from a part or ordinary modifier
changes the solid and modifier volume lists and reslices as before.

## Modifier slices

`SeamPlacer::init()` collects the Precise Seam volumes of each object once:
strong ones in priority order and weak ones reversed. It slices each volume
separately with `PrintObject::slice_single_volume_regions()`, retaining each
region's outer contour and holes as an `ExPolygon`. `prepare_modifier_slices()`
moves these slices into the shared cache, pairing each region with the bounding
box of its exterior. Bounds are computed once before parallel perimeter processing,
not on each extraction call. Empty layers keep their original indices. Both
modifier kinds share this structured cache. Volumes are not merged, so each keeps
its own priority. The result is cached per volume
and indexed by object layer; `Layer::id()` includes raft layers, which are
subtracted. Seam candidates are then gathered in parallel over the layers and
read the cache without locking.

Objects without Precise Seam volumes follow the unchanged seam placement path.
For objects that have them, perimeter extraction also removes consecutive
duplicate points and the repeated closing point of each extrusion loop.
Zero-length edges at path junctions would otherwise prevent point insertion
there. Distinct visits to one point of a self-touching contour are kept.

## Finding perimeter segments

The seam placer works on external perimeter loops, including holes, normalized
to counter-clockwise traversal. Before trying modifiers, it creates one
`PreparedPerimeter`: validation of consecutive/closing duplicates, the perimeter
bounding box, and an open clipping line with the first point repeated at the end.
All strong queries reuse this preparation. If no strong point was inserted, weak
queries reuse it too; weak collects all boundaries before inserting them. A
successful strong insertion ends modifier processing. The preparation borrows the
polygon and must not be used after insertion/refinement changes it. Objects
without Precise Seam modifiers do not perform this preparation.

Extraction clips the prepared perimeter line against each
nearby region of one modifier, keeping its holes attached and applying bounding
box rejection per region using its cached bounds. Geometry and bounds remain
read-only throughout the seam pass. Both strong and weak consume these ready segments.

The fast binding path anchors on the second fragment point and matches interior
vertices exactly, trying either direction and later occurrences of the anchor.
Only cut endpoints need projection. Two-point fragments go straight to the
projection path, which accepts the first matching source edge. Subsequent pairs
must continue on that edge or its neighbor in the established direction. Failure
rolls back and discards only that fragment, with a diagnostic marker. Overlapping
source visits are outside the binding contract. Boundary contacts are accepted
as returned by clipping, without offsets or additional contact rules.

Each segment retains source edge indices and endpoint positions. Neighboring
intervals are joined, including across contour vertex zero. The extractor takes
the modifier mode explicitly. Full containment is detected before strong data
preparation: both consumers skip it, so no lengths or targets are calculated.
Otherwise, strong segments receive their mode point with its source edge index;
arc lengths are measured for Center or comparison of multiple segments. Weak
skips both length and target preparation; its boundary insertion and enforcement
subdivision still perform their own necessary measurements. All modes retain the
complete segment geometry and source bindings, including for debugging.
Insertion snaps points within 1 micrometre of an existing vertex to that vertex.

## Strong modifiers

Strong modifiers are considered in priority order. After full containment is
skipped, each extracted segment already contains its mode point: first point for Left, last point for Right, or half the
arc length for Center. Center temporarily retains edge lengths while measuring
one segment, then uses them and source edge indices to locate the midpoint. This
local array is discarded before processing the next segment; Left and Right do
not fill it. A single Left/Right segment needs no length comparison, so its length
remains zero (unmeasured); its endpoint is still prepared. Selection consumes the
ready points before any perimeter insertion.

The longest segment of that modifier wins, using the sum of Euclidean lengths,
not its chord or vertex count. Exactly equal lengths are resolved by the mode
point: largest bed Y first, then smallest X; a complete tie keeps the first
candidate. Slice coordinates already include instance rotation and have the
bed axes; centering and XY translation do not change this ordering. Nearly equal
lengths are not treated as equal.

Once a modifier yields a segment, no later modifier competes with it. Its chosen
point is inserted with helper points 1 micrometre on either side, after first
and before second to preserve the insertion index. No weak modifiers are then
processed for that perimeter. Multiple extracted segments produce a warning
only for strong modifiers. Full containment retains its separate skip policy.

When candidates are built, the inserted point is the only enforced candidate
and becomes the central enforcer; every other candidate is blocked. The seam
position modes then pick that point: Aligned and Aligned Back prefer the central
enforcer, while Back, Random and Nearest rank enforced candidates above blocked
ones. Alignment and random placement can still move the final position along an
edge. After alignment, `restore_precise_seam_positions()` writes the exact point
and its index back into every perimeter that has a strong seam. Inner walls take
their seam from the external seam as usual, including staggering.

## Weak modifiers

Weak extraction clips the perimeter line against the cached regions of one
modifier, rejecting each region by its bounding box first. It returns all ready
segments with their endpoints and canonical source positions (edge index and
parameter). These positions refer to the original perimeter and are retained as
provenance after insertion/refinement; they are not indices into the modified
polygon. The boundary consumer
does not inspect modifier geometry: it receives the collected segments in
modifier priority order. All boundaries are collected before the perimeter is
modified, so their source edge indices refer to the same contour.

`prepare_weak_modifier_segments()` sorts insertion events by decreasing source
edge index, then decreasing parameter on that edge. This gives arc-length order
without measuring any lengths for sorting. Vertex zero uses the canonical
position `(0, 0)` and is processed last. The segment array keeps modifier priority
order; coincident boundaries continue to share vertices. Each insertion leaves
pending source indices unchanged; a point on the closing edge is appended rather
than inserted at index zero. A helper point is added 1 µm
outside each boundary. Random placement picks a position along the edge that
follows a candidate. These helpers keep that edge 1 µm long at each boundary, so
a zone cannot extend or intrude further than that. Boundaries that coincide
share their helper points.

The zone types are then resolved in priority order, and the edges of enforced
zones are subdivided into steps of at most
`SeamPlacer::enforcer_oversampling_distance` (0.2 mm). The middle candidate of
the longest enforced patch is therefore close to the geometric middle of the
zone. That patch is measured in candidates, across the closing edge, regardless
of where the contour starts; the same rule applies to painted seams.

Candidates first receive their type from seam painting. The weak zones then
overwrite it, lowest priority first. Blocked and Enforced zones therefore take
precedence over painting, and Neutral clears painting inside its zone.

## Unsupported geometry and warnings

- More than one segment for a strong modifier produces a warning; the selected
  longest segment is still used. The count is taken after joining across vertex zero.
- A perimeter fully contained in a modifier is ignored by that modifier.
- An intersection that cannot be bound continuously is discarded with an
  "unable to process intersection" warning and a compact log marker containing
  the object, modifier, layer and failure location. Other segments remain usable.
- Through-body intersections and modifier holes need no separate warnings.
  Multiple weak segments are accepted without a warning.

Atomic warning flags and a failed-fragment counter are shared by all layers and
objects in one `SeamPlacer::init()` call. The first 10 failed fragments receive
detailed log markers; later failures are only counted, without formatting a
message. If the limit is exceeded, one final marker reports the total and omitted
counts. Parallel processing determines which failures are logged first. A new
call starts with a fresh counter; cancellation may omit the final summary.
Callers without shared warning state retain unlimited diagnostic logging.
The limit does not affect discarding fragments or showing the user warning.
After all objects are processed, `SeamPlacer::init()` issues at most one non-critical
warning with the ID `SlicingPreciseSeamWarning`. The warning is a single line
that lists every cause found, because the export warnings dialog shows only the
first line of each warning. Repeated warning events replace this notification
instead of appending text to it.

## User interface

- *Add Precise Seam* in the object menu creates a Center modifier from a
  primitive or a loaded mesh. Text and SVG volumes cannot become Precise Seam
  modifiers: the menu does not offer them, and `ObjectList::set_volume_type()`
  refuses the change.
- *Change Type* has a single *Precise Seam* entry. It converts other volumes to
  Center and keeps the mode of volumes that are already Precise Seam. The
  *Precise Seam Type* submenu appears only when every selected item is a
  Precise Seam volume, including settings rows that resolve to one. It sets the
  chosen mode on all selected volumes.
- Each mode has its own icon in the object list and its own color in the 3D
  view, at 60% opacity: warm oranges for the strong modes, and green, red and
  gray for Enforced, Blocked and Neutral.
- Object list drops map visible rows to volume indices while skipping hidden
  cut connectors, and they refresh the row-to-volume map of the object.
- Precise Seam volumes have no filament, block pasting into SLA, and are exposed
  to Python plugins as `ModelVolumeType` values plus the `is_precise_seam*()`
  methods.

## Implementation and verification

- [PreciseSeam.cpp](../../src/libslic3r/GCode/PreciseSeam.cpp) implements segment
  detection, point insertion, weak-zone resolution and position restoration.
  [SeamPlacer.cpp](../../src/libslic3r/GCode/SeamPlacer.cpp) integrates it into
  candidate gathering and issues the warning.
- [Model.hpp](../../src/libslic3r/Model.hpp) defines the types and their order,
  [PrintApply.cpp](../../src/libslic3r/PrintApply.cpp) handles invalidation, and
  [PrintObjectSlice.cpp](../../src/libslic3r/PrintObjectSlice.cpp) slices the
  modifiers. [bbs_3mf.cpp](../../src/libslic3r/Format/bbs_3mf.cpp) and
  [3mf.cpp](../../src/libslic3r/Format/3mf.cpp) store them.
- [GUI_Factories.cpp](../../src/slic3r/GUI/GUI_Factories.cpp) and
  [GUI_ObjectList.cpp](../../src/slic3r/GUI/GUI_ObjectList.cpp) provide the menus,
  type changes and ordering.
- [Precise Seam tests](../../tests/fff_print/test_precise_seam.cpp) cover the
  strong positions, including a midpoint on an existing vertex or the closing
  edge. They also cover shared and coincident weak boundaries, every warning,
  and the priority order.
- [Seam placer tests](../../tests/fff_print/test_seam_placer.cpp) cover
  enforced-patch selection independent of the contour start, fully painted
  contours, duplicate removal, and `Print::apply()` synchronization through
  type changes and restored model snapshots.
- [3MF tests](../../tests/libslic3r/test_precise_seam_3mf.cpp) cover the round
  trip of every mode and of inactive settings, attribute escaping, and which
  metadata combinations restore a seam mode.
  [Plugin tests](../../tests/slic3rutils/test_precise_seam_plugin.cpp) cover the
  Python bindings.
