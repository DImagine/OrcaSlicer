# Outer wall passes — High Level Design

## Purpose and scope

`split_outer_wall` prints the outer wall of every layer as a stack of `n` thin
loops, the passes, each `h / n` high (`split_outer_wall_count`, layer height
`h`). The inner walls and the infill are printed afterwards at the full layer
height. The outside of the print then shows the fine steps of a thin layer
while the inside prints at the speed of a thick one.

The inner edge of the stack is vertical and common to all the passes. The outer
edge of each pass follows the outline of the object at the middle height of that
pass, so the passes differ in width and reproduce the slope of the surface in
steps of the pass height.

A layer is printed this way only if the widths of all its passes stay within the
range `split_outer_wall_min_width` .. `split_outer_wall_max_width` everywhere.
A layer that does not fit is either split into `n` regular thin layers or printed
the regular way (`split_outer_wall_fallback`).

Only the classic wall generator is supported. The option is off with the spiral
vase mode, a raft, fuzzy skin, Z contouring, objects painted with several
filaments and belt printers; `Print::validate()` warns about it.

## Geometry

For a layer between `lo` and `hi`:

- `S_k` is the outline of the object at `lo + (k + 0.5) * h / n`, the middle of
  pass `k`, and `S_mid` the outline at the middle of the layer. With an odd `n`
  the middle pass is sliced at the middle of the layer.
- `W_nom = (W_min + W_max) / 2` is the nominal width. The wall generator works
  with an outer wall of this width, so the inner edge of the outer wall is where
  a regular wall of the nominal width would end.
- `d_k` is the offset of `S_k` from `S_mid`, positive outwards. Pass `k` is
  `W_nom + d_k` wide and its centerline is shifted by `d_k / 2`.

The layer fits if `-(W_nom - W_min) <= d_k <= W_max - W_nom` holds along the
whole outline for every pass (`outer_wall_passes_fit()`). A part appearing or
vanishing between the passes does not fit. With `D = W_max - W_min` the steepest
slope that fits is `atan((n - 1) * (h / n) / D)` from the horizontal, about 21°
for a 0.2 mm layer in two passes and the default range of a 0.4 mm nozzle.

The outlines are sliced without the XY size compensations. Only the offsets
`d_k` are taken from them and applied to the outer loop the wall generator
produced, which carries the compensations.

## Planning the layers

`PrintObject::plan_outer_wall_passes()` runs in `PrintObject::slice()` between
generating the object layers and slicing the volumes. It slices the outlines of
every candidate layer, applies the criterion and returns the object layers, with
the layers that do not fit split in `n` if configured so. A layer printed in
passes keeps its slices in `Layer::outer_wall_pass_slices`.

The first layer, the layers shrunk by the elephant foot compensation and the
layers whose pass would be thinner than the minimum layer height of the printer
are not candidates.

The layer itself is sliced as usual, at its middle. Everything but the outer
wall is generated from that slice.

## Building the passes

On a layer with passes `LayerRegion::flow()` returns the nominal width for the
outer wall, and the wall generator prints the outer wall before the inner walls
whatever the wall sequence is, as the passes go first. The walls of such a layer
are printed before its infill as well, even with the infill set to go first.

`split_outer_wall_loops()` then replaces every loop of the outer wall
(`OuterWallPasses.cpp`): the loop is sampled, the offsets `d_k` are measured at
the samples as differences of signed distances to `S_k` and `S_mid`, and each
pass is built as a variable width loop. The top pass takes the place of the
original loop in `LayerRegion::perimeters`. The passes below it are kept in
`LayerRegion::outer_wall_lower_passes`, bottom pass first.

A loop is left as it is if it has an overhanging part, if it is held back until
the infill is printed, or if it runs around a hole and `split_outer_wall_holes`
is off. A layer stepping out over the layer below has vertical walls and passes
the criterion, but its loop hangs in the air and must not be thinned.

## G-code

`GCode::extrude_outer_wall_lower_passes()` prints the lower passes right after
the tool change of the filament of the outer wall, before anything else of the
layer. The passes of all the instances are printed lowest first, so the nozzle
never works beside a part of an object printed higher on the same layer. The nozzle
descends on the travel to a pass the way it does at a layer change and returns
to the layer height afterwards. The Z offset of the printer is kept.

All the passes of an outer wall start at the same point: the seam is placed once
and looked up by the first point of the stored loop. A scarf joint is not used on
the passes, as it would start inside the pass below.

The overhang speed and cooling of a pass are measured against what it leans on:
the pass below it, and for the lowest pass the top pass of the layer below, or
that layer as usual if it has no passes
(`ExtrusionQualityEstimator::set_custom_boundaries()`).

## Tool order, supports and the prime tower

On a layer with passes the filament of the split outer walls is printed first
(`LayerTools::outer_wall_passes_region`, the wrapper of `get_custom_seq` in
`ToolOrdering::reorder_extruders_for_minimum_flush_volume()` and
`enforce_outer_wall_passes_order()`). A support interface of another filament
would otherwise stand at the full layer height beside the passes. The walls of
such a layer are not offered for wiping into the object.

The prime tower is printed at the layer height, before the passes or after the
objects. It stands apart from the objects, so it needs no change. Support layers
follow the object layers: the independent support layer height is switched off
with the option on for the whole print and no global setting ruling the passes
out, as it is with a prime tower.

## Restrictions

`Print::validate()` stops with an error for

- a pass thinner than the minimum layer height of the printer,
- a width range that is empty, makes a pass narrower than it is high, or reaches
  beyond the widest line allowed for the nozzle,
- mixed filaments, which are printed in sub layers by a path of their own,
- a skirt around all the objects higher than one layer and a draft shield, which
  are printed at the layer height before the lower passes; a per-object skirt is
  printed with its object, after them,
- support of an independent layer height when the option is set for an object only,
- Organic supports together with layers split into thin layers, as Organic
  supports are built of layers of one height,
- several objects whose layers do not end at the same heights (different layer
  heights, a layer height profile, a raft, the precise Z height on objects of
  different heights), and more than two passes together with layers split into
  thin layers for several objects: layers of another object printed earlier
  would stand above the lower passes. A layer of another object may end no
  higher than the lowest pass of the layer it ends in. Objects printed one by
  one are not restricted,
- outer walls of several filaments among the objects the option is on for.

## Limitations

- A loop with an overhanging part is printed as one regular loop.
- With `split_outer_wall_holes` off the walls of holes are printed in one pass at
  the layer height, but still with the nominal width.
- At corners sharper than about 41° the miter of the offset is limited and the
  outer edge of a pass deviates from the outline by up to half the offset.
- The passes are measured against outlines without the XY size compensations.
- The criterion offsets the middle outline with a miter limit, so a sloped layer
  with a convex corner sharper than about 30° to 40° may be taken for not
  fitting. This only sends the layer to the fallback.

## Tests

`tests/libslic3r/test_outer_wall_passes.cpp` covers the criterion,
`tests/fff_print/test_outer_wall_passes.cpp` the layers, the widths and heights
of the passes and the order of the G-code.
