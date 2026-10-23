# AI seam placement

AI seam exchange replaces the finalized position of external printing loops. Geometry is
exported as one UTF-8 JSON request for the AI and assignments return as one JSON response.
During the experiment an exact `<request>-raw.json` log is written beside every request. The contract
schemas are in `resources/ai_seams`. The runtime API is `Print::export_ai_seams` and
`Print::set_ai_seam_response`; an empty response path disconnects the response. These
methods require an idle slicing worker. No exchange state belongs to PrintConfig,
profiles or 3MF serialization. Loading a different model or clearing Print disconnects it.
Export requires psSlicingFinished, after contouring, simplification, infill and tool
ordering. Attaching a response validates against geometry immediately only at that state;
otherwise generation validates it against the finished registry.

## Geometry and identity

One exchange covers one explicitly selected plate. A source is a closed ExtrusionLoop
inside LayerRegion::perimeters, excluding skirts, with inset_idx == 0 or an external
perimeter path. This includes holes and fully overhanging outer walls. Collections are
traversed recursively. A source ID encodes model index, PrintObject transform group,
layer index, region index and entity-index path. An occurrence appends its PrintInstance
index. Runtime pointers only map emitted loops back to this registry; they never enter
the files. Different copies can have different assignments without changing shared
SeamPlacer Perimeter data or the standard internal-wall seam calculation.
The per-wall-filament splitter clones loops; its pass-scoped aliases explicitly retain
their original source identity and are cleared before the next layer's split storage.

Paths retain XYZ, path boundaries, extrusion roles and z_contoured. Coordinates are in
the PrintObject's centered slicing frame. The raw log keeps integer grid values with
coord_unit_mm = SCALING_FACTOR. The AI request has the same IDs and structure, but every
grid value (points, seams, targets, shift, center_offset) is converted to millimeters and
every millimeter value is rounded to precision_mm = 0.01, so the AI never rescales numbers.
Matrices are unchanged. Paths are already rotated/scaled and centered.
Layout XY = local XY + shift; plate XY = layout XY - plate_origin_mm.xy. Do not transform
paths again using the model matrices.
Path Z is distinct from slice_z_mm and print_z_mm. Markup samples are float-quality
samples from the existing classifier on each real source polygon, not continuous exact
paint intervals or candidates borrowed from a neighboring loop. Types 0/1/2 mean
blocked/neutral/enforced. Overhang ranges identify raw paths.

Strong Precise Seam metadata retains the exact Point before conversion to Vec3f,
the first successful modifier, mode and multiple-intersection flag. Weak metadata
contains the already evaluated zones in application order. A strong match skips weak
evaluation. CCW boundary direction is explicit because printing may run clockwise.
The actual seam and its distance from the strong target are per occurrence. Existing
alignment strings are exported only when all members have reliable source provenance;
these hints precede restoration of strong targets and are not final constraints.

## Digest

`sha256-canonical-json-v1` is lower-case SHA-256 of compact UTF-8 `nlohmann::json::dump()`
with sorted object keys. Arrays retain their order and numbers use the library's default
round-trip serialization. Each object digest covers object_key, model_index, transform,
center_offset, instances (instance_key, shift, model_transform), and layers (layer_index,
slice_z_mm, print_z_mm, loops). Loop geometry includes source_loop_id, kind, loop_role,
inset_idx, paths and derived overhang_ranges. Names, markup, source seams and prompt are
excluded. The digest is computed from the exact integer data, never from the rounded
request. The top-level digest covers digest_algorithm, coordinate_system, coord_unit_mm,
plate, plate_origin_mm, slicer_version and sorted object_digests. export_id is a random
trace ID and is not checked against a surviving export file. Changing serialization
semantics requires a new digest algorithm/version. Structural reorderings can reject an
otherwise identical shape conservatively.

## Application and failure

The response parser rejects invalid format, nonfinite coordinates and duplicate JSON
members/occurrence IDs. Generation compares the digest, then requires exactly one
assignment for every registered occurrence. An XY point in millimeters must be within
0.05 mm of its own exact source loop: request rounding moves path vertices by at most
about 0.007 mm, and the remainder covers the AI's own rounding. Projection uses double
arithmetic on exact segments and rounds XY to the integer grid. The selected
segment/visit survives preflight and supplies the Z interpolation parameter.
Response XY uses the object-local frame of the exported paths. A rejected point that
fits only after removing the instance shift (and plate origin) is reported as layout or
plate coordinates; it is never converted silently.

SeamPlacer receives this effective finalized point before its usual finalized/nearest
selection. Exact splitting bypasses vertex snapping and non-overhang relocation. Scarf
unsupported distance is evaluated at the imported point against the preceding layer,
with explicit first-layer initialization. Gap, scarf, wipe, speeds and subsequent G-code
coordinate quantization use the normal pipeline.
Preflight computes unsupported distance with the region's external flow width and reuses
one previous-layer distance tree per source layer, including across object instances.

A missing, duplicate, unknown or distant assignment fails closed. A geometry mismatch
uses standard seams for the whole generation, with a warning; strict CLI rejects it
before output publication. Every expected occurrence must be visited once. A late
invariant error removes the temporary G-code rather than publishing mixed assignments.

An export performs a separate baseline G-code pass over the current sliced geometry,
so nearest seams reflect nozzle order and instances. It ignores the connected response,
uses a separate processor result and temporary G-code, suppresses working-step warnings
and state completion, and restores Print statistics/options/status. JSON is published
only after all visits and G-code processing succeed, using a sibling temporary file and
rename with no in-place write fallback. The caller's preview/result is not replaced.
An export cannot overwrite its connected response, including an equivalent file path.

## Portable playbook

Both file formats optionally contain `playbook` with format `orca-ai-seams-playbook`,
version 1 and a Markdown `text`. The AI maintains model context, user wishes, visual
landmarks, accepted placement rules and rationale, Precise Seam/overhang policies,
scoped exceptions, unresolved questions and assumptions about nozzle/layer settings.
These notes live inside the existing response file; no third file is required.

The next export copies the connected response's playbook even when its geometry digest
is stale, or takes an explicitly supplied runtime Options playbook. Geometry hashes
exclude it. New layer/nozzle geometry needs fresh IDs, all new projections and a complete
response; old assignments never become valid merely because the notes are reused.
When transferring notes from a stale response through the CLI, combine import and
export without `--strict`: strict preflight rejects changed geometry before exporting.
The prompt requires revisiting invalidated assumptions and excludes local paths,
secrets, tool permissions and old ID/coordinate mappings from these notes. Playbook text
is user-supplied historical context, not executable instructions or renewed tool consent.
The backend validates its structure but does not interpret or execute its Markdown.
There is no text-size limit; keep the notes concise to avoid inflating both files.

## GUI workspace

The global Quality/Seam group contains an Advanced AI Seam checkbox. It expands
session-only working-folder, Apply AI response and Share file paths controls and a
status line. These controls are absent from object/part overrides. They are reset
when the project is replaced and are never saved in profiles or 3MF. File-path sharing
defaults off; disabling Apply preserves the answer and exports a fresh baseline request.

Each plate uses `ai-seam-plateN-request.json` and `ai-seam-plateN-response.json`, with
one-based plate numbers even in a single-plate project. Slice All uses separate pairs.
The folder must already exist. Missing/ignored answers export a request; a matching
answer is applied without rewriting the request. A geometry mismatch exports a new
request carrying the valid answer's playbook and generates standard seams with a warning.
Malformed answers fail the slice. Turning Apply off recovers without reading that answer.

UI setting changes invalidate G-code on every plate. Each explicit Slice with AI Seam
enabled also invalidates G-code before starting the worker, so external file replacements
are observed without scene edits. Slice therefore stays enabled after a valid result. The worker uses one temporary snapshot of the answer
for validation and generation, deletes it on every exit, and publishes the final status
only after processing succeeds. It never edits the public response. No network access,
file watcher or AI execution is introduced. Files are manually exchanged with the AI.

## CLI and automation

`--export-seam-loops file.json` and `--import-seam-points response.json` require an
explicit `--slice N`, N > 0. Cached `--load-slicedata` cannot recover outer-wall inset
provenance and is rejected. Spiral vase mode has no supported normal seam and is rejected.
Both exchange options may be combined: the export remains standard, the actual output
uses the response. Geometry and instances must match between GUI and CLI, including
the internal coordinate scale. `--strict` is recommended for automation.

`--share-seam-paths` opts into local executable/resources/project/exchange paths; it
defaults off. Windows DLL location is informational, launching uses the EXE. The embedded
prompt asks separately for model context, overhang policy, Precise Seam preservation,
Python consent and local-operation consent. Tool permissions cannot be inferred from
sharing paths or answering OK to the model-context question. Automation works on a copy
of a saved project and can use ordinary CLI 3MF export with G-code. The separate response
must be reconnected for future reslicing; the project receives no new seam metadata.
