# Changelog

All notable changes to Cellucid will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.9.1]
<!-- CELLUCID_VERSION -->

Version 0.9.1 is the PyPI submission release.

### Fixed — the two writers now accept exactly the same numeric columns

The float32 range check covered overflow only. Underflow is the quiet end of
the same range: a nonzero value below `2^-149` converts to exactly `0.0`, which
is finite, so it passed. Measured on eight distinct values near `1e-320`, the
export wrote one distinct value — the field published with its every value
replaced by zero. The R writer already refused these inputs, so one export
format had two meanings depending on which package produced it. Python now
refuses them too, at both the embedding and the continuous-obs entry points.
Values `float32` can represent, including genuine subnormals, are unaffected.

### Fixed — a label the browser cannot draw is refused

The rule that keeps a stored label and a drawn label identical covered control
characters, zero-width characters and edge whitespace, but not an unpaired
surrogate. Python holds one, `json.dumps` escapes it, and the browser's
`JSON.parse` yields an unpaired code unit that draws as the replacement
character — the same stored-is-not-drawn defect the rule exists to prevent, and
the only one of its family that reached publication. It is now refused by name,
with the offending code point quoted. R strings cannot hold an unpaired
surrogate, so its writer has nothing to check.

### Fixed — an export cannot publish coordinates the viewer is not told to fetch

`prepare()` reconciled each payload directory against the manifest that
describes it — `obs/`, `var/`, `connectivity/`, `vectors/` — but never
reconciled the export root, and the point payloads live there. They are the only
artifact this format declares by path from the root, in
`dataset_identity.json` under `embeddings.files`, so no axis manifest spoke for
them: the writer captured the path `_write_binary()` actually produced and used
it in a console line only, never comparing it with the name the identity file
declared. Their declared name and their written name were two independent
expressions of the `compression` setting, so a generation that wrote
`points_2d.bin` while declaring `points_2d.bin.gz` published successfully and
then failed in the browser with no coordinates.

The root is now reconciled whole, immediately before the transaction publishes,
which closes the larger sibling as well: an export missing `obs_manifest.json`,
or carrying a stray scratch file or an undeclared directory at its root, used to
publish too. After a successful export the root holds exactly
`dataset_identity.json`, `obs_manifest.json`, the declared point payloads,
whichever of `var_manifest.json` and `connectivity_manifest.json` was written,
and the payload directories the export created — nothing else. Every refusal
rolls the whole generation back, so no partial export reaches the target
directory.

The `compression` setting also decided a payload's name in three different
spellings — `".gz" if compression else ""` in the manifests, `compression and
compression > 0` at the write, and `compression is not None` inside
`_write_binary()`. They agreed only because the validator rejects `0` a thousand
lines away. One helper, `_final_binary_path()`, now names the file for both the
declaration and the write.

`cellucid_prepare()` in the R package had the identical gap and closes it in the
same release, with the same five reconciled surfaces and the same refusals.

### Fixed — a manifest is now checked as bytes, not only as a payload

Reconciliation proved the four manifests *exist* and describe the payloads that
were written. It did not prove their bytes parse to what was intended. Every one
was validated from the in-memory payload and then serialized, and nothing
re-read the result — so any defect in serialization published cleanly and failed
only when the web app opened the export.

`prepare()` now reads each manifest back out of the export transaction, before
publication, and requires it to parse to exactly the payload that was validated:
same kinds, same keys in the same order, same lengths, same values. The parser
is as strict as the reader's, which matters more than it sounds — `json.dumps()`
writes a non-finite float as the bare token `NaN` or `Infinity`, `json.loads()`
reads both back without complaint, and `JSON.parse()` in the browser refuses
them outright. A read-back that used the default parser would have agreed with
the writer and still handed the browser a file it cannot open. A manifest that
does not parse, is not valid UTF-8, or disagrees with its payload on any node
now fails the export and names the node.

`cellucid_prepare()` in the R package gains the same guarantee and the same
message shape in the same release. It also had a real instance of the defect
this closes: `jsonlite` renders an empty unnamed list as `[]`, so an export
carrying no observation field at all — an `obs` with no columns, or an empty
`obs_keys` — published `"_obsSchemas":[]` where this writer published
`"_obsSchemas": {}`. The viewer requires an object there and refused the R
export. Both writers now assert the kind of `_obsSchemas` in the manifest
validator itself, so neither can emit the other's shape.

### Added — a web build can ask what this notebook's protocol accepts

Both servers now answer `GET /_cellucid/protocol` with the wire capabilities the
installed `cellucid` accepts:

```json
{
  "commands": ["clearHighlights", "debug_snapshot", "freeze", "highlight",
               "ping", "requestSessionBundle", "resetCamera", "setColorBy",
               "setVisibility"],
  "events": ["click", "command_error", "debug_snapshot", "hover", "pong",
             "ready", "selection", "session_bundle"]
}
```

The Python package fetches the web UI from `https://www.cellucid.com` without
pinning a build, so a viewer newer than the notebook driving it is the normal
state the day a build deploys. A viewer that emits an event this installation
does not accept does not degrade to silence: `_require_inbound_jupyter_event`
raises, the event endpoint logs a traceback into the kernel and answers
`500 Viewer callback failed`, and the viewer then reports a delivery failure.
So the viewer has to ask first, and nothing already in flight could carry the
answer — the iframe configuration, the health record and every command
acknowledgement are validated in the browser against exact key sets, so a field
added to any of them breaks builds that are already deployed. A route an older
build never requests is the only channel that stays quiet, and because it ships
in the same release as the validator, "the endpoint answers" is exactly
equivalent to "the validator accepts".

`events` is every inbound event type the validator accepts; `commands` is every
command this installation may send, which is also exactly the set of names a
`command_error` may report as having failed, so a viewer whose own command list
has grown ahead of the notebook can avoid naming one it would reject. Both are
lists that grow and the key set stays fixed — a later capability arrives inside
an array rather than as a new key, which would break the very builds the route
exists to protect. There is deliberately no version number: a release string
would be a proxy for the capability rather than evidence of it.

The route is answered from the same declaration the validator reads, so it
cannot announce something the notebook would then reject. It needs no
credentials, exactly like `/_cellucid/health` and `/_cellucid/info`: the body is
a property of the installed version, identical for every viewer and every
dataset, and it creates no state.

No web build fetches this yet, so no `command_error` is emitted and the reader
added below stays idle. This is the precondition for the emitter, not the
emitter.

Also corrects a security page that said events are "routed by `viewerId` but are
not authenticated beyond that". They have been token-authenticated with
`hmac.compare_digest` all along, which the endpoint reference on the neighbouring
page already said.

### Added — a notebook learns when a command it sent was refused

`viewer.set_color_by("CD8A")` for a field the viewer does not have returned
`None` and looked exactly like success. The commands Python sends into the
iframe — `highlight`, `setColorBy`, `setVisibility`, `clearHighlights`,
`resetCamera`, `freeze` — carry no request id, the Python side returned without
waiting, and the protocol defined no inbound error or acknowledgement event, so
there was no channel on which a failure could travel. In the viewer the
rejection reached a `console.error` inside an embedded iframe; in the notebook
it reached nothing at all, and a failed command was indistinguishable from one
that applied and matched no cells.

The protocol now defines a `command_error` event: `command`, which must name
one of the nine command types this protocol defines, and `reason`, one line of
at most 500 characters. Restricting `command` to that closed list is what keeps
browser-supplied text out of the notebook's failure notice.

It reaches a notebook the same way every other event does — `viewer.state`,
`viewer.wait_for_event("command_error")`, `@viewer.on_message` — plus a
dedicated `@viewer.on_command_error` decorator, and
`viewer.debug_connection()` now reports the recent rejections under
`recent_events["command_errors"]`. When no `command_error` hook is registered,
the failure is printed to `sys.stderr` naming the Python call that sent it:

```text
[cellucid] viewer.set_color_by(...) did not take effect in the viewer: ...
```

Registering a hook takes reporting over and silences that notice. No command
method starts raising: the rejection arrives on the server thread after the
cell has finished, so making the call itself fail would mean blocking every
command on a browser round trip, and notebooks that ignore command return
values behave exactly as before.

This release contains the reader only. The web build has to know that the
notebook driving it accepts the event before it may emit one, because
`_require_inbound_jupyter_event` raises on an unknown type, the event endpoint
turns that into a logged traceback and `500 Viewer callback failed`, and the
viewer then reports a delivery failure — louder and less accurate than the
silence being fixed. The Python package fetches the web UI from
`https://www.cellucid.com`, so an older notebook can be driving a newer web
build, and every payload Python sends the viewer is validated there against an
exact key set, which leaves no field free to carry a capability to already
deployed builds. The channel that gates the emitter is `GET /_cellucid/protocol`
above, which ships in this same release; the emitter itself follows, against a
released reader.

### Fixed — `cellucid serve` stops answering a wrong flag with a Python traceback

Forgetting `--dataset-name` printed a traceback through `cli.py` before the line
that said what to do, and so did a port another program was already holding.
Every failure went through one `logger.exception("Unexpected error")`, which
made no distinction between a mistake the person at the terminal can correct and
a defect inside Cellucid. The audience for `cellucid serve` is largely wet-lab,
and a traceback tells that reader nothing they can act on — it also printed the
machine's home directory into whatever they pasted asking for help.

`main()` now classifies the failure. A condition the operator owns is reported
as one `Error:` line naming what failed and what to change, with no traceback:

```text
$ cellucid serve lung_atlas.h5ad --no-browser
cellucid v0.9.1
Error: --dataset-name and --dataset-id are required when serving h5ad or zarr data. An .h5ad file and a Zarr store carry no Cellucid identity of their own, so name the dataset on the command line, for example: --dataset-name "My study" --dataset-id my-study-v1
```

That covers the whole class, not the two reported cases: an inapplicable flag, an
input that is not a dataset, a manifest whose JSON is the wrong shape, a named
`obs` column that is not in the data, a taken port, a host this machine cannot
bind, a port below 1024, an optional package such as `zarr` that is not
installed, and a machine with no browser to open. The operating-system
conditions arrive with an `errno` and no advice, so the CLI supplies the flag
that resolves each one — `--port 0` for a taken port, `--host 127.0.0.1` for an
unbindable host.

Both flags of a missing identity are now named in one message. Reporting only
`--dataset-name` sent someone who supplied neither back to the terminal twice
for one mistake.

A defect keeps its traceback, printed to stderr unconditionally — even under
`--quiet`, and no longer through a logger that `--quiet` left unconfigured — and
the `Error:` line says it is a Cellucid bug and where to report it. `--verbose`
adds the traceback to an operator condition too, so a misclassification is still
reportable. Nothing is discarded.

`stdout` is flushed before anything reaches `stderr`, so the `Error:` line is
genuinely last when the command is piped into a file or a pager; previously it
landed above the banner it should have followed. The message is collapsed to one
line, which is what the documented CLI contract always said it was.

The server tutorial's two failure screenshots showed the removed behaviour and
printed a local home directory into a committed image; they are replaced by the
transcripts the command actually produces, which a contract test now compares
against the command.

### Fixed — `debug_connection()` reports the dataset probes for a published sample state

A prepared export that ships a published default session adds `state_manifest`
and `state_sha256` to its `/_cellucid/datasets` entry — the same pair the
official catalog advertises, and the pair the server itself reads back to decide
which sidecars that dataset may serve. The diagnostics validator still accepted
only `id`, `path` and `name`, so serving such an export made
`viewer.debug_connection()` record

```text
server_datasets_error: "Dataset-list entry 0 must contain exactly id, name, and path"
```

and then skip `dataset_identity_probes` entirely — the part of the report that
proves the server is serving the dataset you think it is, missing exactly for the
datasets that ship a reviewed starting view.

The server was right and the validator's key set was stale: it was written when
an entry had three keys and was never revisited when the published-state pair
was added. The validator now accepts the pair, and only as a complete pair, with
the manifest name matched exactly and the digest required to be 64 lowercase
hexadecimal characters; every other key set is still refused, so no security
boundary was widened to make a diagnostic work. The two declarations of the
entry shape — one in `jupyter/_wire.py`, one in `server/_state.py` — are now
tied together by a test, which is what was missing when they drifted apart.

`server_datasets` carries the pair through into the report, so the diagnostic
also shows which datasets publish a starting state.

### Fixed — every published payload is little-endian, whatever the host is

The export format's dtype names — `float32`, `uint16`, `float64` — carry no
byte-order component, and the web app builds its typed arrays straight from the
bytes it receives, which is host order by definition in JavaScript. There is
therefore no field a reader can check and no error it can raise. `cellucid-r`
already named `endian = "little"` at every `writeBin`; the Python writer emitted
whatever the exporting machine used. On a big-endian host — s390x is the
realistic case — `prepare()` and the AnnData server both produced payloads the
viewer misread as different coordinates, different expression values, and
different edge weights, with no symptom other than a plot that looked entirely
plausible. The defect was unreachable on x86 and ARM.

Byte order is now pinned in `src/cellucid/_byte_order.py` and applied at both
points where this package turns an array into payload bytes: the export writer
(`points_*`, `obs/`, `var/`, `connectivity/`, `vectors/`) and every binary
response the direct AnnData server serves. The conversion is skipped when the
array already holds the published order, so nothing is copied on a
little-endian host and no export produced on one changes by a byte.

The specification now states the rule once and normatively, for every payload
the format defines, instead of only for the connectivity index arrays.

### Fixed — `prepare()` refuses a directory nobody set aside for an export

`prepare()` publishes by renaming the entire existing `out_dir` aside and then
removing it, so `force=True` destroys every file that directory holds. The only
guard was a leaf-name check, which rejects `/` and `.` but accepts the working
directory named by its absolute path and accepts `~`, whose leaf is the user
name. `prepare(out_dir=os.getcwd(), force=True)` therefore deleted unrelated
user files, and `prepare(out_dir="~", force=True)` would have replaced the whole
home directory.

`out_dir` is now validated where it is first resolved, before any path is
created, locked, written, or removed. Refused with a `ValueError` naming the
path and the fix:

- the filesystem root, on any platform;
- the current working directory, whether written `.`, `./`, or in full;
- the home directory, whether written `~` or in full;
- the directory that holds every home, whether written `..` or in full.

A leading `~` and every symbolic link are resolved before the comparison, so an
alias to one of those directories is refused too. `out_dir` must also be a
string or `os.PathLike`, reported as such instead of as a `pathlib` internal.

Naming a dedicated child directory is unchanged, with and without `force`, and
the protected set matches `.validate_output_path()` in the R package so both
writers refuse the same directories.

### Fixed — the web cache refuses the same directories `prepare()` refuses

The web cache had the same defect as `prepare()`, one layer down.
`clear_web_cache(cache_dir=...)` removes the named directory outright and
`ensure_web_ui_cached(cache_dir=..., force=True)` renames it aside and then
removes it, but the only guard refused the filesystem root. Every server and
viewer entry point exposes that argument — `cellucid serve --web-cache-dir`,
`serve(web_cache_dir=)`, `serve_anndata(web_cache_dir=)`,
`show(web_cache_dir=)` — so `clear_web_cache(cache_dir=os.getcwd())` deleted the
working directory and `cache_dir="~"` would have deleted the home directory.

`cache_dir` is now validated in `_require_cache_dir`, the one function every
entry point that accepts one passes through, before any path is inspected,
fetched, staged, renamed, or removed. It refuses the same set as `out_dir`, read
from the same definition so the two cannot drift apart: the filesystem root, the
current working directory (`.`, `./`, or in full), the home directory (`~` or in
full), and the directory that holds every home (`..` or in full).

The comparison now resolves symbolic links instead of only normalizing the path
lexically. A symbolic link *at* the cache path was already refused, but one in
the parent chain was not, so `<alias>/<user>` — a real directory behind an
aliased parent — reached the home directory and removed it. `cache_dir` must
also be a string or `os.PathLike`, reported as such instead of as a `pathlib`
internal.

The default cache location, a dedicated cache directory named relative to the
working directory or through `~`, and a cache reached through a symlinked parent
that is not itself protected all behave exactly as before.

### Fixed — export format: a constant continuous field has an encoding

A gene expressed at one identical level in every exported cell — very often
zero, once an atlas is subset to one lineage — and a continuous `obs` column a
subset flattened are both ordinary scientific data. `prepare()` used to refuse
the whole export when it met one, which cost real features: subsetting Suo to
`LVL0 == "Haematopoeitic_lineage"` leaves four genes detected in none of the
published cells, one of them the nameable `NFIB-AS1`.

compact_v1 now has a named case for them:

- `minValue == maxValue` and every code `0`. The entry keeps its exact shape and
  length — `[index, key, minValue, maxValue]` — so nothing about the manifest
  contract moves.
- The writer takes an explicit branch and never derives a scale, so nothing
  divides by `maxValue - minValue`.
- The reader takes the matching branch and returns `minValue` directly, so the
  constant comes back bit-exact rather than as an approximation of itself.
- Native-double variation finer than float32 resolution collapses to one float32
  value, which is one constant, and is published the same way instead of being
  rejected.

`cellucid_prepare()` in the R package implements the identical case. The only
continuous payload still without an encoding is a categorical field's generated
outlier quantiles when *every* quantile is missing, because no category holds
`centroid_min_points` cells — there is no value to publish at all, and that
failure still names every affected field in one run before any file is written.

### Changed — export format: payload files are named by index

Every scientific payload an export writes is now named by its integer index, or
by a fixed neutral name. No filename anywhere in the tree carries an observation
key, a gene name, or a vector field id:

```text
var/0.values.u8.gz     obs/0.codes.u8.gz     obs/0.outliers.u8.gz
obs/1.values.u8.gz     vectors/0_2d.bin.gz
```

- `_varSchema.pathPattern` and every `_obsSchemas` path pattern now substitute
  `{index}` instead of `{key}`, and each field entry declares its own index as
  element `[0]`:
  - `var` — `[index, name]`, or `[index, name, minValue, maxValue]` quantized;
  - obs continuous — `[index, key]`, or `[index, key, minValue, maxValue]`;
  - obs categorical — `[index, key, categories, dtype, sentinel, centroids]`,
    plus `outlierMinValue`/`outlierMaxValue` when quantized.
  A `var` entry and an obs-continuous entry now share a shape deliberately, so a
  reader must take a field's kind from the array it came from, never from the
  entry's length.
- Within one axis directory the indices are exactly `0 … N-1`, each used once,
  and `obs` shares that one space across `_continuousFields` and
  `_categoricalFields` because both write into `obs/`. `prepare()` asserts this
  against the manifest it has just built, and re-derives every declared payload
  path to compare it with the directory it actually wrote — two fields holding
  one index would write one payload over another and the app would then draw one
  field's values under another field's name, with nothing raised.
- The Jupyter/direct AnnData server serves the same index routes
  (`/var/0.values.f32`, `/obs/0.codes.u8`, `/vectors/0_2d.bin`) and emits
  byte-identical manifest shapes to the batch exporter. Only the exact unpadded
  decimal form resolves; `/var/00.values.f32` is a 404.
- Because an identifier is no longer a path, the filename-portability,
  case-collision, and Windows-reserved-name rules are **removed** from
  observation keys, gene names, and vector field ids. `HLA-DRB1/2`, `CON`,
  `% mito`, `细胞`, and the two distinct fields `Field` and `field` now export
  unchanged. What remains is what an identity is for: non-empty, distinct within
  its axis, and drawable exactly as stored — the same display-text rule every
  string category label already obeyed. `dataset_id` is unaffected: it names a
  real directory and a served URL, so it keeps the full portable-component rule.
- The direct AnnData adapter no longer carries its own copy of that identifier
  validator. It had sanitized unsafe characters into `_` where the batch
  exporter rejected them, so the two producers disagreed about which inputs were
  admissible; both now call the one shared rule.

An index is a position, not a stable name: adding a gene or reordering
`obs_keys` renumbers the payload files. Resolve a payload through the manifest
every time rather than caching or hard-coding a path across exports. Exports are
regenerated from source, so re-export with `force=True`; there is no dual-read
path.

### Changed — the notebook viewers, the session applier, and the data server are packages

Three more modules that had grown past what one file explains are now packages.
Every public name is imported from exactly where it was before:
`cellucid.jupyter` still gives `CellucidViewer`, `AnnDataViewer`, `show`,
`show_anndata`, `HookRegistry`, and `BaseViewer`; `cellucid.anndata_session`
still gives `apply_cellucid_session_to_anndata`; `cellucid.server` still gives
`CellucidServer`, `serve`, and `CORSRequestHandler`.

`src/cellucid/jupyter.py`, 2,175 lines, is now `src/cellucid/jupyter/`. The
exact value checks that guard every message crossing the Python/iframe boundary
are in `_wire.py` and depend on nothing but the standard library; the hook
registry and the latest-event snapshot are in `_hooks.py`; the shared viewer,
the notebook-context detection, the HTTP event-routing bridge, and the
interpreter-exit cleanup are in `_base.py`; and the two concrete viewers are in
`_exported.py` and `_anndata.py`.

`src/cellucid/anndata_session.py`, 1,385 lines, is now
`src/cellucid/anndata_session/`. The closed key sets, chunk profiles, and
chunk-id prefixes that define the session contract are written once in
`_schema.py`, so the chunk validator, the highlight planner, and the field
planner can no longer drift about what a session may contain. The exact value
checks are in `_primitives.py`, the three planners in `_chunks.py`,
`_highlights.py`, and `_fields.py`, and the fingerprint check and the atomic
application in `_apply.py`.

`src/cellucid/server.py`, 1,352 lines, is now `src/cellucid/server/`. The
published-state sidecars, the exported-dataset discovery, and the exact
artifact inventory a dataset is allowed to serve are in `_state.py`,
`_datasets.py`, and `_artifacts.py`; the request handler is in `_handler.py`
and the server lifecycle in `_server.py`.

No behaviour changed and no signature changed: of the 104 top-level definitions
moved, 101 are byte-identical, and the three that are not differ only in the
eight deferred imports whose package depth had to increase by one level.
`import cellucid.server` still loads 191 modules with none of numpy, pandas,
scipy, or tqdm among them, and `import cellucid` still resolves every public
name lazily through `__getattr__`.

### Changed — the export writer is a package, and the servers stopped importing numpy

`src/cellucid/prepare_data.py` had grown to 3,556 lines covering argument
contracts, the manifest format, coordinate normalization, quantization,
centroids, payload writing, the publication transaction, the writer lock, and
the dataset catalog. It is now `src/cellucid/prepare_data/`, one private
submodule per responsibility, with `prepare` and `generate_datasets_manifest`
still imported from `cellucid.prepare_data` exactly as before.

The identity, display-text, and output-directory rules moved out to
`src/cellucid/_contracts.py`, which depends on nothing but the standard
library. `server.py`, `jupyter.py`, `anndata_server.py`, and `web_cache.py`
enforce those same rules and previously reached them through the export writer,
so importing any of them pulled in numpy, pandas, scipy, and tqdm. They no
longer do: `import cellucid.server` drops from 811 modules to 188, and
`web_cache` no longer needs to defer its import to stay light.

No exported byte and no public signature changed. The prepared demo datasets
rebuild byte-identically, and the same rule is now stated once in three places
that had each written it twice: the categorical codes storage a caller
declares, the Windows reparse-point test, and the dead `_select_category_dtype`
duplicate of the adapter's `_categorical_storage`.

### Added

- Exact release-contract validation for package, citation, documentation, and
  downstream recipe metadata.
- Installed wheel and source-distribution gates on Linux, macOS, and Windows
  before PyPI publication.
- Reproducible source-distribution normalization with an exact downstream
  recipe SHA-256 gate.
- Reproducible wheel ZIP timestamps through one canonical
  `SOURCE_DATE_EPOCH`.
- An explicit, strictly validated `prepare(created_at=...)` provenance timestamp
  for byte-identical complete export-directory builds.
- Comprehensive current-contract tests for data identity, embeddings, field
  encoding, weighted connectivity, vector fields, sessions, server lifecycle,
  Jupyter messaging, and hosted viewer assets.
- Session-bundle readers now accept a `cellOrder` record in the manifest's
  `datasetFingerprint`, beside the existing `sourceType`, `datasetId`,
  `cellCount`, and `varCount`. The web app records it so that a dataset
  republished under the same identity with the same cell and gene counts but its
  rows in a different order can no longer silently re-point every saved
  selection at different cells. The field is **accepted, not required**, and its
  digest is never re-derived here: it is taken over exported coordinates, while
  `apply_cellucid_session_to_anndata()` targets an `AnnData` whose row order the
  caller controls, so enforcing that identity belongs to the viewer. Accepting
  both fingerprint shapes is also what keeps existing session files loadable. A
  `cellOrder` that is present is shape-checked strictly: exactly `dimension` and
  `digest`, `dimension` exactly 1, 2, or 3, and `digest` exactly 16 lowercase
  hexadecimal characters.
- An explicit `allowed_hosts=` declaration on `CellucidServer`, `AnnDataServer`,
  `serve`, `serve_anndata`, `CellucidViewer`, `AnnDataViewer`, `show`, and
  `show_anndata`, plus a repeatable `cellucid serve --allowed-host HOST`, naming
  the exact extra `Host` authorities a server answers to. This is what a reverse
  proxy such as `jupyter-server-proxy` requires, because the proxy forwards the
  browser's `Host` verbatim and a proxied request carries no signal that
  separates it from a rebound one. Entries are bare host names — no port, no
  scheme, no wildcard — validated when the server is constructed and matched on
  any port; declaring names also applies the `Host` check to a non-loopback
  bind.

### Changed

- The tested Python range is now 3.11 through 3.14.
- Direct AnnData operation now uses one bounded current runtime:
  AnnData 0.12, Zarr 3, and numcodecs 0.16.
- PyPI publication now uses GitHub OIDC trusted publishing and accepts only a
  pushed tag that exactly matches package version metadata and whose commit is
  contained in `origin/main`.
- Server, Jupyter, cache, session, and prepared-export contracts now reject
  invalid or ambiguous inputs before publishing partial state.
- Every string Cellucid prints verbatim must now read on screen as the value it
  stores. A string category label, `dataset_name`, `dataset_description`,
  `source_name`, `source_url`, and `source_citation` are rejected when they
  carry a control character, one of the zero-width characters `U+200B`,
  `U+2060`, or `U+FEFF`, or leading or trailing whitespace of any kind
  including `U+00A0` NO-BREAK SPACE. An empty category label is rejected, as
  are two labels in one field that a whitespace-collapsing renderer draws
  identically. Nothing is trimmed: trimming would rewrite an annotation the
  caller never asked to change, and would merge `"Liver "` into a separate
  `"Liver"` category and move cells between them. The message names every
  offending label in the field at once and gives the one-line repair. The
  `cellucid` R package enforces the identical rule.
- `dataset_description`, `source_name`, `source_url`, and `source_citation`
  previously accepted padded text in Python while the R writer rejected it;
  the two writers now agree.
- `serve_anndata()` and `show_anndata()` now apply the same category-label rule
  as `prepare()`; the adapter's duplicate label validator was removed in favour
  of the exporter's.
- `generate_datasets_manifest()` now validates a published
  `dataset_identity.json` `description` the same way it validates its `name`.

### Fixed

- `prepare()` no longer rejects an export because of a gene it was never asked
  to export. The two gene-identifier rules that survive the move to
  index-named payloads have two different scopes, and each is now checked at
  its own. Being drawable is a property of a name the viewer shows, so it
  covers the genes `gene_identifiers` selects; a `var` row left out reaches no
  manifest and is not checked, exactly as an `obs` column left out of
  `obs_keys` is not. Distinctness spans the whole `var`, because
  `gene_identifiers` addresses `var` rows by identifier and a repeated one
  names no single row. `cellucid_prepare()` in the R package documents the
  identical pair of scopes.
- Python 3.14 installation on macOS no longer resolves to the old
  `numcodecs<0.16` build path that failed on Apple silicon.
- Dataset preparation and direct AnnData serving now enforce finite,
  float32-representable scientific values and exact cell/gene alignment.
- Prepared generation, replacement, and server shutdown are atomic and retain
  actionable errors at their owning boundary.
- Prepared-directory writers share one persistent exact-target lock with the R
  exporter, reject concurrent publication, and recover after process death
  without leaking descriptors.
- Dataset identifiers, request ranges, and cache paths are validated before
  filesystem or network mutation.
- Prepared exports and direct AnnData responses now use canonical gzip headers,
  eliminating clock and output-path bytes from compressed payloads.

### Security

- Jupyter events, session uploads, web asset inventories, and server artifacts
  use exact authenticated schemas, bounded payloads, and path confinement.
- The data servers validate the `Host` header ahead of routing, file reads, and
  event delivery. A loopback bind stops remote clients but not DNS rebinding:
  once an attacker page's name re-resolves to `127.0.0.1` the browser treats the
  server as same-origin, sends no `Origin`, and runs no CORS check. Only one
  well-formed `Host` naming `localhost` or a loopback IP literal, on the port
  actually bound, is served; every other authority receives `421 Misdirected
  Request` with a body that never echoes the supplied value.

### Documentation

- Reworked the Python, R, and web guides around the current executable
  contracts and added real browser screenshots for primary workflows.
- Corrected embedding guidance: an all-identical embedding is rejected because
  it has no finite normalization range.
- Added a complete standard scVelo Pancreas sample guide covering exact
  catalog selection, 1D/2D/3D navigation, velocity, scientific ownership,
  network payload, provenance, and reproducibility.

## [0.0.9] - 2026-01-01

This release graduates Cellucid Python out of alpha (still pre-1.0).
The GitHub release was tagged `v0.9.0`, while its Python distribution metadata
declared `0.0.9`. Release 0.9.1 restores one exact version across both systems.

### Added

- AnnData-first viewing and serving with `.h5ad` and `.zarr` support (lazy/backed loading where possible).
- Unified CLI serving with auto-detection for AnnData files, Zarr stores, and
  pre-exported dataset directories.
- Jupyter notebook integration (`CellucidViewer`, `AnnDataViewer`, `show()`, `show_anndata()`) with event hooks and session export.
- Session bundle support (`.cellucid-session`) including `CellucidSessionBundle` and `apply_cellucid_session_to_anndata()` for round-tripping highlights and user-defined fields back into AnnData.
- Multi-dimensional embedding exports (1D/2D/3D) and vector-field overlays (RNA velocity / drift) via `prepare()` and `vector_fields` helpers.
- Hosted web UI proxy mode with on-disk caching helpers (`get_web_cache_dir()`, `clear_web_cache()`).

### Changed

- Export format now includes explicit dataset identity metadata (`dataset_identity.json`) for reproducible sharing and session compatibility checks.
- Reduced export size and improved load performance with optimized manifests, improved connectivity edge export, and optional quantization + gzip compression knobs.

### Security

- Session bundles are treated as untrusted input with bounds checks and dataset mismatch policies when applying to AnnData.

### Documentation

- Major Read the Docs expansion and restructuring (Python package + web app guides), plus new publishing and contributing documentation.

## [0.0.1a0] - 2025

### Added

- Initial alpha release
- AnnData visualization with UMAP embeddings (1D, 2D, 3D)
- Gene expression overlays with sparse matrix support
- Cell metadata coloring (categorical and continuous)
- Interactive filtering and cell selection
- KNN connectivity visualization
- Multiple deployment modes:
  - Local demo (browser-only)
  - Browser file picker (h5ad and exported formats)
  - Server CLI (`cellucid serve`)
  - Python API (`serve()`, `serve_anndata()`)
  - Jupyter integration (`show()`, `show_anndata()`)
- Export functionality for web deployment

[0.9.1]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.1
[0.0.9]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.0
[0.0.1a0]: https://github.com/theislab/cellucid-python/releases/tag/v0.0.1a0
