# Decision record: decoders and codecs

This page answers: **"Which decoder decompresses a chunk, why that one, what
was measured, and what was tried and taken out?"** It covers the WebAssembly
decoders (zstd, lz4, c-blosc), the vendored fflate inflater, the platform's
`DecompressionStream`, the JavaScript decoders written for the reader, the
measured choice between them, and the decode pool that runs them. Module
paths are under `cellucid/assets/js/` unless they start with `cellucid/`.
Update this page in the same change as the code it describes; a measurement
quoted here names its engine, machine and date, so a later run that
disagrees replaces the row, it does not sit beside it.

## At a glance

**Audience**
- Developers and coding agents changing a codec, the decode pool, the
  WebAssembly build, or the page's script policy.

**Time**
- 15 minutes.

**Prerequisites**
- {doc}`20_read_path_architecture` (where decoding sits in the read path).

---

## What it does now

1. **The format names the family.** A Zarr codec list (`formats/zarr/codecs/pipeline.js`,
   `decoderFamily`) or an HDF5 filter pipeline (`formats/hdf5/filters.js`,
   `filterFamily`: DEFLATE 1, Blosc 32001, LZ4 32004, bitshuffle 32008,
   Zstandard 32015) names the codec family of each chunk; the read plan
   carries it on the decode task.
2. **The pool takes the chunk** (`formats/io/decode-pool.js`). A chunk that
   decodes to less than `POOL_MIN_BYTES` (32 KiB, measured: a round trip to
   a worker costs about that much decoding) travels with other small chunks
   in one message when together they are worth one, or is decoded in the
   reader thread. Larger chunks go to decode workers. A `background` task (a
   gene index built while no one waits) queues behind every other task.
   Workers are made as work queues (one ahead when a store opens,
   `warmDecodePool`), up to `hardwareConcurrency − 2` and to what the
   memory model's `decodeWorkers` share holds at `DECODE_WORKER_BYTES`
   (40 MiB, measured) each, and end after 10 s idle
   (`workers/idle-retirement.js`).
3. **Each message names its decoder** (`formats/io/decoder-choice.js`).
   Families with more than one decoder are listed in `DECODERS`, the
   reference first: `deflate: native | fflate`, `zstd`, `lz4`,
   `blosc: wasm | js`. The choice is made per family and size class of the
   decoded chunk (`choiceKey`, a power of two). A person's read always goes
   with the key's decided decoder, or, before a decision, with the engine's
   prior (`ENGINE_PRIORS`: Gecko `fflate`/`js`/`js`/`js` for
   DEFLATE/zstd/LZ4/Blosc; Blink `js` for zstd and Blosc, the references
   otherwise; WebKit the references). Never with a decoder untested on this
   machine.
4. **The probe measures only while idle** (`DecoderProbe`). Once no read has
   sent a chunk for `EXPLORE_AFTER_IDLE_MS` (1 s), the pool sends copies of
   chunks it kept (at most one per worker per undecided key) in phases of
   `PHASE_MESSAGES_PER_WORKER` (8) messages per worker, all with one
   decoder, so each decoder is measured under its own contention. The first
   phase of each decoder warms the engine's compiler tiers; then decoders
   alternate, prior first, until each has `PROBE_PHASES` (2) phases after
   its warm-up. A message's time is its worker's busy time, measured in the
   pool's thread. The decoder with the most bytes per busy millisecond over
   all its phases, **the warm-up included**, wins: a person's reads are
   short, so a decoder fast only once optimised is charged its cold cost. A read's chunk arriving mid-phase voids the phase.
   Chunks small enough to travel together are not kept as samples, so their
   keys stay with the prior.
5. **Every non-reference decoder is checked against the reference, but
   never on a read's path** (`formats/io/decode-worker.js`). The first
   measuring message of each such decoder for a key decodes its sample with
   both; a decoder whose bytes differ is dropped for that key. A decoder
   that throws on any chunk, a read's included, where the reference
   decodes is dropped too, and the reference's bytes are used. A read's
   chunk is decoded once: checking it with the reference cost the first
   gene a second decode with a cold WebAssembly module (r5, below). Before a
   key has been measured, a read uses the prior, which is the decoder the
   reader shipped with for that engine.
6. **Decisions are kept** (`formats/io/decoder-memory.js`): IndexedDB
   database `cellucid-decoders`, one record `{version, decided, at}`. The
   version is the WebAssembly module's `?v=`, the `DECODERS` lists, the user
   agent and the core count, so a new module, a browser update or another
   machine measures again. The record is read once when the pool is made;
   until it is read the pool uses the priors and measures nothing. The
   Stored dialog lists it (*Decoder choices measured*,
   `app/ui/modules/kept-stores-dialog.js`) and *Remove all* forgets it.
7. **The WebAssembly module** (`formats/io/wasm-codecs.js`): one file,
   `cellucid/assets/wasm/codecs.wasm` (103,129 bytes, no imports), compiled
   once per thread on first need from its address in
   `cellucid/assets/workers/scripts.js` (`codecsWasm()`, `?v=` the first 12
   hex digits of its SHA-256, so the service worker keeps it like any
   script). One instance per thread, its memory grown to the largest chunk
   so far. An instance that grew past `DECODE_WORKER_BYTES` is let go after
   its call, since WebAssembly memory never shrinks. Each call is given the
   size the format declares and writes no byte past it (`WasmCodecError`
   otherwise).
8. **Blosc frames of zlib streams** (what `prepare()` and
   `cellucid_prepare()` write for shuffled Zarr v3 arrays): c-blosc parses
   the frame, the DEFLATE family's decoder inflates each stream in place
   (`bloscZlib`), and the bytes are unshuffled. No zlib code is in the
   module.
9. **DEFLATE** (`formats/io/inflate.js`): the platform's
   `DecompressionStream`, fed in slices so the declared-size bound holds
   while it decodes, or fflate 0.8.2 (`cellucid/assets/external/fflate-0.8.2/`)
   for a whole chunk of known size (output allocated at the size plus one
   byte, Adler-32 or CRC-32 and size trailer checked). A stream of unknown
   size always goes through the platform's stream.
10. **Script policy**: `cellucid/index.html` grants `'wasm-unsafe-eval'` in
    `script-src` and nothing else that evaluates; no other page grants it.
    Where WebAssembly is refused (Lockdown Mode, a policy), data that needs
    the WebAssembly decoders fails with a message naming the cause and what
    to do (`formats/io/wasm-codecs.js`).

## Decisions and why

| Date | Decision | Why | Rules out |
|---|---|---|---|
| 2026-10-28 | DEFLATE stays on the platform's sliced stream as the reference. | Fastest on one core in every engine (h5py chunks 439 / 275 / 413 MB/s Chromium / Firefox / WebKit) against a whole-chunk stream (135 / 171 / 84) and fflate (213 / 160 / 80). | An own inflater in JavaScript. |
| 2026-10-29 | Zstandard, LZ4 and Blosc decoded by their own libraries (zstd 1.5.6, lz4 1.10.0, c-blosc 1.21.6) compiled to WebAssembly by a repository script. | The libraries are the format references (zstd checksums verified, LZ4's end-of-block rule enforced) and win 1.3 to 4.2× on one core in Chromium and WebKit (zstd 1.7×, LZ4 3.4 to 3.7×, Blosc-LZ4 1.3 to 1.5×, bit shuffle 4.2×). | emscripten builds, a runtime, WASI, committed upstream source. |
| 2026-10-29 | Build pinned and reproducible: zig 0.13.0 (SHA-256 per host from ziglang.org), release tarballs by the SHA-256 their release pages publish (c-blosc publishes none: its tag archive's digest as Homebrew pins it), `wasm32-freestanding`, SIMD, no C library. | `--check` rebuilds byte for byte; the test pins the module and its licences by digest. | Unpinned toolchains, a module nobody can rebuild. |
| 2026-10-29 | fflate 0.8.2 vendored unchanged as the second DEFLATE decoder. | Firefox's `DecompressionStream` contends across threads on small chunks: one worker 292 to 310 MB/s, eight 159 to 160 in total; fflate at eight Firefox workers 882. | Choosing the inflater by engine name. |
| 2026-10-29 | The JavaScript decoders stay as each family's second decoder (`formats/zarr/codecs/javascript/`). | Firefox runs them faster than its WebAssembly, for zstd, LZ4 and Blosc alike (below). | Deleting a decoder that still wins on one engine. |
| 2026-10-29 | `'wasm-unsafe-eval'` added to the app page's `script-src` only. | Compiling WebAssembly needs it; workers inherit the page's policy. JavaScript evaluation stays refused. | `'unsafe-eval'`; granting it on other pages. |
| 2026-10-29 | The decoder is chosen per family and size class by a probe at the pool's width, with every worker busy. | A decoder fast on one core can be slowest at eight (Firefox DEFLATE). Comparing at the same width and dispatch cancels the pool's own costs. | One-core measurements; per-message wall time (see below). |
| 2026-10-30 | Before a decision, a read uses the engine's measured prior; measuring happens only after 1 s with no read; decisions are kept in IndexedDB. | Measuring on a person's read made first genes slower (Firefox zstd first gene 119 → 206 ms in one run). Keeping decisions means only a first visit measures. | Exploring with untested decoders on a read. |
| 2026-10-30 | Chromium starts Zstandard and Blosc from JavaScript; the probe counts its warm-up phase. | Same-tree A/B (r4, 10M CSR Zarr, first gene, JavaScript tree → WebAssembly): Chromium zstd 175–244 → 386–678 ms, bit-shuffled Blosc-LZ4 111–360 → 357–488; each worker compiles the module and a gene's few chunks never leave V8's baseline tier. LZ4 stays WebAssembly (319 → 196). A warmed benchmark (WebAssembly zstd 4,069 against 2,144 MB/s) would keep choosing it, hence the counted warm-up. | Naming the engine in the decision (the prior is only the start; the measurement decides). |
| 2026-10-30 | A read's chunk is never checked against the reference; the measuring messages check each decoder. | r5 (same tree): with Zstandard and Blosc on the same JavaScript decoders as the JavaScript tree, the first gene was still 5 to 10 ms slower: the check decoded the first chunk again with the cold WebAssembly reference. | A second decode on a read's path. |

## Measurements that justify it

Machine: Apple M1 Pro (10 cores). Engines: Chromium 151, Firefox 153,
WebKit 26.5 (Playwright builds). Chunks: the codec corpus's real chunks
(`cellucid/tests/fixtures/codecs/`), byte-identical across decoders.

**The pool, each decoder fixed, then the probe left free (2026-10-30, MB/s
decoded, `scripts/efficiency-floor/cpu-browser.mjs`).** The probe chose the
fastest fixed decoder in 21 of 21 engine × codec cases, with every key
decided while idle.

| Codec | Chromium | Firefox | WebKit |
|---|---|---|---|
| DEFLATE, h5py chunks | **native 1,710**, fflate 1,078 | **fflate 620**, native 115 | **native 992**, fflate 413 |
| zlib, 538 KB chunks | native | **native 992**, fflate 862 to 872 | native |
| Zstandard | **wasm 4,432**, js 2,192 | **js 1,653**, wasm 532 | **wasm 3,306**, js 2,267 |
| LZ4 | **wasm 7,703**, js 3,265 | **js 3,606**, wasm 2,938 | **wasm 11,334**, js 3,606 |
| Blosc-LZ4 | **wasm 6,959**, js 3,707 | **js 1,935**, wasm 902 | wasm |
| Blosc bit shuffle | **wasm 6,198**, js 1,609 | **js 1,150**, wasm 763 | **wasm 6,103**, js 1,526 |
| Blosc-zlib (streams by the DEFLATE decoder) | **native 2,755**, fflate 1,163 | **fflate 786**, native 512 | native |

Cells without numbers: the probe chose the decoder named; that summary kept
no rates for them. Firefox's prior for 538 KB zlib chunks is fflate, and the
probe overruled it with native: the prior is a first guess only.

**One core, the app's WebAssembly against the JavaScript decoder (MB/s,
2026-10-30):** Chromium zstd 662 / 346, LZ4 2,593 / 685, bit shuffle
1,179 / 236; WebKit zstd 620 / 411, LZ4 2,645 / 728; **Firefox zstd 79 / 257,
LZ4 413 / 563**: Firefox's WebAssembly is 3 to 5× slower than its JavaScript
for these decoders, at one core and at the pool's width (zstd 355 against
1,763 MB/s at eight workers, 2026-10-29). Node (V8), one core, 2026-10-29, the JavaScript decoders → the WebAssembly ones:
zstd 345 → 717, LZ4 664 → 2,695, Blosc-LZ4 841 → 1,224, bit shuffle
236 → 1,108, blosclz 710 → 968, Blosc-zlib 614 to 648 → 638 to 685.

**End to end.** Full HLCA (2,282,447 cells) from its CELLxGENE link, first
gene, JavaScript-only tree of 2026-10-29 → WebAssembly tree, one run each:
Chromium 12.8 → 3.9 s, Firefox 12.4 → 5.4 s, WebKit 11.8 → 4.4 s. The gain is
not the decoders' alone: the earlier tree also lacked other read-path work
of the same days. Synthetic 10M-cell CSR Zarr stores per codec over loopback
(three repetitions per engine, 2026-10-30): run-to-run spread is wider than
the differences (Chromium LZ4 first gene 113 to 400 ms on the
JavaScript-only tree, 383 to 663 ms on the WebAssembly tree before the
priors and idle-only measuring landed).

**Same-tree A/B, first gene (r5, 20:00 tree, 10M CSR Zarr over loopback,
three runs, median ms, JavaScript-only copy → this tree; next gene in
brackets).** Chromium zstd 49 → 59 (127 → 109), bit shuffle 49 → 58
(123 → 115), Blosc-LZ4 44 → 42 (111 → 105), LZ4 48 → 50 (110 → 127);
Firefox zstd 32 → 37 (105 → 97), bit shuffle 33 → 39 (95 → 103), Blosc-LZ4
27 → 33 (96 → 101), LZ4 46 → 35 (106 → 101). Zstandard and Blosc ran the same
JavaScript decoders in both trees; the remaining 5 to 10 ms was the read
path's check (decision above, removed; r6 confirms).

## Tried and reverted

| What | When | Result | Would need |
|---|---|---|---|
| numcodecs 0.3.2 (emscripten builds of zstd and c-blosc) | 2026-10-27 | Not adopted: Firefox 0.27× (zstd) and 0.25× (Blosc-LZ4) the app's decoders. `@bokuweb/zstd-wasm` does not load in a browser without a bundler; fzstd 0.1.1 decoded 35 to 72 MB/s (2026-10-28). | Replaced by the lean decode-only build. |
| A decode width per codec, chosen by throughput (`DecodeWidth`) | 2026-10-28 | Reverted: it settled on wrong widths when dispatch, not the decoder, bounded throughput. Chromium `.h5ad` 1M 1,838 → 1,076 MB/s, Zarr 10M 2,093 → 1,149; Firefox unchanged. | Replaced by the decoder probe at a fixed width. |
| WebAssembly as the only decoder for zstd, LZ4 and Blosc | 2026-10-29 | Firefox 3 to 5× slower; the JavaScript decoders came back the same day as candidates. | A family whose WebAssembly wins in every engine at every width may drop its `js` decoder. |
| Probe timing per message (wall time inside the worker) | 2026-10-29 | Chose wrong: Firefox kept native at 292 MB/s where fflate gave 882; WebKit took fflate at 363 where native gave 1,063. Two interleaved messages counted an async decoder's time twice, and too few bytes reached the compiler tiers. | Replaced by phases of one decoder, a warm-up phase, and busy time measured in the pool's thread. |
| Measuring on a person's read | 2026-10-29 | First genes slower in single runs (Firefox zstd 119 → 206 ms; Chromium LZ4 and Blosc 180 to 310 ms slower). Measuring messages also slipped between a read's waves. | Replaced by priors, idle-only measuring after 1 s with no read, and kept decisions. |
| Compiling the module once in the reader thread and posting the compiled `WebAssembly.Module` to each decode worker | 2026-10-30 | Reverted: in Chromium the run ended with no result in 2 of 4 runs with the change, 0 of 4 with it reverted, 0 of 4 on the JavaScript-only tree; an earlier set failed 8 of 12 (0 of 24 in Firefox and WebKit). Cause not found. | The cause understood, then 12 of 12 clean Chromium runs before it lands again. Each worker fetches (from the service worker's cache) and compiles the module itself meanwhile. |
| Cell exports: workers decode into reused buffers and return only the rows an export keeps | 2026-10-30 | Reverted: chunks two runs touch were fetched and decoded twice, nothing was cached, and the peak barely moved (HLCA core `.h5ad` export, Chromium: 754 s before, timed out at 780 s after). | One read per export part listing all its row ranges, each chunk decoded once and cut on the worker: a reader change, see {doc}`25_decisions_exports`. |

## Known limits and open items

- **Firefox and WebAssembly.** Firefox runs these WebAssembly decoders 3 to
  5× slower than its JavaScript decoders (engine behaviour, measured above);
  the probe keeps the JavaScript decoders there.
- **Firefox `DecompressionStream` on small chunks** contends across threads
  (an engine property); fflate is chosen there by measurement.
- **The JavaScript Zstandard decoder does not verify content checksums**
  (`formats/zarr/codecs/javascript/zstd.js`); the WebAssembly one does.
  Where the probe chose JavaScript, the decoded size and the chunk's own
  checksums (CRC-32C in Zarr v3) still guard integrity.
- **Small chunks keep the prior.** Chunks under `POOL_MIN_BYTES` are never
  sampled, so their size classes are decided by the engine's prior, not
  measured.
- **First gene against the JavaScript-only tree** is not yet shown to be no
  slower in every engine × codec; repeated runs are open.
- **Compiled-module sharing** (above) is open until its Chromium failure is
  understood.
- **`.mtx.zst` text matrices** are not recognised by the input-format
  registry (three repositories).
- **bz2 and LZMA** decode at 25 to 66 MB/s (rare in stores; left).
- **R's zstd files**: R 4.5.3 builds here (and conda-forge's) lack zstd, so
  the test files are uncompressed R files compressed with libzstd.

## Where the code and tests live

| Concern | Module(s) | Tests |
|---|---|---|
| WebAssembly build, pins, licences | `cellucid/scripts/build-wasm-codecs.mjs`, `cellucid/scripts/wasm-codecs/`, `cellucid/assets/wasm/codecs.wasm`, `cellucid/assets/wasm/LICENSES/` | `cellucid/tests/wasm-codecs-contract.test.mjs` |
| Loading and calling the module | `formats/io/wasm-codecs.js`, `cellucid/assets/workers/scripts.js` | `wasm-codecs-contract`, `cellucid/tests/worker-bundles-contract.test.mjs` (the module among the workers' assets) |
| Decoder lists, priors, probe | `formats/io/decoder-choice.js` | `cellucid/tests/formats-decode-pool.test.mjs` |
| Kept decisions | `formats/io/decoder-memory.js`, `app/ui/modules/kept-stores-dialog.js` | `formats-decode-pool`, `cellucid/tests/kept-stores-dialog.test.mjs` |
| Pool, batching, priority, idle exploration | `formats/io/decode-pool.js`, `formats/io/decode-worker.js` | `formats-decode-pool`, `cellucid/tests/codec-layer-contract.test.mjs` |
| DEFLATE (native, fflate) | `formats/io/inflate.js`, `cellucid/assets/external/fflate-0.8.2/` | `wasm-codecs-contract`, `codec-layer-contract` |
| Zarr codecs | `formats/zarr/codecs/{zstd,lz4,blosc,pipeline}.js`, `formats/zarr/codecs/javascript/` | `cellucid/tests/zarr-format-codecs.test.mjs`, `wasm-codecs-contract` (JavaScript decoders give the reference's bytes on every corpus chunk) |
| HDF5 filters | `formats/hdf5/filters.js` | `codec-layer-contract` (fixtures `cellucid/tests/fixtures/codecs/hdf5plugin-bitshuffle.h5`, `cellucid/tests/fixtures/codecs/bitshuffle.h5ad`) |
| R zstd streams | `formats/rds/stream.js` | `cellucid/tests/rds-corpus-conformance.test.mjs` |
| Script policy | `cellucid/index.html` | `cellucid/tests/seo-head-contract.test.mjs`, `cellucid/tests/browser/content-security-policy.spec.mjs` |
| Served as `application/wasm` by Python and R | inventory from `cellucid/scripts/generate-web-asset-inventory.mjs` | `cellucid-python/tests/test_web_asset_delivery_contract.py`, `cellucid-r/tests/testthat/test-web-cache.R` |
| Measuring harness | `cellucid/scripts/efficiency-floor/cpu-browser.mjs`, `cellucid/scripts/efficiency-floor/app-run.mjs` | — |

## Related pages

- {doc}`20_read_path_architecture`: where decoding sits between bytes and
  the AnnData reader.
- {doc}`09_data_loading_pipeline_and_caching`: sources, caches and failures.
- {doc}`01_codebase_map_and_entry_points`: the repository map, including
  `assets/wasm/` and `assets/external/`.
- {doc}`03_build_run_and_deployment`: the page's script policy.
- {doc}`04_configuration_env_vars_and_feature_flags`: the IndexedDB
  database `cellucid-decoders` among the browser's stores.
- {doc}`27_decisions_memory_model`: the `decodeWorkers` and
  `decodedBlocks` shares.
- {doc}`23_decisions_partial_reading`: which chunks a share reads at all.
- {doc}`25_decisions_exports`: the export read that would let workers keep
  only the rows an export needs.
