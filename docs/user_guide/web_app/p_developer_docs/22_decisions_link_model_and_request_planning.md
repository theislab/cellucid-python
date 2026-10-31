# Decision record: link model and request planning

This page answers: **"How does the web app measure the link to a file, turn a
reader's byte asks into HTTP requests, and stay gentle with public hosts, and
why is it built that way?"** It records the current design, the decisions
behind it with their dates, the measurements that justify them, and what was
tried and reverted, so nobody re-tries those blind. Update this page in the
same change as the code in `formats/io/` (`link.js`, `byte-reader.js`,
`byte-source.js`, `hosts.js`, `range-cache.js`, `kept-parts.js`) or
`data/host-link-memory.js`. Paths below are under `cellucid/assets/js/`
unless they start with `cellucid/`.

## At a glance

**Audience**
- Developers changing how bytes are fetched, priced,
  retried or kept.

**Time**
- 20 minutes.

**Prerequisites**
- {doc}`20_read_path_architecture` (which module owns each step) and the
  "Bytes" part of {doc}`09_data_loading_pipeline_and_caching`.

---

## What it does now

1. **A byte source per file** (`formats/io/byte-source.js`). A local `File`
   (`blobByteSource`) or a web file (`httpByteSource`: one `Range` header,
   `cache: 'no-store'`, no other header, so no CORS preflight). Each source
   carries a prior, its `ReadPolicy`: `REMOTE_POLICY` (16 KiB metadata
   blocks, 64 KiB merge gap, 8 MiB requests, 6 connections, 100 streams,
   64 MiB cache, no bandwidth assumed) or `LOCAL_POLICY` (64 KiB blocks,
   4 KiB gap, 32 MiB requests, concurrency 4, assumed 1 GB/s disk, read
   ahead from the start). A public S3 object is read from all of its
   equivalent addresses (`formats/io/s3.js`, seven host names on the
   CELLxGENE Census bucket), each with its own six connections
   (`formats/io/mirrors.js`: a request goes to the least-asked address; a
   mirror is dropped only on evidence that it serves another file, or after
   two unreachable answers in a row; a busy answer is no evidence).
2. **One `LinkMeter` per byte reader** (`formats/io/link.js`) measures every
   request (`begin` / `end`) over a sliding window of the latest 16
   (`RECENT_REQUESTS`), as TCP BBR models a path:
   - `latencyMs`: the shortest of the latest requests whose time was mostly
     a round trip (their transfer at most half their time). Used for
     shaping.
   - `loadedLatencyMs`: the median of the latest requests' time less their
     transfer at their share of the bandwidth, counting only requests that
     found a connection already set up. `loadedRoundTripMs(link)` is what
     reads are *priced* at (gene sample, Read for all cells, exports, the
     progressive first view, the analysis run tuner).
   - `loadedRequestMs`: with every slot busy and requests waiting, the time
     between completions less their bytes, over a full window: the host's
     request rate under load.
   - `bytesPerMs`: the fastest of the latest deliveries, per request or over
     a stretch of requests continuously in flight. A delivery bounded by its
     round trip (under `TRANSFER_BOUND`, three round trips) can raise it,
     never lower it (BBR's application-limited samples).
   - `connectionBytesPerMs`: what one connection carries, from requests
     alone on the link and mostly transfer. Under 80 % of the link
     (`CONNECTION_BOUND_SHARE`) the link is per-connection (S3) and
     `connectionRoundTripBytes` is set; the link here is the fourth fastest
     transfer-bound request end ever measured (`LINK_EVIDENCE`), never a
     single delivery.
   - `clockTickMs()`: each request is counted half a clock tick longer
     (`performance.now()` ticks 1 ms in WebKit and Firefox workers, 0.1 ms in
     Chromium), so a same-machine server never measures a round trip of 0.
3. **Protocol and concurrency.** Until the protocol is known, requests are
   dispatched at the host's streams (100) but reads are shaped for six
   (`shapingConcurrency`), and plans are priced under both
   (`concurrencies`). The protocol comes from Resource Timing
   (`watchProtocol`, `nextHopProtocol`) or, where the host hides it (no
   `Timing-Allow-Origin`, as GitHub), from bursts: the seventh completion
   within half a round trip of the first proves multiplexing; two bursts
   whose first six end together and the seventh half a round trip later
   prove HTTP/1.1. A request begun with every known connection busy waited
   in the browser's queue and is not a round-trip sample. Once multiplexed,
   the concurrency doubles while completions per millisecond rise by half
   of what doubling could give, then stays at the last level that paid.
   Background work keeps to `backgroundConcurrency` (one fewer than the
   connections until the host multiplexes, then one fewer than the level
   that paid), so a read a view waits for always finds a connection.
4. **What follows from the link**: `requestCostMs` = the larger of the host
   budget's price (`hosts.js` `requestCostMs`) and `loadedRequestMs`;
   `requestMs` = the larger of `latencyMs / shapingConcurrency` and
   `requestCostMs`; `mergeGapBytes` = `pricingBytesPerMs` x `requestMs`
   (with a local file's prior as its floor); `maxRequestBytes` = bandwidth x
   round trip, never below the prior; `roundTripBytesOf(link)` = the one
   round-trip-in-bytes every cost model uses (`data/sparse-column-index.js`,
   `formats/hdf5/chunk-index.js`, `data/anndata/matrix-access.js`);
   `pricingBytesPerMs` = the bandwidth, doubled while requests wait (so a
   link measured only on small requests still fetches the gaps that would
   pay on a faster one, and measures itself again); `readsAhead` at
   64 MB/s or more (`READ_AHEAD_BYTES_PER_MS`), or always for a local file.
5. **The byte reader** (`formats/io/byte-reader.js`, `createByteReader`):
   reads asked in one task are planned on the next macrotask; parts kept in
   this browser and parts already carried by a request in flight are
   subtracted first (`notInFlight`; a read joins the carrying request);
   `planMergeGap` fetches gaps smallest first, pricing each candidate plan
   at `count x max(latency / concurrency, requestCostMs) + bytes / bandwidth`
   under *every* concurrency the host may have, and takes the plan whose
   worst ratio to the best under each is smallest; no gap is merged across kept bytes or bytes in flight;
   long requests are split by `requestPieceBytes` (down to one connection's
   round trip of bytes where a connection carries only part of the link);
   queued neighbours ride along with a request about to start while round
   trips outweigh bytes (`companionsOf`, coalescing to the link); beyond six
   requests on a host not known to multiplex, more go out only while those
   in flight carry under two round trips of bytes (`slotFree`). Gap bytes no
   read asked for are kept as surplus within `cacheBytes`, pinned while the
   read that fetched them runs, and let go after `SURPLUS_IDLE_MS` (1 s)
   with nothing in flight. Reads marked `cache: true` go through 16 KiB
   aligned metadata blocks in an LRU cache; where a host budget binds,
   neighbouring blocks come along (`metadataReach`). `whenRoom` lets a
   reader gather its own small reads while requests queue.
6. **The host book** (`formats/io/hosts.js`, `realmHosts`, one per realm,
   told to the other realms of the site over a `BroadcastChannel`):
   answer times (seed a new link's stall timeout), deliveries (a request is
   stalled only when no request to its host anywhere in the site received
   anything for a whole timeout), request budgets from `X-RateLimit-*`, the
   IETF `RateLimit-*` fields and their structured form, and `Retry-After`,
   halving the requests at once after a limit answer (`slowDown`: additive
   increase, multiplicative decrease), a learnt pace after a block, and an
   address left unasked after two lost rounds (RFC 6298 back-off to 60 s,
   one probe at a time).
7. **Every request** goes through `fetchOnce` (admission by the host book,
   stall watch from the link's `ResponseClock`: RFC 6298 RTO, 3 s before the
   first answer, never under 1 s, plus the expected bytes' transfer) and
   `withRetries` (waits grown from the RTO, at least three repeats over at
   least 100 s per RFC 1122 R1/R2 on a read its host already delivered to,
   three repeats (R1 alone) for a busy answer before anything was
   delivered; a host-named wait waited exactly and not counted; a host
   that answered before and now refuses this network waited out for as
   long as the read is wanted).
8. **Memory across readers and visits.** `formats/io/range-cache.js` (worker)
   and `data/byte-range-cache.js` (page) share fetched ranges between reader
   workers of the same file version; a page-held piece smaller than the gap
   the link would merge is fetched again with its neighbours instead of
   splitting the request. `formats/io/kept-parts.js` writes every fetched
   byte of a versioned remote file to the origin-private file system as
   journals (manifest first, CRC-framed records, 16 KiB page checks from
   `formats/io/page-checks.js`); a revisit reads them before any request.
   `data/host-link-memory.js` remembers per origin the round trip, round
   trip under load, sustained rate and the bytes they rest on, plus a
   blocked host's pace (`cellucid:host-link:v1`, `cellucid:host-pace:v1`,
   see {doc}`04_configuration_env_vars_and_feature_flags`).
   `pricingLink(snapshot, record)` lifts only `bytesPerMs` and `latencyMs`
   (together, when the record's rate scaled to this path's warm round trip
   is higher) and `loadedLatencyMs` (scaled the same way, when higher), sets
   `measured` and `remembered`, and returns the fresh snapshot unchanged
   once the visit has filled the link (`bandwidthMeasured`) or read as many
   bytes as the record rests on. It does not recompute the snapshot's
   `roundTripBytes`, `mergeGapBytes`, `maxRequestBytes`, `readsAhead` or
   `requestCostMs`: the byte reader's own planning never sees the record.
   A seed whose origin is on this machine (`thisMachine`, a loopback
   address), or a local file's snapshot (`onThisMachine`), is priced as
   this machine at rest (`thisMachinePricing`, marked `thisMachine`):
   `bytesPerMs` lifted to `LOCAL_POLICY.bytesPerMs` whenever it is below
   it, `bandwidthMeasured` or not, `loadedLatencyMs` lowered to
   `latencyMs` and `requestCostMs` to 0; a local file's
   snapshot already carries its priced rate (`_pricedBytesPerMs`, never
   below its prior's).

## Decisions and why

| Date | Decision | Why | Rules out |
|---|---|---|---|
| — | Nothing asks where a file is; every shape follows the measured link, starting from the source's prior | One code path for a disk, `cellucid serve`, a kernel and the internet | Per-host or per-provider tuning tables |
| — | A single `Range` header and `cache: 'no-store'` on every range request | CORS-safelisted (no preflight some CDNs refuse); Chromium lets one request per URL through its cache at a time | Conditional headers on reads; the HTTP cache for ranges |
| — | 16 KiB remote metadata blocks, 512 KiB first head (`HEAD_BYTES`) | anndata keeps each element's header, B-tree node and heap within a few KiB; the head of a Cellucid-ready `.h5ad` holds all of its opening | Larger blocks (fetch mostly data); a smaller head |
| 2026-10-28 | Host request budgets read from headers and price every plan (`requestCostMs`) | A request on Zenodo costs 60 s / 133 = 0.45 s once the budget binds; pricing it at latency / concurrency sent hundreds of tiny reads | Treating requests as free on rate-limited hosts |
| 2026-10-28 | Stall timeout from RFC 6298, 3 s until the first answer; stalls judged per host across realms | Cold first bytes took 0.7-2.4 s and were aborted at 1 s; a request queued behind another reader's is waiting, not stalled | A fixed timeout; per-request silence as a stall |
| 2026-10-28 | Retries from the link's RTO, RFC 1122 R1/R2 (3 repeats, 100 s) | The link's own measure, not a schedule | A fixed retry-delay list (replaced) |
| 2026-10-31 | A busy answer before the host delivered anything to the read: R1 alone (3 repeats), then reported | A catalog answering 503 was asked for 100 s and outlived startup's 60 s step bound, which then said it "did not answer"; the person sees the server's error with Try again | R2's 100 s for requests a person waits on with nothing on screen |
| 2026-10-28 | An address that loses two rounds is left unasked (1 s doubling to 60 s), one probe at a time | A refusing host was asked twice by every read (40 reads, 74 requests) | Every waiting read asking a host that refuses this network |
| 2026-10-28 | Page-held pieces under the merge gap are fetched again with their neighbours (`withRangeCache`) | Holding pieces split one 8 MiB plan into ~190 requests | Splitting a request around every held piece |
| 2026-10-29 | Price reads at the round trip under load (`loadedRoundTripMs`); shape with the shortest | S3 answers in 230 ms idle and 0.35-0.5 s with 42 in flight | Pricing at the idle round trip |
| 2026-10-29 | Kept parts as journals: manifest first, each piece appended with a CRC-framed record | Parts kept in 15 s batches with the manifest last kept nothing from a short visit | Batches written at the end |
| Decision (2026-10-30) | Do not assume a stable connection | The owner's links vary and drop | A link measured once and trusted; whole downloads that restart from byte 0 |
| 2026-10-30 | Windowed link model (16 requests: min round trip, max delivery; application-limited samples never lower bandwidth); `pricingLink` gives way once the visit fills the link | Follows a link that gets faster at once, a slower one within a window, ignores one odd request | Lifetime averages; a remembered rate overriding a link now measured |
| 2026-10-30 | Host request rate under load (`loadedRequestMs`) in `requestCostMs` | CELLxGENE's mirrors complete 70-88 small requests a second with 42 in flight (14 ms a request) where round trip / connections said 6 ms | Pricing a request below what the host completes |
| 2026-10-30 | Split long requests over every connection where one connection carries only part of the link (`requestPieceBytes`), judged only from requests alone on the link | S3 delivers per connection; samples taken while sharing the link split a plan into ever smaller requests | Few multi-MB requests on S3; per-connection rates from shared requests |
| 2026-10-30 | Each byte in flight once: requests registered by their bytes, reads join a carrying request | Plans made apart re-asked bytes on their way (string heaps fetched twice) | Overlapping requests across plans |
| 2026-10-30 | Surplus pinned to the running read, let go 1 s after the reader goes idle | 64 MiB of surplus dropped 1.4-2.6 s before a heap read asked for it; held, a share's scattered reads kept 65 MB after they ended | Unpinned surplus; surplus kept indefinitely |
| 2026-10-30 | A block (no readable answer after most of a minute's requests were answered) sets a pace: half that minute's rate; a `429` sets none | Zenodo refused this network after ~550 range reads in minutes, inside its stated 133 a minute | Trusting stated budgets alone; learning paces from one `429` (gave 60-120 s paces) |
| Decision (2026-10-30) | A pace loosens by 1.25x per visit after a day without refusals (`PACE_RECOVERY_AFTER_MS`); a block at a loosened pace snaps back and the host stays paced | A transient refusal on a fast host would otherwise slow it for every later visit | A pace that never relaxes; one that relaxes within a day |
| 2026-10-31 | A source on this machine (a loopback origin: `cellucid serve`, Python's and R's served stores; a local file) is priced as this machine at rest whatever its link measured: at no less than a disk's rate, at its shortest round trip, no time a request beyond it; a local file's snapshot carries its priced rate (never below its prior's) | An opening's metadata reads measure a fraction of a loopback host (9-45 MB/s in Chromium); priced at that, a 10M-cell store's first gene was read for 32-68 % of its cells where the last pushed commit read every cell in about 6 s. Floored only until a delivery bounded by transfer, a busy machine still sampled: 5 of 6 openings with ten CPU-bound processes beside the browser measured 0.1-13 MB/s from a slowed read or an 863 ms round trip under load, and read the first gene for 19,779-1,947,096 of 10M cells | Sampling genes on this machine from the opening's lower bound; pricing it as measured once a slowed read counts as filling the link (a transient load switched the first gene into a sample); a notebook kernel priced as this machine (it may run across the internet); a tunnel to a remote server (`ssh -L`) told apart from a server here (both are loopback origins: priced as this machine) |

## Measurements that justify it

| What | Conditions | Result |
|---|---|---|
| S3 round trip under load | CELLxGENE Census, us-west-2, from Europe, 2026-10-29 | 230 ms shortest, 0.35-0.5 s per small request with 42 in flight; an HLCA export priced 198 s at the shortest took 271-386 s |
| S3 request rate | same, 2026-10-30 | 70-88 requests/s with 42 in flight (14 ms each) |
| S3 per-connection throughput | full HLCA edges (703.5 MB), 2026-10-30 | ~2,300 requests of 64-256 KB: 172 s; 320-600 requests of several MB: 252-296 s (about 7 of 42 connections busy) |
| Shared-link samples taken as a connection's | 10M-cell export, emulated 30 MB/s HTTP/2, 2026-10-30 | 255 → 7,556 requests, +18 % time, +45 % CPU, same bytes (fixed by `CONNECTION_BOUND_SHARE` and alone-only samples) |
| One fast delivery taken for the link | 10M-cell CSR `.h5ad`, 1 % share, `.h5ad` export, local web host, Chromium, 2026-10-30 | 473 → 3,135 requests, 6.2 → 18.2 s, same bytes, flipping between runs (fixed by `LINK_EVIDENCE`) |
| Zenodo budget pricing | 868 MB `.h5ad`, 2026-10-28 | open + first view 224 s (496 requests, median 16 KB) → 36 s |
| Zenodo budget pacing | live, 2026-10-28 | 150 reads at once from a window with 133 left: all read in 62 s, 0 failed |
| Refusing host | Zenodo refusing this network (HTML 403 without CORS), and an emulated host refusing after 12 answers, 2026-10-28 | Zenodo, 40 reads: 74 → 8 requests; emulated: 23 → 14 requests, the open fails in 1.3 s with the reason |
| First answer time | 2026-10-28 | Dropbox 2.1-2.2 s, S3 us-west-2 0.7-2.4 s, EBI over 1 s (all healthy) |
| Held pieces | dense `.h5ad`, 100 genes, 150 ms / 5 MB/s HTTP/1.1, 2026-10-28 | 208 → 22 requests, 6.09 → 2.92 s (floor 3.01 s) |
| Metadata block size | 1M-cell `.h5ad` over 30 ms/5 MB/s, 80 ms/2 MB/s, 20 ms/50 MB/s | 16 KiB opened in 0.3-0.4x the time of 256 KiB on the slow links (4.2 → 1.4 s; 7.5 → 1.8 MB) |
| Varying link | Node, 10M CSR `.h5ad`, cycling 20 ms/50 MB/s, 300 ms/1 MB/s and a 4 s stall, 2026-10-30 | every cell 42.6 s / 393 MB; 1 % uniform 31.9 s / 393 MB; 10 % blocks 2.4 s / 64 MB; nothing failed |
| Local file prior | 1 % of 40,000 cells, every gene | 4,004 reads in 270 ms → 534 reads in 150 ms with the disk priced at no less than its prior |
| Revisits from kept parts | Chromium, CELLxGENE atlases (HLCA core, full HLCA, 11.4M embryo), 2026-10-30 | open + whole map: 0 requests beyond the version check |
| Public hosts measured | 2026-10-28 to 2026-10-30 | CELLxGENE S3: HTTP/1.1, 7 addresses; GitHub: HTTP/2, 100 streams, no `Timing-Allow-Origin`; Zenodo: 133 requests/min stated; GitHub API: 60/hour signed out; Hugging Face: 500 API requests and 3,000 file resolutions per five minutes stated |
| Loopback opening rate | 10M-cell `.h5ad` CSR / Zarr v3 sharded on a loopback web host (HTTP/1.1), Chromium, 2026-10-31 | the opening (before the gene sample is priced) measured 45 MB/s over 0.9 MB / 9.3 MB/s over 51 KB, 2.5-5 ms round trips, no delivery bounded by transfer; every cell's first gene then read in 5.1-5.2 s / 4.4 s |
| Loopback opening under load | the same stores and host, Chromium, 2026-10-31, ten CPU-bound processes running beside the browser, each opening right after a run of the last pushed commit | 5 of 6 openings priced the first gene far past its time: Zarr: a slowed metadata read counted as filling the link at 100-1,546 B/ms (`bandwidthMeasured`); `.h5ad`: 13,067 B/ms the same way, or a round trip under load of 863 ms (shortest 10.6 ms); first gene for 19,779-1,947,096 of 10M cells. None at rest (6 of 6 the same hour). Priced as this machine at rest: every cell's first gene, sums exact, in 10 of 10 openings alternating with runs of the last pushed commit (4 of them with four CPU-bound processes beside) and in 9 of 10 with ten (5 of which measured 174-2,304 B/ms from a slowed read); the tenth measured even its shortest round trip at 900 ms (doc 23, known issues) |

The time floor every read is judged against (`cellucid/scripts/efficiency-floor/floor.mjs`, corrected
2026-10-28): `max((dependent + ceil(requests / concurrency)) x RTT,
(dependent + 1) x RTT + bytes / bandwidth)`, with concurrency 6 over
HTTP/1.1 and 100 over HTTP/2. The earlier model let bytes arrive before a
request's round trip; 360 of that round's time misses were the model's.

## Tried and reverted

| What | When | Result | What would have to change |
|---|---|---|---|
| Warm the file's own host before cold S3 mirrors | 2026-10-30 | liver 6.35 → 6.0 s, soundlife 8.0-8.3 → 8.9-9.5 s (Node) | Only a browser A/B (browsers keep idle sockets longer than Node's 4 s) |
| Probe the shipped CELLxGENE release beside the release list | 2026-10-30 | the probe ended 0.2 s later (three cold connections instead of one warm) | — |
| Strict "never re-fetch a gap byte fetched before" | 2026-10-30 | duplicates 0, but HDF5 gene test 5 → 8 requests, ~1,000 more requests over HTTP/1.1 at 10M | A merge price that knows a dependent read will want the gap |
| Price every request at the budget interval from the first answer | 2026-10-28 | Zenodo: half the requests, 50 % more bytes, 7 % slower | Knowing the window's demand early |
| Halve a host's concurrency on a lost request too | 2026-10-28 | no fewer requests reached the refusing host; a new file probed in two round trips | — |
| 64 KiB first head | 2026-10-28 | CELLxGENE/Hugging Face within noise; a Cellucid-ready file 1-2 round trips slower (294 → 451/529 ms at 150 ms / 5 MB/s) | — |
| Metadata reach priced by round trips on every link | 2026-10-28 | opening bytes x2, requests unchanged; kept only where a budget binds | — |
| Sampling `loadedRequestMs` from sequential reads | 2026-10-30 | blocks merged and the cache budget test failed; saturation now also requires requests waiting | — |
| Fixed batches, then a union, then streamed batches for HDF5 string heaps | 2026-10-28 to 2026-10-29 | 100M at 20 ms HTTP/2 71.8 → 87.1 s; replaced by the link's room (`whenRoom`): 10M 1 % at 150 ms HTTP/1.1 77.0 s | — |

## Known limits and open items

- `LinkMeter.throttled()` and its `admitted` cap exist, but no path calls
  `throttled()`: the cap in force after a limit answer is the host book's
  (`hosts.js` `slowDown`), and only for `429` or a `403` naming a limit,
  not for `503`. Either wire it or remove it.
- An opening's few requests barely load a link: EBI's gene sample was
  priced at 20 s and took 28.7 s (55 ms measured at open, 135 ms median per
  request on six busy connections; 2026-10-30).
- A range read cut short is asked whole again (whole downloads resume from
  their byte; ranges do not).
- Not built: handing the reader the head's bytes as they arrive (the 512 KiB
  head on a cold connection costs 5-6 round trips of slow start); a probe
  without `HEAD` on hosts that expose `Content-Range`; a GitHub revisit
  without the `HEAD` and 1 KB read.
- Hosts without a readable version (Dropbox: no `ETag` or `Last-Modified`
  exposed) keep no parts: every visit re-reads.
- Zenodo share comparisons need each run in its own budget window.
- Journals accumulate one directory per burst of reads; compaction is not
  built.
- The share-cost residuals that depend on merge pricing (string heaps at 10M
  without the page cache, the share sample's ceiling as the link speeds up)
  are tracked in {doc}`23_decisions_partial_reading`.
- Firefox's Playwright cannot route a reader worker's requests, so
  `cellucid/tests/browser/host-waits.spec.mjs` runs in Chromium and WebKit.

## Where the code and tests live

| Concern | Module(s) | Tests (`cellucid/tests/`) |
|---|---|---|
| Link model, protocol, concurrency, clock | `formats/io/link.js` | `reader-link-measurement-contract.test.mjs`, `formats-io-contract.test.mjs` |
| Planning, merging, in-flight registry, surplus, metadata cache | `formats/io/byte-reader.js` | `formats-io-contract.test.mjs`, `obs-strings-share-contract.test.mjs`, `subsample-share-cost-contract.test.mjs` |
| Sources, priors, `fetchOnce`, retries | `formats/io/byte-source.js` | `http-request-resilience-contract.test.mjs`, `host-budget-contract.test.mjs` |
| Host budgets, stalls across realms, pace, refusals | `formats/io/hosts.js` | `host-budget-contract.test.mjs`, `host-refusal-contract.test.mjs`, `host-refusal-mid-read-contract.test.mjs`, `public-links-contract.test.mjs`, `browser/host-waits.spec.mjs` |
| S3 addresses and mirrors | `formats/io/s3.js`, `formats/io/mirrors.js` | `formats-io-contract.test.mjs`, `zarr-web-host-guards.test.mjs` |
| Bytes shared between readers | `formats/io/range-cache.js`, `data/byte-range-cache.js` | `formats-io-contract.test.mjs` |
| Kept parts and page checks | `formats/io/kept-parts.js`, `formats/io/page-checks.js` | `kept-parts-contract.test.mjs`, `kept-parts-leave-contract.test.mjs`, `page-checks-contract.test.mjs`, `browser/kept-parts.spec.mjs` |
| Remembered link and pace | `data/host-link-memory.js` | `host-link-memory-contract.test.mjs` |
| Read prices from the link | `data/anndata/subsample-plan.js` (`linkReadMs`), `data/anndata/progressive-open.js` | `subsample-share-cost-contract.test.mjs` |
| Opening rounds per source | `data/sources/` | `public-link-opening-rounds-contract.test.mjs` |
| Floor model and harness | `cellucid/scripts/efficiency-floor/floor.mjs`, `cellucid/scripts/efficiency-floor/link.mjs`, `cellucid/scripts/efficiency-floor/public-links.mjs` | `efficiency-floor.test.mjs` |

## Related pages

- {doc}`20_read_path_architecture`: module map of the read path.
- {doc}`09_data_loading_pipeline_and_caching`: sources, caches and failures
  in long form.
- {doc}`19_subsampling_specification`: which cells a subsample keeps.
- {doc}`23_decisions_partial_reading`: share pricing, block ladder, kept
  data and revisits.
- {doc}`21_decisions_decoders_and_codecs`: what happens to the bytes next.
- {doc}`27_decisions_memory_model`: the budgets the byte readers and range
  caches hold to.
- {doc}`26_decisions_served_store_and_widgets`: the same-machine sources.
- {doc}`01_codebase_map_and_entry_points`.
