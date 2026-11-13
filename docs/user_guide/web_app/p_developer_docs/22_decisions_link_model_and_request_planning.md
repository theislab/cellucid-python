# Decision record: link model and request planning

This page answers: **"How does the web app measure the link to a file, turn a
reader's byte asks into HTTP requests, and stay gentle with public hosts, and
why is it built that way?"** Its "What it does now" is the one full
description of the byte layer: the link, the byte reader, requests to a web
host, and the bytes kept across readers and visits;
{doc}`09_data_loading_pipeline_and_caching` is the reference that links here.
It records the decisions behind it with their dates, the measurements that
justify them, and what was tried and reverted, so nobody re-tries those
blind. Update this page in the same change as the code in `formats/io/`
(`link.js`, `byte-reader.js`, `byte-source.js`, `hosts.js`, `range-cache.js`,
`kept-parts.js`) or `data/host-link-memory.js`. Paths below are under
`cellucid/assets/js/` unless they start with `cellucid/`.

## At a glance

**Audience**
- Developers changing how bytes are fetched, priced,
  retried or kept.

**Time**
- 30 minutes.

**Prerequisites**
- {doc}`20_read_path_architecture` (which module owns each step) and the
  "Bytes" layer of {doc}`09_data_loading_pipeline_and_caching`.

:::{note}
A date written "by 2026-11-12" is the day the entry was first written into
these records, where the day of the decision itself was not recorded.
:::

---

## What it does now

### The link and the byte reader

1. **A byte source per file** (`formats/io/byte-source.js`). A local `File`
   (`blobByteSource`) or a web file (`httpByteSource`: one `Range` header,
   `cache: 'no-store'`, no other header, so no CORS preflight). Each source
   carries a prior, its `ReadPolicy`: `REMOTE_POLICY` (16 KiB metadata
   blocks, 64 KiB merge gap, 8 MiB requests, 6 connections, 100 streams,
   64 MiB cache, no bandwidth assumed) or `LOCAL_POLICY` (64 KiB blocks,
   4 KiB gap, 32 MiB requests, concurrency 4, assumed 1 GB/s disk, read
   ahead from the start). Until it is measured, a local file is assumed to be
   a disk and a web host nothing. File sizes come from `HEAD` or the
   provider, because many hosts do not expose `Content-Range`. An `.h5ad`
   opens one round trip after its link is followed: the `HEAD` that asks its
   size goes out with a read of its first 512 KiB (`HEAD_BYTES`,
   `probeHttpResource`; GitHub's provider sends it beside its own `HEAD`),
   and the reader starts from those bytes (`createByteReader`'s `head`, which
   also measures the link); a Cellucid-ready `.h5ad` keeps everything
   opening reads there. A global heap collection read through bytes the
   reader holds is read only as far as they go (`cachedEnd`): reading on for
   a round trip's worth would cost the round trip it saves. A host that
   content-codes a file even in a range (`Content-Encoding`, hidden from the
   page; NCBI does so for GEO's `.rds`) is told by the probe (the ranged read
   fails while the whole file answers), and the file is downloaded whole.
2. **S3 mirrors.** A public Amazon S3 object (an unsigned URL of a regional
   S3 endpoint: the CELLxGENE Census bucket every CELLxGENE link reads, or
   any S3 link) is read through every address S3 answers it at
   (`formats/io/s3.js`: virtual-hosted and path-style, each also on the
   dual-stack endpoint and on the region's dashed endpoint `s3-<region>`,
   and virtual-hosted on the global endpoint; one `ETag`, CORS open on each;
   seven host names on the Census bucket in us-west-2), each with its own
   six connections, so a read of many small requests runs on 42 instead of
   six (`httpByteSource` `mirrors`, and the same for a Zarr folder on S3,
   `httpStore` `mirrors`). Each request goes to the address with the fewest
   requests asked of it (`formats/io/mirrors.js`, shared by both); a mirror
   that fails, or answers another status, size or `ETag` than the file's own
   would, is dropped for the rest of the read; one that stalls or cannot be
   reached is dropped at its second failure in a row (a Region without a
   dashed endpoint does not resolve it); one that says it is busy (a 503
   every address of a bucket gives alike) is kept. Either way the request is
   asked of the file's own URL at once, so a mirror adds connections and
   never a failure.
3. **One `LinkMeter` per byte reader** (`formats/io/link.js`; and the request
   layer of a Zarr folder on a web host) measures every request (`begin` /
   `end`) over a sliding window of the latest 16 (`RECENT_REQUESTS`), since a
   connection changes while it is read, as TCP BBR models a path by a
   windowed minimum and maximum. Nothing asks where a file is, except that a
   picked Zarr folder is this machine's disk. A link that gets faster is
   followed at its first faster request, one that gets slower once a window
   of requests saw nothing faster, and one odd request (a stall) moves
   neither.
   - `latencyMs`: the shortest of the latest requests whose time was mostly
     a round trip (their transfer at most half their time). Used for
     shaping (merge gaps, request sizes).
   - `loadedLatencyMs`: the median of the latest requests' time less their
     transfer at their share of the bandwidth, counting only requests that
     found a connection already set up (one that opened its connection
     waited for its handshakes too). `loadedRoundTripMs(link)` is what
     reads are *priced* at (gene sample, Genes read for, exports, the
     progressive first view, the analysis run tuner's read-ahead).
   - `loadedRequestMs`: with every slot busy and requests waiting, the time
     between completions less their bytes' transfer at the fastest delivery
     seen, over a full window: the host's request rate under load. Part of
     `requestCostMs`, so of `requestMs`, every merge gap and `linkReadMs`: a
     read's floor is the longer of its requests at that rate and its bytes at
     the bandwidth, and plans merge through the gaps that pay at that price.
   - `bytesPerMs`: the fastest of the latest deliveries, per request or over
     a stretch of requests continuously in flight. A delivery bounded by its
     round trip (under `TRANSFER_BOUND`, three round trips: small reads, a
     subsample's kept rows) can raise it, never lower it (BBR's
     application-limited samples).
   - `connectionBytesPerMs`: what one connection carries, measured only by a
     transfer-bound request that was alone on the link from start to end.
     Under 80% of the link (`CONNECTION_BOUND_SHARE`) the link is
     per-connection (S3) and `connectionRoundTripBytes` is set; the link here
     is the fourth fastest transfer-bound request end ever measured
     (`LINK_EVIDENCE`), never a single delivery. On a host whose one
     connection carries the link (most web hosts and CDNs) nothing is split.
   - `clockTickMs()`: each request is counted half a clock tick longer
     (`performance.now()` ticks 1 ms in WebKit and Firefox workers, 0.1 ms in
     Chromium), so the shortest of many lies in the middle of where the true
     one does and a same-machine server, ending on the tick it began on,
     never measures a round trip of 0.
4. **Protocol and concurrency.** Concurrency starts at the host's HTTP/2
   stream limit: 100, the minimum RFC 9113 recommends and the value GitHub
   announces. A browser keeps only six connections to an HTTP/1.1 host, and
   the page often cannot see which protocol a host speaks, so until the
   protocol is known requests are dispatched at the host's streams but reads
   are shaped for six (`shapingConcurrency`), and plans are priced under both
   (`concurrencies`): the plan taken is the one nearest the best under
   either. The protocol comes from Resource Timing (`watchProtocol`,
   `nextHopProtocol`) or, where the host hides it (no `Timing-Allow-Origin`,
   as GitHub), from bursts: over HTTP/1.1 the seventh request of a burst
   completes a whole round trip after the first, so the seventh completion
   within half a round trip of the first proves multiplexing, and two bursts
   whose first six end together followed by half a round trip of silence
   prove HTTP/1.1 (the reader drops to six). A request belongs to the burst
   it began in even when its round trip is timed from later (a read that
   first asked the page's byte cache, `range-cache.js`). A request begun with
   every known connection busy waited in the browser's queue and is not a
   round-trip sample. Once multiplexed, the concurrency doubles while
   completions per millisecond rise by half of what doubling could give,
   then stays at the last level that paid. Until a host is proven to
   multiplex, requests beyond six start only while those in flight carry
   less than two round trips' worth of bytes (`slotFree`). Background work
   keeps to `backgroundConcurrency` (one fewer than the connections until the
   host multiplexes, then one fewer than the level that paid: never more
   than five on an HTTP/1.1 host), so a read a view waits for always finds a
   connection. A local file needs none of this: it is read ahead from the
   start, whatever one small read measures.
5. **What follows from the link**: `requestCostMs` = the larger of the host
   budget's price (`hosts.js` `requestCostMs`) and `loadedRequestMs`;
   `requestMs` = the larger of `latencyMs / shapingConcurrency` and
   `requestCostMs`; `mergeGapBytes` = `pricingBytesPerMs` x `requestMs`
   (with a local file's prior as its floor); `maxRequestBytes` = bandwidth x
   round trip, never below the prior: the bandwidth measured is shared by
   every request in flight, so a request that long already keeps the link
   busy for its own round trip, and a longer one only serializes its bytes
   on one connection; `roundTripBytesOf(link)` = the one round-trip-in-bytes
   every cost model uses (`data/sparse-column-index.js`,
   `formats/hdf5/chunk-index.js`, `data/anndata/matrix-access.js`; the sparse
   cost model prices a request at one measured round trip);
   `pricingBytesPerMs` (`LinkMeter.pricingBytesPerMs`) = the bandwidth,
   doubled while requests wait (so a link measured only on small requests
   still fetches the gaps that would pay on a faster one, and measures itself
   again); `readsAhead` at 64 MB/s or more (`READ_AHEAD_BYTES_PER_MS`, about
   the rate at which the reader indexes stored bytes, measured), or always
   for a local file: such a link is not what a read waits for, so the file is
   read as a local one (the gene index starts when the file opens, and every
   value is read with it). Below that, and until the link is measured,
   nothing is read before a view needs it.
6. **The byte reader** (`formats/io/byte-reader.js`, `createByteReader`, one
   per file): reads asked in one task are planned on the next macrotask;
   parts kept in this browser and parts already carried by a request in
   flight are subtracted first (`notInFlight`; a read joins the carrying
   request, and no plan fetches a gap across bytes in flight);
   `planMergeGap` fetches gaps smallest first, pricing each candidate plan
   at `count x max(latency / concurrency, requestCostMs) + bytes / bandwidth`
   under *every* concurrency the host may have, and takes the plan whose
   worst ratio to the best under each is smallest; no gap is merged across
   kept bytes or bytes in flight; long requests are split by
   `requestPieceBytes` (the plan's bytes over `shapingConcurrency`, never
   below one connection's round trip of bytes,
   `LinkMeter.connectionRoundTripBytes`, nor below `mergeGapBytes` or a
   block; at most `maxRequestBytes`) where a connection carries only part of
   the link; queued neighbors ride along with a request about to start while
   round trips outweigh bytes and a gap is no wider than the waiting
   requests' own plan would fetch, whichever reads asked for them
   (`companionsOf`, coalescing to the link). Gap bytes no read asked for are
   kept as surplus within `cacheBytes`, pinned while the read that fetched
   them runs, and let go after `SURPLUS_IDLE_MS` (1 s) with nothing in
   flight, so the next wave of a scattered read finds its bytes. Reads marked
   `cache: true` go through 16 KiB aligned metadata blocks in an LRU cache;
   where a host budget binds, neighboring blocks come along
   (`metadataReach`). `whenRoom` lets a reader gather its own small reads
   while requests queue. The reads of a gene index built over a subsample's
   kept rows (background work the first gene waits for) are merged across
   their gaps by the same price as any read.
7. **HDF5 chunk indexes and Zarr shard indexes.** A chunk B-tree of an HDF5
   file (`.h5ad`, `.h5mu`, `.loom`, a 10x `.h5`) is read from its root as
   reads need it: the levels below the root are read with the table only
   while a level is one wave of the link (its nodes transfer in less than a
   round trip and go out at once), every other node -- internal or leaf --
   when a read first needs the chunks below it, a level a round trip. What a
   dataset stores is estimated from the nodes read so far (every key holds
   the size of the first chunk below it), so pricing a read costs no request.
   Once a tree's leaves are all known, the rest of them follow in the
   background when they transfer in less than a round trip, so a later gene
   is one round trip, not two. The reads asked for in one task (the genes of
   one call, the arrays of one sparse gene) wait for the leaves of all of
   them together (`chunk-index.js`, `indexBurst`), so their chunks go out in
   one task and the byte reader plans them as one set. A read that goes on
   where the last one of the same array ended (rows in file order: a
   cell-major matrix read slab by slab, a subsample's kept rows wave by wave)
   has the leaves after it read ahead in the background, as many as it
   needed, asked for as its own leaves land, in the same task as its chunks
   (a leaf sits among the chunks it indexes, so the next one comes in their
   request when the gap is worth it); a read elsewhere reads nothing ahead.
   Background work waits for a leaf already being read ahead instead of
   asking for it again. What an array stores (`storedBytes`, which cost
   models and progress read before a read) is known without reading a leaf:
   a leaf's parent holds the key of the leaf's first chunk, its stored size
   included, so a chunk whose leaf is not read yet counts that size (exact
   for unfiltered chunks; within a few per cent for gzip ones), and a chunk
   whose leaf is read counts its own. A sharded Zarr v3 array does the same
   with its shard indexes (`formats/zarr/array.js`): a read asks for every
   shard index it needs at once and for the inner chunks of all its shards
   once they are in; on a web host the shard indexes in hand count exactly,
   one more is read, and the other shards' inner chunks count the mean
   stored size of the indexes in hand; a folder or ZIP, which knows its
   shards' sizes, reads every index. Of an unfiltered HDF5 chunk (anndata's
   default sparse `.h5ad`), or a Zarr chunk stored as its values' bytes (no
   compressor or filter, C order: what `cellucid serve` and `cellucid_serve`
   answer), only the run of bytes a selection needs is fetched
   (`partialChunks`; a few columns of some rows of such a chunk read those
   rows' bytes). An `.h5ad` read fetches its chunks in groups of neighbors,
   with as many bytes in flight as {doc}`21_decisions_decoders_and_codecs`
   says keeps decoding busy.
8. **One `abort` listener per `AbortSignal`** (`utils/abort.js`), however
   many reads wait on it: thousands of reads under one signal (the rows of a
   large selection) register in constant time each, and an abort still
   rejects every one at once.

(link-model-shared-reads)=
### Shared reads and the page's requests

The reader client sends the gene reads asked for in one task as one request
per matrix and form (dense or sparse), each gene read once however many
callers asked. Aborting a read rejects the caller at once and sends the
worker `{cancel: id}`, whose reads then stop (a request the network has not
sent is never sent; a gene index other reads still need goes on). A batched
gene request is canceled once every caller in it has canceled. The same rule
holds below it: **a shared read is canceled only when every reader waiting
on it withdrew** (`utils/abort.js`). Reads merged into one request
(neighboring ranges of a file, a Zarr store's chunks), a read two openings of
the same file both wait on, and a read the page's byte cache hands to another
reader each stop for one reader at once -- that reader's wait is rejected --
and the request itself is aborted only when no reader is left waiting. The
worker's open reply and its transfer events say how many reads the source
answers at once (unbounded for a local file, the link's background
concurrency for a web host), which the page's download scheduler follows for
long runs (`setSourceRequestConcurrency`). Reads of many genes or cells at
once (`getMatrixColumns`, `readMatrixRows`, or many single reads issued
together, as analysis does) go out as one set of coalesced, parallel
requests. `readMatrixRows` (a cell export) reads every cell's row in one
call: the rows' pointers once, then parts of about a slab of entries in the
cells' order, each part asked for before the one before it is handed out
(`partial` messages from the reader worker), so the link stays busy across
parts and the caller holds one part at a time; a read that arrived frees its
request's share of the in-flight budget at once (`matrix-access.js`
`pipelined`), so a slow read at the head never holds the reads behind it
below the link's concurrency. An index built ahead of need is background
work: its requests wait while any read a view waits for is queued and never
take the last request slot, so coloring by a field or drawing an embedding is
not held up by it; once a gene waits for the index, the index is no longer
background.

(link-model-zarr-folders)=
### Zarr folders and archives

Over an HTTP/1.1 host a Zarr folder stays where six requests at a time put
it: each chunk is its own object, and nothing merges two objects into one
request. A picked or dropped `.zarr` folder is read file by file
(`formats/zarr/stores/directory.js`), each chunk a
`File.slice().arrayBuffer()`:

- its link is this machine's disk, a local file's policy, and is **never
  measured**: a folder read is one file each, whose time is the file
  system's cost of opening a file plus however long the reader thread is
  busy before it takes the bytes, not the disk's speed;
- 16 file reads are in flight at once, a view's ahead of background work's
  (an index built ahead of need never takes the last slot): enough that
  chunk files arrive faster than the decode pool consumes them; more only
  buffers files ahead of the decoders;
- the inner chunks of a shard (Zarr v3 `sharding_indexed`) a read needs are
  read in one slice per run of neighbors, up to 1 MiB: merging saves a
  read's fixed cost, a longer slice only delays decoding its first chunk;
- copying a chunk into a read walks the fewest axes the layout allows (axes
  both sides lay out as one are merged, axes of extent 1 dropped): an
  embedding chunked by rows is one block move per chunk instead of one per
  row;
- the queue of waiting reads (and of waiting chunk decodes) is kept per
  activity, so starting the next one never scans every waiting read.

An anndata store's matrix is CSR, so its gene index is built in the
background from the moment the folder opens, as for a local `.h5ad`. A
`.zarr.zip` is read through the byte reader as one file. Opening reads its
last 256 KiB (the footer and central directory), and every entry lying in
those bytes is served from them: an archive no longer than that is read
once, in one request, locally or over HTTP. Every entry is checked against
its CRC-32, computed sixteen bytes at a time from 32-bit words.

(link-model-web-host-requests)=
### Requests to a web host

Every request to a web host (a byte range of an `.h5ad` or `.zarr.zip`, a
Zarr store's objects, the size probe, catalog files and provider APIs) goes
through `fetchOnce` and `withRetries` in `formats/io/byte-source.js`
(`fetchWhenNotBusy` in `data/data-source.js` for catalog files and provider
APIs). A request asked once (a person's "Try again" of a catalog, a Cellucid
server's health, info and catalog, whose failures are terminal) goes through
`fetchOnce` too (`data-source.js` `fetchAnswer`). A raw `fetch` is left only
where it has its own semantics, each named with its reason in
`tests/one-implementation-contract.test.mjs` (`RAW_FETCH`): probes whose
answer must not count in the host book (the no-cors reachability probe,
Drive's one-byte held-back probe), OAuth token requests and keepalive
revocations, the model APIs, and the page's own files. A raw `fetch` that
fails rejects with a bare `TypeError`, which the app reads as a problem
inside Cellucid (`isInternalFailure`); a failure that got no answer is told
by `isUnanswered`, never by an error's message. None is left waiting for
ever:

- **Transient failures are asked again**, with waits that are the link's
  own measure, as TCP retransmits: the first is the link's retransmission
  timeout (`ResponseClock.timeoutMs`, below: at least one second), each next
  one twice the last (RFC 6298 5.5), each stretched by up to a half by a
  draw from the request's own URL and range (`retryWaitMs`, `spreadOf`:
  requests a busy host turned away together come back apart, and the same
  request always waits the same), never shorter than a `Retry-After` the
  page may read. A request is reported as failed once it was asked again at
  least three times and its waits add up to 100 seconds (`RETRY_AT_LEAST`,
  `RETRY_FOR_MS`: RFC 1122 4.2.3.5's R1 and R2, the least a TCP connection
  is given before it is dropped). That patience is for a read its host
  already delivered to (`withRetries` `delivered`: a file whose ranges, a
  folder whose objects arrived, a download resumed part way): a busy or
  failing answer to a request whose host has delivered nothing to it yet
  (a catalog file, a provider's API, a file's or a folder's opening) is
  reported once asked again three times (R1 alone, a few seconds), so a
  server answering with errors is said to the person, with **Try again**,
  instead of holding a startup step past its bound; a request that gets no
  answer at all keeps the whole patience either way. Transient: HTTP 408,
  429 and 5xx (but 501 and 505), a 403 that says it is a rate limit
  (Google's APIs name `rateLimitExceeded` in its JSON; GitHub's API sets
  `X-RateLimit-Remaining: 0`, `isLimitAnswer`), and a connection that fails
  or drops. GitHub Pages answers a burst of requests with an occasional 503
  that carries no CORS headers, which every browser reports to the page as a
  failed connection, so a failed connection is transient too.
- **A wait the host names is waited for exactly, and is no failure.** A
  `Retry-After` of a second or more that the page may read (the host lists
  it in `Access-Control-Expose-Headers`), or a request budget with nothing
  left and its reset stated (below), is the host saying when: the request is
  asked again then, however often the host says so, and it counts toward
  neither the repeats nor the 100 seconds. A `503` with `Retry-After` is the
  server preparing what was asked for (`cellucid serve` laying a matrix out
  gene-major answers so, with how far it is in a `text/plain` line, `Retry-
  After` one second until it has a pace): it holds that request only, never
  the server's other keys, and the read's activity carries the server's line
  (`ReadActivity.waiting`), which the gene's notification shows ("Reading
  gene CD3E: X is being laid out gene-major for the viewer (37% done, about
  12 s left).") and the opening progress shows as **Waiting for the host**
  (the transfer log's `waiting`), never as a stall or a retry.
- **A host's request budget paces every request to it** (`formats/io/hosts.js`,
  `realmHosts`, one book of hosts per realm). Where a host states its budget
  on its answers -- `X-RateLimit-*` (`X-RateLimit-Limit/Remaining/Reset`: Zenodo: 133 a
  minute, every range read counted, the reset in epoch seconds; GitHub's
  API: 60 an hour without a sign-in), the IETF `RateLimit-*` fields, or
  their structured form (`RateLimit: "api";r=499;t=227` with
  `RateLimit-Policy`, Hugging Face's) -- requests are sent while what remains
  exceeds the requests on their way, and wait for the next answer, or with
  none on its way for the reset, when it does not. A `429` with a readable
  `Retry-After` holds every request to the host until then, and the host is
  asked fewer at once after it (halved from what was on its way, one more
  per round it answers). A block -- a host that answered this site (most of
  the minute's requests) and then sends nothing a page may read to two
  rounds in a row (Zenodo's and Google's "unusual traffic" pages) -- sets
  the host's **pace** (`learnPace`; a `429` is an answer, handled as above,
  and sets none): requests to it are spaced so a minute holds at most half
  the requests sent in the minute before the block, a further block halving
  again, whatever its headers allow; `requestCostMs` prices each request at
  no less, so plans read through larger gaps in fewer requests. The pace is
  told to the other realms, kept per origin by the page
  (`host-link-memory.js` `rememberHostPace`, `cellucid:host-pace:v1`) and
  seeded into the page's and the reader worker's books when the host is
  opened again (`seedPace`), so a later visit never sends the burst that met
  the block. It recovers slowly: a host that blocked once and has refused
  nothing for a day (`PACE_RECOVERY_AFTER_MS`) starts each later visit a
  step looser (`visitPace`: a quarter more requests a minute, once a visit);
  a block met at a loosened pace snaps the kept pace back to the one the
  first block taught and marks the host as one that blocks again, which then
  stays paced, as does one that blocked twice. A host on this machine
  (`utils/loopback.js`) learns no pace and keeps none
  (`rememberedHostPace`). What one realm learns of a host (its budget, a
  wait it named) is told to the others of this site -- the page, its reader
  workers, other tabs -- on a `BroadcastChannel`, since the host counts them
  all. A request whose answer is optional (GitHub's commit of a file,
  without which the file is read at its branch) is never held for a spent
  budget.
- **A host's budget prices every request to it** (`hosts.js`
  `requestCostMs`, `link.js` `requestMs`): its interval, the window over the
  limit (Zenodo: 60 s / 133, 0.45 s; the window learned from the longest
  reset stated with a whole window), in the share of the window already
  spent -- nearly nothing while most of it remains, so early reads stay
  exact, and the whole interval once nothing remains for what is asked or a
  request has waited for the budget -- however fast the host answers. The
  byte reader's plans (`planMergeGap`), the requests a waiting one takes
  along, the link's merge gap and the cost models (`linkReadMs`: auto-fit,
  the progressive plan) price a request at no less, so as the budget runs
  down gaps worth that long at the link's rate are read through; and the
  metadata blocks around a cached read (an HDF5 header, heap or chunk-index
  node) are read with it, half of what the link moves in that interval each
  side, so the file's structure comes in a few large requests. A host that
  states no budget (or no window) is priced at the link alone.
- **A host that says it is busy is asked less at once.** Each 429 (or a
  403 that names a rate limit, `isLimitAnswer`) halves the requests sent to
  that host at once, from those on their way, for every reader and store of
  this site that asks it (`formats/io/hosts.js`, `slowDown`), and they grow
  back by one for every round of requests it answers (additive increase,
  multiplicative decrease). A 5xx is retried, not counted as a limit.
- **A stalled request is aborted and asked again at once**, except one to
  a notebook's kernel (`formats/io/host-fetch.js`): the kernel answers its
  requests one after another, and none while a cell runs, so asking again
  would only queue a copy; its silence is reported ("still waiting") and
  the request waited for. The link's `ResponseClock` (`formats/io/link.js`)
  applies TCP's retransmission timeout (RFC 6298) to the time requests wait
  for their answers: the timeout is the smoothed wait plus four times its
  variation, never under one second, and three seconds before the first
  answer, which waits for the connection to be made too
  (`FIRST_ANSWER_TIMEOUT_MS`, RFC 6298 5.7 and RFC 1122 4.2.3.1), doubled by
  every stall and computed afresh by the next answer. A request is stalled
  when neither it nor any other request of the link, nor any request to its
  host anywhere in this site, has received anything (an answer, or bytes of
  a body) for a whole timeout: a browser keeps six connections to an
  HTTP/1.1 host for every reader and tab of the site, so a request queued
  behind another reader's (a catalog listing several files of one host,
  GEO's three MatrixMarket files, an old and a new reader worker of one
  file) is waiting, not stalled, while the host delivers (`hosts.js` counts
  each host's deliveries, and tells the other realms of the site at most
  every half second, half the shortest timeout). A range of known length is
  also given the time its bytes take at the measured bandwidth, since a host
  that buffers a range (a proxy, a CDN fetching the object from its origin)
  answers only once it holds all of it. Bodies are read as streams, so a
  body that stops arriving is caught as well as an answer that never comes.
  The file's head, read beside the probe of its size, seeds the timeout with
  how long the host took to answer it (and a file whose size its provider
  gave without a probe, from Zenodo or Hugging Face, from how long the host
  took to answer the provider's own request, `hosts.js` `rememberAnswer`),
  so a host slower than a second to answer is not taken for stalled on every
  first read.
- **A file that changed between its description and its first read is
  refused**: the `ETag` or, from a host that sends none, the `Last-Modified`
  the probe saw (`describeRemoteFile` carries both; with the size they are
  the file's version, `writeVersion`, unless its provider names its content)
  must be what every answer says, a whole download's included
  (`downloadRemoteFile`). The probe's `HEAD` passes the HTTP cache by
  (`no-store`), so the version is the host's at this opening; it goes out
  beside the head read, and a file whose two answers name two versions
  (their `ETag`s, weak marks aside, or `Last-Modified` dates, where both
  state one) was replaced between them: it is probed again, and refused as
  `changed` when the next two disagree too. Nothing keeps a
  description past the opening that made it: a link opened again asks its
  provider again (`link-source.js`), and a dataset catalog describes a file
  at each opening that reads it.
- **An answer of other bytes is refused, never read.** Where the host lets
  the page see `Content-Range` (`cellucid serve`, `cellucid_serve()`, most
  object stores), a 206 must name exactly the range asked for (another range:
  `range-ignored`) of a file of the size it was opened at (another size:
  `changed`); and an answer whose `ETag` differs from the file's, or, from a
  host that sends no `ETag`, whose `Last-Modified` differs, means the file
  changed while it was open (`changed`: "open it again to read the new
  version"), whether the answer came from the file's URL, a redirect or a
  mirror. A Zarr folder's ranges (a shard's index and inner chunks) are held
  to the same `Content-Range` check.
- **A failure that will stay is reported at once** (`worthAskingAgain`).
  The signs that a host will answer later keep the whole schedule: a
  `Retry-After`, a status that says the host or its gateway is busy (408,
  429, 502, 503, 504: S3's `503 SlowDown` has the same body and no
  `Retry-After` every time), a connection dropped in the middle of a body, a
  stall. Any other error status without `Retry-After` that repeats exactly
  (the same status and `Content-Length`: Zenodo answers a missing record
  with the same 500 every time) is reported after its second answer, and a
  request that got no answer at all (no connection, or a host that does not
  let this site read it) after it was asked once more; the momentary 503
  GitHub Pages sends without CORS headers is answered on that second
  request.
- **A host that stops answering is left alone, not asked by every read**
  (`hosts.js`). A request that gets no answer a page may read is a lost
  request. A file's range reads count their losses by address (its URL;
  Drive's file and its details are two addresses, as Google holds back only
  the file), requests on their way together as one round: a second lost
  round in a row (the retries of the first) leaves the address unasked for
  RFC 6298's retransmission timeout, one second, doubled for every further
  lost round up to its 60 s maximum. Meanwhile its reads fail at once,
  unsent (`HostNotAnswering`, `sent: false`), and then one request asks
  again; the first readable answer ends it. A read of a file whose ranges
  (or a folder whose objects) were delivered before waits such a refusal
  out instead of failing (`withRetries` `delivered`, `refusesForNow`):
  asked again after a back-off from the link's retransmission timeout
  doubled to 60 s, uncounted, its activity told "`<host>` stopped answering
  this site, which it answered before … Asking again in N s; what was read
  is kept.", until it is answered or no longer wanted; a long read (a
  graph's edges, genes for all cells) then goes on where it was. While this
  computer says it is offline (`navigator.onLine`, looked at every 0.5 s: a
  dedicated worker is not reliably sent the `online` event), the read waits
  for the connection instead, told "This computer is offline: reading goes
  on as soon as it is back; what was read is kept.", and asks again the
  moment it is back: a loss streak of offline rounds holds its address quiet
  only until then, and no such loss (nor one failing sooner than half the
  host's quickest answer) teaches the host a pace (`hosts.js`
  `lostOnThisSide`). A file that never delivered a range fails at once (CORS,
  a host that cannot serve ranges). A provider's own requests (a record, a
  file's details) are asked as the person asked, each twice. The failure
  says what the page knows: from a host that answered this site before
  (`answeredBefore`), it is not CORS, so "`<host>` stopped answering this
  site, which it answered before: it may be refusing requests from this
  network for now (hosts do so for a while after many requests from one
  network), or failing, or this computer lost its connection. Try again
  later."; a provider's first request without an answer says
  "`<provider>` (`<host>`) did not answer this site. It lets other websites
  read it, so this computer's connection may be down, `<provider>` may be
  failing, or it may be refusing requests from this network for now … Try
  again later." "Every host failing at once" is this computer's connection
  (`connection`, `hosts.js` `connectionDown`): at least two unrelated sites
  (`siteOf`) that had answered this site before lost a request within one
  maximum retransmission timeout, and no host answered since.
- **A failure that outlasts its repeats names the request**, with the
  host's own words when it sent a `text/plain` line:
  `https://…/suo.zarr/obsm/X_umap_3d/3.0: the host could not answer (HTTP
  503). Asked 8 times.` A failure that is not transient (404 on a file, a
  short or changed answer, a range ignored) is reported at once.
- **A host that does not give what was asked for is said by whose problem
  it is** (`HOST_FAULTS`, `hostFaultSentence` in `formats/io/byte-source.js`;
  `hostFaultFailure`, `hostFaultOf`, `hostFaultOfResponse`,
  `describeHostFault` in `data/data-source.js`). A host-side status
  (`isHostSideStatus`: any 5xx, 408, 429) that outlasts its repeats is one
  sentence for the person (`hostFaultSentence`, `withRetries`): "`<host>` is
  not answering (HTTP 503, service unavailable) for
  `suo.zarr/obsm/X_umap_3d/3.0`: a problem on `<host>`'s side, not with the
  link or with Cellucid. It usually passes; try again later." (429: "is
  limiting requests from this network … Try again in 30 s, as `<host>`
  asks"), the requests above under its **Details** (its cause). A request no
  answer arrived for keeps its own words (they name its address, for
  **Details**) and carries the fault where it is known (`unansweredFault`):
  this computer `offline`; every host asked lately silent (`connection`); or
  a host that answered this site before now `silent` ("`<host>` stopped
  answering this site, which it answered before
  (`garcia.zarr/.zmetadata`): it may be refusing requests from this network
  for now …, or failing, or this computer lost its connection. Try again
  later."), so a sentence names it without its address. A file's probe and
  whole download (`remote-file.js`) say every class the same way
  (`missing`, `refused`, `refusing`, `offline`), every provider's own
  requests do (`hostFaultFailure`, `hostFaultOfResponse`), and a link names
  the provider the host serves (`link-source.js` `providerFault`: "NCBI GEO
  (`ftp.ncbi.nlm.nih.gov`) is not answering …"), turning a 401/403 of a
  provider whose files are all public (`capabilities.public`: GEO, the
  Census) into `refusing`. Faults that may pass are never settled: the panel
  over the empty viewer offers **Try again** first. A provider whose listing
  a failing host still answers asks before it says a file is missing: a GEO
  series' `filelist.txt` answering a host-side status, or `geo/README.txt`
  (one byte, asked only then) doing so, is said as NCBI failing
  (`geoFailing`, `providers/geo.js`).
- **A sign-in never follows a redirect to another host.** A redirected URL
  is reused until it fails or its stated lifetime ends; reused on another
  host than the one asked, it is asked with `Range` alone, so a token given
  for Google Drive or Microsoft Graph never reaches a CDN it redirected to
  (a browser drops `Authorization` on a cross-origin redirect it follows
  itself).
- **A whole download** (a host that serves no part of a file, or a
  compressed copy read from its start) is asked as any request is, its body
  watched for silence as it is written to storage (`fetchOnce`'s `stream`):
  a busy answer is asked again, a host that stops sending fails it as
  stalled instead of holding the opening, and a body cut short goes on from
  the byte it stopped at (`Range: bytes=<received>-`, under the `ETag` and
  date it started with) where the host serves ranges, so a drop on an
  unsteady link never starts a large compressed copy over; only from a host
  that ignores ranges or content-codes the file is it downloaded again from
  its start, for as long as a request would be asked again. Each repeat is
  said in the Transfers log.
- Each repeat is reported to the transfer log while it happens (`retrying`,
  naming the URL; the console's `[cellucid:transfer]` warning; the Dataset
  box's one "Retrying a read" notification, closed when bytes arrive again),
  and bytes of a request still under way are reported as they arrive
  (`arriving`, at most every 250 ms), so a read of one large request shows
  progress before it completes. The bytes also count into the read's
  activity as they arrive (the `[cellucid:read]` progress, its rate and time
  left), and an attempt that fails has its bytes taken back, so what is
  fetched never passes what was planned before it has arrived.
- Cancellation: the caller's signal aborts the request in flight and any
  wait to ask again at once, and nothing is asked after it. A range read's
  planned bytes are counted once, however often its request was made.

(link-model-kept-bytes)=
### Bytes kept across readers and visits

1. **The page's byte cache.** The page keeps the bytes of remote files its
   reader workers fetched (`data/byte-range-cache.js`, bounded by its
   `fetchedRanges` memory share, least recently used ranges dropped first),
   and every reader worker reaches it through a `MessagePort` posted when
   the worker starts (`formats/io/range-cache.js`), closed when the reader
   closes or the page is left (`AnnDataReaderClient`); the cache itself
   lives as long as the page. A remote file's reads ask the page for the
   bytes it holds first and fetch only the rest, which the page then keeps
   (except bytes only a cell export reads, `ReadActivity.readOnce`: read once
   for a file written out, they would only push out what other reads use, at
   a copy and a message each), and finds by binary search, under a key
   naming the file's bytes exactly (its identity's key and version, or its
   URL, size and `ETag`; a file known by no validator is not kept). A reader
   asks the page once a file which ranges it holds of it, adds the ranges it
   hands the page itself, and asks for a read's bytes only when the read
   touches them: a read of bytes the page does not hold goes to the host at
   once. A read that touches them waits for the page's answer no longer than
   one round trip of its link (before the link is measured, no longer than
   the host is given to answer), since the page answers from its main
   thread, which a long task may hold; then it fetches, and the late answer
   is dropped. A canceled read stops waiting at once, and the request's
   round trip is timed from when it went to the host. So a dataset opened
   again in a new worker (a worker the page ended and started again, the
   same file in another reader) fetches nothing the old one fetched. A read
   answered in part from the page measures no link. Pieces the page holds
   split a read into the parts it does not, each a request of its own; a
   held piece smaller than the gap the link makes worth fetching to save a
   request (`planMergeGap`, the byte reader's own pricing, once the link is
   measured) is fetched again with its neighbors instead (`withRangeCache`).
2. **A Zarr folder's objects in the page.** A Zarr folder on a web host
   (every object a file of its own: the sample datasets, GitHub, Hugging
   Face, S3 and GCS folders, a Python or R server's store) keeps its objects
   in the same page memory, for the tab, under a prefix naming the store's
   version: its address (the link's identity, never its credentials) and the
   version its provider names its content by (a GitHub commit, a Hugging Face
   git tree: shared even where the host shows no validator, as GitHub's raw
   host does not), else the validator its host gave its consolidated
   metadata (`range-cache.js` `storeObjectsKey`); a folder named by neither
   keeps nothing. Its metadata is read from the host at each opening (from
   this browser's kept copy only under a version that names its content,
   below), so a rewritten store is read under a new prefix, and its old
   objects only age out. A reader asks the page once a store which objects
   it holds (not once an object), then asks for an object only when the page
   holds it, within the same one round trip; whole objects, ranges of shards
   and shard indexes are kept apart (`storeObjects`). So an embedding or
   field read again, the first chunk the counts probe read, and every
   opening of the store in another reader (a share raised or lowered) fetch
   nothing the tab fetched before.
3. **Downloaded parts, kept across visits.** Beneath the page's memory, what
   a reader fetches of a remote file is kept in this browser across visits:
   its **downloaded parts** (`formats/io/kept-parts.js`, `withKeptParts`, in
   the reader worker). Every range fetched is kept as it arrives, whatever
   the visit's length: each reader keeps a **journal**
   (`transcode/space.js` `openKeptJournal`), one directory whose manifest is
   written first, so it is kept from its first byte; each piece is appended
   to its `d.bin` through a synchronous access handle, then a record of
   where it is of the file (`j.bin`: framed by its length and CRC-32, with
   the piece's page checksums). A read does not wait for its write, and a
   write is the operating system's once it returns, so a tab closed the
   moment its reads ended, a worker ended from outside or a browser that
   crashes loses at most the pieces still waiting to be written, which the
   next visit fetches again; a machine that loses power may lose what the
   disk had not yet stored, whose records then do not read whole or do not
   match their checksums and are not used (`parseKeptJournal` stops at the
   first). The page asks its reader workers to have the disk hold what they
   wrote when it is hidden or left (`flushKeptParts` on `visibilitychange`
   and `pagehide`; a reader closed, or a page left, asks its worker to leave,
   `{leave: true}`, which stops its reads, flushes and ends itself). A
   journal is closed, never appended to again, once its reader's fetches go
   quiet for the link's retransmission timeout, when the page is hidden or
   left and when the reader closes (the next piece begins a new one): a
   journal still being written (its writer holds its writing lock: another
   tab, or the previous generation of the dataset in another worker during
   a read) is left to the next opening, because a file another worker still
   appends to cannot be read in place. A file whose reads name no version
   (Google Drive's public reads) records a confirmation once its provider
   confirms the version it opened as still holds (asked once its reads go
   quiet for the link's retransmission timeout and when the page is hidden,
   `sources/access.js` `version`), and a later opening uses only the pieces a
   confirmation follows. Every piece has a CRC-32 of every 16 KiB page
   (`page-checks.js`); every journal is filed under `{parts: identity,
   size}`: the file's identity as this opening established it (its version,
   {doc}`09_data_loading_pipeline_and_caching`, *When a dataset changes*) and
   its exact size.
4. **The order a read is answered in.** A read is answered from the reader's
   own memory (cached metadata blocks, the surplus of merged requests: within
   its budget, plus the requests a read still running made, held until that
   read ends, `ReadActivity.finished`, and counted under `fetchedRanges`, so
   a heap collection fetched inside a reference chunk's request is never
   fetched again before its references are decoded; what all of them hold is
   bounded by the memory model's `fetchedRanges` budget, the least recently
   used going first, `configureByteReaders`) first, then from the kept parts
   (a read of this machine's disk, waited for only once, before the reader's
   first plan: `kept.ready`; the kept ranges are held in order,
   `SortedRanges`, so finding them and every lookup cost a sort and a search
   whatever their number), then from the page's memory, and only what none
   holds is fetched. The byte reader's plans treat kept bytes as present
   (`byte-reader.js`): no request asks for them, two ranges are never merged
   across kept bytes, a partly kept read fetches only its missing pieces, and
   a metadata read whose block holds kept bytes is read exactly (a whole
   block would fetch them again). `ByteReader.bytesToFetch(start, end)` says
   what a range still costs (neither kept nor cached; `keptReady` resolves
   once the kept parts are known), and `storedBytes(selection, {toFetch:
   true})` (`formats/store.js`) what an array read would fetch: the cost
   models price with them (the first gene's notice, auto-fit by gene time, a
   progressive opening's first view).
5. **What keeps no parts.** A file a transcoder reads whole (an R object, a
   MatrixMarket file) keeps no parts: its converted store is what is kept of
   it. Nothing is kept as parts of a source on this machine or a notebook
   kernel (`kept-parts.js` `sameMachine`: a loopback address — `localhost`,
   `127.x`, `::1`, `cellucid serve`, a tunnel — or the widget's host,
   `EMBEDDED_HOST_ORIGIN`), nor of a store a Python or R server serves
   (`served-anndata.js` `keptParts: false`; a store served from an object in
   memory is named anew every run): a copy spares no download.
6. **A Zarr folder's objects across visits.** A Zarr store folder on a web
   host keeps its objects (`formats/zarr/stores/http.js` `keeping`: whole
   objects and a shard's ranges, each with the validator the host answered
   it with; under such a version nothing of it changes, so its consolidated
   metadata is kept too, and the objects its host answered as absent, which
   a later visit does not ask for again: its only request is the provider's
   version) only under a version that names its content: a Hugging Face
   folder, named by its git tree (`hf-tree:<tree id>`, the folder's `oid` the
   Hub's `paths-info` answers: a hash of every file below it by its own
   content name, so a commit elsewhere in the repository leaves it),
   threaded from the provider through `link-source.js` and
   `loadRemoteFolder`; a GitHub folder, read at the commit that last wrote
   anything below it (`commit:<sha>`, one request to the commits API) and
   filed under its link's address (`key`: the ref the link names, so another
   commit's kept objects are stale). A folder on a plain host (GitHub Pages,
   S3, GitHub by branch when its API refuses) keeps nothing across openings
   (identity `{kind: 'url'}`, its opening's only): its consolidated
   metadata's validators say when the metadata was written, not that every
   object is as it was (`rclone`/`gsutil rsync` and in-place writes rewrite
   objects under the same metadata).
7. **Witnesses of a store's bytes.** While a folder is open, an object
   answered with another validator than before (in this opening, or kept) is
   a question its bytes answer (`formats/io/read-witness.js`): GitHub Pages
   gives every file of a site new validators at every deployment although
   their bytes are the same. The same range read before is compared with a
   checksum of what it held; a range or an object never read has nothing of
   its own, so the consolidated metadata and one read of each array read
   before (at most 64 KiB each, one `HEAD` when the deployment is known
   already) are asked for again and compared, also when an object never read
   says it was written after the store was opened. The same bytes read on,
   and the store is still the one it was opened as (`unchanged`). Any other
   bytes, another size or an object gone are refused as changed, and so is
   every later read of that opening, so nothing mixes two versions; a retry
   cannot help then (`isChangedOnHost`, `data/data-source.js`), so a field's
   or a growth step's failure says to open the store again. What this cannot
   see, and no client can: a rewrite that keeps every byte read so far and
   changes only bytes never read. A file on a web host is the same
   (`byte-source.js`): its first, middle and last reads are its witnesses,
   the head read beside its probe among them, and with nothing read at all a
   new validator is refused.
8. **A probe when something is kept.** When anything is kept of a file at
   its address (`space.js` `keptAtAddress`, from the directory names alone),
   its probe asks only for the version (a `HEAD`), not its first 512 KiB,
   which then come from what is kept when the version is the same
   (`data/sources/remote-file.js`). Google Drive and OneDrive read their
   file's head beside their own API request (it proves the file can be read
   without a sign-in, and resolves the download link): one byte of it when
   something is kept (`headBytesFor`). Figshare reads its head after its item
   names the file, one byte of it when something is kept; a link naming the
   file by its id reads the head beside the item, whole (512 KiB transfer in
   less than the round trip a read after the item would cost).

(link-model-remembered-links)=
### Remembered links

An opening sizes its gene sample (by gene time, `subsample-plan.js`
`geneTimeFit`) and a progressive first view at the speed of its link. Its
first requests are small, a round trip each, so what it measures is a lower
bound. So what each host delivered is kept across visits
(`data/host-link-memory.js`, `localStorage` `cellucid:host-link:v1`, see
{doc}`04_configuration_env_vars_and_feature_flags`): by origin
(`scheme://host[:port]`, never a path, a dataset, a query or a token), the
round trip, the round trip under load (`LinkMeter.loadedLatencyMs`: a
revisit's first reads are priced at it, scaled to the opening's round trip
and never below its own, until the opening rests on as many bytes), the
rate, the bytes they rest on (`LinkMeter.deliveredBytes`) and when they were
measured, plus a blocked host's pace (`cellucid:host-pace:v1`). The rate is
what the reads sustained: bytes delivered over the time any request was in
flight (`LinkMeter.busyMs`), round trips and all, never the link's fastest
moment (`bytesPerMs`), which overstates a long read. A visit's measurement is
blended into the record, each weighed by its bytes, the record at half its
weight (`rememberHostLink`); at most `HOST_LINK_RECORDS` (32) hosts are kept,
the least recently measured dropped. The page records it from the reader's
link measurements and hands the host's record to each opening (`hostLink`);
the reader prices with `pricingLink(snapshot, record)`: the opening's own
measurement once it has measured its bandwidth with a delivery bounded by
transfer (`bandwidthMeasured`: a slow day is priced as slow within a few
requests, never at a remembered fast rate), or once it rests on as many
bytes as the record; before that the record's rate scaled by the ratio of
the two round trips (a longer path delivers less) and never below the
opening's own rate. The opening's round trip is one measured on a connection
already made (`LinkMeter.warmLatencyMs`: a request begun after another
completed), never its first requests', whose time holds the connection's
set-up; until it has one, the record's rate stands as it was measured.
`pricingLink` lifts only `bytesPerMs` and `latencyMs` (together, when the
record's rate scaled to this path's warm round trip is higher) and
`loadedLatencyMs` (scaled the same way, when higher), sets `measured` and
`remembered`, and returns the fresh snapshot unchanged once the visit has
filled the link or read as many bytes as the record rests on. It does not
recompute the snapshot's `roundTripBytes`, `mergeGapBytes`,
`maxRequestBytes`, `readsAhead` or `requestCostMs`: the byte reader's own
planning never sees the record. `rememberedHosts` lists the hosts, and
`forgetHostLinks` removes them.

A host on this machine (`onThisMachine`: a loopback origin -- `localhost`,
`127.0.0.0/8`, `[::1]`: `cellucid serve`, Python's and R's served stores, a
web server run here) is marked in the seed (`hostLinkSeed` `thisMachine`),
and `pricingLink` prices it as this machine at rest, whatever its link
measured (`thisMachinePricing`, marked `thisMachine`): at no less than a
local file's prior (a disk, `LOCAL_POLICY.bytesPerMs`, 1 GB/s; lifted
whenever below it, `bandwidthMeasured` or not), at its shortest round trip
(`latencyMs`, never its round trip under load: `loadedLatencyMs` lowered to
`latencyMs`) and at no time a request beyond that (`requestCostMs` 0). What
its reads measure past that is the machine's own load, never a network. A
local file's own link is priced the same way: its rate by the meter
(`LinkMeter` `_pricedBytesPerMs`, carried by `snapshot().bytesPerMs`, never
below its prior's), the rest by `pricingLink`, which its snapshot asks for
(`snapshot().onThisMachine`: its prior is a disk). A tunnel to a remote
server (`ssh -L` to a loopback port) is priced so too. A notebook kernel's
host (`EMBEDDED_HOST_ORIGIN`) is not marked: a kernel may run across the
internet (Colab, a JupyterHub).

## Decisions and why

### The link and request planning

| Date | Decision | Why | Rules out |
|---|---|---|---|
| — | Nothing asks where a file is; every shape follows the measured link, starting from the source's prior | One code path for a disk, `cellucid serve`, a kernel and the internet | Per-host or per-provider tuning tables |
| — | A single `Range` header and `cache: 'no-store'` on every range request | CORS-safelisted (no preflight some CDNs refuse); Chromium lets one request per URL through its cache at a time, so ranges of one file asked for together went out one round trip after another (Measurements) | Conditional headers on reads; the HTTP cache for ranges |
| — | 16 KiB remote metadata blocks, 512 KiB first head (`HEAD_BYTES`) | anndata keeps each element's header, B-tree node and heap within a few KiB beside its data, so larger blocks fetched mostly data; the head of a Cellucid-ready `.h5ad` holds all of its opening | Larger blocks (fetch mostly data); a smaller head |
| by 2026-11-12 | Plans are priced under every concurrency the host may have until its protocol is known, and the plan nearest the best under each is taken | Merging for six fetched 2.6× the bytes over HTTP/2, merging for a hundred costs round trips over HTTP/1.1 (Measurements) | Assuming HTTP/1.1 or HTTP/2 before the host shows it |
| by 2026-11-12 | The request size is one round trip's worth of the whole link | 24 connections at 50 MB/s and 150 ms asked 180 MB requests before; 7.5 MB now: a longer request only serializes its bytes on one connection | Requests sized per connection |
| by 2026-11-12 | A request belongs to the burst it began in, even when its round trip is timed from later | Timing those from the host's answer had put every completion outside its burst, so no host was ever taken for HTTP/1.1 by bursts | Bursts made by answer time |
| by 2026-11-12 | An HDF5 chunk B-tree is read from its root as reads need it, levels below the root with the table only while a level is one wave; what an array stores is counted from the keys in hand | Every leaf's parent read with the table cost a CELLxGENE atlas 1,099 requests before its first view, and reading every leaf before the first gene 11,706 requests in 10 minutes (Measurements) | Reading a chunk index whole at open |
| by 2026-11-12 | Reads of one task wait for the leaves of all of them together (`indexBurst`) | Their chunks then go out in one task, planned as one set instead of wave by wave as each read's leaves land | Each read planned when its own leaves land |
| by 2026-11-12 | While requests wait, gaps are priced at twice the measured bandwidth | A link measured only while its connections were full of small requests measured its round trips, not its bandwidth, and priced at that no gap would ever pay | Pricing gaps at a bandwidth measured on round-trip-bound requests |
| by 2026-11-12 | A public S3 object is read through every address S3 answers it at | S3 speaks only HTTP/1.1, and a browser keeps six connections to one host, so a read of many small requests waited on six round trips at a time; over seven hosts it runs on 42 | Six connections for every S3 read |
| 2026-10-28 | Host request budgets read from headers and price every plan (`requestCostMs`) | A request on Zenodo costs 60 s / 133 = 0.45 s once the budget binds; pricing it at latency / concurrency sent hundreds of tiny reads | Treating requests as free on rate-limited hosts |
| 2026-10-28 | Page-held pieces under the merge gap are fetched again with their neighbors (`withRangeCache`) | Holding pieces split one 8 MiB plan into ~190 requests | Splitting a request around every held piece |
| 2026-10-29 | Price reads at the round trip under load (`loadedRoundTripMs`); shape with the shortest | S3 answers in 230 ms idle and 0.35-0.5 s with 42 in flight; an export priced at the shortest said half its time | Pricing at the idle round trip |
| 2026-10-30 | Windowed link model (16 requests: min round trip, max delivery; application-limited samples never lower bandwidth); `pricingLink` gives way once the visit fills the link | Follows a link that gets faster at once, a slower one within a window, ignores one odd request | Lifetime averages; a remembered rate overriding a link now measured |
| 2026-10-30 | Host request rate under load (`loadedRequestMs`) in `requestCostMs` | CELLxGENE's mirrors complete 70-88 small requests a second with 42 in flight (14 ms a request) where round trip / connections said 6 ms | Pricing a request below what the host completes |
| 2026-10-30 | Split long requests over every connection where one connection carries only part of the link (`requestPieceBytes`), judged only from requests alone on the link, and only on four transfer-bound request ends (`LINK_EVIDENCE`) | S3 delivers per connection; samples taken while sharing the link split a plan into ever smaller requests, and one delivery timed over a stretch little longer than the clock's tick can be any rate | Few multi-MB requests on S3; per-connection rates from shared requests; one delivery taken for the link |
| 2026-10-30 | Each byte in flight once: requests registered by their bytes, reads join a carrying request | Plans made apart re-asked bytes on their way (string heaps fetched twice) | Overlapping requests across plans |
| 2026-10-30 | Surplus pinned to the running read, let go 1 s after the reader goes idle | 64 MiB of surplus dropped 1.4-2.6 s before a heap read asked for it; held, a share's scattered reads kept 65 MB after they ended | Unpinned surplus; surplus kept indefinitely |
| 2026-10-31 | A source on this machine (a loopback origin: `cellucid serve`, Python's and R's served stores; a local file) is priced as this machine at rest whatever its link measured: at no less than a disk's rate, at its shortest round trip, no time a request beyond it; a local file's snapshot carries its priced rate (never below its prior's) | An opening's metadata reads measure a fraction of a loopback host (9-45 MB/s in Chromium); priced at that, a 10M-cell store's first gene was read for 32-68% of its cells where the code before the change read every cell in about 6 s. Floored only until a delivery bounded by transfer, a busy machine still sampled: 5 of 6 openings with ten CPU-bound processes beside the browser measured 0.1-13 MB/s from a slowed read or an 863 ms round trip under load, and read the first gene for 19,779-1,947,096 of 10M cells | Sampling genes on this machine from the opening's lower bound; pricing it as measured once a slowed read counts as filling the link (a transient load switched the first gene into a sample); a notebook kernel priced as this machine (it may run across the internet); a tunnel to a remote server (`ssh -L`) told apart from a server here (both are loopback origins: priced as this machine) |
| by 2026-11-12 | A remembered rate is what the reads sustained (`busyMs`), not the fastest delivery | HLCA's fastest delivery was 6.2 MB/s while its gene reads sustained 1.6 MB/s, and sizing from the fastest gave a later opening's first gene three times its target | Remembering `bytesPerMs` |
| — | Every read through an open store runs within the store's lifetime: `close()` aborts it, the reader closes its store on leave, and a failed opening cancels the requests it started | A source disconnected while the reader kept a built gene index on disk had that background read (`X/indptr/0` of a Zarr folder over HTTP) go on with no signal of its own, retrying a host that had stopped answering for as long as the worker lived; found on 2026-11-01 under the slow-delivery preload, when the contract suite's process never ended | A background read with no owner to stop it |
| — | On a clock that stands still the link is never measured, and the source's prior shapes the queue | On the real clock a slow runner priced the gaps between far-apart reads (the entries of a ZIP archive) as worth fetching, because the link was measured from how long the machine took to read the archive's end, and merged the queued reads | Measuring the link from a clock that cannot tell a read's time |
| — | A page's answer about a held range is waited for one round trip of the link, never a fixed short time | On a link whose round trip is 1 ms, a page slower than 1 ms (a Windows runner, 2026-11-01) was given no time and the whole range was fetched | A fixed wait for the page |

### Requests to a web host

| Date | Decision | Why | Rules out |
|---|---|---|---|
| 2026-10-28 | Stall timeout from RFC 6298, 3 s until the first answer; stalls judged per host across realms | Cold first bytes took 0.7-2.4 s and were aborted at 1 s; a request queued behind another reader's is waiting, not stalled | A fixed timeout; per-request silence as a stall |
| by 2026-11-12 | The head read beside the probe seeds the stall timeout with the host's answer time | S3's colder storage classes take 1–3 s for a first byte (measured on the CELLxGENE bucket): every first read was taken for stalled and asked twice | One first-answer timeout for every host |
| 2026-10-28 | Retries from the link's RTO, RFC 1122 R1/R2 (3 repeats, 100 s) | The link's own measure, not a schedule | A fixed retry-delay list (replaced) |
| 2026-10-31 | A busy answer before the host delivered anything to the read: R1 alone (3 repeats), then reported | A catalog answering 503 was asked for 100 s and outlived startup's 60 s step bound, which then said it "did not answer"; the person sees the server's error with Try again | R2's 100 s for requests a person waits on with nothing on screen |
| 2026-10-28 | An address that loses two rounds is left unasked (1 s doubling to 60 s), one probe at a time | A refusing host was asked twice by every read (40 reads, 74 requests) | Every waiting read asking a host that refuses this network |
| 2026-10-30 | A read that was delivered to before waits a host's refusal out, as long as it is wanted | CELLxGENE's S3 refused this machine part way through HLCA's graph twice on 2026-10-30, and the whole read failed | Failing a long read on a refusal for now |
| Decision (2026-10-30) | Do not assume a stable connection | The owner's links vary and drop | A link measured once and trusted; whole downloads that restart from byte 0 |
| 2026-10-30 | A block (no readable answer after most of a minute's requests were answered) sets a pace: half that minute's rate; a `429` sets none | Zenodo refused this network after ~550 range reads in minutes, inside its stated 133 a minute | Trusting stated budgets alone; learning paces from one `429` (gave 60-120 s paces) |
| Decision (2026-10-30) | A pace loosens by 1.25x per visit after a day without refusals (`PACE_RECOVERY_AFTER_MS`); a block at a loosened pace snaps back and the host stays paced | A transient refusal on a fast host would otherwise slow it for every later visit | A pace that never relaxes; one that relaxes within a day |
| Decision (2026-11-10) | A host on this machine (`utils/loopback.js`: localhost, 127.0.0.0/8, ::1, 0.0.0.0) learns no pace and keeps none (`hosts.js` `learnPace`, `host-link-memory.js` `rememberedHostPace`) | It blocks no network: its silence is a server restarting or down. Measured: lost answers from the local catalog server taught a pace that, kept for a day, read Suo at 5–60 KB/s for 95 s on every visit (3.1 s unpaced, measured the same day) | Pacing it as any host |
| Decision (2026-11-11) | Requests lost on this side of the link teach a host no pace (`hosts.js` `lostOnThisSide`): this computer offline (`navigator.onLine`), or a request failing in less than half the quickest answer the host gave — a host's refusal (an "unusual traffic" page read as no answer for its missing CORS headers) is an answer and comes no sooner than its answers | A link that drops for a moment is no block. Measured: a 3 s drop while Suo opened from GitHub Pages (requests failing at once, `ERR_INTERNET_DISCONNECTED`) taught `theislab.github.io` a 2.4 s pace; the rest of the opening grew at 70 KB/s over a fast link (two steps: 24 s and 17 s), and the pace was kept for a day. Without it: 0.2–0.3 s a step; and an 18 s outage on 3G no longer adds up to 20 s of back-off before Suo resumes | Pacing every lost round; a pace that relaxes on its own (it could meet a real block again) |
| Decision (2026-10-31) | A host that does not give what was asked for is said by **whose** problem it is, at one layer every source goes through (`HOST_FAULTS`, `hostFaultSentence` in `formats/io/byte-source.js`; `hostFaultFailure`, `hostFaultOf`, `hostFaultOfResponse` in `data/data-source.js`): `failing` (5xx, 408), `limiting` (429, a rate limit; the wait the host names is said), `refusing` (401/403 from a provider whose files are all public, `capabilities.public`: GEO, the Census, or from a host that answered this site before), `refused` (the link's permissions), `missing` (404, 410), `offline`; the link source names the provider (`providerFault`, "NCBI GEO (ftp.ncbi.nlm.nih.gov)"); a read that outlasts its repeats on a host-side status is said so by `withRetries` itself, so an opening, a gene and an analysis say the same; the faults that may pass are never settled (Try again first) | NCBI's GEO file service answered 500 (and 403) for every file under `/geo/` on 2026-10-31 (from 16:47 CET, for hours) while its folder listings answered with no files, or some: the live app said "GEO GSE327505 has no supplementary file Cellucid opens" (**Open without this link** first) and "`<host>` failed to give `<file>`. Try again later." — read as a wrong link or Cellucid's failure. Each layer wrote its own sentence, most without the status or the provider | A sentence per provider; a status-only message; a listing taken at its word while its host fails (GEO asks `filelist.txt`, else one byte of `geo/README.txt`, before saying a file is missing) |
| Decision (2026-10-31) | Every surface says a host's answer the same way (`describeHostFault`, `data/data-source.js`): the link source and every read, the sample catalog's notice, a served source's (`cellucid serve`) and a notebook or R viewer's connection (their requests now carry `url` and `status`), the panel over the empty viewer, and a link's notification, which offers **Try again** for a fault that may pass | One sentence per class, wherever the person meets it; a status in the sentence is what tells a person the host is at fault | Per-surface wording ("the server reported a problem of its own" without the host or status) |
| Decision (2026-10-31) | "Every host failing at once" is this computer's connection (`connection`): at least two unrelated sites (`siteOf`: the last two labels of a name, so one operator's addresses count once) that had answered this site before lost a request within one maximum retransmission timeout (60 s, RFC 6298 2.5), and no host answered since the earliest of those losses (`hosts.js` `connectionDown`) | A request is lost by TCP's own measure once its longest retransmission timeout passes; one host's outage leaves the others answering, a dropped connection silences them all within that time. A host that never answered this site is no witness: its silence may be CORS or a mistyped name | A count of failed requests; any silence within an arbitrary window; treating a new host's CORS refusal as the connection |
| 2026-10-23 | Every request into a Zarr store on a web host bypasses the HTTP cache (`cache: 'no-store'`), and the suite asserts the cache mode rather than reproducing the hang | WebKit (Safari, and Playwright's WebKit) never finished reading the bodies of chunks it served from its HTTP cache to the AnnData reader worker: Suo opened from the sample catalog loaded once, then every reload, with the store unchanged, waited for those bodies until the startup step's deadline. A small store (the browser fixtures) does not trigger it | The default cache mode or `no-cache` revalidation for store objects; a query string as a store's version |
| — | A host whose answers no cache may reuse (Google Drive: `private, max-age=0, must-revalidate`, no validators) gives each request in flight its own URL, `<url>&lane=<n>` with `n` the lowest number no request in flight holds, read under `cache: 'default'` (`byte-source.js` `lane`, `data/sources/providers/gdrive.js`) | Chromium's HTTP cache lets one request of a URL through at a time, so ranges of one URL went out one round trip after another; with `no-store` each is preflighted instead, and with lanes each lane's preflight is asked for once and kept (Measurements). No other host names a lane, so nothing is ever answered from the cache | `no-store` ranges on a host that preflights each; a lane on a host whose answers a cache may reuse |
| — | Offline (`navigator.onLine` false), a read waits for the connection, not a back-off, and is asked again the moment it is back; the wait looks every `ONLINE_POLL_MS` (500 ms) | A back-off outlasts the outage: a read waited up to 20 s more after an 18 s outage ended. A dedicated worker is not reliably told by an `online` event in Chromium (measured 2026-11-11), so it looks | A back-off while the browser reports no connection; waiting on the `online` event |
| 2026-11-10 | A sample that stopped answering is said by the host-fault sentence for an unanswered request | It was said with the generic "Try again", its URL-led sentence being unshowable | The generic retry sentence for an unanswered host |
| — | A sample listing that fails to read (a busy host, a file taken down) while the person works on the open dataset leaves the open dataset shown and the picker recoverable | The samples are listed from their files once the open dataset has its first view; one failing used to empty the whole picker with "The dataset list could not be loaded", as if the open dataset were gone, when only that sample's listing failed (`tests/browser/data-source-failure-modes.spec.mjs`) | One failed listing emptying the picker |
| — | A catalog entry's own name is read from the catalog, not only the dataset's id | The picker and the opening say "Opening Garcia", not "Opening garcia", before the file is read (`data/sources/anndata-catalog.js`) | An id shown as the name until the file is read |
| — | A host that sends half a range and stops is reported as stalling, not read as content-coding its files | Read as content-coded, the whole file's first bytes came, so the opening downloaded the file whole, which stalled again and said nothing of a stall | Reading a half-sent range as a content coding |
| — | A server's 409 for a file changed on disk is refused as changed, in the server's words (`changedWhereServed`); the store reads nothing more of either version | A file written again while `cellucid serve` served it had every chunk answered 500 with an HDF5 error, said as the server failing ("try again later"); the server now answers 409 with what to do | Reading it as the host failing |
| — | A failure the person can act on is one sentence naming the file and what to do: a Zarr host answering a missing object with its single-page site (200 `text/html`) is read as absent, as a 404 is; a catalog file answered with a page is "not found", not a network failure; a Loom or 10x file missing the part its layout needs, an HDF5 or Zarr store of no format, and an AnnData store without `obs` are said so | They were read as AnnData and refused as "has no cell embedding" (`tests/reader-resilience-contract.test.mjs`) | A refusal worded by the last check that failed |
| — | A host with a public and a private folder is not said to refuse "this network for now, which it gives everyone" for the private file; a host whose public folder sends CORS headers and whose other folder does not is told apart from one that stopped answering; a Zarr folder on a host without CORS headers is said as a file there is; a dataset switched to after going offline is not said to be refused "for now, as hosts do ... after many requests"; a host that stops sending is said to stall, not to refuse CORS | No host blocks a network with a 401. The store's error crosses the reader worker as a plain `Error` named `ByteSourceError`, which `isUnanswered` did not know, so a missing-CORS folder was said as a host silent. The offline case was measured in a browser test | One sentence for every refusal |
| 2026-10-25 | A Microsoft sign-in window cut off from the opener (Cross-Origin-Opener-Policy) is reported as closed by Chromium and Firefox, not WebKit: the step says so and offers the button again, and a sign-in finished in that window still counts; a work account whose organization lets only an administrator approve Cellucid (seen on one research institute's tenant) is said so | The opener cannot see the window; reading the closed window as a refusal left the person with no way on | Reading a closed window as a failed sign-in |

### Bytes kept across readers and visits

| Date | Decision | Why | Rules out |
|---|---|---|---|
| 2026-10-29 | Kept parts as journals: manifest first, each piece appended with a CRC-framed record | Parts kept in 15 s batches with the manifest last kept nothing from a short visit | Batches written at the end |
| by 2026-11-12 | A reader asks the page which ranges it holds once a file, and asks for a read's bytes only when the read touches them | On HLCA, 99 of the first 245 reads had waited a round trip each for a page busy opening the file | Asking the page before every read |
| by 2026-11-12 | Bytes only a cell export reads are not kept in the page's cache (`ReadActivity.readOnce`) | Read once for a file written out, they would only push out what other reads use, at a copy and a message each: 812 MB for 21,285 cells of HLCA | Keeping every fetched byte in the page |
| by 2026-11-12 | What the readers hold of fetched requests is bounded by the `fetchedRanges` budget, least recently used first (`configureByteReaders`) | A long read such as a cell export asks for each region once and held every request it made until it ended: 210 of 216 fetched MB for a 21,285-cell export of a 10M-cell file | Holding a running read's requests unbounded |
| by 2026-11-12 | A store folder's objects across visits are kept only under a version that names its content (a Hugging Face git tree, a GitHub commit) | A plain host's metadata validators say when the metadata was written, not that every object is as it was (`rclone`/`gsutil rsync` and in-place writes rewrite objects under the same metadata) | Keeping a plain host's objects by its metadata's validators |
| by 2026-11-12 | An object answered with a new validator is checked by its bytes (`read-witness.js`) before it is refused | GitHub Pages gives every file of a site new validators at every deployment although their bytes are the same | Refusing every new validator |

## Measurements that justify it

Machine: Apple M1 Pro, 16 GB. "Simulated" or "emulated" means
`scripts/benchmark-anndata.mjs` or the efficiency-floor harness over a link
of the stated round trip and rate (`--protocol http/1.1|h2`: a host that
answers six requests at once or all of them; in Node, decoding shares the
reader thread). Suo is 562k cells × 8,192 genes.

### The link and planning

| What | Conditions | Result |
|---|---|---|
| S3 round trip under load | CELLxGENE Census, us-west-2, from Europe, 2026-10-29 | 230 ms shortest, 0.35-0.5 s per small request with 42 in flight; an HLCA export priced 198 s at the shortest took 271-386 s; an opening's first requests priced an export's round trip at 883 ms against 365 ms as it ran, each also setting up its connection (so `loadedLatencyMs` counts only requests made on a connection already set up) |
| S3 request rate | same, 2026-10-30 | 70-88 requests/s with 42 in flight (14 ms each, where a round trip over the connections said 6 ms) |
| S3 per-connection throughput | full HLCA edges (703.5 MB), 2026-10-30 | ~2,300 requests of 64-256 KB: 172 s; 320-600 requests of several MB: 252-296 s (about 7 of 42 connections busy): S3 delivers about 0.5 MB/s a connection and tens of MB/s in all, which is why one connection's bound is tested against the link's best (`CONNECTION_BOUND_SHARE`) |
| Shared-link samples taken as a connection's | 10M-cell export, emulated 30 MB/s HTTP/2, 2026-10-30 | 255 → 7,556 requests of 16–40 KB, +18% time, +45% CPU, same bytes (fixed by `CONNECTION_BOUND_SHARE` and alone-only samples) |
| One fast delivery taken for the link | 10M-cell CSR `.h5ad`, 1% share, `.h5ad` export, local web host, Chromium, 2026-10-30 | 473 → 3,135 requests, 6.2 → 18.2 s, same bytes, flipping between runs (fixed by `LINK_EVIDENCE`) |
| Ranges through the HTTP cache | five 100-byte ranges of one URL | five round trips in Chromium, one in Firefox and WebKit; one in all three with `no-store` |
| The head beside the probe | Suo `.h5ad` over 80 ms at 50 MB/s; the app in Chromium | probe and opening 327 → 98 ms (7 requests in three dependent waves → the `HEAD` and the head in one); first embedding request 0.81 → 0.22 s after the link is followed, the dataset ready 1.9 → 1.3 s |
| Metadata block size | 1M-cell `.h5ad` over 30 ms/5 MB/s, 80 ms/2 MB/s, 20 ms/50 MB/s | 16 KiB opened in 0.3-0.4x the time of 256 KiB on the slow links (4.2 → 1.4 s; 7.5 → 1.8 MB) |
| Plans priced under both protocols | 99 genes of Suo's `.h5ad` at 80 ms / 50 MB/s, the protocol not yet known | merging for six fetched 2.6× the bytes over HTTP/2 (9.3 MB for 3.6); over HTTP/2 316 → 186 ms (11.3 → 4.6 MB; one round trip and the bytes take 151 ms), over HTTP/1.1 688 → 780 ms (the quickest plan of single-range requests, simulated, takes 621 ms: six connections, not bytes, set that floor) |
| Learning the protocol | Zarr folder on a GitHub-like host (simulated, 80 ms / 50 MB/s) | open 307 → 203 ms, 57 obs fields 475 → 199 ms, 99 genes 332 → 196 ms; `.h5ad` 1,000 cell rows 5.98 → 1.51 s |
| HDF5 tree at open | a CELLxGENE `.h5ad` | about half of its objects load at open (HLCA's 259) |
| Chunk B-tree read as needed | 11.4M-cell CELLxGENE atlas (`X/indices` in 2.26M chunks), S3 | every leaf's parent read with the table: 1,099 requests and 13.5 s before the first view → 58 index and structure requests (5.2 s); 1,281 → 240 requests for the whole opening |
| Leaves of one task together (`indexBurst`) | 100 Loom genes over 150 ms / 5 MB/s | 11.1 → 7.8 s for every cell, 2.1 → 1.6 s at 1% |
| Stored bytes without the leaves | the 2.3M-cell lung atlas on CELLxGENE (`X` CSR, 470,000 gzip chunks of `indices`, 235,000 of `data`), Chromium | all 12,382 leaves of both, 2 KB each and 850 KB apart, one request each, before its first gene could start: 11,706 requests in the first 10 minutes → none |
| Shard indexes on a web host | anndata's automatic shards (about 11.5M values) | a matrix of 4 billion values has some 350 shards: one request, not 350 |
| Background coalescing | atlas-shaped file over 150 ms / 50 MB/s (24 connections), 20% blocks | first gene 10.6 → 6.9 s, 1,617 → 1,013 requests (more bytes, fewer round trips) |
| S3 mirrors | `.h5ad` export of 21,285 scattered cells of CELLxGENE's HLCA core (Node, the app's reader, six connections a host) | 4 hosts 493 s, 7 hosts 383 s; read as one stream of parts (`readMatrixRows`) 271 s, every connection busy throughout |
| Mirrors and stored bytes together | the lung atlas's first gene over about 2 MB/s and 0.3–1 s to S3 (Node, six connections a host) | in its first 90 s, 34 → 135 MB with mirrors alone; with both changes, in its first 120 s, 2,500 requests of index nodes (12 MB, none of the matrix) → 940 requests and 200 MB of the matrix, the link's full rate |
| Varying link | Node, 10M CSR `.h5ad`, cycling 20 ms/50 MB/s, 300 ms/1 MB/s and a 4 s stall, 2026-10-30 | every cell 42.6 s / 393 MB; 1% uniform 31.9 s / 393 MB; 10% blocks 2.4 s / 64 MB; nothing failed |
| HDF5 string heaps loaded as parts arrive (`formats/hdf5/nodes.js`) | a CSR `.h5ad`'s obs index of strings, 1% uniform | each part's collections loaded as they arrived 85.8 s at 10M over 150 ms HTTP/1.1 (2,067 requests) and 71.8 s at 100M over 20 ms HTTP/2; one union after every part 68.8 s and 87.7 s (a barrier: the link idle while every reference arrived); batches of a round of the link's concurrent requests 80.2 s and 87.1 s; with the room rule (`whenRoom`: load while the link has a slot, gather what arrives while requests queue) 77.0 s at 10M against 84.6 s without the gathering (the Tried and reverted rows below) |
| Local file prior | 1% of 40,000 cells, every gene | 4,004 reads in 270 ms → 534 reads in 150 ms with the disk priced at no less than its prior |
| First head read (`HEAD_BYTES`, 512 KiB) | Suo `.h5ad` (5,103 genes), every link measured | its opening reads 416 KB of metadata, all inside the head; 512 KiB transfer in less than the round trip another read would cost (the slowest link, 150 ms at 5 MB/s, carries 750 KB per round trip) |
| Gene index from a compressed row-major matrix (`READ_AHEAD_BYTES_PER_MS`, 64 KB/ms) | the reader building the index of 562k cells × 8k genes (52M entries, gzip), Apple M1 Pro | the column indices (80 MB stored) in 1.2 s, about 67 MB/s; the indices and the values about 90 MB/s: a link delivering at least this is not what the index waits for |
| Small reads of a local file | a gene's counted entries (a few hundred bytes); a local `.zarr.zip`'s opening reads; Node | they measured a disk at 141 MB/s and the archive's opening reads at 62.7 KB/ms, under the read-ahead line; a 2,000-cell share was sampled at 100-200 cells before the local prior was floored (the same case as the 4,004 → 534 reads above) |
| Surplus of a heap collection | a collection read in merged requests of 17-40 MB of reference chunks, dropped from a 64 MiB surplus 1.4-2.6 s before its references were decoded | 86 MB of a 10M-cell file's names fetched twice; a 1% share of a 1.2M-cell category column read 126.3 MB against 124.5 MB for every cell (why a running read pins its requests) |
| A file on seven S3 addresses (six connections each) | HTTP/1.1, 235 ms; requests dispatched at the streams (100) wait in the browser for one of 42 connections | measured from dispatch they looked like 0.6-1.3 s round trips and a multiplexing host, and Genes read for at 100% was priced about 10× too high; the link meter measures at the connection |

**The reader against an earlier revision** (`scripts/benchmark-anndata.mjs`
on Suo, and the same calls against `cellucid serve` on the same machine),
before → now:

| Store and link | Step | Before | Now |
|---|---|---|---|
| CSR `.h5ad`, `cellucid serve` on localhost, a gene chosen 3 s after opening | first gene / 100 genes | 1,322 / 574 ms | 3 / 22 ms |
| CSR `.h5ad`, 20 ms / 100 MB/s | 99 genes / 1,000 cells | 682 / 4,448 ms | 13 / 1,902 ms |
| CSR `.h5ad`, 80 ms / 50 MB/s | 1,000 cells | 16.7 s | 7.1 s |
| Cellucid-ready Zarr folder, 80 ms / 50 MB/s, HTTP/2 host | 99 genes / 1,000 cells | 1.13 / 21.7 s | 0.29 / 3.7 s |
| Cellucid-ready Zarr folder, 150 ms / 5 MB/s, HTTP/2 host | 99 genes / 1,000 cells | 2.3 / 43.6 s | 1.2 / 25.5 s |
| Cellucid-ready `.h5ad`, 80 ms / 50 MB/s, HTTP/2 host | 99 genes | 1.08 s | 0.33 s |

Over an HTTP/1.1 host a Zarr folder stays where six requests at a time put
it (1,000 cells: 22 s at 80 ms).

### Zarr folders and archives

| What | Conditions | Result |
|---|---|---|
| A folder's link measured | from its first small files (a 21 KB `.zmetadata` in 1 ms) | a disk delivering gigabytes a second fell below the 64 MB/s read-ahead line in most opens, and a CSC matrix's cell rows were indexed without their values (1,000 rows of Suo, Node: 2.3–2.7 s instead of 1.4–1.6 s): hence never measured |
| Fewest axes copied | 24-column `scanvi` embedding, 562k cells | 69 → 59 ms in Chromium; 127 → 109 ms in Node |
| Queue per activity | 25,000 queued background reads | 2.4 s of queue work before (100,000 took over a minute: the pass looked at every waiting task a start, now at each waiting activity), 14 ms now |
| `.zarr.zip` tail | the 39 KB demo archive `synthetic-cell-types-2d.zarr.zip`, opened and every array read | 30 requests and 56 KB before, now 1 request and 39 KB |
| Listing a dropped Zarr folder with calls in flight at once (`data/local-folder.js`) | 22,000-file and 6,272-file Zarr folders | WebKit answers each file on a background thread, so its listing takes about 0.73 of the time; Chromium's file system answers one call at a time however many are asked (its directory handles are the faster way, read where a drop gives one: {doc}`23_decisions_partial_reading`); Firefox answers a file in microseconds, with nothing to overlap |
| Node against a browser | `File` from `fs.openAsBlob`, which opens and reads a file on the reading thread (about 0.35 ms per file) | a folder's 1,580 chunk files of one cell-row read cost about half a second more than the same bytes in one `.zarr.zip` (1,000 cells: 1.34 s against 0.72 s); in Chromium the browser reads files off the worker's thread |

**Zarr stores in a browser** (a Chromium worker, Playwright, headless; the
same calls as `scripts/benchmark-anndata.mjs`, best of three, Apple M1 Pro,
files in the page cache), on Suo written by `cellucid.prepare()` (Zarr v2
zlib 6, v3 gzip 6, v3 sharded, and v2 with 4 KB chunks, 25,260 files) and by
anndata (Zarr v2 Blosc-LZ4, v3 Zstandard, CSR float64), whole runs in
milliseconds, before → now:

| Store | Total | 1,000 cells | Graph | Embedding (24 columns) |
|---|---|---|---|---|
| `.h5ad` (for comparison) | 1,424 → 1,195 | 490 → 485 | 657 → 425 | 59 → 63 |
| Cellucid-ready `.zarr.zip` | 1,576 → 1,165 | 596 → 526 | 611 → 343 | 138 → 78 |
| Cellucid-ready folder, Zarr v2 | 1,236 → 1,082 | 551 → 581 | 440 → 255 | 69 → 59 |
| Cellucid-ready folder, Zarr v3 | 1,197 → 1,027 | 530 → 514 | 410 → 261 | 71 → 52 |
| Cellucid-ready folder, Zarr v3 sharded | 1,308 → 996 | 590 → 488 | 420 → 286 | 71 → 51 |
| Cellucid-ready folder, 4 KB chunks | 4,366 → 3,955 | 3,528 → 3,349 | 482 → 272 | 76 → 57 |
| anndata `.zarr.zip` | 2,220 → 1,780 | 246 → 174 | 551 → 282 | 174 → 141 |
| anndata folder, Zarr v2 | 1,894 → 1,670 | 138 → 138 | 474 → 257 | 85 → 105 |
| anndata folder, Zarr v3 | 2,460 → 2,149 | 229 → 232 | 555 → 357 | 110 → 123 |

(The graph and the `.h5ad` row also move with changes to the reader outside
`formats/zarr` made in the same release.) A Cellucid-ready folder is faster
in total than the same data as `.h5ad` or `.zarr.zip` in a browser, and
within 1.2× of the faster of the two at every step that takes more than a
few milliseconds, except the first obs field (39 ms against 19 ms from the
`.h5ad`; about 34 ms from the `.zarr.zip`, mostly starting the decode
workers). A folder of 25,260 files of 4 KB stays bound by the file count:
1,000 cells take 3.3 s against 0.5 s, the cost of opening every file plus
decoding chunks too small to send to a decode worker on the reader thread.
An anndata store's gene index, built in the background from the moment the
folder opens (it was not, whenever the folder's link had been mis-measured),
makes the first obs field read meanwhile take 39–55 ms instead of 19 ms,
while the index competes for the reader thread.

### Requests to a web host

| What | Conditions | Result |
|---|---|---|
| Zenodo budget pricing | 868 MB `.h5ad`, 2026-10-28 | open + first view 224 s (496 requests, median 16 KB) → 36 s |
| Zenodo budget pacing | live, 2026-10-28 | 150 range reads of one file at once from a window with 133 left: all read in 62 s, the 35 refused (sent before the first answer told the budget) waited for the reset, 0 failed |
| Refusing host | Zenodo refusing this network (an HTML "unusual traffic" page, HTTP 403 without CORS headers, from 2026-10-28 on, after many requests from one network), and an emulated host refusing after 12 answers, 2026-10-28 | Zenodo, 40 reads of one file: 74 → 8 requests to the refusing host; emulated: 23 → 14 requests, the open fails in 1.3 s with the reason |
| First answer time | 2026-10-28 | Dropbox share link 2.1-2.2 s, S3 us-west-2 0.7-2.4 s, EBI over 1 s (all healthy); S3's colder storage classes 1–3 s |
| Injected faults | Suo served with injected faults in WebKit, Chromium and Firefox (the app's own requests: a 503 without CORS headers, a 429 with `Retry-After: 1`, a request never answered, a body that stops after its first bytes, a dropped connection) | each recovers and the dataset opens; a request that is never answered on any attempt ends the startup step with the error that names it |
| Public hosts measured | 2026-10-28 to 2026-10-30 | CELLxGENE S3: HTTP/1.1, 7 addresses; GitHub: HTTP/2, 100 streams, no `Timing-Allow-Origin`; Zenodo: 133 requests/min stated; GitHub API: 60/hour signed out; Hugging Face: 500 API requests and 3,000 file resolutions per five minutes stated |
| Bytes read to open from a link (the table of {doc}`../b_data_loading/14_open_from_link`) | Chromium on www.cellucid.com, 2026-10-31; the OneDrive and Google Drive rows 2026-10-25 (reads the opening still had on their way are not counted) | Tabula Sapiens thymus (a 1.6 GB CELLxGENE file) showed its embedding after reading 2.6 MB; the 21.9 GB Human Lung Cell Atlas all 2,282,447 of its cells after 19 MB, in 17 s |
| Google Drive's answer time | the Drive API with Cellucid's public key, 2026-10-25 | each read answered in about 0.65 s (0.56-0.98 s; a few kilobytes and a few megabytes alike), at about 4 MB/s a request; Cellucid opened each file with the same requests and bytes as from a host answering ranges after the same 0.65 s |
| OneDrive and SharePoint without a sign-in | Microsoft Graph and the OneDrive API, a share link shared with "Anyone with the link", checked 2026-10-25 | HTTP 401: such a link opens only after Sign in with Microsoft |
| Lanes for a Drive-shaped host | 650 ms a request, 200 ms a preflight, Chromium | a signed first gene took 7.0 s with a preflight per request (`no-store`) and 5.3 s with lanes, as a request needing no preflight does |
| A round trip and a request on public hosts | from Europe, read off the real hosts | a round trip is 0.15 to 0.7 s, and a request spends one of a rate-limited host's budget (Zenodo: 133 a minute); the floor of dependent round trips and requests per source is pinned by `public-link-opening-rounds-contract.test.mjs` |
| Zenodo's refusal and record | 2026-10-28; 2026-10-31 | an HTML "unusual traffic" page (HTTP 403 without CORS headers) for a while after many requests, so a page saw only failed requests; on 2026-10-31 the record answered HTTP 200 with `Retry-After: 60` |
| Google Drive's hold-back | seen three times on 2026-10-25 | HTTP 403 without CORS headers, each after some 200-400 MB from one address within minutes, for about 50 minutes; file details still answered and a request carrying a sign-in still reached the API |
| NCBI's file service failing | 2026-10-31 (16:47-17:07 CET, and for hours after) | every file under `ftp.ncbi.nlm.nih.gov/geo/` (`geo/README.txt` and every `filelist.txt` included) answered HTTP 500, and 403 for some, while folder listings answered 200 with no files or only some: hence a listing that lacks what the link names asks `geo/README.txt` for one byte before it is believed. GEO and the Census serve only public files, so a 403 there is the host refusing this network |
| Figshare's `ndownloader` bot challenge | AWS WAF, HTTP 202, 2026-10-23 | Chromium's page passes and its workers do not (a worker's request carries no client hints), WebKit passes from both, headless Chromium never passes: hence the page asks `ndownloader` for the signed link and the reader reads S3 only |
| OneDrive's `$select` | a personal OneDrive's `/shares/{id}/driveItem`, 2026-10-25 | `$select=…,@microsoft.graph.downloadUrl` answered every other field it named but not the annotation: hence the whole item is asked for |
| Hugging Face's repository listing | a large dataset repository (a Zarr v2 folder adds one entry per chunk) | the listing names every file of every folder: 6.9 MB and 4.6 s; hence only the path the link names is asked for, through `paths-info` |
| Resuming a body cut short | GEO's `.h5ad.gz`, 616 MB, read at 1.5-8 MB/s | a download that starts over on every drop never finishes on an unsteady link; `Range: bytes=<received>-` resumes it |
| Public datasets the link tests are shaped after | GEO and Zenodo files | GSE336826 (Cell Ranger ARC, 32,058 genes and 239,160 peaks) failed after its whole download because the MuData `var` index joined each modality's features with one spread call, past the call stack; GSE131907's `raw_UMI_matrix.rds.gz` (a 604 MB data frame) was refused after 194 s and its whole download, though its first word says it is a list; a MatrixMarket size line read the first 12 MB of a 7.6 GB GEO matrix (three 4 MB requests) and its decompressor and reads went on after it; GSE158055's (1.46M cells) `…_covid19_counts.mtx.gz` was "not a file Cellucid opens from a link" |
| WebKit and the HTTP cache | Suo from the sample catalog, a persistent WebKit profile, the store served with GitHub Pages' `Cache-Control: max-age=600`, 2026-10-23 | the default cache mode and `no-cache` revalidation both hung on every reload; `no-store` opened it three reloads in a row. A reload reads the store again (Suo's first view: 5.5 MB) |

### Kept bytes and remembered links

| What | Conditions | Result |
|---|---|---|
| Held pieces | dense `.h5ad`, 100 genes read after a few others, 150 ms / 5 MB/s HTTP/1.1, 2026-10-28 | 208 → 22 requests, 6.09 → 2.92 s (floor 3.01 s), fetching 4% of the bytes again |
| A dataset opened again in a new worker | atlas-shaped CSR `.h5ad` over 150 ms / 50 MB/s, re-opening at 10% after 5% blocks | 1.1 s and 12 requests before, none now; its first gene 6.0 s / 81 MB before, 2.7 s / 34 MB now |
| Kept parts lost on close | Firefox under load, a tab closed the moment its reads ended | one range in about one close in twelve is fetched again next visit |
| Kept ranges in order | a first visit whose analyses read 200,000 ranges | opens again in about a second |
| Revisits from kept parts | Chromium, CELLxGENE atlases (HLCA core, full HLCA, 11.4M embryo), 2026-10-30 | open + whole map: 0 requests beyond the version check |
| An opening's measure | HLCA | the opening measured 244 ms and about 640 KB/s (small reads, mostly a round trip each); its gene reads then ran at 1.5–3.4 MB/s, so the gene-time auto-fit sized from the opening's lower bound opened two to three times fewer cells than the link reads within its target |
| A first visit's link | an 11.4M-cell CELLxGENE atlas, S3 | 1 MB/s at its first resolve and 11.5 MB/s on its first gene: until the link is measured the page's own downlink (`pageDownlinkBytesPerMs`) stands in for the rate |
| A revisit priced from the first requests | GEO | a revisit priced from one cold `HEAD` and one `GET` (376 ms with the connection's set-up in it) read GEO at a third of its rate; round trips are compared like with like (`warmLatencyMs`) |
| The identity of a folder kept as a handle | the SHA-256 of a line per file (path, size, modification time), Node | about 0.4 s for 10^6 files |
| Chromium's off-the-record storage | a private window or a test browser's new context, once the site's share is spent, 2026-10-29 | a write reports every byte and keeps part: a 40,000,000-byte payload kept as 29,060,298 with `navigator.storage.estimate()` still saying 3 GB; a 10M-cell Seurat object's UMAP read back as 94,031,122 of its 160,000,000 bytes ("Offset is outside the bounds of the DataView"): hence every write is checked against the file's size and a short read is refused as a full storage |
| Loopback opening rate | 10M-cell `.h5ad` CSR / Zarr v3 sharded on a loopback web host (HTTP/1.1), Chromium, 2026-10-31 | the opening (before the gene sample is priced) measured 45 MB/s over 0.9 MB / 9.3 MB/s over 51 KB, 2.5-5 ms round trips (1.5-7 ms at rest elsewhere), no delivery bounded by transfer; priced at that, every cell's first gene (5-7 s in fact) was priced at 15-37 s, so the first gene was read for a sample of 22-68% of the cells; priced as this machine, every cell's first gene read in 5.1-5.2 s / 4.4 s |
| Loopback opening under load | the same stores and host, Chromium, 2026-10-31, ten CPU-bound processes running beside the browser, each opening right after a run of the code before the change | 5 of 6 openings priced the first gene far past its time: Zarr: a slowed metadata read counted as filling the link at 100-1,546 B/ms (`bandwidthMeasured`); `.h5ad`: 13,067 B/ms the same way, or a round trip under load of 863 ms (shortest 10.6 ms); first gene for 19,779-1,947,096 of 10M cells. None at rest (6 of 6 the same hour). Priced as this machine at rest: every cell's first gene, sums exact, in 10 of 10 openings alternating with runs of the code before the change (4 of them with four CPU-bound processes beside) and in 9 of 10 with ten (5 of which measured 174-2,304 B/ms from a slowed read); the tenth measured even its shortest round trip at 900 ms ({doc}`23_decisions_partial_reading`, Known limits) |

The time floor every read is judged against (`cellucid/scripts/efficiency-floor/floor.mjs`, corrected
2026-10-28): `max((dependent + ceil(requests / concurrency)) x RTT,
(dependent + 1) x RTT + bytes / bandwidth)`, with concurrency 6 over
HTTP/1.1 and 100 over HTTP/2. The earlier model let bytes arrive before a
request's round trip; 360 of the time misses then counted were the model's.

## Tried and reverted

| What | When | Result | What would have to change |
|---|---|---|---|
| Warm the file's own host before cold S3 mirrors | 2026-10-30 | liver 6.35 → 6.0 s, soundlife 8.0-8.3 → 8.9-9.5 s (Node) | Only a browser A/B (browsers keep idle sockets longer than Node's 4 s) |
| Probe the shipped CELLxGENE release beside the release list | 2026-10-30 | the probe ended 0.2 s later (three cold connections instead of one warm) | — |
| Strict "never re-fetch a gap byte fetched before" | 2026-10-30 | duplicates 0, but HDF5 gene test 5 → 8 requests, ~1,000 more requests over HTTP/1.1 at 10M | A merge price that knows a dependent read will want the gap |
| Price every request at the budget interval from the first answer | 2026-10-28 | Zenodo: half the requests, 50% more bytes, 7% slower | Knowing the window's demand early |
| Halve a host's concurrency on a lost request too | 2026-10-28 | no fewer requests reached the refusing host; a new file probed in two round trips | — |
| 64 KiB first head | 2026-10-28 | CELLxGENE/Hugging Face within noise; a Cellucid-ready file 1-2 round trips slower (294 → 451/529 ms at 150 ms / 5 MB/s) | — |
| Metadata reach priced by round trips on every link | 2026-10-28 | opening bytes x2, requests unchanged; kept only where a budget binds | — |
| Sampling `loadedRequestMs` from sequential reads | 2026-10-30 | blocks merged and the cache budget test failed; saturation now also requires requests waiting | — |
| Fixed batches, then a union, then streamed batches for HDF5 string heaps | 2026-10-28 to 2026-10-29 | 100M at 20 ms HTTP/2 71.8 → 87.1 s; replaced by the link's room (`whenRoom`): 10M 1% at 150 ms HTTP/1.1 77.0 s | — |

## Known limits and open items

- `LinkMeter.throttled()` and its `admitted` cap exist, but no path calls
  `throttled()`: the cap in force after a limit answer is the host book's
  (`hosts.js` `slowDown`), and only for `429` or a `403` naming a limit,
  not for `503`. Either wire it or remove it.
- An opening's few requests barely load a link: EBI's gene sample was
  priced at 20 s and took 28.7 s (55 ms measured at open, 135 ms median per
  request on six busy connections; 2026-10-30). On Census S3 the opening's
  resolve measured 0.96–1.1 MB/s where the first gene then measured 11 MB/s.
- A range read cut short is asked whole again (whole downloads resume from
  their byte; ranges do not).
- Not built: handing the reader the head's bytes as they arrive (the 512 KiB
  head on a cold connection costs 5-6 round trips of slow start); a probe
  without `HEAD` on hosts that expose `Content-Range`; a GitHub revisit
  without the `HEAD` and 1 KB read.
- Hosts without a readable version (Dropbox: no `ETag` or `Last-Modified`
  exposed) keep no parts: every visit re-reads.
- A host naming only a tag made from the bytes shows a rewrite of a store
  only when its consolidated metadata's bytes change: one keeping every shape
  and attribute and changing only values is not seen
  ({doc}`23_decisions_partial_reading`, reading again).
- For up to ten minutes after a store is rebuilt, GitHub Pages' CDN may still
  answer some of its objects from before the rebuild; a chunk that does not
  fit the new metadata is refused naming its URL.
- Zenodo share comparisons need each run in its own budget window.
- Journals accumulate one directory per burst of reads; compaction is not
  built.
- The share-cost residuals that depend on merge pricing (string heaps at 10M
  without the page cache, the share sample's ceiling as the link speeds up)
  are tracked in {doc}`23_decisions_partial_reading`.
- `connection` needs two sites that answered this site before: the first
  opening of a visit with the connection already down is said as its one
  host not answering (or `offline`, when the browser knows).
- Firefox's Playwright cannot route a reader worker's requests, so
  `cellucid/tests/browser/host-waits.spec.mjs` runs in Chromium and WebKit.

## Where the code and tests live

| Concern | Module(s) | Tests (`cellucid/tests/`) |
|---|---|---|
| Link model, protocol, concurrency, clock | `formats/io/link.js` | `reader-link-measurement-contract.test.mjs`, `formats-io-contract.test.mjs` |
| Planning, merging, in-flight registry, surplus, metadata cache | `formats/io/byte-reader.js` | `formats-io-contract.test.mjs`, `obs-strings-share-contract.test.mjs`, `subsample-share-cost-contract.test.mjs` |
| HDF5 chunk indexes, Zarr shard indexes | `formats/hdf5/chunk-index.js`, `formats/zarr/array.js` | `hdf5-reader-structures.test.mjs`, `zarr-format-conformance.test.mjs` |
| Sources, priors, `fetchOnce`, retries | `formats/io/byte-source.js` | `http-request-resilience-contract.test.mjs`, `host-budget-contract.test.mjs` |
| Who is at fault when a host does not give (host faults) | `formats/io/byte-source.js` (`hostFaultSentence`, `withRetries`), `data/data-source.js` (`hostFaultFailure`), `data/link-source.js` (`providerFault`), `data/sources/remote-file.js`, `data/sources/providers/geo.js` (`geoFailing`) | `host-fault-contract.test.mjs`, `host-connection-contract.test.mjs`, `host-refusal-contract.test.mjs`, `startup-data-source-failure-messages.test.mjs`, `browser/host-faults.spec.mjs` |
| Host budgets, stalls across realms, pace, refusals | `formats/io/hosts.js` | `host-budget-contract.test.mjs`, `host-refusal-contract.test.mjs`, `host-refusal-mid-read-contract.test.mjs`, `public-links-contract.test.mjs`, `browser/host-waits.spec.mjs` |
| S3 addresses and mirrors | `formats/io/s3.js`, `formats/io/mirrors.js` | `formats-io-contract.test.mjs`, `zarr-web-host-guards.test.mjs` |
| Zarr folders on disk and on web hosts, witnesses of a store's bytes | `formats/zarr/stores/directory.js`, `formats/zarr/stores/http.js`, `formats/io/read-witness.js` | `formats-io-contract.test.mjs`, `zarr-folder-reads.test.mjs`, `zarr-web-host-guards.test.mjs`, `kept-parts-contract.test.mjs`, `host-redeploy-content-contract.test.mjs` |
| Bytes shared between readers | `formats/io/range-cache.js`, `data/byte-range-cache.js` | `formats-io-contract.test.mjs` |
| Kept parts and page checks | `formats/io/kept-parts.js`, `formats/io/page-checks.js` | `kept-parts-contract.test.mjs`, `kept-parts-leave-contract.test.mjs`, `page-checks-contract.test.mjs`, `browser/kept-parts.spec.mjs` |
| Remembered link and pace | `data/host-link-memory.js` | `host-link-memory-contract.test.mjs` |
| Read prices from the link | `data/anndata/subsample-plan.js` (`linkReadMs`), `data/anndata/progressive-open.js` | `subsample-share-cost-contract.test.mjs` |
| Opening rounds per source | `data/sources/` | `public-link-opening-rounds-contract.test.mjs` |
| Floor model and harness | `cellucid/scripts/efficiency-floor/floor.mjs`, `cellucid/scripts/efficiency-floor/link.mjs`, `cellucid/scripts/efficiency-floor/public-links.mjs` | `efficiency-floor.test.mjs` |

## Related pages

- {doc}`20_read_path_architecture`: module map of the read path.
- {doc}`09_data_loading_pipeline_and_caching`: the reference: sources, the
  start-up path, the caches and failures.
- {doc}`19_subsampling_specification`: which cells a subsample keeps.
- {doc}`23_decisions_partial_reading`: share pricing, block ladder, kept
  data and revisits.
- {doc}`21_decisions_decoders_and_codecs`: what happens to the bytes next.
- {doc}`27_decisions_memory_model`: the budgets the byte readers and range
  caches hold to.
- {doc}`26_decisions_served_store_and_widgets`: the same-machine sources.
- {doc}`01_codebase_map_and_entry_points`.
