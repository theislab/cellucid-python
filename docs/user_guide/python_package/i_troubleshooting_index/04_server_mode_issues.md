# Server mode issues

## Address already in use

Another process owns the requested port. Identify that process or choose a
different explicit port (`--port 9000`, or `--port 0` for any free port and the
Viewer URL it prints); do not terminate an unknown process. A notebook viewer
created without `port=` never fails this way: it takes 8765, or the next free
port when 8765 is taken, and prints a note saying so.

## The server starts but another machine cannot connect

`127.0.0.1` is reachable only from the same machine. Remote access requires an
intentional bind address plus an appropriate firewall, SSH tunnel, or reverse
proxy. Do not expose sensitive datasets to an untrusted network.

## The banner shows a URL I cannot open

A server bound to a wildcard address accepts connections on every interface:
`--host 0.0.0.0`, `--host ::`, and an empty host all mean that. A wildcard is a
statement about which interfaces accept connections, not an address a client can
connect *to*, so it is never rendered inside a URL. The banner reports the
loopback origin as the Local URL and the Viewer URL, and adds the addresses
another machine can use when this machine has a name that resolves off itself.
The banner transcript is in
{doc}`../d_viewing_apis/12_remote_servers_ssh_tunneling_and_cloud`.

Read the banner as follows:

- Through an SSH tunnel, open the loopback Viewer URL, query included. The
  tunnel terminates on your own machine's loopback interface, so that is the
  address the browser needs.
- The `From another machine` block is absent when the machine's own hostname
  resolves to loopback, because such a name addresses nothing outside this
  machine. A loopback bind never prints the block at all.
- An IPv6 bind is bracketed as a URL requires, `http://[::1]:8765`, and a
  loopback origin derived from a `::` bind is `[::1]`, not `127.0.0.1`.

`--host 0.0.0.0` is what a compute node needs when the tunnel terminates on a
different machine, such as a cluster login node; `cellucid serve --help` says
so. It also has no authentication: every user who can reach the port reads the
dataset.

## The browser reports an HTTP error

Record the exact URL, method, status, and response body. The status decides
what to do:

- **404** — the route or the name in it does not exist: an unknown path, a
  route under `/_cellucid/anndata/` other than the one `/_cellucid/datasets`
  names, or a Zarr store key the store does not hold (a Zarr reader takes that
  404 as "absent", so an unwritten chunk is not an error). Verify that you
  supplied a file or folder Cellucid reads (an `.h5ad`, an AnnData or MuData
  Zarr store, `.h5mu`, `.loom`, a 10x `.h5` or `outs/` folder, a MatrixMarket
  folder, an `.rds`/`.rda`), then test the same URL from the browser host.
- **416** — a `Range` request the file cannot satisfy, or more than one range
  in one request. The viewer's reader sends exactly one satisfiable range; a
  416 in the browser points at a proxy that rewrote the request.
- **500** — a genuine internal fault. Read the server log, which carries the
  traceback.

A gene or column the viewer cannot draw (one holding an infinity) is not an
HTTP error: the viewer reads the values itself and says
so in a notification — see
{doc}`../d_viewing_apis/15_troubleshooting_viewing`.

## The process does not stop cleanly

Keep the returned server/viewer handle and use its documented lifecycle method.
In a CLI session, use the terminal interrupt described by the command banner.

## Serving says nothing for minutes

Silence here is work, not a stall. The direct-AnnData server runs four numbered
steps — `[1/4] Detecting format`, `[2/4] Reading metadata`,
`[3/4] Checking launch choices`, `[4/4] Starting server` — and reads only the
store's metadata before it binds, so the last line on screen names the phase
you are waiting on. The full startup transcript is in
{doc}`../../web_app/b_data_loading/04_server_tutorial`.

Where the time goes:

- **Step 1** plans what it serves: an in-memory AnnData, and an AnnData file
  not in priority order, as a Zarr store in priority order whose chunks are
  made as the browser asks for them (`Serving: a Zarr store in priority order,
  made as the browser reads it`; a store holding arrays the browser cannot
  decode is decoded on the kernel this way too, after `Browser cannot decode: …`
  names them); it writes a temporary `.h5mu` for an in-memory MuData. Any other
  format holding such arrays is refused with that reason.
- **Step 4** downloads and verifies the current viewer build from its
  configured source, at every start, before binding; that is where a machine
  with no outbound access stalls. Serve the data endpoints alone with
  `--no-web-ui` on such a machine.

## Cellucid says a package is missing that I know is installed

`cellucid serve` classifies the import condition instead of advising `pip` for
every one of them, so read which of the three lines you got:

```text
Error: Cellucid needs the 'anndata' package to read this input and it is not installed in this environment. Install it and run the command again: pip install anndata
```

The distribution is genuinely absent from the interpreter running Cellucid.
Confirm which interpreter that is before installing anything: `which cellucid`,
then `python -m pip show cellucid` and `python -c "import sys; print(sys.executable)"`
in the environment you believe is active. When a module's import name differs
from the name that installs it, the advice names the installing one: advising
the import name would name a package that does not exist on PyPI.

```text
Error: Cellucid could not read a name it needs from the installed 'anndata' at /path/to/env/lib/python3.11/site-packages/anndata/__init__.py: cannot import name 'read_zarr' from 'anndata'. The installed version of that package is not one this Cellucid can use. Upgrade it in this environment and run the command again: pip install --upgrade anndata
```

The module imported and a name inside it did not. The package is present but too
old to hold the name Cellucid asked of it, so installing the same version again
repeats the failure. The line names the file it read, which identifies the
environment that owns it.

```text
Error: Cellucid needs the standard-library module 'tomllib', and this Python build cannot import it. Python 3.10.13 at /path/to/python3.10 is running Cellucid, which requires Python <3.15,>=3.11. Run Cellucid on a Python in that range. No package installs this module: the interpreter itself provides it.
```

This one is not a `pip` condition at all, and no `pip install` is offered. The
module belongs to the running interpreter, so the message names that
interpreter's version, its executable, and the range this Cellucid supports. Run
Cellucid on a Python inside that range; see {doc}`01_installation_and_dependency_issues`.

## Dataset-shape symptoms

A missing edge overlay, an embedding Cellucid will not read, and an `X_umap` that
used to be refused are covered in
{doc}`../d_viewing_apis/15_troubleshooting_viewing`.
