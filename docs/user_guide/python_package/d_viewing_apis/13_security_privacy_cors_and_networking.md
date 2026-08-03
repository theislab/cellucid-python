# Security, privacy, CORS, and networking

This page documents the security/privacy model of Cellucid viewing workflows, especially:
- what is (and is not) sent over the network,
- why `--host 0.0.0.0` is risky,
- what a wildcard bind publishes on a shared machine, and how to bound it,
- and how CORS behaves.

## At a glance

**Audience**
- Anyone working with sensitive data (patient, unpublished, internal)
- Anyone serving on shared networks
- Anyone serving from a multi-user HPC compute node

## What network traffic happens in the recommended workflow?

### Your dataset

In the recommended workflows (`cellucid serve`, `show`, `show_anndata`):
- your **dataset bytes are served by your local Python server**,
- and your **browser downloads the dataset from that local server**.

Cellucid does not “upload your dataset to cellucid.com” as part of this workflow.

### The verified viewer UI generation

The Python server downloads the **viewer UI assets** from:

```text
https://www.cellucid.com
```

and publishes them as one verified local generation. This lets the UI and
dataset API use the same origin.

Implications:
- every viewer-serving startup requires access to the configured source,
- the complete declared generation is verified before the server starts,
- source failure is reported instead of substituting a previous generation,
- the local files are UI assets, not your dataset.

Configure the local generation directory with `--web-cache-dir PATH` or the
`web_cache_dir=` Python argument.

## The biggest risk: exposing the dataset server

### Default behavior (safe)

By default, the server binds to:

```text
127.0.0.1
```

Meaning: only your machine can access it.

A loopback bind is still not the whole boundary. Any page you have open can
re-point its own DNS name at `127.0.0.1` (DNS rebinding), after which the
browser treats this server as same-origin and applies no CORS check. Cellucid
closes that by validating `Host` on every request before routing: only one
well-formed authority naming `localhost` or a loopback IP literal on the bound
port is served, and everything else gets `421 Misdirected Request`.

Behind a reverse proxy (`jupyter-server-proxy`, a Colab proxy port), the
forwarded `Host` is the proxy's, so it must be declared:
`--allowed-host hub.example.org` or `allowed_hosts=["hub.example.org"]`. One
bare host name per entry, no port, no scheme, no wildcard. Details:
{doc}`../b_concepts_mental_models/06_privacy_security_and_offline_vs_online`.

### Risky behavior: `--host 0.0.0.0`

Binding to `0.0.0.0` makes the server reachable from other machines on the network (subject to firewall rules).

This can expose:
- embeddings and cell metadata (obs fields),
- gene expression values,
- and any other exported artifacts the server can serve.

```{warning}
CORS does not “secure” a server. CORS is a browser policy; any machine that can reach your server can still download data using non-browser clients (curl, Python, etc.).
Treat the server as a dataset-serving HTTP endpoint.
```

Cellucid has **no authentication of any kind** — no password, no token, no
login, and nothing to configure that would add one. Unlike a notebook server,
there is no token in the URL: the address Cellucid prints is the whole thing a
client needs, and it is not a secret. Reads are unauthenticated whatever the
bind address; on `127.0.0.1` only your own machine can perform them, and on
`0.0.0.0` every machine that can route to the port can.

Recommended alternative for remote access: {doc}`12_remote_servers_ssh_tunneling_and_cloud` (SSH tunnel).

### On a shared machine, “the network” is every other user

On a multi-user HPC compute node, `--host 0.0.0.0` publishes the dataset to
every user who can route to that port — not only to you, and not only to the
tunnel you opened. Anyone with a shell on that node, and usually anyone with a
shell on any node of the cluster, can fetch:

- the embeddings you served,
- the obs columns you served,
- and every expression value the viewer is able to request.

They need the port number and nothing else. A single `curl` against
`http://compute-node:8765/` is enough — no browser, no Cellucid installation, no
account on your project, and no interaction with you. Batch schedulers routinely
place other users’ jobs on the same node, so a shared node is the ordinary case
rather than the unlucky one.

If the data is patient data, unpublished, or otherwise governed, that is the
decision to weigh before typing the flag — not the tunnel command afterwards.

### When the exposure is nonetheless the right trade

Bind the wildcard when the tunnel cannot terminate on the machine that holds the
data. On most clusters only the login node accepts your SSH connection, so the
forward is `ssh -N -L 8765:compute-node:8765 you@login-node`. The login node
resolves `compute-node` itself and opens an ordinary TCP connection to it; that
connection arrives from another machine, and a socket bound to `127.0.0.1` will
not accept it. No loopback-only arrangement works in that shape, so the
alternative to the wildcard bind is not a safer bind — it is no viewer at all.

This is a property of the site, not of Cellucid. The Helmholtz Munich / ICB
cluster and many others accept SSH on login nodes only, sometimes admitting a
compute-node connection just while a job of yours is running there. Whether you
ever need the wildcard is therefore decided by your cluster's policy: where the
client can reach the server directly, keep the loopback bind and expose nothing.

Note also that the wildcard changes the bind, not the address you browse to. The
viewer is still opened at `http://127.0.0.1:<port>/…` through the tunnel, so
seeing a loopback URL in the banner is not evidence that the port is private.

The full recipe (finding the node, Slurm job scripts, port collisions) is in
{doc}`17_hpc_slurm_and_compute_node_serving`; the decision itself is in
{doc}`12_remote_servers_ssh_tunneling_and_cloud`.

### How to bound a wildcard bind

- **Serve only while you are looking.** The exposure lasts exactly as long as
  the process. Start it when you sit down, `Ctrl+C` when you stand up, and do
  not leave it running inside a multi-day job.
- **Pick an unusual port.** `8765` is the documented default and the first
  number anyone would try; `--port 47231` is not secrecy, but it takes the
  server out of the path of casual scans and of a colleague who happens to open
  the default port on the node you share.
- **Prefer a loopback bind plus a direct tunnel whenever the client can reach
  the machine directly.** A lab server, a cloud VM, a container host, or a
  cluster that admits `ssh -J you@login-node … you@compute-node` all support
  `ssh -N -L 8765:127.0.0.1:8765 you@remote-host` against the default bind,
  which publishes nothing to anyone.
- **Serve less.** `--obs-key` restricts the served obs columns to the ones you
  name, so an exposed server cannot hand out the clinical columns you never
  asked it to read.
- **Add firewall rules** where the site allows them, and treat a site policy
  against non-loopback binds as authoritative over this page.

Together, on a compute node:

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset \
    --host 0.0.0.0 --no-browser --port 47231 \
    --obs-key cell_type --obs-key total_counts
```

### A wildcard bind turns the `Host` check off

The loopback rule above is enforced whenever the bound address is loopback. A
deliberate non-loopback bind — `--host 0.0.0.0`, `--host ::` — is an explicit
request for network exposure, and the legitimate names of an exposed deployment
(the node’s own name, a site alias, an IP literal) cannot be known from inside
the process, so the check is not applied there by default. Requests naming any
authority are routed.

Declaring names turns it back on, because a declaration is exactly that missing
knowledge:

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset \
    --host 0.0.0.0 --allowed-host compute-node
```

The server then answers `compute-node` on any port, plus loopback authorities on
the bound port, and refuses everything else with `421 Misdirected Request`. The
tunnelled loopback URL keeps working either way, so the declaration never breaks
the URL Cellucid printed.

### What the banner prints when a wildcard is bound

`0.0.0.0` states which interfaces accept connections; it is not an address
anything can connect *to*, so Cellucid never renders it inside a URL. The banner
reports the loopback origin — correct for a browser on the serving machine and
for the far end of any tunnel — and then names the machine itself for clients
elsewhere on the network:

A wildcard bind prints the loopback origin plus the machine's own name; see {doc}`12_remote_servers_ssh_tunneling_and_cloud`.

Read that second block as the audit line it is: those are the addresses your
neighbours on the cluster can open. When the machine has no name that resolves
outside itself, the block is absent — there is no travelling address to print,
and a name that resolves to loopback would be a URL that cannot leave the host.
An IPv6 bind is bracketed as a URL requires, so `--host ::` reports
`http://[::1]:8765` and not a bare `::1`.

## CORS behavior (what it is and why it exists)

Cellucid’s viewer UI may run on a different origin than your dataset server in some workflows.
To make browser requests work, the server sends CORS headers.

High-level behavior:
- The server only allows a limited set of origins (loopback + Cellucid’s hosted origin).
- Requests from random websites should not be able to read your dataset from a browser tab.

## Notebook-specific security notes

Notebook mode uses two channels:

1) Python → frontend commands via `postMessage` (includes a per-viewer `viewerToken`)
2) Frontend → Python events via `POST /_cellucid/events` (identified by `viewerId`)

Practical guidance:
- keep the server bound to localhost unless you have a strong reason,
- don’t embed the viewer inside untrusted web pages,
- and treat session bundles as potentially sensitive (they can encode selections, highlights, and analysis state).

## Troubleshooting security/networking issues

- Remote access (SSH tunnels): {doc}`12_remote_servers_ssh_tunneling_and_cloud`
- Viewer source/generation issues: {doc}`15_troubleshooting_viewing`
