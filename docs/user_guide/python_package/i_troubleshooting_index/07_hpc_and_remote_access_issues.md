# HPC and remote access issues

Symptom-first triage for serving a dataset that lives on a cluster or another
machine. The full recipe is
{doc}`../d_viewing_apis/17_hpc_slurm_and_compute_node_serving`; this page is for
when something has already gone wrong.

## SSH will not reach the machine with the data

### `ssh: connect to host <node> port 22: Operation timed out`

Nothing answered. The node's address is almost certainly private to the cluster
network and not routable from where you are. No SSH option fixes a missing
route.

**Do this instead:** forward to the node *from* the login node, which is inside
that network:

```bash
ssh -N -L 8765:compute-node:8765 you@login-node
```

and serve with `--host 0.0.0.0` so the connection the login node opens is
accepted.

### `ssh: connect to host <node> port 22: Connection refused`

Something answered and rejected the connection: packets reach the machine, but
nothing is listening on that port. The route exists. If the host normally runs
sshd on a non-standard port, use `-p`; on a cluster, prefer the login-node
forward above.

### `Permission denied (publickey,keyboard-interactive)` on a compute node

The node answered and refused your credentials. Common and expected: many
clusters accept SSH on login nodes only, or admit a compute-node connection only
while you hold a running job on that node.

`ssh -J login-node compute-node` does **not** help. `ProxyJump` still opens a
real SSH session to the compute node and uses the login node only as a relay, so
the same refusal applies.

Use the login-node port forward, which needs no SSH session to the compute node
at all.

### The password is not accepted

Check which host is asking. Run the command with `-v` and read the
`Authenticating to <host>` line: a two-hop command can prompt for the far host
while you are typing the near host's credential. Cluster pools also often use
different credentials from the general institute account.

## The tunnel connects but the viewer does not load

### Browser says the site cannot be reached

The forward exists but nothing is serving behind it. Check, in order:

1. Is the server still running? An interactive job's server dies with the SSH
   session that started it -- run it inside `tmux`.
2. Is it on the port you forwarded? `--port 0` picks an arbitrary free port.
3. Is it on the node you named? Re-check with `squeue -u "$USER"`; a new
   allocation usually means a new node.

### `channel 2: open failed: connect failed: Connection refused`

The tunnel works; the far end refused. Either nothing is listening on that port
on that node, or the server is bound to `127.0.0.1` while your forward names the
node by name. A loopback socket accepts connections from its own machine only,
and the login node is a different machine. Add `--host 0.0.0.0`.

### I bound `0.0.0.0` — so which address do I open?

`http://127.0.0.1:<port>/…`, unchanged. This is the single most common confusion
on a cluster, and the two settings answer different questions:

| | Decides | Value on a compute node |
| --- | --- | --- |
| `--host` | which interfaces the **server accepts on** | `0.0.0.0`, so the login node's connection is accepted |
| Browser address | where the **tunnel ends** | `127.0.0.1`, because that end is on your own machine |

Never type `0.0.0.0` into a browser: it means "every interface" and resolves
nowhere. Do not substitute the compute node's name either, unless you are already
inside the cluster network and not using a tunnel at all.

### The banner shows `http://0.0.0.0:8765` in an older transcript

That was a defect, now fixed: the wildcard was printed where a URL belonged.
Current versions print the loopback origin plus the machine's own name. Through a
tunnel, always use the loopback Viewer URL including its query —
`http://127.0.0.1:8765/?anndata=true` for direct AnnData, `/?source=remote` for a
prepared export.

### `421 Misdirected Request`

The `Host` header does not name this server. This is the reverse-proxy case, not
the tunnel case: declare the proxy's host with `--allowed-host hub.example.org`.
A tunnel sends `Host: 127.0.0.1:<port>`, which a loopback bind already accepts.

## The server itself

### `Port <n> is already in use`

Compute nodes are shared, and popular ports collide. Choose another
(`--port 12765`), or use `--port 0` and read the port from the banner -- but
then start the tunnel *after* you know the number.

### Startup prints nothing for minutes

It is working, not hung. A large object spends that time inside step 2 of five,
building the adapter. The steps report as they go:

```text
[2/5] Loading AnnData...
      Mode: read-only backed h5ad
      Obs columns: classifying 34
      Embeddings: resolving obsm keys
      Embeddings: 2D from obsm['X_umap']
```

If you asked for connectivity on a large graph, the `Connectivity: reading
obsp['connectivities']` line is followed by the longest wait of the run, and a
graph above 50,000,000 stored neighbors logs a warning saying so.

### Startup aborts while establishing the web build

Every viewer-serving start verifies the current web build from its configured
source, and a node with no outbound network access cannot. Serve the data
endpoints alone with `--no-web-ui`, or run on a node with access. See
{doc}`../b_concepts_mental_models/06_privacy_security_and_offline_vs_online`.

### Cellucid says a package is missing that is definitely installed

Read the whole line; three different conditions produce three different
messages:

- *"is not installed in this environment"* -- the distribution is genuinely
  absent from the interpreter running Cellucid. Check you are in the environment
  you think you are: `which cellucid` and `cellucid --version`. A shell that
  cached an older path is a common cause; `hash -r` clears it.
- *"could not read a name it needs from the installed …"* -- the package is
  present but too old to hold the name Cellucid asked of it. Upgrade it.
- *"Cellucid needs the standard-library module …"* -- the module belongs to the
  interpreter, and no package supplies it. The message names the running Python
  version, its executable, and the range Cellucid supports. This is what an
  out-of-range interpreter looks like: run Cellucid on a supported Python.

## The dataset opens but something is missing

### The neighbor-graph overlay is not there

Connectivity is served only when asked for. The dataset report says which state
you are in:

```text
      Connectivity: not served (obsp['connectivities'] is present; pass serve_connectivity=True, or --connectivity, to draw it)
```

Add `--connectivity` on the command line, or `serve_connectivity=True` in
Python, and expect startup to take longer: the whole graph is read and validated
before the server binds.

### `Connectivity was asked for, and adata.obsp has no 'connectivities' matrix`

You asked for a graph this object does not carry. Compute one with
`sc.pp.neighbors(adata)` before serving, or serve without asking for it.

### The UMAP is not found, but `obsm['X_umap']` exists

A bare `X_umap` of 1, 2, or 3 columns is read at the dimension its column count
states, so this should work. If it is refused, the message names the shape it
found: an array of some other width -- a latent space named after a plot -- names
a dimension no viewer draws. Assign the columns you mean to draw under the key
for their own count.

## Related pages

- {doc}`../d_viewing_apis/17_hpc_slurm_and_compute_node_serving`
- {doc}`../d_viewing_apis/12_remote_servers_ssh_tunneling_and_cloud`
- {doc}`04_server_mode_issues`
- {doc}`06_performance_and_memory_issues`
