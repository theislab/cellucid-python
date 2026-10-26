# Serve from a Cluster or over SSH

**Audience:** analysts whose objects live on a computing cluster or a remote
server
**What you'll do:** serve an object for a browser on your own computer.

```r
library(cellucid)
viewer <- cellucid_serve(seu)
#> Serving seu at http://127.0.0.1:8765
#> Open http://127.0.0.1:8765/?viewer=true&viewerId=…&viewerToken=…&anndata=true in a browser (from another machine: ssh -L 8765:127.0.0.1:8765 <this host>, then open the address there).
#> Stop it with cellucid_stop().
```

`cellucid_serve()` is `cellucid_view()` without a pane. On your computer, open
the tunnel the message names and then the address:

```sh
ssh -L 8765:127.0.0.1:8765 me@cluster.example.org
```

The server listens on the cluster node's loopback address only
(`host = "127.0.0.1"`), so nobody else on the network reaches it; the tunnel is
how you do. Keep the same port on both ends: Google sign-in (community
annotation) accepts only `localhost:8765`. When 8765 is taken on the node, the
next free port is used and named in the message; tunnel that one.

## Behind RStudio Server or Posit Workbench

These proxies forward requests under their own host name, which the server
refuses unless it is declared (the check that stops other web pages from
reading your data):

```r
viewer <- cellucid_view(seu, allowed_hosts = "rstudio.example.org")
```

## Several viewers, stopping, restarting

Each viewer has its own port; `cellucid_viewers()` lists the running ones,
`cellucid_stop(viewer)` stops one (the page shows it was stopped) and frees its
port, `cellucid_stop_all()` stops every one. A port you name with `port =` is
used or refused with what to do.

## Offline

The web app is cached after the first download and is used as it is when
www.cellucid.com cannot be reached. To copy it to a machine with no internet,
copy the folder `file.path(tools::R_user_dir("cellucid", "cache"), "web")` and pass
`web_cache_dir =` and `web_source_url = NULL`.
