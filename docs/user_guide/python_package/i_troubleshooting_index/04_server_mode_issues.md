# Server mode issues

## Address already in use

Another process owns the requested port. Identify that process or choose a
different explicit port; do not terminate an unknown process.

## The server starts but another machine cannot connect

`127.0.0.1` is reachable only from the same machine. Remote access requires an
intentional bind address plus an appropriate firewall, SSH tunnel, or reverse
proxy. Do not expose sensitive datasets to an untrusted network.

## The browser reports an HTTP error

Record the exact URL, method, status, and response body. Verify that you supplied
the export root or supported AnnData input expected by the command, then test the
same URL from the browser host.

## The process does not stop cleanly

Keep the returned server/viewer handle and use its documented lifecycle method.
In a CLI session, use the terminal interrupt described by the command banner.

See {doc}`../d_viewing_apis/04_cli_cellucid_serve_quickstart`,
{doc}`../d_viewing_apis/09_server_mode_advanced`,
{doc}`../d_viewing_apis/12_remote_servers_ssh_tunneling_and_cloud`, and
{doc}`../h_developer_docs/09_server_mode_architecture_endpoints_and_security`.
