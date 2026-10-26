# Google Drive, OneDrive and private sharing

**Audience:** anyone who keeps data on Google Drive or OneDrive, and labs that
share unpublished data between collaborators.
**Time:** 5 minutes.

Cellucid opens an `.h5ad` file, a Zarr ZIP (`.zarr.zip`) or any other single
file it reads straight from Google Drive or OneDrive, read in place like every
other link: only what each view needs is read. A file you share privately with
a collaborator opens for them after one sign-in. The file and the sign-in go
only between the browser and Google or Microsoft: no Cellucid server is
involved.

:::{important}
Drive and OneDrive work only on a site registered with Google and Microsoft
(an OAuth client, and for Google an API key). www.cellucid.com and a local
`cellucid serve` on port 8765 are registered. A copy of the app hosted
anywhere else answers a Drive link with "Google Drive is not configured on
this site yet", and a sign-in with "Signing in to Google / Microsoft is not
configured on this site yet"; open the link on www.cellucid.com instead.
:::

## At a glance

| | Google Drive | OneDrive and SharePoint |
|---|---|---|
| Shared with anyone who has the link | Opens with no sign-in | Opens after a Microsoft sign-in (OneDrive serves nothing anonymously) |
| Shared privately with you | Opens after **Sign in to Google** and one click on **Select** in Google's file picker | Opens after **Sign in with Microsoft**; a colleague's file in an organisation may need your administrator's approval once |
| Accounts | Any Google account | Personal Microsoft accounts, and work or school accounts |
| What Cellucid can open | One file: an `.h5ad`, a `.zarr.zip`, or another file Cellucid reads | One file: an `.h5ad`, a `.zarr.zip`, or another file Cellucid reads |
| Folders | Refused ([why](#why-folders-are-refused)) | Refused ([why](#why-folders-are-refused)) |
| How long the sign-in lasts | This tab, until it is closed; a Google sign-in needs one click to renew after about an hour | This tab, until it is closed; renewed without a click for 24 hours |

Paste the file's share link into the **Link or server** bar of the Session
panel, or open `https://www.cellucid.com/?data=<link>` ({doc}`14_open_from_link`).

## Google Drive

**Public files.** A file shared as "Anyone with the link" opens at once:
Cellucid reads it through Google's Drive API with its own public key, the only
way Google lets another website read a Drive file's bytes. Its first view
starts one round trip after the paste, and like every other link only what
each view needs is read ({doc}`14_open_from_link` lists public examples with
what each one read).

**When Google holds back downloads.** After heavy use without a Google
sign-in (many files or genes read in a few minutes), Google may hold back
every download without a sign-in from your network address for up to about
an hour. On a shared lab or university network, other people's downloads
count too. Signing in to Google always works: requests that carry a sign-in
are not held back, and read the same way (a 1.7 GB `.h5ad` opened signed in on
9.3 MB in 84 requests; a private file shared with a second account, 3.3 MB in
20). For a dataset you publish to everyone, prefer Zenodo, Figshare, Hugging
Face or GitHub, which need no sign-in and hold nothing back
({doc}`14_open_from_link`); Drive and OneDrive suit sharing with
collaborators, above all privately. While the hold lasts, Google answers each download with an
"unusual traffic" page (after about seven seconds) instead of the file,
though it still answers for the file's name and size. A browser cannot read
that page from another site, so Cellucid tells it from a failed read, and
after those seconds says:

> Google Drive sends this file only with a Google sign-in for now.
> **[Sign in to Google]** ⓘ

Signing in and choosing the file in Google's picker (as for a private file,
below) opens it with your sign-in. A read held back after the file opened
(a gene, another embedding) is told the same way, as soon as Google's answer
shows it (about nine seconds): Cellucid asks Drive for the file's details,
which answer, and for one byte of it, which does not. It then says "Google
Drive is not sending this file to this browser without a Google sign-in for
now … Open the file again to read it with your Google account, or try again
in a few minutes": opened again, the file asks you to sign in, or, when you
are signed in already, opens Google's picker on it. A read that fails
because the connection dropped is asked again as usual and says so.
Google does not publish its limit. In our tests (2026-10-25) the hold began
three times, each after a few hundred megabytes had been downloaded from one
address within minutes (about 200–400 MB in 3–10 minutes, whether in 30
requests or 400), and lasted about 50 minutes. An address released from a
hold was held again much sooner: the same day, after some 30 MB a quarter of
an hour after its hold ended. Reading in place keeps a
first view to a few megabytes, but a large R object is read whole, and the
first gene of a large `.h5ad` stored by cell reads the matrix's index once
(125 MB for a 257 MB file). A file you
granted to Cellucid in the picker is always read with your sign-in while you
are signed in, public or not.

**When you are signed in already.** Cellucid's sign-in reaches only the
files you chose for it in Google's picker. A Drive link to any other file,
public or private, opens Google's picker on that file at once, with no
sign-in step: one click on **Select** grants it, and it is read with your
sign-in, which Google never holds back. Close the picker instead and a public
file is read without your sign-in, unless Google is holding such downloads
back from your network just then.

Read with your sign-in, a file takes the same requests, and the same time,
as read without it. A request that carries a sign-in makes the browser ask
Google first whether this site may send it (a CORS preflight); each request in
flight goes to an address of its own (`…&lane=<n>`, which Google ignores),
which the browser asks about once and then remembers, instead of asking before
every request.

**Private files.** A file shared only with you shows one line under the address
bar, and the same as a notification:

> Private file on Google Drive. **[Sign in to Google]** ⓘ

```{figure} ../../../_static/screenshots/data_loading/gdrive-private-sign-in.png
:alt: The Link or server bar holding a Google Drive file link, with one line under it reading "Private file on Google Drive.", a Sign in to Google button and an (i).
:width: 496px

A private Drive link pasted into the bar: one line under it and **Sign in to
Google**. Nothing is asked of Google until you click.
```

1. Click **Sign in to Google** and choose your account in Google's window. The
   first time, Google asks you to allow Cellucid access to "the specific Google
   Drive files you use with this app".
2. Google's file picker opens, already showing the file the link names. Click
   **Select**.
3. The file opens. Next time, in this tab or after closing it, the same file
   opens without the picker: you granted it once.

Cellucid asks Google for the `drive.file` permission only: it covers the files
you choose for Cellucid in the picker and nothing else in your Drive. Google
words it as "see, edit, create and delete" because that permission has no
read-only form; Cellucid only ever reads.

When the Google sign-in expires (about an hour), the reads that need it pause
and the line says **Google sign-in expired. [Sign in again]**. One click
continues; Google's window usually closes by itself. However long Google's
window, the consent screen or the picker stays open, the file waits for you:
time spent there never counts as a stalled download.

### Google Workspace accounts: administrator approval

A Google Workspace organisation can restrict which apps its accounts may
sign in to. Google's window then says **"Access blocked"**, or Cellucid says
"Your Google Workspace administrator does not allow Cellucid for this
account", offers **Use another account** (a personal Google account works)
and links here. For administrators: in the Google Admin console, **Security
→ Access and data control → API controls → Manage Third-Party App Access →
Add app → OAuth App Name Or Client ID**, search for Cellucid's client ID
`209459713812-l2nm220ahvies0dto43b6v6korefveeq.apps.googleusercontent.com`,
and set its access to **Trusted** (or **Limited** with the `drive.file`
scope). Cellucid asks only for `drive.file` (the files a person chooses for
Cellucid) and has no server that receives data ([What Cellucid sees and
keeps](#what-cellucid-sees-and-keeps)).

```{figure} ../../../_static/screenshots/data_loading/gdrive-admin-blocked.png
:alt: The Link or server bar after Sign in to Google, for an account whose Google Workspace administrator does not allow Cellucid: the line under it says so and what to do, with a link "How an administrator allows Cellucid" and the button Use another account.
:width: 496px

What a Google Workspace account sees when its administrator does not allow
Cellucid: the reason on a row of its own, the link to the administrator steps
above, and **Use another account**, which opens Google's account chooser.
```

## OneDrive and SharePoint

OneDrive no longer serves files to anonymous visitors, not even files shared
with "Anyone with the link", so every OneDrive file opens after you sign in
with a Microsoft account:

> OneDrive file: sign in with Microsoft to open it. **[Sign in with Microsoft]** ⓘ

```{figure} ../../../_static/screenshots/data_loading/onedrive-sign-in.png
:alt: The Link or server bar holding a OneDrive share link, with one line under it reading "OneDrive file: sign in with Microsoft to open it.", a Sign in with Microsoft button and an (i).
:width: 496px

A OneDrive share link: every OneDrive file opens after a Microsoft sign-in,
and the line says so before anything is asked of Microsoft.
```

Click it and sign in in Microsoft's window. Cellucid asks for permission to read
your files and the files shared with you (`Files.Read`). The first time you open
a share link, Microsoft adds the file to your **Shared with me**, exactly as
opening the link in a browser does; afterwards the link opens by the file's own
id, and the address bar and saved sessions show `onedrive:<drive>/<item>`,
never the share link.

The sign-in is renewed without a click for 24 hours (Microsoft's limit for a
web page), then the line says **Sign in again**.

### Work and school accounts: administrator approval

With a work or school account, `Files.Read` covers your own files. A file a
colleague shared with you needs one more permission, `Files.Read.All` ("read all
files that you can access"). Cellucid asks for it only then:

> This file may belong to someone else in your organisation. … **[Ask for access]**

Microsoft's default policy lets only an administrator grant that permission, so
in many organisations Microsoft's window says that **approval is required**,
and may offer to send your administrator a request. Some organisations let
only an administrator approve any app, or any app from an unverified
publisher: their accounts see **"Need admin approval"** (`AADSTS90094`) at the
first sign-in already, with Cellucid named "unverified" (Cellucid's publisher
is not verified with Microsoft yet; once it is, organisations that let users
approve verified publishers' apps with low-risk permissions let their users
approve it themselves, which depends on each organisation's policy). Cellucid
then says:

> OneDrive file: sign in with Microsoft to open it. Your organisation
> (*its name*) lets only an administrator approve Cellucid. …
> **How an administrator approves Cellucid** **[Use another account]**

**Use another account** opens Microsoft's account chooser: a personal
Microsoft account works. "Return to the application without granting consent"
on Microsoft's page comes back to Cellucid exactly as declining does, so then
Cellucid says both: sign in again and accept, or, if Microsoft said "Need admin
approval", ask your administrator. Once an administrator approves Cellucid, it
works for everyone in the organisation.

For administrators: Cellucid is a single-page application; for OneDrive it
has no server, stores nothing and receives no data (see [What Cellucid sees and
keeps](#what-cellucid-sees-and-keeps)). It asks for the delegated Microsoft
Graph permissions `Files.Read`, `openid` and `offline_access`, and
`Files.Read.All` when a person opens a colleague's file. To approve it for the
organisation, open the Microsoft Entra admin center, **Enterprise applications
→ Cellucid → Permissions → Grant admin consent** (Cellucid's application ID is
`18aa2a4e-6431-4ef0-858d-fb6fdcab528e`; it is listed there once someone in the
organisation has tried to sign in). Not tested yet with real work or school
accounts (the one tried, at Helmholtz Zentrum München, requires an
administrator's approval, `AADSTS90094`): a colleague's file, sharing across
organisations (a guest opening another organisation's file), and the
SharePoint download host.

Checked with real personal accounts on 2026-10-25: a second personal account
opened a file shared only with it ("People you choose", view: 2.8 MB read in
30 requests) and someone else's "Anyone" link (a 5.7 MB Seurat `.rds`, read
whole in 2 requests), signed in with `Files.Read` alone: Microsoft never asked
for more.

## Share private data with a collaborator

1. Write the dataset as one file: an `.h5ad`, or a Zarr store zipped
   (`.zarr.zip`). `cellucid.prepare()` (Python) and `cellucid_prepare()` (R)
   write either; any other single file Cellucid reads (an `.h5mu`, a `.loom`,
   an `.rds`) works too.
2. Upload it to Google Drive or OneDrive.
3. Share it with your collaborator's Google or Microsoft address (Drive: **Share**
   → add their address, role Viewer; OneDrive: **Share** → specific people).
4. Send them the file's link, or a Cellucid link: open the file in Cellucid
   yourself and use **Share** (the link is `www.cellucid.com/?data=gdrive:<id>`
   or `?data=onedrive:<drive>/<item>`).
5. They open it, click **Sign in to Google** (and **Select** in the picker) or
   **Sign in with Microsoft**, and see the data.

Anyone who opens the link without access to the file sees "Your Google account
cannot open this file" or "Your Microsoft account cannot open this file", and
nothing else about it.

## What Cellucid sees and keeps

- **Where the data goes:** your browser reads the file straight from Google or
  Microsoft. Neither the file nor the sign-in passes through a Cellucid
  server.
- **What Cellucid can read:** on Google Drive, the files you chose for Cellucid
  in Google's picker; on OneDrive, the files your account can read, and only
  when you open one.
- **What is kept:** the sign-in (an access token, and for Microsoft a renewal
  token) is kept in this browser tab's session storage. It survives a reload of
  the tab and is gone when the tab closes. It is never saved to disk by
  Cellucid, never put in a link, a session file, a figure, the console or
  analytics, and never sent anywhere but to Google or Microsoft.
- **What a session file or a share link records:** the provider and the file's
  id (`gdrive:<id>`, `onedrive:<drive>/<item>`), never a token or a download
  address.

The privacy policy is at [www.cellucid.com/privacy.html](https://www.cellucid.com/privacy.html)
and the terms of use at [www.cellucid.com/terms.html](https://www.cellucid.com/terms.html),
linked as **Privacy** and **Terms** at the foot of the sidebar.

## Why folders are refused

A Cellucid-ready Zarr store is hundreds to thousands of small files. On Google
Drive each file is its own address, and a signed-in read of each costs two
round trips; Drive also does not grant a folder's files to Cellucid when you
pick the folder. OneDrive could list a folder, at a round trip per directory.
The same store zipped as one `.zarr.zip` reads in place at the speed of the
folder, with none of that cost, so a folder link is refused with what to share
instead:

> This link is a Google Drive folder. Cellucid opens single files from Google
> Drive: share the file itself (… a list of the file kinds Cellucid opens …)
> and paste its link.

## Revoke access

While you are signed in, the line under the **Link or server** bar in the
**Session** panel says so ("Signed in to Google") with **Sign out**; the line
goes when the sign-in expires (Google's after about an hour, Microsoft's
after 24 hours). Sign out
forgets the sign-in in this tab (and revokes Google's token at Google), and
any open file read with it asks to sign in again on its next read. To remove
Cellucid's permission from your account as well:

- **Google:** [myaccount.google.com/permissions](https://myaccount.google.com/permissions)
  → Cellucid → **Delete all connections**. Closing the tab already ends the
  sign-in.
- **Microsoft, personal account:** [account.live.com/consent/Manage](https://account.live.com/consent/Manage)
  → Cellucid → **Remove these permissions**.
- **Microsoft, work or school account:** [myapps.microsoft.com](https://myapps.microsoft.com)
  → Cellucid → **Manage your application** → **Revoke permissions** (or ask your
  administrator).

Cellucid does not sign you out of Google or Microsoft themselves.

## Where Drive and OneDrive work

Google and Microsoft answer only the addresses Cellucid is registered with.
Google matches an address exactly, port included, for public files (Cellucid's
API key) and for sign-ins alike. Microsoft ignores the port of a `localhost`
address, and accepts plain `http://` only for this computer.

| Cellucid opened at | Google Drive | OneDrive |
|---|---|---|
| www.cellucid.com | Yes | Yes |
| `cellucid serve` on its default port (`http://127.0.0.1:8765`, `http://localhost:8765`) | Yes | Yes |
| `cellucid serve --port <another port>` | No: Cellucid says to start it without `--port` | Yes |
| A notebook viewer (`show_anndata`), on 8765 unless another server holds it | On 8765, public files (private files are not tested yet); on another port Cellucid says to free 8765 | Only when the notebook page is on the same host as the viewer (Jupyter at `127.0.0.1`, not tested yet). From `localhost`, VS Code or any other site the browser keeps Microsoft's window from reaching the embedded viewer, and Cellucid says to open Cellucid in its own tab |
| Another computer's address (`cellucid serve --host 0.0.0.0`, `http://192.168.…`) | No | No: Cellucid says it needs `https://` or this computer's own address |
| JupyterHub, Colab, a remote server, any other site | No | No |

Where a row says no, open the file on www.cellucid.com.

Checked against Google and Microsoft on 2026-10-25 from `http://127.0.0.1:8765`:
public Drive files open with no sign-in, and Microsoft's sign-in page accepts
all three addresses Cellucid returns to (`www.cellucid.com`, `127.0.0.1` and
`localhost`, on any port). With real accounts: a Google sign-in and Google's
picker (a second Google account opened a file shared only with it: 3.3 MB in
20 requests after **Select**; a public 1.7 GB `.h5ad` opened on 9.3 MB in 84
requests after the picker), and personal Microsoft accounts on OneDrive (a
file of the owner's own, one shared only with another account, and another
person's "Anyone" link). Work and school accounts are not tested yet.

## Related pages

- {doc}`14_open_from_link`: every link Cellucid opens.
- {doc}`08_troubleshooting_data_loading`: blocked pop-ups, sign-in errors.
- {doc}`../p_developer_docs/04_configuration_env_vars_and_feature_flags`: where the tokens are stored.
