# Google Drive, OneDrive and private sharing

**Audience:** anyone who keeps data on Google Drive or OneDrive, and labs that
share unpublished data between collaborators.
**Time:** 5 minutes.

Cellucid opens an `.h5ad` file or a Zarr ZIP (`.zarr.zip`) straight from Google
Drive or OneDrive, read in place like every other link: only what each view
needs is read. A file you share privately with a collaborator opens for them
after one sign-in. The file and the sign-in go only between the browser and
Google or Microsoft: no Cellucid server is involved.

## At a glance

| | Google Drive | OneDrive and SharePoint |
|---|---|---|
| Shared with anyone who has the link | Opens with no sign-in | Opens after a Microsoft sign-in (OneDrive serves nothing anonymously) |
| Shared privately with you | Opens after **Sign in to Google** and one click on **Select** in Google's file picker | Opens after **Sign in with Microsoft**; a colleague's file in an organisation may need your administrator's approval once |
| Accounts | Any Google account | Personal Microsoft accounts, and work or school accounts |
| What Cellucid can open | One `.h5ad` or `.zarr.zip` file | One `.h5ad` or `.zarr.zip` file |
| Folders | Refused ([why](#why-folders-are-refused)) | Refused ([why](#why-folders-are-refused)) |
| How long the sign-in lasts | This tab, until it is closed; a Google sign-in needs one click to renew after about an hour | This tab, until it is closed; renewed without a click for 24 hours |

Paste the file's share link into the **Link or server** bar of the Session
panel, or open `https://www.cellucid.com/?data=<link>` ({doc}`14_open_from_link`).

## Google Drive

**Public files.** A file shared as "Anyone with the link" opens at once:
Cellucid reads it through Google's Drive API with its own public key, the only
way Google lets another website read a Drive file's bytes. Its first view
starts one round trip after the paste.

**Private files.** A file shared only with you shows one line under the address
bar, and the same as a notification:

> Private file on Google Drive. **[Sign in to Google]** ⓘ

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
continues; Google's window usually closes by itself.

## OneDrive and SharePoint

OneDrive no longer serves files to anonymous visitors, not even files shared
with "Anyone with the link", so every OneDrive file opens after you sign in
with a Microsoft account:

> OneDrive file: sign in with Microsoft to open it. **[Sign in with Microsoft]** ⓘ

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
and may offer to send your administrator a request. Once an administrator
approves Cellucid, it works for everyone in the organisation.

For administrators: Cellucid is a single-page application; for OneDrive it
has no server, stores nothing and receives no data (see [What Cellucid sees and
keeps](#what-cellucid-sees-and-keeps)). To approve it for the organisation, open
the Microsoft Entra admin center, **Enterprise applications → Cellucid →
Permissions → Grant admin consent**. Sharing across organisations (a guest
opening another organisation's file) is not tested yet.

## Share private data with a collaborator

1. Write the dataset as one file: an `.h5ad`, or a Zarr store zipped
   (`.zarr.zip`). `cellucid.prepare()` (Python) and `cellucid_prepare()` (R)
   write either.
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

The privacy policy is at [www.cellucid.com/privacy.html](https://www.cellucid.com/privacy.html).

## Why folders are refused

A Cellucid-ready Zarr store is hundreds to thousands of small files. On Google
Drive each file is its own address, and a signed-in read of each costs two
round trips; Drive also does not grant a folder's files to Cellucid when you
pick the folder. OneDrive could list a folder, at a round trip per directory.
The same store zipped as one `.zarr.zip` reads in place at the speed of the
folder, with none of that cost, so a folder link is refused with what to share
instead:

> This link is a Google Drive folder. Cellucid opens single files from Google
> Drive: share the .h5ad or .zarr.zip file itself and paste its link.

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

## Signing in from a local server

On `cellucid serve`, Microsoft sign-in works on any port. Google matches the
page's address exactly, so Google sign-in works on www.cellucid.com and on
`cellucid serve`'s default port only (`http://127.0.0.1:8765` or
`http://localhost:8765`); on another port, Cellucid says so. Public Drive
files open on any port.

## Related pages

- {doc}`14_open_from_link`: every link Cellucid opens.
- {doc}`08_troubleshooting_data_loading`: blocked pop-ups, sign-in errors.
- {doc}`../p_developer_docs/04_configuration_env_vars_and_feature_flags`: where the tokens are stored.
