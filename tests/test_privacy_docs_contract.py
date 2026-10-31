"""The documentation says what the code does about privacy, in every page.

Each sentence below was true once and is not any more: what is read is also
kept in this browser's storage (not "memory only"); sign-ins are short-lived
tokens that survive a reload and are revoked on Sign out (not "session-only"
or "gone when the tab closes"); a link's access key goes nowhere but the
address bar; community annotation signs in with GitHub on www.cellucid.com
(not Google on port 8765); a server on ``127.0.0.1`` is readable by every
user of that machine; Python's web cache is per user; the session format is
3; every filter keeps applying. A page that says one of them again fails here.
"""

from __future__ import annotations

import re
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DOCS = REPOSITORY_ROOT / "docs" / "user_guide"

STALE = {
    r"stays in this tab's memory": "what is read is also kept in this browser's storage",
    r"nothing Cellucid reads is stored": "what is read is also kept in this browser's storage",
    r"memory only: nothing is saved to disk": "what is read is also kept in this browser's storage",
    r"closing the tab signs you out \(token is gone\)": "sign-ins survive a reload and are revoked on Sign out",
    r"tokens are session-only": "GitHub tokens expire after 8 hours",
    r"token \(session-only\)": "GitHub tokens expire after 8 hours",
    r"Even though tokens are session-only": "GitHub tokens expire after 8 hours",
    r"including any access key it carries": "a link's access key goes nowhere but the address bar",
    r"Google sign-in \(community annotation\)": "community annotation signs in with GitHub on www.cellucid.com",
    r"\(Google Drive data, community annotation\) works only": "community annotation signs in with GitHub on www.cellucid.com",
    r"This means \*\*only your machine\*\* can access it": "every user of a shared machine can read it",
    r"process-independent temporary directory": "the web cache is per user",
    r"CELLUCID_SESSION/2\\n": "the session format is 3",
    r"session file of format 2 …": "the session format is 3",
    r"currently 2 — and restores": "the session format is 3",
    r"active-field scoped": "every filter keeps applying until removed",
    r"\(\*Representation\* ": "the Details row is Categories",
    r"Not covered by this review: the app is www.cellucid.com": "a notebook widget is a viewer page and sends no analytics",
    r"which is the case in this release\) says \"not configured": "the hosted app is registered with Google and Microsoft",
    r"Outside the analytics beacon, nothing is sent that you did not initiate": "the hosted app reads its sample catalog at startup",
    r"The server only allows a limited set of origins \(loopback \+ Cellucid’s hosted origin\)": "CORS names the hosted app's origins exactly",
    r"https://\*\.cellucid\.com": "CORS names the hosted app's origins exactly",
}


def test_no_documentation_page_repeats_a_privacy_claim_the_code_contradicts() -> None:
    found: list[str] = []
    for path in sorted(DOCS.rglob("*.md")):
        if "_build" in path.parts:
            continue
        text = path.read_text(encoding="utf-8")
        for pattern, truth in STALE.items():
            for match in re.finditer(pattern, text):
                line = text.count("\n", 0, match.start()) + 1
                found.append(f"{path.relative_to(REPOSITORY_ROOT)}:{line}: {match.group(0)!r} (now: {truth})")
    assert not found, "\n".join(found)


def test_the_privacy_model_states_the_current_storage_sign_in_and_link_rules() -> None:
    page = (DOCS / "web_app" / "o_accessibility_privacy_security" / "02_privacy_model.md").read_text(encoding="utf-8")
    for required in (
        "(link-access-keys)=",
        "(shared-machines)=",
        "expires after 8 hours",
        "is held in the page's memory",
        "oauth2.googleapis.com/revoke",
        "a viewer page never loads analytics",
        "user-provided data collection",
        "`cellucid:opening-cost:v2`",
        "`cellucid:analysis-cost:v2`",
        "`cellucid_antialias`",
    ):
        assert required in page, required
