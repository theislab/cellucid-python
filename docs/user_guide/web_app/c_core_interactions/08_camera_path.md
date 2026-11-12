# Camera path

**Audience:** anyone presenting a dataset: a talk, a lab meeting, a tour for a
collaborator  
**Time:** 5 minutes  
**What you’ll learn:**
- How to save camera positions as keyframes and play them as one smooth film
- What each keyframe row and each Playback setting does
- What happens with reduced motion, in sessions and history, and what an
  agent can do with a path

A camera path is a list of saved camera positions (**keyframes**) that
Cellucid flies through in order, like the shots of a documentary. Nothing
moves until you press **Play**.

---

## Make and play a path

1) Open the **Camera path** section of the sidebar (it starts collapsed).
2) Frame the first view you want (drag, zoom, or switch dimension) and click
   **Save position**. The keyframe appears in the list as **KF 1**.
3) Move to the next view and click **Save position** again. Repeat for every
   shot.
4) With two keyframes or more, a playback bar appears over the top edge of
   the plot: **Stop and return to the start**, **Play** and a timeline you can
   drag. Press **Play**.

The playback bar hides itself after ten idle seconds; move the pointer over
the plot, or press {kbd}`Tab`, to bring it back. It stays visible while the
path plays.

```{note}
There is no video export. To keep the film as a video, record the screen
with your operating system's screen recorder while the path plays. For a
still picture of any frame, drag the timeline to it and use
{doc}`../k_figure_export/index`.
```

---

## The Keyframes list

Each keyframe is one row:

- **Its label** (**KF 1**, **KF 2**, …): click it to rename the keyframe, to
  name the shot ("Beta cells", "Overview"). Labels are unique.
- **Move up** / **Move down** change its place in the film, **Go to position**
  moves the camera there, **Delete keyframe** removes it.
- Under it, its **hold**: the seconds the camera rests there before moving on
  (empty for none). The last keyframe's hold is the rest at the end.
- Then **→** and the time to the next keyframe, in seconds. Leave it empty
  (**auto**) and **Default speed** paces it by the distance the camera
  travels.

Below the list:

- **End with a return to the first keyframe** (two keyframes or more) adds a
  last shot, **↩ Return to start**, that flies back to the first keyframe:
  the path ends where it began, ready to loop.
- **Set all** … **s** with **Apply** gives every transition the same time
  (leave it empty and apply to put every one back to auto).
- **Clear all** deletes every keyframe.

---

## Playback settings

| Setting | What it does | Starts at |
|---|---|---|
| **Default speed** | How fast the camera travels on a transition left to auto. Keyframes with a set time keep it. | Slow |
| **Position interpolation** | The camera's route: **Catmull-Rom spline** (a smooth curve through every keyframe), **Linear** (straight segments) or **Bezier** | Catmull-Rom spline |
| **Rotation interpolation** | How the camera turns: **SLERP (smooth)** or **Linear (Euler)** | SLERP (smooth) |
| **Easing** | How the whole path starts and ends: **Linear**, **Ease in-out**, **Ease in**, **Ease out** | Linear |
| **Loop playback** | Play the path again from the start when it ends | Off |
| **Autoplay** | Start a ready path at once, and start it when a session saved with Autoplay on is loaded | Off |
| **Captions** | While the path plays, show each keyframe you named as a caption over the plot | On |

### Captions

A keyframe's label is also its caption, a documentary's title: while the
path plays, the label of the keyframe the film is at appears low over the
plot, from halfway through the move into that keyframe to halfway through
the move out of it, so it is up while the camera approaches, rests on and
leaves its subject. It fades in and out; when your system asks for reduced
motion it is shown and hidden at once (and a path then goes to its end at
once, so no caption is up). A keyframe still called by its automatic name
(**KF 3**) and the return to start have no caption: name a keyframe ("Beta
cells", "Overview") to caption it. Pausing keeps the caption up; stopping
takes it away. **Captions** is saved with the path in a session. A picture
of the screen (an agent's `view_picture` with source "screen") shows the
caption as drawn, and a figure, which leaves captions out, says which one
was up.

---

## How a path plays

- **Each keyframe is reached at its time.** A transition set to 5 s takes
  5 s, whatever the easing; an empty duration is paced by **Default speed**
  from the distance the camera travels. The marks on the playback bar are
  where the camera is at each keyframe, and a keyframe that holds has a
  wider band after its mark: the time the camera rests there.
- **Easing is how the path starts and ends**: **Ease in-out** starts from
  rest and comes to rest at the last keyframe, **Ease in** and **Ease out**
  one of the two, **Linear** moves at once and stops at once. It never
  stretches one transition at another's expense. (Before 0.9.8 easing
  warped the whole timeline: in an eased path the middle transitions played
  at about half their time and the first and last at about twice.)
- **Catmull-Rom** (the default) is the smooth one: the camera's speed does
  not jolt at a keyframe, even when the transitions on either side have
  different times, and no coordinate overshoots the keyframes on either
  side. The camera never comes closer to the cells than the keyframes it
  moves between, and never tilts past a keyframe over the top (where the
  view would flip).
- **A seamless loop**: with **Loop playback**, **End with a return to the
  first keyframe** and **Linear** easing, the last transition flows into the
  first with no change of speed or direction at the seam (a loop behind a
  talk). With an eased loop the camera rests at the seam each time round.
- **A hold**: the row under each keyframe starts with its hold, the seconds
  the camera rests there before moving on (empty for none; the last
  keyframe's holds the end). The camera eases into it, rests, and eases out:
  the documentary pause on a subject. A hold adds to the path's length and
  leaves every transition its own time; keyframes are saved with their holds
  in sessions.
- **A close-up shows every cell in it.** A keyframe's camera orbits a point
  (the subject an agent's **Frame cells with the camera** centers). Cells
  between the camera and that point are drawn whole and opaque, larger the
  nearer they are (up to 160 pixels), so a shot through another population
  shows that population in front: choose a direction with a clear view,
  frame less tightly, or lower **Point size** or **Perspective size
  scaling** ({doc}`03_render_modes_points_vs_volumetric_smoke`).
- The playback bar's **Pause** holds the camera where the film is (**Play**
  goes on from there); its **Stop** ends the film and returns the camera to
  the first keyframe, as its tooltip says. Every other stop leaves the camera
  where the film stopped, and the bar shows that place: a film played to its
  end, an agent's **Stop the camera path**, **Stop** on the agent's turn. One
  rule: **Play** goes on from where the bar shows, and from the end it starts
  again. Changing the path (a keyframe, a timing) puts the bar back at the
  start of the new path. An agent's **Play the camera path** always starts
  from the first keyframe.
- An agent's **Stop the camera path** stops a playing path where the camera
  is, at once, also while the agent's **Play the camera path** still waits
  for it: that play then answers that it was stopped. **Stop** on the
  agent's turn does the same: the camera stays where the film stopped, and
  the play's one step holds that camera (Undo puts the camera back where it
  was before the play). A path the agent set looping belongs to its turn
  too: **Stop** on the turn stops the loop where it is; a loop the turn left
  playing when it ended plays on until you pause or stop it.
- A path made by an agent in the chat ({doc}`../g_agents/01_chat_bar`) is
  the same path: its keyframes are listed here, and Play, the timings and
  every setting work on it as on yours. On a narrow window, where the chat
  stands over the foot of the plot, it folds to its button while the path
  plays and comes back when it stops.

---

## Reduced motion

When your system asks for reduced motion (macOS *Reduce motion*, Windows
*Animation effects* off), **Play** does not fly: the camera goes straight to
the last keyframe, as the film would end. The timeline still reaches every
keyframe and every moment between them: drag it to the shot you want. The
same holds for a path an agent plays.

---

## Sessions, history and datasets

- Every change to the path (a keyframe saved, renamed, moved or deleted, a
  time, a hold, a setting) is a step of the history: **Undo** takes it back
  ({doc}`../m_state_history/index`).
- A session file keeps the keyframes, their times and holds and every
  Playback setting ({doc}`../l_sessions_sharing/index`). A path saved with
  **Autoplay** on starts playing once the session has loaded; otherwise it
  waits for **Play**.
- A path belongs to its dataset: opening another dataset clears it, so a
  path made on one dataset never moves the camera over another.

---

## Camera paths made by an agent

An agent in the chat ({doc}`../g_agents/01_chat_bar`) can build and play a
path for you, as a director would: frame groups of cells, save them as named
keyframes with holds, set the timing, turn on the return to the start, set
the Playback options, play the path, stop it, and read it back. Its
keyframes are listed here like yours, and every row and setting works on
them. Each of its changes is a step you can undo.

---

## Troubleshooting

### “I don’t see the playback bar”

- The path needs two keyframes or more.
- The bar hides after ten idle seconds: move the pointer over the plot or
  press {kbd}`Tab`.

### “Play jumps straight to the end”

Reduced motion is on in your system settings (above).

### “The path started on its own”

**Autoplay** is on, in this session or in the session you loaded. Turn it
off under **Playback**.

### “My keyframes are gone”

Another dataset was opened (a path belongs to its dataset), or the steps were
undone. Load the session you saved, or redo.

---

## Interface reference

```{figure} ../../../_static/screenshots/web_app/camera-path-configured-panel.png
:alt: The Camera path panel after two positions were saved: keyframes KF 1 and KF 2 with the KF 1 → KF 2 duration on auto, End with a return to the first keyframe, Set all with Apply, and Playback with Default speed, Position interpolation Catmull-Rom spline, Rotation interpolation SLERP (smooth), Easing Linear, Loop playback and Autoplay.
:width: 532px

Two valid keyframes expose camera-path timing and interpolation settings.
Autoplay remains off until the user explicitly enables it.
```
