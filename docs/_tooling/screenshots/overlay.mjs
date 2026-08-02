// The pointer and annotation layer a documentation screenshot needs.
//
// A browser screenshot never contains the operating system's cursor, so an
// image of "click this button" shows a button and nothing else. The pointer
// here is drawn into the page immediately before the capture, anchored to the
// bounding rectangle of the element being described, so it lands on the control
// the prose names and moves with it when the layout changes.
//
// Two deliberate choices:
//
//   * The pointer is drawn *in addition to* a real `page.hover()`. The hover is
//     what gives the control its true `:hover` styling; the drawing is what
//     tells the reader where the pointer is. Drawing without hovering produces
//     a picture of an arrow floating over an idle button, which is worse than
//     no arrow at all.
//   * The arrow is a filled white glyph with a dark outline and a soft shadow.
//     A flat black arrow disappears over dark UI and a flat white one
//     disappears over light UI; the outlined form reads on both, and the shadow
//     makes it legible as a deliberate annotation rather than a rendering
//     artefact.
//
// Everything is `position: fixed`, `pointer-events: none`, and painted above
// the application, so nothing here can reflow the page or change what the
// application does. Every function in this module is serialised into the page
// by Playwright and therefore closes over nothing.

/** Element id of the injected layer, used to remove it again. */
export const OVERLAY_ROOT_ID = 'cellucid-docs-capture-overlay';

/**
 * Draw the pointer and annotations. Runs inside the page.
 *
 * @param {{
 *   rootId: string,
 *   cursor: {x: number, y: number, state: 'hover'|'press'}|null,
 *   rings: Array<{x: number, y: number, width: number, height: number,
 *                 label: string|null}>,
 * }} plan
 */
export function paintOverlay(plan) {
  const previous = document.getElementById(plan.rootId);
  if (previous !== null) previous.remove();

  const root = document.createElement('div');
  root.id = plan.rootId;
  root.setAttribute('aria-hidden', 'true');
  root.style.cssText = [
    'position:fixed',
    'inset:0',
    'pointer-events:none',
    'z-index:2147483647',
    'contain:strict',
  ].join(';');

  const ACCENT = '#f2711c';

  for (const ring of plan.rings) {
    const box = document.createElement('div');
    const inset = 4;
    box.style.cssText = [
      'position:absolute',
      `left:${ring.x - inset}px`,
      `top:${ring.y - inset}px`,
      `width:${ring.width + inset * 2}px`,
      `height:${ring.height + inset * 2}px`,
      'border-radius:8px',
      `border:2px solid ${ACCENT}`,
      `box-shadow:0 0 0 4px rgba(242,113,28,0.22)`,
      'box-sizing:border-box',
    ].join(';');
    root.appendChild(box);

    if (ring.label !== null) {
      const badge = document.createElement('div');
      badge.textContent = ring.label;
      badge.style.cssText = [
        'position:absolute',
        `left:${ring.x - inset - 11}px`,
        `top:${ring.y - inset - 11}px`,
        'width:22px',
        'height:22px',
        'border-radius:11px',
        `background:${ACCENT}`,
        'color:#ffffff',
        'font:700 13px/22px ui-sans-serif,-apple-system,Helvetica,Arial,sans-serif',
        'text-align:center',
        'box-shadow:0 1px 3px rgba(0,0,0,0.45)',
      ].join(';');
      root.appendChild(badge);
    }
  }

  if (plan.cursor !== null) {
    // viewBox units are the classic arrow outline; the hotspot is (1,1), so the
    // glyph is offset by the same fraction of its rendered size to put the tip
    // exactly on the requested point.
    const WIDTH = 22;
    const HEIGHT = 32;
    const HOTSPOT_X = (1 / 17) * WIDTH;
    const HOTSPOT_Y = (1 / 25) * HEIGHT;

    const pressed = plan.cursor.state === 'press';
    const holder = document.createElement('div');
    holder.style.cssText = [
      'position:absolute',
      `left:${plan.cursor.x - HOTSPOT_X}px`,
      `top:${plan.cursor.y - HOTSPOT_Y}px`,
      `width:${WIDTH}px`,
      `height:${HEIGHT}px`,
      'filter:drop-shadow(0 1px 2px rgba(0,0,0,0.55))',
    ].join(';');
    holder.innerHTML =
      `<svg width="${WIDTH}" height="${HEIGHT}" viewBox="0 0 17 25" ` +
      'xmlns="http://www.w3.org/2000/svg">' +
      '<path d="M1 1 L1 22.4 L6.2 17.5 L9.5 24.8 L12.9 23.2 L9.6 16.1 L16.4 16.1 Z" ' +
      `fill="${pressed ? '#e6e6e6' : '#ffffff'}" stroke="#141414" ` +
      'stroke-width="1.3" stroke-linejoin="round"/></svg>';
    root.appendChild(holder);

    if (pressed) {
      // A click is an event, not a state, so it needs a mark of its own. Two
      // concentric rings at the hotspot read as "activated here" without
      // implying an animation the still cannot show.
      const pulse = document.createElement('div');
      const size = 30;
      pulse.style.cssText = [
        'position:absolute',
        `left:${plan.cursor.x - size / 2}px`,
        `top:${plan.cursor.y - size / 2}px`,
        `width:${size}px`,
        `height:${size}px`,
        `border-radius:${size / 2}px`,
        'border:2px solid rgba(242,113,28,0.95)',
        'box-shadow:0 0 0 5px rgba(242,113,28,0.20)',
        'box-sizing:border-box',
      ].join(';');
      root.insertBefore(pulse, holder);
    }
  }

  document.body.appendChild(root);
}

/**
 * Remove the overlay. Runs inside the page.
 *
 * @param {string} rootId
 */
export function clearOverlay(rootId) {
  document.getElementById(rootId)?.remove();
}
