/**
 * Cross-component UI state (Svelte 5 runes).
 *
 * `analysisMaximized` lives here, not inside ContextPanel, so the root App
 * keydown handler can see it. Escape must close the fullscreen analysis overlay
 * WITHOUT also firing App's "Escape clears the active vessel" shortcut — both
 * used to run off the same `window` keydown, and `stopPropagation()` does not
 * stop sibling listeners on the same target, so closing the overlay wiped the
 * FAI result. A single owner of the flag keeps that decision in one place.
 */

let analysisMaximized = $state(false);

export const uiStore = {
  get analysisMaximized(): boolean {
    return analysisMaximized;
  },
  set analysisMaximized(v: boolean) {
    analysisMaximized = v;
  },
};
