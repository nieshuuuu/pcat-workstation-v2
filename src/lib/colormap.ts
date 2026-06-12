/**
 * Shared colormap + CT-window constants for the material-decomposition views.
 *
 * Single source of truth so the 2D cross-section overlay (SnakeEditor), the
 * whole-volume Water/Lipid map (WaterLipidView) and their colorbars all render
 * the same fraction/density as the same color, and the grayscale CT underlay the
 * same way. Previously jet() and the CT window were copy-pasted into each view.
 */

/** Window for the grayscale CT underlay, in Hounsfield units. */
export const CT_LO = -160;
export const CT_HI = 240;

/** Classic jet (matplotlib): t∈[0,1] → blue→cyan→green→yellow→red. Returns
 *  8-bit RGB. Input is clamped to [0, 1]. */
export function jet(t: number): [number, number, number] {
  const u = Math.max(0, Math.min(1, t));
  const r = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * u - 3)));
  const g = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * u - 2)));
  const b = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * u - 1)));
  return [Math.round(r * 255), Math.round(g * 255), Math.round(b * 255)];
}
