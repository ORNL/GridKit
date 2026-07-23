// Entry point for the latkit browser bundle used by the interactive oneline
// diagrams (see docs/_ext/oneline.py and docs/_static/js/oneline.js).
// `npm run bundle` compiles this into docs/_static/js/latkit.bundle.js,
// which is generated at docs build time and never committed.
export { createNetwork } from "@latkit/network";
export { requestDevice, GpuUnavailableError } from "@latkit/gpu";
export { colormap } from "@latkit/colormaps";
