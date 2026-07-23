// Renders interactive oneline diagrams inside the `.oneline` containers
// produced by the `oneline` Sphinx directive (docs/_ext/oneline.py).
//
// Progressive enhancement: diagrams initialize lazily when scrolled into
// view, all diagrams on a page share one GPU device, and on any failure
// (no WebGPU, compatibility-mode device, lost device) the container hides
// itself so the page reads exactly as it does today.

const hosts = document.querySelectorAll(".oneline");

// Containers are display:none in CSS; only reveal them when WebGPU can at
// least be attempted. IntersectionObserver never fires on hidden elements,
// so activation must happen before observing.
if (hosts.length > 0 && globalThis.navigator?.gpu) {
  for (const host of hosts) host.classList.add("oneline--active");
  // Promise of [latkit module, GPUDevice], shared by every diagram on the
  // page. A rejection is cached too: later diagrams fail fast.
  let shared = null;

  const acquire = () => {
    shared ??= import(new URL("./latkit.bundle.js", import.meta.url).href).then(
      async (latkit) => [latkit, await latkit.requestDevice()],
    );
    return shared;
  };

  const states = new WeakMap();
  const observer = new IntersectionObserver(onIntersect, { rootMargin: "200px" });
  for (const host of hosts) observer.observe(host);

  function onIntersect(entries) {
    for (const entry of entries) {
      const state = states.get(entry.target);
      if (!state) {
        if (!entry.isIntersecting) continue;
        const created = { net: null, attempts: 0 };
        states.set(entry.target, created);
        mount(entry.target, created).catch(() => unavailable(entry.target));
      } else if (state.net) {
        // Keep rendering only while on screen; latkit already handles
        // document-level visibility itself.
        if (entry.isIntersecting) state.net.resume();
        else state.net.pause();
      }
    }
  }

  function unavailable(host) {
    host.classList.remove("oneline--active");
  }

  async function mount(host, state) {
    const data = JSON.parse(
      host.querySelector('script[type="application/json"]').textContent,
    );
    const acquired = acquire();
    const [latkit, device] = await acquired;

    const canvas = document.createElement("canvas");
    canvas.setAttribute("aria-hidden", "true");
    host.append(canvas);

    let palette;
    try {
      palette = latkit.colormap(host.dataset.colormap || "viridis");
      // Unknown names return a closure that only throws once called;
      // evaluate here so the fallback can catch it.
      palette(0);
    } catch {
      palette = latkit.colormap("viridis");
    }

    const net = await latkit.createNetwork(device, canvas, {
      colormap: palette,
      daylight: false,
      earthAxis: false,
    });
    state.net = net;

    const edgeCount = data.edges.length / 2;
    net.load({
      vertexCount: data.vertexCount,
      vertexCoords: data.coords ? new Float32Array(data.coords) : undefined,
      edges: new Uint32Array(data.edges),
      polylineStart: new Uint32Array(edgeCount + 1),
    });
    if (data.vm) {
      // Per-unit voltage magnitude from the case's initial condition.
      net.setChannel("vertexColor", new Float32Array(data.vm), [0.9, 1.1]);
    }
    net.fit(true);
    net.fadeIn();

    const caption = document.createElement("div");
    caption.className = "oneline__caption";
    const idle =
      `${data.name || "network"} — ${data.vertexCount} buses, ` +
      `${edgeCount} branches · drag to pan, scroll to zoom`;
    caption.textContent = idle;
    host.append(caption);

    net.on("hover", (kind, index) => {
      caption.textContent = describe(data, kind, index) ?? idle;
    });

    net.on("deviceLost", () => {
      net.destroy();
      canvas.remove();
      caption.remove();
      // Every diagram shares one device, so one loss fires every handler;
      // only the first should discard the cached promise or later handlers
      // would clobber the replacement and split diagrams across devices.
      if (shared === acquired) shared = null;
      state.net = null;
      if (state.attempts++ < 2) {
        mount(host, state).catch(() => unavailable(host));
      } else {
        unavailable(host);
      }
    });
  }

  function describe(data, kind, index) {
    if (index === null) return null;
    if (kind === "vertex") {
      const kv = data.kv[index];
      return `Bus ${data.names[index]}${kv ? ` — ${kv} kV` : ""}`;
    }
    if (kind === "edge") {
      const from = data.names[data.edges[index * 2]];
      const to = data.names[data.edges[index * 2 + 1]];
      return `Branch ${from} — ${to}`;
    }
    return null;
  }
}
