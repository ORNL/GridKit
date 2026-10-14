import { readFile } from 'node:fs/promises';
import { createServer } from 'node:http';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const require = createRequire(import.meta.url);
const { build } = require('esbuild');
const { chromium } = require('playwright');
const directory = path.dirname(fileURLToPath(import.meta.url));
const readJson = async name => JSON.parse(await readFile(new URL(name, import.meta.url), 'utf8'));
const model = await readJson('../../../../cases/PhasorDynamics/WECC240/WECC240.case.json');
const partition = await readJson('./WECC240.partition.json');
const width = 1600, height = 1280;
const panels = {
  Northwest: { x: 80, y: 100, w: 600, h: 450 },
  Southwest: { x: 920, y: 100, w: 600, h: 450 },
  NorCal: { x: 80, y: 760, w: 600, h: 450 },
  SoCal: { x: 920, y: 760, w: 600, h: 450 },
};
const sets = partition.partitions.map(p => new Set([...p.internal, ...p.external]));
const caseBuses = new Map(model.buses.map(b => [b.number, b]));
const shared = new Map();
for (const p of partition.partitions) {
  for (const b of p.external) {
    if (!shared.has(b)) shared.set(b, []);
    shared.get(b).push(p.name);
  }
}
const pair = (a, b) => [...shared].filter(([, names]) => names.includes(a) && names.includes(b)).map(([bus]) => bus).sort((a, b) => a - b);
const ports = new Map();
const key = (name, bus) => `${name}:${bus}`;
const spread = (n, a, b) => Array.from({ length: n }, (_, i) => n === 1 ? (a + b) / 2 : a + (b - a) * i / (n - 1));
function join(a, b, positions) {
  const numbers = pair(a, b);
  numbers.forEach((bus, i) => {
    const [first, second] = positions(i, numbers.length);
    ports.set(key(a, bus), first);
    ports.set(key(b, bus), second);
  });
}
join('Northwest', 'Southwest', (i, n) => {
  const y = spread(n, 160, 450)[i];
  return [{ x: 680, y, side: 'right' }, { x: 920, y, side: 'left' }];
});
join('Northwest', 'NorCal', (i, n) => {
  const x = spread(n, 180, 420)[i];
  return [{ x, y: 550, side: 'bottom' }, { x, y: 760, side: 'top' }];
});
join('Southwest', 'SoCal', (i, n) => {
  const x = spread(n, 1025, 1460)[i];
  return [{ x, y: 550, side: 'bottom' }, { x, y: 760, side: 'top' }];
});
join('NorCal', 'SoCal', (i, n) => {
  const y = spread(n, 900, 1100)[i];
  return [{ x: 680, y, side: 'right' }, { x: 920, y, side: 'left' }];
});
join('NorCal', 'Southwest', (i, n) => [
  { x: spread(n, 570, 630)[i], y: 760, side: 'top' },
  { x: 920, y: spread(n, 500, 540)[i], side: 'left' },
]);

const buses = [], locations = new Map();
for (const p of partition.partitions) {
  const panel = panels[p.name];
  const originals = p.internal.map(number => caseBuses.get(number));
  const xs = originals.map(b => b.extension.longitude), ys = originals.map(b => b.extension.latitude);
  const minX = Math.min(...xs), maxX = Math.max(...xs), minY = Math.min(...ys), maxY = Math.max(...ys);
  const occupied = [];
  for (const number of p.internal) {
    const geo = caseBuses.get(number).extension;
    let x = panel.x + 70 + (geo.longitude - minX) / (maxX - minX) * (panel.w - 140);
    let y = panel.y + 65 + (maxY - geo.latitude) / (maxY - minY) * (panel.h - 130);
    // Keep co-located substation buses distinct without changing the case data.
    const origin = { x, y };
    for (let i = 0; occupied.some(p => Math.hypot(p.x - x, p.y - y) < 18); ++i) {
      const r = 8 * Math.sqrt(i + 1), angle = i * 2.399963229728653;
      x = Math.max(panel.x + 45, Math.min(panel.x + panel.w - 45, origin.x + r * Math.cos(angle)));
      y = Math.max(panel.y + 45, Math.min(panel.y + panel.h - 45, origin.y + r * Math.sin(angle)));
      if (i > 10000) throw new Error('Unable to separate bus markers');
    }
    occupied.push({ x, y });
    locations.set(key(p.name, number), buses.length);
    buses.push({ number, partition: p.name, x, y, external: false });
  }
  for (const number of p.external) {
    const position = ports.get(key(p.name, number));
    if (!position) throw new Error(`Missing border port: ${p.name}:${number}`);
    locations.set(key(p.name, number), buses.length);
    buses.push({ number, partition: p.name, ...position, external: true });
  }
}
const branches = model.devices.filter(d => d.class === 'Branch').map(d => {
  const a = d.ports.bus1, b = d.ports.bus2;
  const i = sets.findIndex(s => s.has(a) && s.has(b));
  if (i < 0) throw new Error(`Unassigned branch ${d.id}`);
  const name = partition.partitions[i].name;
  return [locations.get(key(name, a)), locations.get(key(name, b))];
});
const links = [];
for (const [bus, names] of shared) {
  for (let i = 1; i < names.length; ++i) links.push([locations.get(key(names[i - 1], bus)), locations.get(key(names[i], bus))]);
}
const data = { buses, branches, links, width, height };
const escape = text => String(text).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&apos;' }[c]));
const text = (x, y, value, attrs = '') => `<text x="${x}" y="${y}" ${attrs}>${escape(value)}</text>`;
const frames = `<rect x="35" y="60" width="1530" height="1185" rx="8" stroke-dasharray="8 6"/>` +
  Object.values(panels).map(p => `<rect x="${p.x}" y="${p.y}" width="${p.w}" height="${p.h}" rx="8"/>`).join('');
let labels = text(55, 60, 'WECC240', 'class="heading"');
for (const p of partition.partitions) {
  const f = panels[p.name];
  labels += text(f.x + 16, f.y, p.name, 'class="heading"');
  labels += text(f.x + 16, f.y + 28, `${p.internal.length} internal, ${p.external.length} external`, 'class="count"');
}
for (const b of buses.filter(b => b.external)) {
  const dx = b.side === 'right' ? -10 : 10;
  const anchor = b.side === 'right' ? 'end' : 'start';
  labels += text(b.x + dx, b.y, b.number, `class="port" text-anchor="${anchor}"`);
}
const bundle = await build({ entryPoints: [path.join(directory, 'network.js')], bundle: true, format: 'iife', write: false,
  nodePaths: process.env.NODE_PATH?.split(path.delimiter) ?? [] });
const html = `<!doctype html><html><head><meta charset="utf-8"><title>WECC240 partitions</title><style>
*{box-sizing:border-box}body{margin:0;background:white;font-family:Arial,sans-serif;color:#171717}
#figure{width:1600px;padding-top:24px}h1{font-size:28px;margin:0 35px 12px}.legend{margin:0 35px;font-size:16px;display:flex;gap:30px;align-items:center}
.dot{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:8px}.red{background:#d64545}.blue{background:#007aff}
.line{display:inline-block;width:30px;border-top:2px dashed;margin-right:8px}.branch{border-top-style:solid;border-color:#9b9fa5}
#scene{position:relative;width:${width}px;height:${height}px}#scene>svg,#network{position:absolute;inset:0;width:100%;height:100%}
svg{overflow:visible}.frames{fill:none;stroke:#202020;stroke-width:1.5}.heading{font-size:19px;dominant-baseline:middle;paint-order:stroke;stroke:white;stroke-width:9px;stroke-linejoin:round}
.port{font-size:13px;fill:#007aff;dominant-baseline:middle;paint-order:stroke;stroke:white;stroke-width:4px}
.count{font-size:13px;fill:#666}
</style></head><body><main id="figure"><h1>WECC240 · four-region partition</h1>
<div class="legend"><span><i class="dot red"></i>Internal bus</span><span><i class="dot blue"></i>External bus copy</span><span><i class="line"></i>Same bus across partitions</span><span><i class="line branch"></i>Physical branch</span><span>243 buses · 447 branches</span></div>
<div id="scene"><svg class="frames" viewBox="0 0 ${width} ${height}">${frames}</svg><img id="network" alt="WECC240 network rendered by Latkit"><svg viewBox="0 0 ${width} ${height}">${labels}</svg></div></main>
<script>window.figureData=${JSON.stringify(data)};</script><script src="/network.js"></script></body></html>`;
const server = createServer((req, res) => {
  res.setHeader('Content-Type', req.url === '/network.js' ? 'text/javascript' : 'text/html');
  res.end(req.url === '/network.js' ? bundle.outputFiles[0].contents : html);
});
await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
let browser;
try {
  browser = await chromium.launch({ headless: true, args: ['--enable-unsafe-webgpu', '--use-angle=swiftshader', '--no-sandbox'] });
  const page = await browser.newPage({ viewport: { width, height: 1450 }, deviceScaleFactor: 2 });
  const errors = [];
  page.on('pageerror', error => errors.push(String(error)));
  page.on('console', message => { if (message.type() === 'error') errors.push(message.text()); });
  await page.goto(`http://127.0.0.1:${server.address().port}`);
  await page.waitForFunction(() => window.renderReady || window.renderError, undefined, { timeout: 120000 });
  const result = await page.evaluate(() => ({ error: window.renderError, stats: window.renderStats }));
  if (result.error || errors.length) throw new Error(result.error ?? errors.join('\n'));
  if (result.stats.vertices !== buses.length || result.stats.edges !== branches.length + links.length) throw new Error('Latkit rendered unexpected graph counts');
  await page.locator('#figure').screenshot({ path: path.join(directory, 'WECC240.png') });
  console.log(JSON.stringify({ buses: buses.length, physicalBranches: branches.length, sharedBusLinks: links.length, latkit: result.stats }, null, 2));
} finally {
  await browser?.close();
  await new Promise(resolve => server.close(resolve));
}
