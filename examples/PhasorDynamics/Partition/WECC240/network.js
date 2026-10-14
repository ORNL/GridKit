import { createData, textColumn } from '@latkit/model';
import { createGpu, createColormap } from '@latkit/gpu';
import { createNetwork } from '@latkit/network';

const red = [214 / 255, 69 / 255, 69 / 255, 1];
const blue = [0, 122 / 255, 1, 1];
const index = type => ({ source: 'WECC240-partitions', type, version: '1' });
const numeric = values => ({ kind: 'numeric', offset: 0, length: values.length, values: Float64Array.from(values) });
const reference = values => ({ kind: 'reference', offset: 0, length: values.length, values: Uint32Array.from(values), index: index('Bus') });
const rows = (type, count, columns) => ({ kind: 'rows', index: index(type), rows: { kind: 'range', offset: 0, count }, columns });

async function render() {
  const { buses, branches, links, width, height } = window.figureData;
  const edgeFields = { from: { type: { kind: 'reference', to: 'Bus' } }, to: { type: { kind: 'reference', to: 'Bus' } } };
  const source = createData({ types: {
    Bus: { fields: { x: { type: 'float64' }, y: { type: 'float64' }, external: { type: 'float64' }, label: { type: 'text' } } },
    Branch: { fields: edgeFields },
    SharedBus: { fields: edgeFields },
  } }, [
    rows('Bus', buses.length, {
      x: numeric(buses.map(b => b.x)), y: numeric(buses.map(b => -b.y)),
      external: numeric(buses.map(b => Number(b.external))),
      label: textColumn(buses.map(b => b.external ? '' : String(b.number))),
    }),
    ...[['Branch', branches], ['SharedBus', links]].map(([type, edges]) => rows(type, edges.length, {
      from: reference(edges.map(e => e[0])), to: reference(edges.map(e => e[1])),
    })),
  ]);
  const gpu = await createGpu({ validate: true });
  const network = createNetwork(gpu, {
    source,
    vertices: { Bus: {
      x: 'x', y: 'y',
      color: { field: 'external', domain: [0, 1], colormap: createColormap({ kind: 'categorical', colors: [red, blue] }) },
      radiusPx: { field: 'external', domain: [0, 1], range: [3.5, 5] },
      labels: { field: 'label', color: red, fontSizePx: 11, maxCount: buses.length },
    } },
    edges: {
      Branch: { ends: ['from', 'to'], color: [0.3, 0.33, 0.38, 0.38], widthPx: 1 },
      SharedBus: { ends: ['from', 'to'], color: [0, 0, 0, 1], widthPx: 1.5, dash: true },
    },
    camera: { projection: 'flat', fit: false, center: [width / 2, -height / 2], scale: 1 },
    background: [1, 1, 1, 0], grid: false, shadows: false, input: 'none',
    edgeSpacingPx: 2, vertexMinRadiusPx: 3.5, labelHaloPx: 2,
  });
  try {
    const blob = await network.image({ width, height, pixelRatio: 1 });
    const image = document.getElementById('network');
    image.src = URL.createObjectURL(blob);
    await image.decode();
    window.renderStats = network.stats();
    window.renderReady = true;
  } finally {
    network.destroy();
    gpu.destroy();
  }
}

render().catch(error => { window.renderError = String(error.stack ?? error); });
