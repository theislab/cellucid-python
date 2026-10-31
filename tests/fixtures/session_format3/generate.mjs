// Writes `session.cellucid-session` and `expected.json` beside this file with
// the web app's own session encoders (`cellucid/assets/js/app/session/**`), so
// the Python reader is tested against bytes the viewer writes.
//
//   node tests/fixtures/session_format3/generate.mjs   (from cellucid-python;
//   the cellucid repository must be a sibling checkout)
import { writeFileSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const web = path.resolve(here, '../../../../cellucid/assets/js/app/session');
const load = async relative => import(path.join(web, relative));
const { writeBundle } = await load('bundle/writer.js');
const { gzipCompress } = await load('codecs/gzip.js');
const { digestCellSet } = await load('codecs/cell-set.js');
const highlightsCells = await load('contributors/highlights-cells.js');
const highlightsMeta = await load('contributors/highlights-meta.js');
const userDefinedCodes = await load('contributors/user-defined-codes.js');
const analysisWindows = await load('contributors/analysis-windows.js');
const analysisDocked = await load('contributors/analysis-docked.js');
const { DE_OPTIONS } = await import(path.resolve(web, '../analysis/options/de-options.js'));
const { validateDataMapping } = await import(path.resolve(web, '../../data/data-mapping.js'));

const N = 1024;
// obs "ctype": A for cells 0-99, B for 100-899, C for 900-1023.
const categories = ['A', 'B', 'C'];
const codes = new Uint8Array(N).map((_, i) => (i < 100 ? 0 : i < 900 ? 1 : 2));
const ctype = { key: 'ctype', kind: 'category', categories, codes };
const range = (a, b) => Array.from({ length: b - a }, (_, i) => a + i);
const subsetB = range(100, 900).filter(i => i % 4 !== 1);
const groups = [
  { id: 'highlight_1', type: 'lasso', label: 'Two cells', cells: [3, 50] },
  { id: 'highlight_2', type: 'lasso', label: 'Every other', cells: range(0, N).filter(i => i % 2 === 0) },
  { id: 'highlight_3', type: 'lasso', label: 'A block', cells: range(10, 411) },
  { id: 'highlight_4', type: 'category', label: 'ctype: A', cells: range(0, 100), category: 0 },
  { id: 'highlight_5', type: 'category', label: 'ctype: B (visible)', cells: subsetB, category: 1 },
];
const liveGroups = groups.map(g => {
  const base = { id: g.id, type: g.type, label: g.label, enabled: true, cellCount: g.cells.length, cellIndices: Uint32Array.from(g.cells) };
  if (g.type !== 'category') return base;
  return { ...base, fieldKey: 'ctype', fieldIndex: 0, fieldSource: 'obs', categoryIndex: g.category, categoryName: categories[g.category] };
});
const userField = { _userDefinedId: 'user_cat_1', kind: 'category', key: 'regions', codes: new Uint8Array(N).map((_, i) => (i * 7919) % 3) };
const state = {
  pointCount: N,
  obsData: { fields: [ctype] },
  varData: { fields: [{ key: 'gene_0', kind: 'continuous' }] },
  activeFieldSource: null,
  activeFieldIndex: -1,
  activeVarFieldIndex: -1,
  viewContexts: new Map(),
  getHighlightPages: () => [{ id: 'page_1', name: 'Page 1', color: '#2563eb', highlightedGroups: liveGroups }],
  getActivePageId: () => 'page_1',
  getUserDefinedFieldsRegistry: () => ({ getAllFields: () => [userField] }),
};
const ctx = { state };

const overlays = {
  renames: { fields: {}, categories: {} },
  deletedFields: { deleted: [], purged: [] },
  userDefinedFields: [{
    id: 'user_cat_1', source: 'obs', kind: 'category', key: 'regions', categories: ['R1', 'R2', 'rest'],
    isDeleted: false, isPurged: false, codesLength: N, codesType: 'Uint8Array', sourceField: null, operation: null,
    sourcePages: [], overlapStrategy: 'first', overlapLabel: null, intersectionLabels: null, uncoveredLabel: null, createdAt: 1,
  }],
};
// The analysis chunks go through the app's own contributors, with records in
// the shape a form-based analysis exports (`{options, selectedPages, view}`)
// and the window manager's remembered sizes.
const deSettings = {
  options: DE_OPTIONS.withDefaults({ pageA: 'page_1', pageB: 'page_2' }),
  selectedPages: ['page_1', 'page_2'],
  view: { disclosures: {} },
};
const windowsManager = {
  exportSessionWindows: () => [{
    modeId: 'differential', open: true, rect: { left: 1, top: 2, width: 300, height: 200 }, z: 1,
    settings: structuredClone(deSettings), result: { genes: ['gene_0'], log2fc: new Float32Array([1.5]) }, headerTitle: 'DE',
  }],
  exportViewSizes: () => ({ window: { differential: { width: 300, height: 200 } }, expandedView: {} }),
};
const dockedModule = {
  exportSessionState: () => ({
    openMode: 'differential',
    modes: { differential: { settings: structuredClone(deSettings), result: { genes: ['gene_0'], pValue: new Float64Array([0.001]) } } },
  }),
};
const [windowsChunk] = analysisWindows.capture({ analysisWindowManager: windowsManager });
const [dockedChunk] = analysisDocked.capture({ comparisonModule: dockedModule });
const json = value => new TextEncoder().encode(JSON.stringify(value));
const chunks = [
  { id: 'core/field-overlays', contributorId: 'field-overlays', priority: 'eager', kind: 'json', codec: 'gzip', label: 'Field overlays', datasetDependent: true, bytes: json(overlays) },
  ...userDefinedCodes.capture(ctx).map(c => ({ ...c, bytes: c.payload })),
  { id: 'cinematic/camera', contributorId: 'cinematic-camera', priority: 'eager', kind: 'json', codec: 'gzip', label: 'Cinematic camera path', datasetDependent: true, bytes: json({ notMaterializedByPython: true }) },
  { id: 'core/state', contributorId: 'core-state', priority: 'eager', kind: 'json', codec: 'gzip', label: 'Core state', datasetDependent: true, bytes: json({ notMaterializedByPython: true }) },
  { id: 'ui/dockable-layout', contributorId: 'dockable-layout', priority: 'eager', kind: 'json', codec: 'gzip', label: 'Floating panels', datasetDependent: false, bytes: json({ panels: [] }) },
  { id: 'ui/figure-export', contributorId: 'figure-export', priority: 'eager', kind: 'json', codec: 'gzip', label: 'Figure export', datasetDependent: false, bytes: json({ version: 1 }) },
  { ...windowsChunk, bytes: windowsChunk.payload },
  { ...dockedChunk, bytes: dockedChunk.payload },
  ...highlightsMeta.capture(ctx).map(c => ({ ...c, bytes: json(c.payload) })),
  ...highlightsCells.capture(ctx).map(c => ({ ...c, bytes: c.payload })),
];
// Eager chunks first, each group in contributor order (user-defined codes of
// a field no view shows are lazy and go before the highlight cells).
const ordered = [...chunks.filter(c => c.priority === 'eager'), ...chunks.filter(c => c.priority === 'lazy')];
const stored = [];
const manifestChunks = [];
for (const c of ordered) {
  const gz = await gzipCompress(c.bytes, { signal: null });
  stored.push(gz);
  manifestChunks.push({ id: c.id, contributorId: c.contributorId, priority: c.priority, kind: c.kind, codec: c.codec, label: c.label, datasetDependent: c.datasetDependent, storedBytes: gz.byteLength, uncompressedBytes: c.bytes.byteLength });
}
const manifest = {
  createdAt: '2026-10-22T00:00:00.000Z',
  datasetFingerprint: {
    datasetId: 'format3', name: 'Format 3 fixture', cellCount: N, varCount: 1,
    cellOrder: { embedding: { key: 'X_umap', columns: [0, 1] }, dimension: 2, digest: '0123456789abcdef' },
    // Checked by the app's own validator, so the record is in its current shape.
    dataMapping: validateDataMapping({
      kind: 'anndata',
      embeddings: { enabled: [{ key: 'X_umap', columns: [0, 1] }], defaults: { '1d': null, '2d': { key: 'X_umap', columns: [0, 1] }, '3d': null } },
      knn: null, geneNames: 'index', expression: 'X', log1p: {}, highlyVariable: null, latent: null, vectorFields: [],
    }),
  },
  source: { kind: 'file', name: 'format3.h5ad', dataset: 'format3' },
  requires: { obs: ['ctype'], var: [] },
  chunks: manifestChunks,
};
const blob = writeBundle({ manifest, chunks: stored });
writeFileSync(path.join(here, 'session.cellucid-session'), new Uint8Array(await blob.arrayBuffer()));
const encodings = Object.fromEntries(highlightsCells.capture(ctx).map(c => [c.id, c.payload[0]]));
writeFileSync(path.join(here, 'expected.json'), JSON.stringify({
  cellCount: N,
  ctype: Array.from(codes, c => categories[c]),
  groups: Object.fromEntries(groups.map(g => [g.id, g.cells])),
  encodings,
  regionsCodes: Array.from(userField.codes),
  regionsEncoding: userDefinedCodes.capture(ctx)[0].payload[0],
  digests: { A: digestCellSet(Uint32Array.from(range(0, 100))), B: digestCellSet(Uint32Array.from(range(100, 900))), large: digestCellSet(Uint32Array.from([0, 1, 65535, 65536, 561946, 4294967295])) },
}) + '\n');
console.log('wrote', Object.values(encodings).join(','));
