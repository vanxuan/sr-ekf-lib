# Technology Stack

## Language

- **TypeScript 5.8** — the library is authored in modern TypeScript, compiled to ES modules (`"type": "module"`), shipping both JS output and `.d.ts` type declarations.
- Strict, fixed-size, allocation-free coding style targeting mobile JS runtimes.

## Dependencies

- **Runtime dependencies: none.** `sr-ekf` is a single-file, zero-runtime-dependency library (only `src/sr-ekf.ts` is shipped as the implementation).
- **Dev dependencies:** `typescript` (build) and `vitest` (testing) only.

## Build

- **Compiler:** `tsc` (`npm run build`) emitting `dist/sr-ekf.js` + `dist/sr-ekf.d.ts`.
- **Package exports:** dual `import`/`types` entries via `exports` map in `package.json`.
- **Publish hygiene:** `prepublishOnly` runs the build; `files` whitelist (`dist`, `src`, `AGENTS.md`) keeps the package lean.

## Testing

- **Vitest 3.x** — fast, native-ESM test runner.
- **Suites:** `tests/qr-verification.test.ts` (9 numerical QR/Cholesky checks) and `tests/sr-ekf.test.ts` (69 unit/behavior tests). All deterministic simulations.
- **Commands:** `npm test` (run once), `npm run test:watch`, and `npm run benchmark` for performance baselines.

## CI / Tooling

- **GitHub Actions** (`.github/workflows`) for automated test runs.
- **Config:** `tsconfig.json` (build), `vitest.config.ts` (tests), `vitest.benchmark.ts` (benchmarks).
- **Docs:** `AGENTS.md` is the authoritative algorithm/API reference; `README.md` is the public-facing guide.

## Platform Targets

- **Mobile:** modern iOS/Android WebViews and React Native/Hermes JS engines.
- **Runtime constraints:** Float64Array-backed matrices, fixed `O(N³)` operations, zero allocations on the `predict()` hot path.
