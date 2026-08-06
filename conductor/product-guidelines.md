# Product Guidelines

## Product Communication & Positioning

- **Lead with the problem:** describe `sr-ekf` in terms of the navigation problem it solves (smooth, drift-resistant fusion) before any implementation detail.
- **Be concrete and quantitative:** prefer numbers over adjectives ("zero allocations on the predict hot path", "78 tests", "2D CTRA 8-state model") so engineers can assess fit immediately.
- **Match the established voice:** `README.md` and `AGENTS.md` use terse, precise, specification-like prose. New documentation must match that tone.

## Documentation Standards

- **AGENTS.md is the single source of truth** for the algorithm and API. Any behavior change must update it in the same commit.
- **Every public API member is documented** in AGENTS.md (signature, semantics, units, defaults).
- **Rationale over recipes:** document *why* a design exists (e.g., "removed a_bias_y — unobservable from forward acceleration alone") alongside *what* it does.
- **Symptom-driven entries:** when a fix addresses a reported behavior, name the symptom it prevents (e.g., "car icon moves on a line", "heading stuck after standing still").

## Code & API Quality

- **Determinism is a feature:** never introduce time-of-day reads; all timestamps must be explicit parameters.
- **Mobile-first performance:** the `predict()` hot path must remain allocation-free; new features must fit the fixed-size matrix architecture.
- **Numerical stability:** any covariance/cross-covariance manipulation must preserve the Cholesky (square-root) invariant; QR diagonals must stay positive.
- **Minimal API surface:** expose only what applications need; keep internals private to preserve future flexibility.
- **Semantic versioning:** breaking API changes require a major version; additive config/behavior defaults must not silently change existing callers.

## Testing Standards

- **Every behavior change ships with a regression test** that names the symptom it guards against.
- **Deterministic simulations only** — no randomness, no wall-clock dependence.
- **Unit + numerical verification:** keep both the QR-verification suite and behavior tests green before any merge.
- **Document the verification command** in the validation section (`npm test`, `npm run build`).

## UX Principles (for integrators)

- **Fail safe:** corrupt covariance must trigger safeguards (coast/recovery), never silent NaN propagation.
- **Observable state:** expose diagnostics and debug snapshots so integrators can diagnose sensor issues (misaligned compass, stale GPS) without internals access.
- **Predictable noise behavior:** gates and inflation caps should degrade gracefully in urban multipath rather than hard-rejecting.
