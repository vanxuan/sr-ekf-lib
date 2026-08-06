# Project Tracks

This file tracks all major tracks for the project. Each track has its own detailed plan in its respective folder.

---

- [x] **Track: Separate device stillness from motion stillness: introduce a fused motionStillness metric and route velocity-domain consumers (ZUPT, coasting velocity damping, GPS stationaryWeight, diagnostics) to it, keeping device stillness for device-domain logic (mag adaptive noise, rest-rotation omegaScale)**
  *Link: [./tracks/motion_stillness_20260806/](./tracks/motion_stillness_20260806/)* — complete (3 phases, 10 tasks, 11 commits, 90 tests, 3 checkpoints)

---

- [~] **Track: Restructure code — remove 1-file rule, split into modules by concern, keep sr-ekf.ts entry point, keep zero dependencies**
  *Link: [./tracks/code_restructure_20260806/](./tracks/code_restructure_20260806/)*
