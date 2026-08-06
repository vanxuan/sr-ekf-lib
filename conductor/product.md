# Product Guide

## # Initial Concept

Production-grade Square-Root Extended Kalman Filter for real-time IMU + GPS sensor fusion on mobile devices. Fuses accelerometer, gyroscope, GPS, and magnetometer data into a smooth, drift-resistant navigation solution using a 2D CTRA (Constant Turn Rate and Acceleration) motion model.

## Product Vision

`sr-ekf` is a modular zero-runtime-dependency TypeScript library (single entry point `sr-ekf.ts`) that turns noisy, low-rate GPS plus high-rate IMU (accelerometer/gyroscope) into a smooth, accurate, drift-resistant navigation solution on resource-constrained mobile devices. It is designed to be embedded directly into applications — no build step, no dependencies, no server.

## Target Audience

- **Mobile navigation / telematics apps** that need smooth position and heading between GPS fixes (tunnels, urban canyons, parking garages).
- **Handheld & wearable devices** (phones, dashcams) where the device orientation relative to the vehicle is arbitrary.
- **Engineers** integrating sensor fusion into TypeScript/JavaScript applications who value deterministic behavior, numerical stability, and predictable mobile performance.

## Core Value Proposition

- **Drift-resistant dead reckoning:** IMU prediction with on-line bias compensation (accelerometer, gyroscope) sustains accurate position/heading during GPS outages.
- **Robust GPS fusion:** Mahalanobis-gated updates, direction-aware outlier rejection, robust M-estimation, and anisotropic noise models keep the filter stable in urban multipath environments.
- **Sensor-rich heading:** Magnetometer heading with auto-calibrating magnetic declination as an EKF state; GPS velocity direction authority at speed.
- **Zero-dependency, modular:** `npm install`, import, use. Internally organized by concern (math, config, CTRA, diagnostics). Fully deterministic — no clock reads, all timestamps explicit.

## Key Features & Capabilities

- **Square-Root EKF** — Cholesky-factored covariance (`P = S·Sᵀ`) with QR-based prediction/update for numerical stability.
- **8-state 2D CTRA model** — position, speed, heading, sideslip, accelerometer/gyroscope biases, magnetic declination.
- **GPS latency compensation** — rewind/replay circular buffer for delayed GPS timestamps.
- **Zero-Velocity Update (ZUPT)** and **Zero Angular Rate Update (ZARU)** — bias calibration at rest.
- **Magnetometer fusion** — magnetic declination auto-calibration; speed-gated heading authority.
- **Robust M-estimation** (Cauchy/Huber) and **adaptive noise scaling** against GPS outliers.
- **Anisotropic & adaptive process/measurement noise** — speed-scaled, IMU-energy-scaled, direction-aware.
- **Mobile-optimized** — zero allocations on the `predict()` hot path, fixed-size matrix math, `Float64Array`-backed buffers.
- **Barometric ramp-speed estimation** — constrains forward speed on grades via pitch when orientation is set.

## Design Principles

1. **Deterministic:** fully reproducible behavior — no `Date.now()`; every timestamp is an explicit parameter.
2. **Production-grade:** NaN/Inf safeguards, covariance health checks, coasting recovery, auto-divergence detection.
3. **Mobile-first performance:** allocation-free hot paths, fixed-size `O(N³)` matrix operations.
4. **Single source of truth:** one `AGENTS.md` documenting the full algorithm and API surface.
5. **Observable:** diagnostics, debug snapshot, and zero-allocation state reads for render loops.

## Non-Goals

- Not a complete GNSS/INS suite (no 3D attitude estimation, barometer altitude state, or GNSS raw-observation processing).
- Not a sensor abstraction/hardware layer — consumes already-orientated ENU-frame IMU data (with optional device-orientation alignment).
- Not server-side processing — designed for on-device, real-time, low-latency fusion.
