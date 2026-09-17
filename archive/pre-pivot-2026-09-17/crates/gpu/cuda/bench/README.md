# S12 GPU-spike benchmarks (measured, RTX 4070 Ti SUPER)

Standalone CUDA micro-benchmarks that produced the S12 throughput + determinism
numbers (PLAN S12; ruling #9). Not part of any build — compile with `nvcc` on the
GPU box and run directly. The `crucible-gpu` crate (`../..`) carries the Rust→CUDA
binding + the CPU↔GPU tolerance cross-check against `crucible_solvers`.

**Box of record:** RTX 4070 Ti SUPER (16 GB, 256-bit, ~672 GB/s advertised peak),
driver 591.86, CUDA 13.3, WSL2 Ubuntu 24.04. *(Ruling #9 said "RTX 4080"; the
actual card is a 4070 Ti SUPER — same Ada gen, same 16 GB, f64 at 1/64 FP32,
~6% less bandwidth. Corrected in ruling #9 + §3.)*

Build + run (on the box, over ssh):
```
export PATH=/usr/local/cuda/bin:$PATH
nvcc -O3 -arch=sm_89 stream.cu      -o stream      && ./stream
nvcc -O3 -arch=sm_89 sweep_flux.cu  -o sweep_flux  && ./sweep_flux
nvcc -O3 -arch=sm_89 eos_project.cu -o eos_project && ./eos_project
```

## Measured (f64)

| Benchmark | Result | Determinism |
|---|---|---|
| `stream.cu` — STREAM triad (the bandwidth ceiling) | **401 GB/s** sustained (60% of 672 peak) | — |
| `sweep_flux.cu` — PPM + HLLC-Batten z-sweep, NCOMP=7 | **~0.96×10⁸ cell-updates/s** | run-twice BIT-IDENTICAL |
| `eos_project.cu` — Illinois (p,h,Z) projection | **~4.9×10⁸ cell-projections/s** (root err 2.8e-15) | run-twice BIT-IDENTICAL |

## The bottleneck finding (why the sweep is the governing number)

The sweep's ~10⁸ cups sits far below **both** ceilings — ~10 GB/s effective DRAM
traffic vs the 401 GB/s bandwidth ceiling, and the f32 twin is only **~2×** faster
(1.78×10⁸), not the 20–60× a f64-ALU wall would give. `ptxas` reports **146
registers/thread ⇒ ~25% occupancy** (register-limited): the kernel is
occupancy/latency-bound, not bandwidth- or f64-bound. So ruling #9's
"bandwidth-bound" premise is *reachable* — the path is register reduction (staged
kernels, fewer live stencil arrays) in the S13 residency + S14 profiling, not
blocked by the consumer f64 rate. The EOS projection is already ~5× the sweep and
is not the bottleneck. Class-D diffusion (CG) is unspiked (S13).

## Simplifications (spike, not residency — S13 generalizes)

- Sweep is z-direction on a uniform (r,z) field, N_θ=1: z-face annular areas cancel
  in the divergence (metric-exact for z); the r-sweep's `face_radius` weighting +
  cut apertures + mixed N_θ are S13.
- The EOS surface is a synthetic **multilinear** (p,h,Z) grid (exercises the exact
  trilinear gather + Illinois control flow); binding the real HDF5 equilibrium
  surface + the multi-root scan guard (SOLV-4 0.4.8) is S13.
- Fused single-kernel sweep (not the staged predictor/correct SDC pipeline).
