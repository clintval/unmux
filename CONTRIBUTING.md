# Contributing

## Toolchain

Use [`pixi`](https://pixi.sh/) for a reproducible environment including a Rust toolchain and the external tools used for  tests.

## Building

The matcher (`sassy`) uses SIMD and requires AVX2 (x86) or NEON (ARM).
On x86 you must tell the compiler the target supports it:

```bash
# Native build (uses every instruction the build host has):
RUSTFLAGS="-C target-cpu=native" cargo build --release

# Portable AVX2 build (runs on any x86-64-v3 host):
RUSTFLAGS="-C target-cpu=x86-64-v3" cargo build --release
```

On hosts without AVX2/NEON, build the scalar fallback; it runs correctly at reduced throughput and needs no `target-cpu` flag:

```bash
cargo build --features scalar
```

The `dist` profile (`cargo build --profile dist`) simply creates a release optimized binary.
But for distribution, one can run `pixi run dist-multivers` (`cargo multivers --profile dist`) which packs one slice per CPU tier (`x86-64-v3` and `x86-64-v4`, per `[package.metadata.multivers.x86_64]`) behind a CPUID dispatcher chosen at runtime, so a single binary runs natively on both AVX2-only and AVX-512 hosts and a non-AVX-512 host never executes AVX-512 code.
Both slices run sassy's 4-lane matcher; the `x86-64-v4` slice gets AVX-512-aware codegen for it.
Sassy's wider 8-lane `avx512` feature (about 2x on the batch path) is not yet forwarded: it fails to compile in sassy 0.2.3 (`search.rs` hardcodes `wide::u64x4` while the feature widens the lane type to `u64x8`), so enabling the 8-lane path is blocked on an upstream fix.

## Mandatory Pre-submission Checks

All must pass before a change is submitted:

```bash
cargo build --release
cargo test --all --no-fail-fast --verbose
cargo fmt -- --check
cargo clippy --all-targets --all-features -- -D warnings
```

Equivalent pixi tasks:

```bash
pixi run build
pixi run test
pixi run fmt-check
pixi run lint
```
