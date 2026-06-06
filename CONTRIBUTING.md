# Contributing

## Toolchain

Use [`pixi`](https://pixi.sh/) for a reproducible environment including a Rust toolchain and the external tools used for differential tests.
A recent system Rust toolchain also works for the core crate.

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
