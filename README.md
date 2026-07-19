[![Build](https://github.com/fg-labs/unum/actions/workflows/check.yml/badge.svg)](https://github.com/fg-labs/unum/actions/workflows/check.yml)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fg-labs/unum/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21445800.svg)](https://doi.org/10.5281/zenodo.21445800)

# unum

`unum` is a pure-Rust HLA/KIR genotyper — a port of [T1K](https://github.com/mourisl/T1K) that rebuilds its kmer-based genotyping pipeline in Rust. It started as a byte-identical port, validated by differential-testing against the original C++ implementation, but has since diverged from T1K with deliberate bug fixes and improvements. Owning the port lets us correct T1K's latent bugs and improve behavior where `unum` can do better; every such divergence is tracked in [docs/DIVERGENCES.md](docs/DIVERGENCES.md). Much of the original byte-for-byte parity remains frozen into self-contained golden/unit tests, and `unum` has no C++ dependency.

<p>
<a href="https://fulcrumgenomics.com">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/fg-labs/unum/main/.github/logos/fulcrumgenomics-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/fg-labs/unum/main/.github/logos/fulcrumgenomics-light.svg">
  <img alt="Fulcrum Genomics" src="https://raw.githubusercontent.com/fg-labs/unum/main/.github/logos/fulcrumgenomics-light.svg" height="100">
</picture>
</a>
</p>

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your bioinformatics with `unum` and beyond.

> Note: the badge/logo URLs above point at `fg-labs/unum`; adjust them if the repository lands under a different path.

## Subcommands

- `unum build` — build reference FASTAs from an IPD-IMGT/HLA or IPD-KIR `.dat` database and a genome GTF.
- `unum extract` — extract candidate reads from FASTQ or BAM/CRAM input against a reference.
- `unum genotype` — call a genotype from candidate reads against a reference.

## Building

```bash
cargo build
```

The build is pure Rust with no C++ dependency — `unum-core` and `unum` never invoke a C++ compiler.

## Development

CI runs three gates, each wired to a `cargo` alias in [`.cargo/config.toml`](.cargo/config.toml) so they are identical locally and in CI (both under the toolchain pinned in [`rust-toolchain.toml`](rust-toolchain.toml)):

```bash
cargo ci-fmt    # rustfmt --check
cargo ci-lint   # clippy with -D warnings and the project's pedantic set
cargo ci-test   # the nextest suite
```

To catch `format`/`lint` failures before they reach CI, enable the shared pre-push hook once per clone:

```bash
git config core.hooksPath .githooks
```

This runs `cargo ci-fmt` and `cargo ci-lint` on every `git push` (bypass a single push with `git push --no-verify`). Worktrees inherit the setting, but `core.hooksPath` is per-clone, so run the command again in any freshly cloned checkout.

## Divergences from T1K

`unum` began as a byte-identical port but intentionally diverges from T1K where it can be more correct or robust. Each divergence — what changed, why, and the T1K source it departs from — is recorded in [docs/DIVERGENCES.md](docs/DIVERGENCES.md).

## License

MIT — see [LICENSE](LICENSE).

This project is a Rust port of, and derives from, [T1K](https://github.com/mourisl/T1K) (MIT, © Li Song, Bo Li, Heng Li).
