[![Build](https://github.com/fg-labs/unum/actions/workflows/check.yml/badge.svg)](https://github.com/fg-labs/unum/actions/workflows/check.yml)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fg-labs/unum/blob/main/LICENSE)

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

## Divergences from T1K

`unum` began as a byte-identical port but intentionally diverges from T1K where it can be more correct or robust. Each divergence — what changed, why, and the T1K source it departs from — is recorded in [docs/DIVERGENCES.md](docs/DIVERGENCES.md).

## License

MIT — see [LICENSE](LICENSE).

This project is a Rust port of, and derives from, [T1K](https://github.com/mourisl/T1K) (MIT, © Li Song, Bo Li, Heng Li).
