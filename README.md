[![Build](https://github.com/fulcrumgenomics/fg-t1k/actions/workflows/check.yml/badge.svg)](https://github.com/fulcrumgenomics/fg-t1k/actions/workflows/check.yml)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fulcrumgenomics/fg-t1k/blob/main/LICENSE)

# fg-t1k

Fulcrum-owned Rust HLA/KIR genotyper — a strangler-fig port of [T1K](https://github.com/mourisl/T1K), rebuilding its kmer-based genotyping pipeline in Rust while retaining the original C++ implementation as a differential-testing oracle during the migration.

<p>
<a href="https://fulcrumgenomics.com">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/fg-t1k/main/.github/logos/fulcrumgenomics-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/fg-t1k/main/.github/logos/fulcrumgenomics-light.svg">
  <img alt="Fulcrum Genomics" src="https://raw.githubusercontent.com/fulcrumgenomics/fg-t1k/main/.github/logos/fulcrumgenomics-light.svg" height="100">
</picture>
</a>
</p>

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your bioinformatics with fg-t1k and beyond.

## Status

This repo is under active development and is not yet published. We are in Phase 0/1 of the port (workspace scaffold + vendored T1K oracle for differential testing). See the roadmap in [issue #1](https://github.com/fulcrumgenomics/fg-t1k/issues/1) for the full plan and current progress.

See `docs/superpowers/specs/` for the design and `docs/superpowers/plans/` for the implementation plan.

## Building

```bash
cargo build
```

The default build is pure Rust with no C++ dependency — `fg-t1k-core` and `fg-t1k` never invoke a C++ compiler.

The `t1k-sys` cargo feature builds the vendored T1K C++ oracle (the bundled `samtools-0.1.19` plus the original `genotyper`/`analyzer`/`fastq-extractor`/`bam-extractor` binaries) for differential testing against the Rust port. It requires a C++ toolchain and zlib:

```bash
cargo build --features t1k-sys
cargo test --features t1k-sys
```

## License

MIT — see [LICENSE](LICENSE).

This project vendors and derives from [T1K](https://github.com/mourisl/T1K) (MIT, © Li Song, Bo Li, Heng Li) under `vendor/t1k/`, with the upstream license preserved and provenance tracked in [`vendor/t1k/PROVENANCE.md`](vendor/t1k/PROVENANCE.md).
