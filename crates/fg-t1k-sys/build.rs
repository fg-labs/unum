use std::path::{Path, PathBuf};
use std::process::Command;

/// Probes whether `-l{lib_name}` resolves for the given compiler by attempting to
/// link a trivial `int main(){}` translation unit against it. Used to build a
/// portable system-lib set for the `bam-extractor` oracle: htslib 1.19.1's static
/// `libhts.a` pulls in a different set of system libs depending on how it was
/// built (CRAM/curl backends, libdeflate vs bundled deflate, etc.), and that set
/// differs between macOS (which ships bz2/curl/lzma via the system) and a bare
/// Linux runner (where only what's `apt-get install`ed is present). Rather than
/// hardcode a single platform's lib set, probe for each candidate and link only
/// the ones that actually resolve.
fn probe_lib(cxx: &str, out_dir: &Path, lib_name: &str) -> bool {
    let probe_src = out_dir.join(format!("probe_{lib_name}.cpp"));
    std::fs::write(&probe_src, "int main() { return 0; }\n").expect("write probe source");
    let probe_bin = out_dir.join(format!("probe_{lib_name}"));
    let status = Command::new(cxx)
        .arg(&probe_src)
        .arg("-o")
        .arg(&probe_bin)
        .arg(format!("-l{lib_name}"))
        .status();
    status.map(|s| s.success()).unwrap_or(false)
}

fn main() {
    if std::env::var("CARGO_FEATURE_T1K_SYS").is_err() {
        return;
    }
    let vendor = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../vendor/t1k")
        .canonicalize()
        .expect("vendor/t1k");
    let out = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    let cxx = std::env::var("CXX").unwrap_or_else(|_| "c++".into());

    // `hts-sys` (a *direct* dependency, see Cargo.toml) forwards its `links =
    // "hts"` metadata as `DEP_HTS_*` env vars: the include dir containing
    // `htslib/sam.h`, and the dir containing `libhts.a`. This is the SAME
    // htslib 1.19.1 that `rust-htslib` links into the Rust side -- building
    // the oracle against it (rather than the old bundled samtools-0.1.19)
    // means exactly one htslib copy exists in any process that links both
    // the oracle's `libbam`-equivalent code and rust-htslib, eliminating the
    // duplicate-symbol collision the old two-htslib setup had.
    let hts_include = PathBuf::from(
        std::env::var("DEP_HTS_INCLUDE").expect("DEP_HTS_INCLUDE not set by hts-sys"),
    );
    let hts_libdir =
        PathBuf::from(std::env::var("DEP_HTS_LIBDIR").expect("DEP_HTS_LIBDIR not set by hts-sys"));
    assert!(
        hts_include.join("htslib/sam.h").exists(),
        "expected {}/htslib/sam.h to exist (from hts-sys)",
        hts_include.display()
    );

    // T1K's vendored alignments.hpp hardcodes `#include
    // "htslib-1.15.1/htslib/sam.h"` under `#ifdef HTSLIB` (vendored source,
    // not ours to edit). Build an `-I` alias directory so that literal path
    // resolves to hts-sys's headers instead: OUT_DIR/htsalias/htslib-1.15.1/
    // is a directory containing a symlink named `htslib` pointing at
    // `{DEP_HTS_INCLUDE}/htslib`.
    let htsalias = out.join("htsalias");
    let htsalias_versioned = htsalias.join("htslib-1.15.1");
    std::fs::create_dir_all(&htsalias_versioned).expect("create htsalias dir");
    let alias_link = htsalias_versioned.join("htslib");
    if alias_link.symlink_metadata().is_err() {
        std::os::unix::fs::symlink(hts_include.join("htslib"), &alias_link)
            .expect("symlink htslib alias");
    }

    // (binary name, source, needs libhts)
    for (name, src, needs_hts) in [
        ("genotyper", "Genotyper.cpp", false),
        ("analyzer", "Analyzer.cpp", false),
        ("fastq-extractor", "FastqExtractor.cpp", false),
        ("bam-extractor", "BamExtractor.cpp", true),
    ] {
        let mut c = Command::new(&cxx);
        c.args(["-O3", "-g", "-o"]).arg(out.join(name)).arg(vendor.join(src));
        if needs_hts {
            // Only BamExtractor.cpp touches BAM I/O (via alignments.hpp);
            // the other three oracle binaries use ReadFiles.hpp/kseq.h and
            // need nothing beyond zlib + pthread.
            //
            // hts-sys's static `libhts.a` (htslib 1.19.1) pulls in a system-lib
            // set that varies by platform/build (CRAM + curl backends, and
            // htslib 1.19 prefers libdeflate when present). Hardcoding one set
            // works only on whichever platform provided all of them incidentally
            // (macOS does; a minimal `apt-get install build-essential zlib1g-dev`
            // Linux runner does not, and fails with `cannot find -lbz2`/`-lcurl`).
            // Instead: always link the libs `libhts.a` unconditionally needs
            // (pthread, hts itself, z), then probe for the optional/backend libs
            // and link only those that are actually present. Order matters for
            // static linking: `-lhts` must come BEFORE the system libs it
            // depends on, since `ld` resolves undefined symbols in static
            // archives left-to-right.
            c.arg("-DHTSLIB")
                .arg(format!("-I{}", htsalias.display()))
                .arg(format!("-L{}", hts_libdir.display()))
                .args(["-lpthread", "-lhts", "-lz"]);
            for lib in ["bz2", "lzma", "curl", "deflate", "m"] {
                if probe_lib(&cxx, &out, lib) {
                    c.arg(format!("-l{lib}"));
                }
            }
        } else {
            c.args(["-lpthread", "-lz"]);
        }
        assert!(
            c.status().unwrap_or_else(|e| panic!("compile {name}: {e}")).success(),
            "compile {name} failed"
        );
    }
    println!("cargo:rustc-env=FG_T1K_ORACLE_DIR={}", out.display());
    println!("cargo:rerun-if-changed={}", vendor.display());

    // Header-only shim: includes T1K headers (e.g. KmerCode.hpp) but links zero
    // T1K .cpp files. Defines the extern nucToNum/numToNuc tables itself since
    // they are declared (not defined) in the headers.
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    cc::Build::new()
        .cpp(true)
        .file(manifest_dir.join("shim/shim.cpp"))
        .include(&vendor) // T1K headers
        .include(&htsalias) // resolves alignments.hpp's "htslib-1.15.1/htslib/sam.h" to hts-sys
        .include(manifest_dir.join("shim"))
        .define("HTSLIB", None) // alignments.hpp: build against real htslib, matching the oracle above
        .flag_if_supported("-std=c++11")
        .opt_level(3)
        .compile("t1k_shim");
    println!("cargo:rerun-if-changed={}", manifest_dir.join("shim").display());

    // SeqSet.hpp (pulled in by shim.cpp for the SeqSet opaque handle) drags
    // in ReadFiles.hpp (zlib-backed FASTA/FASTQ reading via kseq.h) and uses
    // pthread_mutex_t/pthread_mutex_init/pthread_mutex_destroy directly --
    // link the same two system libs for the shim's final link step.
    println!("cargo:rustc-link-lib=z");
    println!("cargo:rustc-link-lib=pthread");

    // alignments.hpp (pulled in by shim.cpp for the Alignments opaque handle)
    // is now compiled WITH -DHTSLIB (see the `.define("HTSLIB", None)` call
    // above), matching the bam-extractor oracle binary. Unlike that oracle
    // binary -- a standalone `cc`-invoked executable that must link `-lhts`
    // itself -- the shim is compiled into a static lib (`libt1k_shim.a`)
    // that gets linked into Rust *test* binaries alongside `rust-htslib`.
    // Those test binaries already pull in hts-sys's `libhts` transitively
    // through rust-htslib, so we deliberately do NOT emit a second
    // `cargo:rustc-link-lib=hts` here: adding one would either be a no-op
    // (same static lib, still one copy) or, if ever this crate stopped
    // depending on rust-htslib in a test binary, a hard link error -- either
    // way it is not this crate's job to decide how htslib gets linked into
    // its own reverse-dependents. This is precisely the "exactly one htslib
    // project-wide" invariant this task establishes.
}
