use std::path::PathBuf;
use std::process::Command;

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

    // Build bundled samtools (libbam.a) via its own Makefile.
    let sam = vendor.join("samtools-0.1.19");
    if !sam.join("libbam.a").exists() {
        assert!(
            Command::new("make")
                .current_dir(&sam)
                .status()
                .expect("make")
                .success(),
            "samtools build failed"
        );
    }

    // (binary name, source, needs -lbam)
    for (name, src, bam) in [
        ("genotyper", "Genotyper.cpp", false),
        ("analyzer", "Analyzer.cpp", false),
        ("fastq-extractor", "FastqExtractor.cpp", false),
        ("bam-extractor", "BamExtractor.cpp", true),
    ] {
        let mut c = Command::new(&cxx);
        c.args(["-O3", "-g", "-o"])
            .arg(out.join(name))
            .arg(vendor.join(src))
            .arg(format!("-I{}", sam.display()))
            .arg(format!("-L{}", sam.display()))
            .args(["-lpthread", "-lz"]);
        if bam {
            c.arg("-lbam");
        }
        assert!(
            c.status()
                .unwrap_or_else(|e| panic!("compile {name}: {e}"))
                .success(),
            "compile {name} failed"
        );
    }
    println!("cargo:rustc-env=FG_T1K_ORACLE_DIR={}", out.display());
    println!("cargo:rerun-if-changed={}", vendor.display());
}
