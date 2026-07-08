// `unum-core` stays `#![forbid(unsafe_code)]` -- see that crate's
// `apply_reference` doc comment. This binary crate is downgraded to
// `#![deny(unsafe_code)]` (still -D warnings under `cargo ci-lint`, so any
// new `unsafe` must be an explicit, reviewed `#[allow(unsafe_code)]`) solely
// to permit the single `set_var` call in `main` below, which is the only
// process-wide place htslib's default (network-capable) CRAM reference chain
// can be neutralized before any BAM/CRAM reader -- including a stdin CRAM,
// whose container format can't be sniffed from a non-seekable pipe -- is
// opened.
#![deny(unsafe_code)]

mod cli;
mod stages;

use anyhow::Context as _;
use clap::Parser;
use cli::{Cli, Commands};

/// Use mimalloc as the global allocator. The genotyper's hot read-assignment
/// loop allocates many short-lived per-read/per-overlap buffers across rayon
/// workers; the system allocator (macOS libmalloc) serializes these on
/// per-size-class locks, capping multi-thread scaling. mimalloc's per-thread
/// heaps remove that contention. Allocator choice does not change any output
/// (byte-identical). The `#[global_allocator]` static is a safe declaration --
/// mimalloc's `GlobalAlloc` impl (with its own unsafe) lives inside the crate,
/// so `#![deny(unsafe_code)]` on this crate still holds (no `unsafe` block is
/// needed here).
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// How many freshly timestamped names [`create_private_ref_path_dir`] will
/// try before giving up. `std::fs::DirBuilder::create` (non-recursive) fails
/// with `AlreadyExists` if the candidate path is already occupied --
/// including by an attacker who pre-created it -- so a collision is a signal
/// to retry under a new name, never to reuse what's there. A handful of
/// retries turns a pathological run of collisions into a clean startup
/// error instead of an infinite loop; it is not a security boundary (the
/// exclusivity of `create_dir`/`DirBuilder::create` is).
const REF_PATH_DIR_CREATE_ATTEMPTS: u32 = 8;

/// Creates a fresh, empty, `0700` (owner-only) directory under the OS temp
/// namespace and returns its path, for use as `REF_PATH`'s local
/// reference-search directory (see
/// [`neutralize_cram_ref_path_network_fallback`]).
///
/// This replaces an earlier, simpler design that pointed `REF_PATH` at a
/// single fixed path (`/tmp/unum-no-network-ref-path-sentinel`) computed at
/// compile time. That fixed path is itself an attack surface on a
/// multi-tenant host: `/tmp` is world-writable, and htslib treats a `REF_PATH`
/// directory as a lookup root (`<dir>/<md5[:2]>/<md5[2:]>`) -- so a co-resident
/// attacker who pre-creates `<fixed-path>/xx/yyy...` before this process
/// starts can make a `-r`-less CRAM decode silently succeed against
/// attacker-supplied bytes instead of failing. Making the directory (a)
/// per-process (a fresh name every run, unpredictable before the process
/// exists), (b) created with `std::fs::DirBuilder::create` in
/// non-recursive mode -- which is exclusive: it fails with `AlreadyExists`
/// if anything (file, directory, or an attacker's pre-created entry) already
/// occupies that path, so a collision is *detected*, never silently trusted
/// -- and (c) `0700` (owner-only, so even a correct guess of the name can't
/// be written into after the fact) removes that surface while preserving
/// the exact same no-network guarantee: the directory is created empty and
/// never populated, so htslib's `<dir>/<md5>` lookup always misses and a
/// `-r`-less CRAM still fails LOCALLY.
///
/// An even simpler alternative -- `REF_PATH=""` (present but empty, as
/// opposed to unset) -- was considered and empirically ruled out: htslib's
/// `cram_populate_ref` (`cram/cram_io.c`, checked against this workspace's
/// pinned `hts-sys` 2.2.0) tests `!ref_path || *ref_path == '\0'`, treating
/// an empty string identically to an unset variable and substituting the
/// built-in EBI URL either way. A throwaway probe confirmed this at
/// runtime: with `REF_PATH=""` and `REF_CACHE` unset, decoding an
/// externally-referenced CRAM (no local reference, `-r`-less) produced
/// htslib's own `cram_populate_ref` info log --
/// `Populating local cache: .../hts-ref/%2s/%2s/%s` -- followed by an
/// actual outbound attempt to
/// `https://www.ebi.ac.uk/ena/cram/md5/<checksum>` (a libcurl connection
/// error, not a local "no such path" error). So `REF_PATH=""` reaches the
/// network fallback exactly like unset `REF_PATH` does, and cannot be used
/// here.
///
/// The returned directory is deliberately never deleted (kept for the full
/// process lifetime, per [`neutralize_cram_ref_path_network_fallback`]'s
/// contract) and its emptiness is never populated by this process -- an
/// empty, unused `0700` directory left behind in the OS temp namespace after
/// exit is a negligible, harmless byproduct, unlike the fixed-path design it
/// replaces.
///
/// # Errors
///
/// Returns an error if [`REF_PATH_DIR_CREATE_ATTEMPTS`] consecutive
/// timestamped names all collide with an existing path, or if directory
/// creation fails for any other reason (e.g. an unwritable temp namespace).
fn create_private_ref_path_dir() -> anyhow::Result<std::path::PathBuf> {
    use std::os::unix::fs::DirBuilderExt as _;
    use std::os::unix::fs::PermissionsExt as _;

    let base = std::env::temp_dir();
    let pid = std::process::id();
    for attempt in 0..REF_PATH_DIR_CREATE_ATTEMPTS {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos());
        let candidate = base.join(format!("unum-no-network-ref-path-{pid}-{nanos}-{attempt}"));
        match std::fs::DirBuilder::new().mode(0o700).create(&candidate) {
            Ok(()) => {
                // Belt-and-braces: force the mode explicitly rather than
                // relying solely on `DirBuilderExt::mode`, which is still
                // subject to the process umask (POSIX `mkdir(2)` semantics).
                // A `0o700` request can only be *narrowed* by a umask that
                // masks owner bits, which no realistic umask does, but
                // setting it again here removes any doubt.
                std::fs::set_permissions(&candidate, std::fs::Permissions::from_mode(0o700))
                    .with_context(|| {
                        format!("setting 0700 permissions on {}", candidate.display())
                    })?;
                return Ok(candidate);
            }
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(e) => {
                return Err(e).with_context(|| {
                    format!("creating private REF_PATH directory {}", candidate.display())
                });
            }
        }
    }
    anyhow::bail!(
        "failed to create a private REF_PATH directory under {} after \
         {REF_PATH_DIR_CREATE_ATTEMPTS} attempts (unexpected repeated name collisions)",
        base.display()
    );
}

/// Neutralizes htslib's default (network-capable) CRAM reference chain for
/// the lifetime of this process, by pointing the `REF_PATH` environment
/// variable at a fresh, private, `0700` directory (see
/// [`create_private_ref_path_dir`]) instead of leaving it unset (which would
/// fall back to htslib's built-in EBI URL) or empty (which, per
/// [`create_private_ref_path_dir`]'s doc comment, htslib treats the same as
/// unset). `REF_PATH` is a pure fallback: it is only ever consulted for a
/// CRAM reader that has NOT been given an explicit reference via
/// `set_reference` (i.e. no `-r` was threaded through), so this has no
/// effect on any decode that supplies `-r` (BAM, or CRAM with `-r` -> local
/// `.fai`-backed decode, unaffected). It IS the only backstop for CRAM
/// arriving on stdin, whose container format can't be sniffed from a
/// non-seekable pipe -- see `stages::extract::run_bam_no_alignment`'s stdin
/// arm.
///
/// # Errors
///
/// Returns an error if [`create_private_ref_path_dir`] cannot create the
/// directory `REF_PATH` is pointed at.
///
/// # Safety
///
/// Must run as the very first statement of `main`, before `Cli::parse()` and
/// before anything else reads or writes the process environment, and while
/// the process is still single-threaded (no other thread has been spawned
/// yet -- `rayon`'s thread pool is built lazily on first use, well after
/// this call). Both are edition-2024 `set_var` soundness preconditions:
/// nothing else can be racing this mutation or observing a torn read. The
/// directory creation itself is ordinary, safe I/O and runs before the
/// `unsafe` block; only the `set_var` call needs (and gets) the `unsafe`
/// scope.
#[allow(unsafe_code)]
fn neutralize_cram_ref_path_network_fallback() -> anyhow::Result<()> {
    let ref_path_dir = create_private_ref_path_dir()?;

    // SAFETY: called as the first statement of `main`, before `Cli::parse()`
    // and before any other code (this crate's or a dependency's) reads or
    // writes the environment, and before any additional thread exists --
    // satisfying edition 2024's `set_var` precondition that no concurrent
    // access can race this write.
    unsafe {
        std::env::set_var("REF_PATH", &ref_path_dir);
    }
    Ok(())
}

fn main() -> anyhow::Result<()> {
    neutralize_cram_ref_path_network_fallback()?;

    let cli = Cli::parse();
    match cli.command {
        Commands::Build(args) => stages::build::run(&args),
        Commands::Extract(args) => stages::extract::run(&args),
        Commands::Genotype(args) => stages::genotype::run(&args),
        Commands::Analyze(args) => stages::analyze::run(&args),
        Commands::Run(args) => stages::run(&args),
    }
}
