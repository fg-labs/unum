#![forbid(unsafe_code)]

// Scaffold placeholder: the CLI entry point is fallible once command dispatch lands,
// so the `Result` return type is intentional even though it's an unconditional `Ok`
// today. Remove this allow once real command handling is wired up.
#[allow(clippy::unnecessary_wraps)]
fn main() -> anyhow::Result<()> {
    Ok(())
}
