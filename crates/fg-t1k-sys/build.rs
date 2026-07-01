fn main() {
    // No-op unless t1k-sys is enabled (keeps default builds C++-free).
    if std::env::var("CARGO_FEATURE_T1K_SYS").is_err() {
        return;
    }
    // Populated in Task 0.3 / 0.4.
}
