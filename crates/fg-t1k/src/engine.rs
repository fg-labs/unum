//! Per-stage strangler router: each pipeline stage runs the C++ oracle or the Rust port.
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Engine {
    Cpp,
    Rust,
}

#[derive(Default)]
pub struct EngineOverrides(HashMap<String, Engine>);

impl EngineOverrides {
    pub fn parse(specs: &[String]) -> anyhow::Result<Self> {
        let mut m = HashMap::new();
        for spec in specs {
            let (stage, eng) =
                spec.split_once('=').ok_or_else(|| anyhow::anyhow!("bad --engine spec: {spec}"))?;
            let eng = match eng {
                "cpp" => Engine::Cpp,
                "rust" => Engine::Rust,
                other => anyhow::bail!("unknown engine: {other}"),
            };
            m.insert(stage.to_string(), eng);
        }
        Ok(Self(m))
    }
    pub fn engine_for(&self, stage: &str) -> Engine {
        self.0.get(stage).copied().unwrap_or(Engine::Cpp) // default: C++ oracle
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_valid_spec() {
        let overrides = EngineOverrides::parse(&["extract=rust".to_string()]).unwrap();
        assert_eq!(overrides.engine_for("extract"), Engine::Rust);
    }

    #[test]
    fn parse_bad_spec_missing_equals() {
        match EngineOverrides::parse(&["extractrust".to_string()]) {
            Err(err) => assert!(err.to_string().contains("bad --engine spec")),
            Ok(_) => panic!("expected an error for a spec without '='"),
        }
    }

    #[test]
    fn parse_unknown_engine_name() {
        match EngineOverrides::parse(&["extract=python".to_string()]) {
            Err(err) => assert!(err.to_string().contains("unknown engine")),
            Ok(_) => panic!("expected an error for an unknown engine name"),
        }
    }

    #[test]
    fn engine_for_defaults_to_cpp() {
        let overrides = EngineOverrides::default();
        assert_eq!(overrides.engine_for("genotype"), Engine::Cpp);
    }

    #[test]
    fn engine_for_returns_override() {
        let overrides = EngineOverrides::parse(&["genotype=rust".to_string()]).unwrap();
        assert_eq!(overrides.engine_for("genotype"), Engine::Rust);
        assert_eq!(overrides.engine_for("analyze"), Engine::Cpp);
    }
}
