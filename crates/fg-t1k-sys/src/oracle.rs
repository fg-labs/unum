use std::path::PathBuf;

#[derive(Clone, Copy, Debug)]
pub enum OracleStage {
    Genotyper,
    Analyzer,
    FastqExtractor,
    BamExtractor,
}

impl OracleStage {
    fn binary_name(self) -> &'static str {
        match self {
            OracleStage::Genotyper => "genotyper",
            OracleStage::Analyzer => "analyzer",
            OracleStage::FastqExtractor => "fastq-extractor",
            OracleStage::BamExtractor => "bam-extractor",
        }
    }
}

/// Absolute path to the compiled T1K oracle binary for `stage`.
pub fn binary_path(stage: OracleStage) -> PathBuf {
    PathBuf::from(env!("FG_T1K_ORACLE_DIR")).join(stage.binary_name())
}
