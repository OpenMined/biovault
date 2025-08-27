use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Error)]
pub enum Error {
    #[error("Project folder does not exist: {0}")]
    ProjectFolderMissing(String),

    #[error("Project configuration not found: {0}")]
    ProjectConfigMissing(String),

    #[error("Workflow file not found: {0}")]
    WorkflowMissing(String),

    #[error("Patient file does not exist: {0}")]
    PatientFileMissing(String),

    #[error("Patient not found: {0}")]
    PatientNotFound(String),

    #[error("Patient not found in patient file: {0}")]
    PatientNotFoundInFile(String),

    #[error("No patients specified")]
    NoPatientsSpecified,

    #[error("Templates not found. Please run 'bv init' first")]
    TemplatesNotFound,

    #[error("File not found: {file}: {details}")]
    FileNotFound { file: String, details: String },

    #[error("Checksum verification failed for {0}")]
    ChecksumFailed(String),

    #[error("HTTP request failed with status: {0}")]
    HttpRequestFailed(String),

    #[error("{0} patient(s) failed processing")]
    ProcessingFailed(usize),

    #[error(transparent)]
    Anyhow(#[from] anyhow::Error),

    #[error(transparent)]
    Io(#[from] std::io::Error),

    #[error(transparent)]
    Yaml(#[from] serde_yaml::Error),
}
