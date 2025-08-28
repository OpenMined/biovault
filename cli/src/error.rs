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

    #[error("Participant file does not exist: {0}")]
    ParticipantFileMissing(String),

    #[error("Participant not found: {0}")]
    ParticipantNotFound(String),

    #[error("Participant not found in participant file: {0}")]
    ParticipantNotFoundInFile(String),

    #[error("No participants specified")]
    NoParticipantsSpecified,

    #[error("Templates not found. Please run 'bv init' first")]
    TemplatesNotFound,

    #[error("File not found: {file}: {details}")]
    FileNotFound { file: String, details: String },

    #[error("Checksum verification failed for {0}")]
    ChecksumFailed(String),

    #[error("HTTP request failed with status: {0}")]
    HttpRequestFailed(String),

    #[error("{0} participant(s) failed processing")]
    ProcessingFailed(usize),

    #[error("BioVault not initialized. Please run 'bv init' first")]
    NotInitialized,

    #[error("SyftBox config file not found: {0}")]
    SyftBoxConfigMissing(String),

    #[error("Datasites directory not found: {0}")]
    DatasitesDirMissing(String),

    #[error(transparent)]
    Anyhow(#[from] anyhow::Error),

    #[error(transparent)]
    Io(#[from] std::io::Error),

    #[error(transparent)]
    Yaml(#[from] serde_yaml::Error),

    #[error(transparent)]
    Json(#[from] serde_json::Error),
}
