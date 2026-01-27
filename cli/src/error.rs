use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Error)]
pub enum Error {
    #[error("Module folder does not exist: {0}")]
    ModuleFolderMissing(String),

    #[error("Module configuration not found: {0}")]
    ModuleConfigMissing(String),

    #[error("Workflow file not found: {0}")]
    WorkflowMissing(String),

    #[error("Participant not found: {0}")]
    ParticipantNotFound(String),

    #[error("No participants specified")]
    NoParticipantsSpecified,

    #[error("Templates not found. Please run 'bv init' first")]
    TemplatesNotFound,

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
