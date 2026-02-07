-- BioVault Unified Schema
-- Version: 2.0.0
-- Migration from: messages.db (v1) → biovault.db (v2)

-- Schema versioning
CREATE TABLE IF NOT EXISTS schema_version (
    version TEXT PRIMARY KEY,
    applied_at DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- Existing: Messages (from messages.db)
CREATE TABLE IF NOT EXISTS messages (
    id TEXT PRIMARY KEY,
    thread_id TEXT,
    parent_id TEXT,
    from_address TEXT NOT NULL,
    to_address TEXT NOT NULL,
    subject TEXT,
    body TEXT NOT NULL,
    message_type TEXT DEFAULT 'text',
    metadata TEXT,
    status TEXT NOT NULL,
    sync_status TEXT NOT NULL,
    created_at TEXT NOT NULL,
    sent_at TEXT,
    received_at TEXT,
    read_at TEXT,
    rpc_request_id TEXT,
    rpc_ack_status INTEGER,
    rpc_ack_at TEXT
);

CREATE INDEX IF NOT EXISTS idx_thread_id ON messages(thread_id);
CREATE INDEX IF NOT EXISTS idx_from_address ON messages(from_address);
CREATE INDEX IF NOT EXISTS idx_to_address ON messages(to_address);
CREATE INDEX IF NOT EXISTS idx_status ON messages(status);
CREATE INDEX IF NOT EXISTS idx_created_at ON messages(created_at);
CREATE INDEX IF NOT EXISTS idx_message_type ON messages(message_type);

-- NEW: Participants (shared CLI/Desktop)
CREATE TABLE IF NOT EXISTS participants (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    participant_id TEXT UNIQUE NOT NULL,
    inferred_sex TEXT,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_participant_id ON participants(participant_id);

-- NEW: Files (shared CLI/Desktop)
-- participant_id is NULLABLE - files can exist without participant assignment
CREATE TABLE IF NOT EXISTS files (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    participant_id INTEGER,
    file_path TEXT UNIQUE NOT NULL,
    file_hash TEXT NOT NULL,
    file_type TEXT,
    file_size INTEGER,
    data_type TEXT DEFAULT 'Unknown',
    metadata TEXT,
    status TEXT DEFAULT 'complete',
    processing_error TEXT,
    queue_added_at DATETIME,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (participant_id) REFERENCES participants(id) ON DELETE SET NULL
);

CREATE INDEX IF NOT EXISTS idx_files_participant_id ON files(participant_id);
CREATE INDEX IF NOT EXISTS idx_files_file_type ON files(file_type);
CREATE INDEX IF NOT EXISTS idx_files_hash ON files(file_hash);
-- idx_files_data_type and idx_files_status are created by migration after adding columns

-- NEW: Datasets
CREATE TABLE IF NOT EXISTS datasets (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL UNIQUE,
    version TEXT NOT NULL DEFAULT '1.0.0',
    author TEXT NOT NULL,
    description TEXT,
    schema TEXT NOT NULL,
    public_url TEXT,
    private_url TEXT,
    http_relay_servers TEXT, -- JSON array
    extra TEXT,             -- JSON for future fields
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS dataset_assets (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    dataset_id INTEGER NOT NULL,
    asset_key TEXT NOT NULL,
    asset_uuid TEXT NOT NULL,
    kind TEXT NOT NULL,
    url TEXT NOT NULL,
    private_ref TEXT,
    mock_ref TEXT,
    extra TEXT,
    private_file_id INTEGER,
    mock_file_id INTEGER,
    private_path TEXT,
    mock_path TEXT,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (dataset_id) REFERENCES datasets(id) ON DELETE CASCADE,
    UNIQUE(dataset_id, asset_key)
);

CREATE INDEX IF NOT EXISTS idx_dataset_assets_dataset_id ON dataset_assets(dataset_id);

-- NEW: Genotype Metadata (for files where data_type = 'Genotype')
CREATE TABLE IF NOT EXISTS genotype_metadata (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    file_id INTEGER UNIQUE NOT NULL,
    source TEXT,
    grch_version TEXT,
    row_count INTEGER,
    chromosome_count INTEGER,
    inferred_sex TEXT,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (file_id) REFERENCES files(id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_genotype_file_id ON genotype_metadata(file_id);

-- NEW: Variant Metadata (for files where data_type = 'Variants')
CREATE TABLE IF NOT EXISTS variant_metadata (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    file_id INTEGER UNIQUE NOT NULL,
    source TEXT,
    grch_version TEXT,
    format TEXT,
    reference_file_id INTEGER,
    reference_index_file_id INTEGER,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (file_id) REFERENCES files(id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_variant_file_id ON variant_metadata(file_id);

-- NEW: Aligned Metadata (for files where data_type = 'Aligned')
CREATE TABLE IF NOT EXISTS aligned_metadata (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    file_id INTEGER UNIQUE NOT NULL,
    source TEXT,
    grch_version TEXT,
    format TEXT,
    reference_file_id INTEGER,
    reference_index_file_id INTEGER,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (file_id) REFERENCES files(id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_aligned_file_id ON aligned_metadata(file_id);

-- NEW: Reference Metadata (for files where data_type = 'Reference')
CREATE TABLE IF NOT EXISTS reference_metadata (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    file_id INTEGER UNIQUE NOT NULL,
    source TEXT,
    grch_version TEXT,
    format TEXT,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (file_id) REFERENCES files(id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_reference_file_id ON reference_metadata(file_id);

-- NEW: Database Metadata (for files where data_type = 'VariantDatabase')
-- Used for variant annotation databases like ClinVar, dbSNP, gnomAD, etc.
CREATE TABLE IF NOT EXISTS database_metadata (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    file_id INTEGER UNIQUE NOT NULL,
    source TEXT,                          -- Source name (e.g., "ClinVar", "dbSNP", "gnomAD")
    grch_version TEXT,                    -- Reference genome version (e.g., "GRCh38", "GRCh37")
    format TEXT,                          -- File format (e.g., "vcf.gz")
    version TEXT,                         -- Database version/release date
    index_file_id INTEGER,                -- Associated index file (.tbi)
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (file_id) REFERENCES files(id) ON DELETE CASCADE,
    FOREIGN KEY (index_file_id) REFERENCES files(id) ON DELETE SET NULL
);

CREATE INDEX IF NOT EXISTS idx_database_file_id ON database_metadata(file_id);
CREATE INDEX IF NOT EXISTS idx_database_source ON database_metadata(source);

-- Add columns to existing files table if they don't exist (migration)
-- SQLite doesn't have IF NOT EXISTS for ALTER TABLE, so we check first
-- This is handled by separate migration code if needed

-- NEW: Modules (Desktop-only but in shared DB)
CREATE TABLE IF NOT EXISTS modules (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL,
    version TEXT NOT NULL DEFAULT '1.0.0',
    author TEXT NOT NULL,
    workflow TEXT NOT NULL,
    template TEXT NOT NULL,
    module_path TEXT UNIQUE NOT NULL,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(name, version)
);

-- NEW: Development Environments (virtualenvs for Jupyter, etc.)
CREATE TABLE IF NOT EXISTS dev_envs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    module_path TEXT UNIQUE NOT NULL,
    python_version TEXT NOT NULL,
    env_type TEXT DEFAULT 'jupyter',
    jupyter_installed INTEGER DEFAULT 0,
    jupyter_port INTEGER,
    jupyter_pid INTEGER,
    jupyter_url TEXT,
    jupyter_token TEXT,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    last_used_at DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_dev_envs_module_path ON dev_envs(module_path);
CREATE INDEX IF NOT EXISTS idx_dev_envs_env_type ON dev_envs(env_type);

-- NEW: Flows (shared CLI/Desktop)
CREATE TABLE IF NOT EXISTS flows (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL UNIQUE,
    flow_path TEXT NOT NULL,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- NEW: Unified Flow Runs Table (handles both flow and standalone module runs)
CREATE TABLE IF NOT EXISTS flow_runs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    flow_id INTEGER,                    -- NULL for standalone module runs
    module_id INTEGER,                  -- NULL for flow runs (references modules table)
    status TEXT NOT NULL,
    work_dir TEXT NOT NULL,
    results_dir TEXT,
    participant_count INTEGER,          -- Only for standalone module runs
    metadata TEXT,                      -- JSON: { "input_overrides": {...}, "parameter_overrides": {...} }
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    completed_at DATETIME,
    FOREIGN KEY (flow_id) REFERENCES flows(id) ON DELETE CASCADE,
    FOREIGN KEY (module_id) REFERENCES modules(id) ON DELETE SET NULL,
    CHECK ((flow_id IS NULL) != (module_id IS NULL))  -- Exactly one must be set
);

CREATE INDEX IF NOT EXISTS idx_flow_runs_flow_id ON flow_runs(flow_id);
CREATE INDEX IF NOT EXISTS idx_flow_runs_module_id ON flow_runs(module_id);
CREATE INDEX IF NOT EXISTS idx_flow_runs_status ON flow_runs(status);

-- NEW: Run Participants (for standalone module runs only)
CREATE TABLE IF NOT EXISTS flow_run_participants (
    run_id INTEGER NOT NULL,
    participant_id INTEGER NOT NULL,
    FOREIGN KEY (run_id) REFERENCES flow_runs(id) ON DELETE CASCADE,
    FOREIGN KEY (participant_id) REFERENCES participants(id) ON DELETE CASCADE,
    PRIMARY KEY (run_id, participant_id)
);

-- NEW: Run Configurations (saved flow input configurations)
CREATE TABLE IF NOT EXISTS flow_run_configs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    flow_id INTEGER NOT NULL,
    name TEXT NOT NULL,
    config_data TEXT NOT NULL,  -- JSON: { "inputs": {...}, "parameters": {...} }
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (flow_id) REFERENCES flows(id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_flow_run_configs_flow_id ON flow_run_configs(flow_id);

-- NEW: Sessions (for collaborative Jupyter environments with beaver)
-- Sessions allow data scientists to create collaborative Jupyter environments
-- with optional peer invitations for data sharing via beaver
CREATE TABLE IF NOT EXISTS sessions (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    session_id TEXT NOT NULL UNIQUE,     -- 12-char hex hash, matches beaver session format
    name TEXT NOT NULL,                  -- User-friendly name
    description TEXT,
    session_path TEXT UNIQUE NOT NULL,   -- Path to session folder in BioVault Home
    owner TEXT NOT NULL,                 -- Owner email
    peer TEXT,                           -- Invited peer email (optional)
    role TEXT DEFAULT 'owner',           -- 'owner' or 'collaborator'
    status TEXT DEFAULT 'active',        -- 'pending', 'active', 'closed'
    jupyter_port INTEGER,
    jupyter_pid INTEGER,
    jupyter_url TEXT,
    jupyter_token TEXT,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_sessions_session_id ON sessions(session_id);
CREATE INDEX IF NOT EXISTS idx_sessions_name ON sessions(name);
CREATE INDEX IF NOT EXISTS idx_sessions_peer ON sessions(peer);
CREATE INDEX IF NOT EXISTS idx_sessions_status ON sessions(status);

-- NEW: Session Messages (chat within a session)
CREATE TABLE IF NOT EXISTS session_messages (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    session_id INTEGER NOT NULL,
    sender TEXT NOT NULL,
    body TEXT NOT NULL,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (session_id) REFERENCES sessions(id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_session_messages_session_id ON session_messages(session_id);

-- NEW: Session Datasets (datasets associated with a session)
-- Links sessions to datasets for collaborative analysis
CREATE TABLE IF NOT EXISTS session_datasets (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    session_id TEXT NOT NULL,             -- References sessions.session_id
    dataset_public_url TEXT NOT NULL,     -- syft:// URL to the dataset manifest
    dataset_owner TEXT NOT NULL,          -- Owner email of the dataset
    dataset_name TEXT NOT NULL,           -- Name of the dataset
    role TEXT DEFAULT 'shared',           -- 'shared' (using shared data) or 'yours' (your own data)
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (session_id) REFERENCES sessions(session_id) ON DELETE CASCADE,
    UNIQUE(session_id, dataset_public_url)
);

CREATE INDEX IF NOT EXISTS idx_session_datasets_session_id ON session_datasets(session_id);
CREATE INDEX IF NOT EXISTS idx_session_datasets_owner ON session_datasets(dataset_owner);
