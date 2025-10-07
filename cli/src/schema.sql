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
    metadata TEXT,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (participant_id) REFERENCES participants(id) ON DELETE SET NULL
);

CREATE INDEX IF NOT EXISTS idx_files_participant_id ON files(participant_id);
CREATE INDEX IF NOT EXISTS idx_files_file_type ON files(file_type);
CREATE INDEX IF NOT EXISTS idx_files_hash ON files(file_hash);

-- NEW: Projects (Desktop-only but in shared DB)
CREATE TABLE IF NOT EXISTS projects (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT UNIQUE NOT NULL,
    author TEXT NOT NULL,
    workflow TEXT NOT NULL,
    template TEXT NOT NULL,
    project_path TEXT NOT NULL,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- NEW: Runs (Desktop-only)
CREATE TABLE IF NOT EXISTS runs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    project_id INTEGER NOT NULL,
    work_dir TEXT NOT NULL,
    participant_count INTEGER NOT NULL,
    status TEXT NOT NULL,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (project_id) REFERENCES projects(id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_runs_project_id ON runs(project_id);
CREATE INDEX IF NOT EXISTS idx_runs_status ON runs(status);

-- NEW: Run Participants (Desktop-only, but can join with files)
CREATE TABLE IF NOT EXISTS run_participants (
    run_id INTEGER NOT NULL,
    participant_id INTEGER NOT NULL,
    FOREIGN KEY (run_id) REFERENCES runs(id) ON DELETE CASCADE,
    FOREIGN KEY (participant_id) REFERENCES participants(id) ON DELETE CASCADE,
    PRIMARY KEY (run_id, participant_id)
);
