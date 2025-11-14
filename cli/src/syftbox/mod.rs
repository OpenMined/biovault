pub mod app;
pub mod control;
pub mod endpoint;
pub mod rpc;
pub mod storage;
pub mod types;

pub use app::SyftBoxApp;
pub use control::{
    detect_mode, is_syftbox_running, start_syftbox, state as syftbox_state, stop_syftbox,
    SyftBoxMode, SyftBoxState,
};
pub use endpoint::Endpoint;
pub use rpc::{check_requests, send_response};
pub use storage::SyftBoxStorage;
pub use types::{RpcHeaders, RpcRequest, RpcResponse};

use anyhow::Result;
use std::path::Path;

/// Initialize a SyftBox app with the given name in the specified data directory
pub fn init_app(data_dir: &Path, email: &str, app_name: &str) -> Result<SyftBoxApp> {
    SyftBoxApp::new(data_dir, email, app_name)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn init_app_creates_structure() {
        let td = TempDir::new().unwrap();
        let app = init_app(td.path(), "u@example.com", "bv").unwrap();
        assert!(app.app_data_dir.exists());
        assert!(app.rpc_dir.exists());
    }
}
