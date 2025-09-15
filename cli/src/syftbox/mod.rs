pub mod app;
pub mod endpoint;
pub mod rpc;
pub mod types;

pub use app::SyftBoxApp;
pub use endpoint::Endpoint;
pub use rpc::{check_requests, send_response};
pub use types::{RpcHeaders, RpcRequest, RpcResponse};

use anyhow::Result;
use std::path::Path;

/// Initialize a SyftBox app with the given name in the specified data directory
pub fn init_app(data_dir: &Path, email: &str, app_name: &str) -> Result<SyftBoxApp> {
    SyftBoxApp::new(data_dir, email, app_name)
}
