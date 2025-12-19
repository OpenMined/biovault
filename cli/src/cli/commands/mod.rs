pub mod biobank;
pub mod check;
pub mod config_cmd;
pub mod contacts;
pub mod daemon;
pub mod datasets;
pub mod fastq;
pub mod files;
pub mod hard_reset;
pub mod inbox;
pub mod info;
pub mod init;
pub mod jupyter;
pub mod key;
pub mod messages;
pub mod network;
pub mod participant;
pub mod participants;
pub mod pipeline;
pub mod project;
pub mod project_management;
pub mod python;
pub mod run;
pub mod run_dynamic;
pub mod sample_data;
pub mod samplesheet;
pub mod sessions;
pub mod setup;
pub mod sql;
pub mod submit;
pub mod syc;
pub mod syftbox;
pub mod syftboxd;
pub mod update;

/// Configure a child process command to avoid flashing a console window when the CLI library
/// is used from the Desktop GUI on Windows.
#[cfg(target_os = "windows")]
pub(crate) fn configure_child_process(cmd: &mut std::process::Command) {
    let hide = std::env::var_os("BIOVAULT_HIDE_CONSOLE").is_some()
        || std::env::var_os("BIOVAULT_DESKTOP").is_some();

    if !hide {
        return;
    }

    use std::os::windows::process::CommandExt;
    const CREATE_NO_WINDOW: u32 = 0x08000000;
    cmd.creation_flags(CREATE_NO_WINDOW);
}

#[cfg(not(target_os = "windows"))]
pub(crate) fn configure_child_process(_cmd: &mut std::process::Command) {}
