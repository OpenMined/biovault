use crate::config;
use crate::Result;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
#[cfg(target_family = "unix")]
use std::os::unix::fs::MetadataExt;
use sysinfo::{Disks, System};

#[derive(Debug, Serialize, Deserialize)]
pub struct SystemInfo {
    pub os: String,
    pub os_version: String,
    pub cpu_arch: String,
    pub cpu_count: usize,
    pub ram_gb: f64,
    pub disk_free_gb: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub syftbox: Option<SyftBoxInfo>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SyftBoxInfo {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub binary: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub version: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub data_dir: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub email: Option<String>,
}

pub async fn execute(json: bool) -> Result<()> {
    let mut sys = System::new_all();
    sys.refresh_all();

    // Collect OS info
    let os_name = System::name().unwrap_or_else(|| "Unknown".to_string());
    let os_version = System::os_version().unwrap_or_default();
    let cpu_arch = System::cpu_arch().unwrap_or_else(|| "Unknown".to_string());
    let cpu_count = sys.cpus().len();

    // RAM
    let total_memory = sys.total_memory();
    let total_memory_gb = total_memory as f64 / (1024.0 * 1024.0 * 1024.0);

    // Disk Free Space (deduplicated by filesystem device and excluding virtual FS)
    let disks = Disks::new_with_refreshed_list();

    // Filesystems to ignore when reporting free space
    let virtual_fs: [&str; 11] = [
        "tmpfs", "devtmpfs", "proc", "sysfs", "efivarfs", "cgroup", "cgroup2", "overlay",
        "squashfs", "ramfs", "zram",
    ];

    #[cfg(target_family = "unix")]
    let mut device_max_free: HashMap<u64, u64> = HashMap::new();
    #[cfg(not(target_family = "unix"))]
    let mut device_max_free: HashMap<String, u64> = HashMap::new();

    for disk in disks.iter() {
        // Skip virtual filesystems that shouldn't count toward usable disk
        let fs_type = disk.file_system().to_string_lossy().to_lowercase();
        if virtual_fs.iter().any(|v| fs_type == *v) {
            continue;
        }

        let avail = disk.available_space();

        #[cfg(target_family = "unix")]
        {
            let mount_point = disk.mount_point();
            if let Ok(meta) = std::fs::metadata(mount_point) {
                let dev_id = meta.dev();
                device_max_free
                    .entry(dev_id)
                    .and_modify(|v| {
                        if avail > *v {
                            *v = avail;
                        }
                    })
                    .or_insert(avail);
            }
        }

        #[cfg(not(target_family = "unix"))]
        {
            // On non-Unix, fall back to grouping by disk name
            let key = disk.name().to_string_lossy().to_string();
            device_max_free
                .entry(key)
                .and_modify(|v| {
                    if avail > *v {
                        *v = avail;
                    }
                })
                .or_insert(avail);
        }
    }

    let total_free_space: u128 = device_max_free.values().map(|v| *v as u128).sum();
    let free_space_gb = total_free_space as f64 / (1024.0 * 1024.0 * 1024.0);

    // Collect SyftBox environment info if available
    let syftbox_info = if config::is_syftbox_env() {
        Some(SyftBoxInfo {
            binary: config::get_syftbox_binary(),
            version: config::get_syftbox_version(),
            data_dir: std::env::var("SYFTBOX_DATA_DIR").ok(),
            email: std::env::var("SYFTBOX_EMAIL").ok(),
        })
    } else {
        None
    };

    let info = SystemInfo {
        os: os_name.clone(),
        os_version: os_version.clone(),
        cpu_arch: cpu_arch.clone(),
        cpu_count,
        ram_gb: total_memory_gb,
        disk_free_gb: free_space_gb,
        syftbox: syftbox_info.clone(),
    };

    if json {
        println!("{}", serde_json::to_string_pretty(&info)?);
    } else {
        println!("System Information");
        println!("==================");
        println!("OS: {} {}", os_name, os_version);
        println!("CPU Arch: {}", cpu_arch);
        println!("CPUs: {}", cpu_count);
        println!("RAM: {:.2} GB", total_memory_gb);
        println!("DISK FREE: {:.2} GB", free_space_gb);

        if let Some(syftbox) = syftbox_info {
            println!("\nSyftBox Environment");
            println!("===================");

            if let Some(binary) = syftbox.binary {
                println!("SyftBox Binary: {}", binary);
            }

            if let Some(version) = syftbox.version {
                println!("SyftBox Version: {}", version);
            }

            if let Some(data_dir) = syftbox.data_dir {
                println!("SyftBox Data Dir: {}", data_dir);
            }

            if let Some(email) = syftbox.email {
                println!("SyftBox Email: {}", email);
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn info_execute_runs() {
        // Just ensure the function executes without error; it queries local sysinfo only
        execute(false).await.unwrap();
    }
}
