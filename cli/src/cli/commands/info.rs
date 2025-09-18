use crate::Result;
use std::collections::HashMap;
#[cfg(target_family = "unix")]
use std::os::unix::fs::MetadataExt;
use sysinfo::{Disks, System};

pub async fn execute() -> Result<()> {
    let mut sys = System::new_all();
    sys.refresh_all();

    println!("System Information");
    println!("==================");

    // OS
    println!(
        "OS: {} {}",
        System::name().unwrap_or_else(|| "Unknown".to_string()),
        System::os_version().unwrap_or_default()
    );

    // CPU Architecture
    println!(
        "CPU Arch: {}",
        System::cpu_arch().unwrap_or_else(|| "Unknown".to_string())
    );

    // CPU count
    let cpu_count = sys.cpus().len();
    println!("CPUs: {}", cpu_count);

    // RAM
    let total_memory = sys.total_memory();
    let total_memory_gb = total_memory as f64 / (1024.0 * 1024.0 * 1024.0);
    println!("RAM: {:.2} GB", total_memory_gb);

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
    println!("DISK FREE: {:.2} GB", free_space_gb);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn info_execute_runs() {
        // Just ensure the function executes without error; it queries local sysinfo only
        execute().await.unwrap();
    }
}
