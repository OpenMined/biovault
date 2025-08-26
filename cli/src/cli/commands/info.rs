use crate::Result;
use sysinfo::{Disks, System};

pub async fn execute() -> Result<()> {
    let mut sys = System::new_all();
    sys.refresh_all();
    
    println!("System Information");
    println!("==================");
    
    // OS
    println!("OS: {} {}", 
        System::name().unwrap_or_else(|| "Unknown".to_string()),
        System::os_version().unwrap_or_else(|| "".to_string())
    );
    
    // CPU Architecture  
    println!("CPU Arch: {}", 
        System::cpu_arch().unwrap_or_else(|| "Unknown".to_string())
    );
    
    // CPU count
    let cpu_count = sys.cpus().len();
    println!("CPUs: {}", cpu_count);
    
    // RAM
    let total_memory = sys.total_memory();
    let total_memory_gb = total_memory as f64 / (1024.0 * 1024.0 * 1024.0);
    println!("RAM: {:.2} GB", total_memory_gb);
    
    // Disk Free Space
    let disks = Disks::new_with_refreshed_list();
    let mut total_free_space: u64 = 0;
    for disk in disks.iter() {
        total_free_space += disk.available_space();
    }
    let free_space_gb = total_free_space as f64 / (1024.0 * 1024.0 * 1024.0);
    println!("DISK FREE: {:.2} GB", free_space_gb);
    
    Ok(())
}