use anyhow::{anyhow, Result};
use blake3::Hasher;
use indicatif::{ProgressBar, ProgressStyle};
use std::time::Instant;
use regex::Regex;
use std::collections::{BTreeMap, HashSet, HashMap};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use chrono::{DateTime, FixedOffset};

#[derive(Debug, Clone)]
struct FastqFile {
    path: PathBuf,
    size: u64,
    compression: CompressionType,
    sequence_number: Option<u32>,
    base_name: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum CompressionType {
    None,
    Gzip,
    Bzip2,
    Xz,
    Zip,
    Zstd,
}

impl CompressionType {
    fn from_path(path: &Path) -> Self {
        let path_str = path.to_string_lossy().to_lowercase();
        
        if path_str.ends_with(".gz") {
            CompressionType::Gzip
        } else if path_str.ends_with(".bz2") {
            CompressionType::Bzip2
        } else if path_str.ends_with(".xz") {
            CompressionType::Xz
        } else if path_str.ends_with(".zip") {
            CompressionType::Zip
        } else if path_str.ends_with(".zst") {
            CompressionType::Zstd
        } else {
            CompressionType::None
        }
    }
}

#[derive(Debug)]
struct FastqStats {
    num_seqs: u64,
    total_len: u64,
    min_len: u64,
    avg_len: f64,
    max_len: u64,
    avg_qual: f64,
    gc_content: f64,
}

#[derive(Debug, Clone)]
struct NanoporeMetadata {
    runid: String,
    flow_cell_id: String,
    start_time: DateTime<FixedOffset>,
    protocol_group_id: String,
    sample_id: String,
    basecall_model: String,
}

pub async fn combine(
    input_folder: String, 
    output_file: String,
    validate: Option<bool>,
) -> Result<()> {
    let input_path = Path::new(&input_folder);
    
    // Auto-generate output filename if output_file is a directory
    let output_path = if output_file.ends_with('/') || (Path::new(&output_file).exists() && Path::new(&output_file).is_dir()) {
        generate_output_filename(input_path, &output_file)?
    } else {
        PathBuf::from(&output_file)
    };
    
    if !input_path.exists() {
        return Err(anyhow!("Input folder does not exist: {}", input_folder));
    }
    
    if !input_path.is_dir() {
        return Err(anyhow!("Input path is not a directory: {}", input_folder));
    }
    
    // Create output directory if it doesn't exist
    if let Some(parent) = output_path.parent() {
        fs::create_dir_all(parent)?;
    }
    
    println!("Scanning for FASTQ files in: {}", input_folder);
    
    // Find all FASTQ files
    let mut fastq_files = find_fastq_files(input_path)?;
    
    if fastq_files.is_empty() {
        return Err(anyhow!("No FASTQ files found in {}", input_folder));
    }
    
    // Check for and remove duplicate files
    let duplicates = check_for_duplicates(&fastq_files);
    if !duplicates.is_empty() {
        println!("\n⚠️  Duplicate files detected (will be skipped):");
        for dup in &duplicates {
            println!("  - {}", dup.display());
        }
        fastq_files.retain(|f| !duplicates.contains(&f.path));
    }
    
    // Group files by base name
    let grouped_files = group_files_by_base(&fastq_files);
    
    // Display file groups
    for (base_name, files) in &grouped_files {
        println!("\n📁 File group: {}", base_name);
        
        // Check for missing sequence numbers
        let sequence_numbers: Vec<u32> = files.iter()
            .filter_map(|f| f.sequence_number)
            .collect();
        
        if !sequence_numbers.is_empty() {
            let min_seq = *sequence_numbers.iter().min().unwrap();
            let max_seq = *sequence_numbers.iter().max().unwrap();
            
            let mut missing = Vec::new();
            for i in min_seq..=max_seq {
                if !sequence_numbers.contains(&i) {
                    missing.push(i);
                }
            }
            
            println!("  Files: {} (sequences {}-{})", files.len(), min_seq, max_seq);
            
            if !missing.is_empty() {
                println!("  ⚠️  Missing sequences: {:?}", missing);
            }
        } else {
            println!("  Files: {} (no sequence numbers detected)", files.len());
        }
        
        // Calculate total size
        let total_size: u64 = files.iter().map(|f| f.size).sum();
        println!("  Total size: {}", format_size(total_size));
        
        // Check compression types
        let compression_types: HashSet<CompressionType> = files.iter()
            .map(|f| f.compression.clone())
            .collect();
        
        if compression_types.len() > 1 {
            println!("  ⚠️  Multiple compression types detected:");
            let mut compression_files: HashMap<CompressionType, Vec<String>> = HashMap::new();
            for file in files {
                let file_name = file.path.file_name().unwrap().to_string_lossy().to_string();
                compression_files.entry(file.compression.clone())
                    .or_insert_with(Vec::new)
                    .push(file_name);
            }
            for (ct, file_list) in &compression_files {
                println!("    - {:?}:", ct);
                for (i, file) in file_list.iter().enumerate() {
                    if i < 3 {
                        println!("        {}", file);
                    } else if i == 3 {
                        println!("        ... and {} more", file_list.len() - 3);
                        break;
                    }
                }
            }
            return Err(anyhow!("Cannot combine files with different compression types"));
        }
    }
    
    // Select which group to combine (for now, we'll combine all files)
    let all_files: Vec<FastqFile> = grouped_files.into_iter()
        .flat_map(|(_, files)| files)
        .collect();
    
    // Parse Nanopore metadata from first file of each group
    println!("\n🔬 Checking for Oxford Nanopore metadata...");
    let metadata_map = parse_nanopore_metadata(&all_files)?;
    
    if !metadata_map.is_empty() {
        display_nanopore_summary(&metadata_map)?;
        
        // Check for multiple sample IDs
        let sample_ids: HashSet<String> = metadata_map.values()
            .filter_map(|m| m.as_ref())
            .map(|m| m.sample_id.clone())
            .collect();
        
        if sample_ids.len() > 1 {
            println!("\n⚠️  Warning: Multiple sample IDs detected:");
            for id in &sample_ids {
                println!("    - {}", id);
            }
            print!("❓ Do you want to continue combining these files? (y/n): ");
            io::stdout().flush()?;
            
            let mut input = String::new();
            io::stdin().read_line(&mut input)?;
            if input.trim().to_lowercase() != "y" {
                return Err(anyhow!("Aborted: Multiple sample IDs detected"));
            }
        }
    }
    
    println!("\n📊 Summary:");
    println!("  Total files: {}", all_files.len());
    let total_size: u64 = all_files.iter().map(|f| f.size).sum();
    println!("  Total size: {}", format_size(total_size));
    
    // Check if user wants to validate
    let should_validate = if let Some(v) = validate {
        v
    } else {
        print!("\n❓ Do you want to validate all files before combining? (y/n): ");
        io::stdout().flush()?;
        
        let mut input = String::new();
        io::stdin().read_line(&mut input)?;
        input.trim().to_lowercase() == "y"
    };
    
    let mut stats_map = BTreeMap::new();
    
    if should_validate {
        // Check if seqkit is available
        if !check_seqkit_available() {
            return Err(anyhow!("seqkit is not installed. Please install it first.\nInstallation: https://bioinf.shenwei.me/seqkit/download/"));
        }
        
        println!("\n🔍 Validating FASTQ files...");
        
        let pb = ProgressBar::new(all_files.len() as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{bar:40.cyan/blue}] {pos}/{len} {msg}")
                .unwrap()
                .progress_chars("#>-")
        );
        
        let mut validation_errors = Vec::new();
        
        for (idx, file) in all_files.iter().enumerate() {
            pb.set_position(idx as u64);
            pb.set_message(format!("Validating {}", file.path.file_name().unwrap().to_string_lossy()));
            
            match validate_fastq_file(&file.path) {
                Ok(stats) => {
                    stats_map.insert(file.path.clone(), stats);
                }
                Err(e) => {
                    validation_errors.push((file.path.clone(), e));
                }
            }
        }
        
        pb.finish_with_message("Validation complete");
        
        if !validation_errors.is_empty() {
            println!("\n❌ Validation errors found:");
            for (path, error) in &validation_errors {
                println!("  {}: {}", path.display(), error);
            }
            return Err(anyhow!("Validation failed for {} files", validation_errors.len()));
        }
        
        println!("✅ All files validated successfully!");
        
        // Save validation stats to file
        let stats_file = output_path.with_extension("pre_combine_stats.txt");
        let file_hashes = generate_file_hashes(&all_files)?;
        save_stats_to_file(&stats_map, &stats_file, Some(&file_hashes))?;
        println!("📄 Validation stats saved to: {}", stats_file.display());
    }
    
    // Combine files
    println!("\n🔗 Combining FASTQ files...");
    combine_fastq_files(&all_files, &output_path)?;
    
    // Generate Blake3 hash
    println!("\n🔐 Generating Blake3 hash...");
    let hash = generate_blake3_hash(&output_path)?;
    let hash_file = output_path.with_extension("blake3");
    fs::write(&hash_file, hash)?;
    println!("📄 Blake3 hash saved to: {}", hash_file.display());
    
    // Validate combined file if validation was requested
    if should_validate {
        println!("\n🔍 Validating combined file...");
        
        let combined_stats = validate_fastq_file(&output_path)?;
        
        // Compare stats
        let total_expected_seqs: u64 = stats_map.values().map(|s| s.num_seqs).sum();
        
        if combined_stats.num_seqs != total_expected_seqs {
            println!("⚠️  Warning: Sequence count mismatch!");
            println!("  Expected: {} sequences", total_expected_seqs);
            println!("  Got: {} sequences", combined_stats.num_seqs);
        } else {
            println!("✅ Combined file validation successful!");
            println!("  Total sequences: {}", combined_stats.num_seqs);
            println!("  Average length: {:.2}", combined_stats.avg_len);
            println!("  GC content: {:.2}%", combined_stats.gc_content);
        }
        
        // Save final stats
        let final_stats_file = output_path.with_extension("stats.txt");
        let mut final_stats = BTreeMap::new();
        final_stats.insert(output_path.to_path_buf(), combined_stats);
        let combined_hash = generate_blake3_hash(&output_path)?;
        let mut final_hashes = BTreeMap::new();
        final_hashes.insert(output_path.to_path_buf(), combined_hash.clone());
        save_stats_to_file(&final_stats, &final_stats_file, Some(&final_hashes))?;
        println!("📄 Final stats saved to: {}", final_stats_file.display());
    }
    
    println!("\n✨ FASTQ files successfully combined to: {}", output_file);
    
    Ok(())
}

fn find_fastq_files(dir: &Path) -> Result<Vec<FastqFile>> {
    let mut files = Vec::new();
    
    // Define patterns for FASTQ files
    let patterns = vec![
        r"\.fastq$",
        r"\.fq$",
        r"\.fastq\.gz$",
        r"\.fq\.gz$",
        r"\.fastq\.bz2$",
        r"\.fq\.bz2$",
        r"\.fastq\.xz$",
        r"\.fq\.xz$",
        r"\.fastq\.zip$",
        r"\.fastq\.zst$",
    ];
    
    let combined_pattern = format!("({})", patterns.join("|"));
    let re = Regex::new(&combined_pattern)?;
    
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        
        if path.is_file() {
            let file_name = path.file_name().unwrap().to_string_lossy();
            
            if re.is_match(&file_name.to_lowercase()) {
                let metadata = entry.metadata()?;
                let compression = CompressionType::from_path(&path);
                
                // Extract base name and sequence number
                let (base_name, sequence_number) = extract_base_and_sequence(&file_name);
                
                files.push(FastqFile {
                    path,
                    size: metadata.len(),
                    compression,
                    sequence_number,
                    base_name,
                });
            }
        }
    }
    
    // Sort files by base name and sequence number
    files.sort_by(|a, b| {
        match a.base_name.cmp(&b.base_name) {
            std::cmp::Ordering::Equal => {
                a.sequence_number.cmp(&b.sequence_number)
            }
            other => other,
        }
    });
    
    Ok(files)
}

fn extract_base_and_sequence(file_name: &str) -> (String, Option<u32>) {
    // Try to extract sequence number from patterns like:
    // - file_001.fastq.gz
    // - file_1.fastq.gz
    // - file_0001.fastq.gz
    
    let patterns = vec![
        // Pattern with underscore separator
        (Regex::new(r"^(.+?)_(\d+)\.").unwrap(), 2),
        // Pattern with dash separator
        (Regex::new(r"^(.+?)-(\d+)\.").unwrap(), 2),
        // Pattern with dot separator before number
        (Regex::new(r"^(.+?)\.(\d+)\.").unwrap(), 2),
    ];
    
    for (re, seq_group) in patterns {
        if let Some(caps) = re.captures(file_name) {
            let base_name = caps.get(1).unwrap().as_str().to_string();
            let seq_str = caps.get(seq_group).unwrap().as_str();
            
            if let Ok(seq_num) = seq_str.parse::<u32>() {
                return (base_name, Some(seq_num));
            }
        }
    }
    
    // If no pattern matches, use the whole name without extension as base
    let base_name = file_name.split('.').next().unwrap_or(file_name).to_string();
    (base_name, None)
}

fn group_files_by_base(files: &[FastqFile]) -> BTreeMap<String, Vec<FastqFile>> {
    let mut groups = BTreeMap::new();
    
    for file in files {
        groups.entry(file.base_name.clone())
            .or_insert_with(Vec::new)
            .push(file.clone());
    }
    
    groups
}

fn check_seqkit_available() -> bool {
    Command::new("seqkit")
        .arg("version")
        .output()
        .map(|output| output.status.success())
        .unwrap_or(false)
}

fn validate_fastq_file(path: &Path) -> Result<FastqStats> {
    let output = Command::new("seqkit")
        .arg("stats")
        .arg("-T")  // Tab-separated output
        .arg("-a")  // All statistics
        .arg(path)
        .output()?;
    
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(anyhow!("seqkit validation failed: {}", stderr));
    }
    
    let stdout = String::from_utf8_lossy(&output.stdout);
    parse_seqkit_stats(&stdout)
}

fn parse_seqkit_stats(output: &str) -> Result<FastqStats> {
    // Parse seqkit stats output
    // Format is typically:
    // file  format  type  num_seqs  sum_len  min_len  avg_len  max_len  Q1  Q2  Q3  sum_gap  N50  Q20(%)  Q30(%)  GC(%)
    
    let lines: Vec<&str> = output.lines().collect();
    if lines.len() < 2 {
        return Err(anyhow!("Invalid seqkit stats output"));
    }
    
    // Skip header and parse data line
    let data_line = lines[1];
    let fields: Vec<&str> = data_line.split('\t').collect();
    
    if fields.len() < 15 {
        return Err(anyhow!("Insufficient fields in seqkit stats output"));
    }
    
    Ok(FastqStats {
        num_seqs: fields[3].replace(',', "").parse().unwrap_or(0),
        total_len: fields[4].replace(',', "").parse().unwrap_or(0),
        min_len: fields[5].replace(',', "").parse().unwrap_or(0),
        avg_len: fields[6].replace(',', "").parse().unwrap_or(0.0),
        max_len: fields[7].replace(',', "").parse().unwrap_or(0),
        avg_qual: 0.0, // Not directly provided by seqkit stats
        gc_content: fields[14].parse().unwrap_or(0.0),
    })
}

fn save_stats_to_file(stats: &BTreeMap<PathBuf, FastqStats>, path: &Path, hashes: Option<&BTreeMap<PathBuf, String>>) -> Result<()> {
    let mut file = File::create(path)?;
    
    writeln!(file, "FASTQ Validation Statistics")?;
    writeln!(file, "{}", "=".repeat(50))?;
    writeln!(file)?;
    
    let mut total_seqs = 0u64;
    let mut total_len = 0u64;
    
    for (file_path, stats) in stats {
        writeln!(file, "File: {}", file_path.display())?;
        writeln!(file, "  Sequences: {}", stats.num_seqs)?;
        writeln!(file, "  Total length: {}", stats.total_len)?;
        writeln!(file, "  Min length: {}", stats.min_len)?;
        writeln!(file, "  Avg length: {:.2}", stats.avg_len)?;
        writeln!(file, "  Max length: {}", stats.max_len)?;
        writeln!(file, "  GC content: {:.2}%", stats.gc_content)?;
        if let Some(h) = hashes {
            if let Some(hash) = h.get(file_path) {
                writeln!(file, "  Blake3 hash: {}", hash)?;
            }
        }
        writeln!(file)?;
        
        total_seqs += stats.num_seqs;
        total_len += stats.total_len;
    }
    
    writeln!(file, "{}", "=".repeat(50))?;
    writeln!(file, "Total sequences: {}", total_seqs)?;
    writeln!(file, "Total length: {}", total_len)?;
    
    Ok(())
}

fn combine_fastq_files(files: &[FastqFile], output_path: &Path) -> Result<()> {
    // Check that all files have the same compression type
    let compression_types: HashSet<CompressionType> = files.iter()
        .map(|f| f.compression.clone())
        .collect();
    
    if compression_types.len() > 1 {
        return Err(anyhow!("Cannot combine files with different compression types"));
    }
    
    let _compression = files.first().unwrap().compression.clone();
    
    // Calculate total bytes to process
    let total_bytes: u64 = files.iter().map(|f| f.size).sum();
    
    let mut output_file = File::create(output_path)?;
    
    let pb = ProgressBar::new(total_bytes);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{bar:40.cyan/blue}] {bytes}/{total_bytes} {bytes_per_sec} | ETA: {eta} | {msg}")
            .unwrap()
            .progress_chars("#>-")
    );
    
    let start_time = Instant::now();
    let mut bytes_processed = 0u64;
    
    for file in files {
        let file_name = file.path.file_name().unwrap().to_string_lossy();
        pb.set_message(format!("Processing: {}", file_name));
        
        let mut input_file = File::open(&file.path)?;
        let mut buffer = [0; 8192];
        
        loop {
            let bytes_read = input_file.read(&mut buffer)?;
            if bytes_read == 0 {
                break;
            }
            
            output_file.write_all(&buffer[..bytes_read])?;
            bytes_processed += bytes_read as u64;
            pb.set_position(bytes_processed);
        }
    }
    
    let elapsed = start_time.elapsed();
    let throughput = if elapsed.as_secs() > 0 {
        format_size(bytes_processed / elapsed.as_secs())
    } else {
        format_size(bytes_processed)
    };
    
    pb.finish_with_message(format!("Files combined - {} in {:.1}s ({})", 
        format_size(bytes_processed), 
        elapsed.as_secs_f32(),
        throughput + "/s"
    ));
    
    Ok(())
}

fn generate_blake3_hash(path: &Path) -> Result<String> {
    let mut hasher = Hasher::new();
    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(8192, file);
    
    loop {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() {
            break;
        }
        hasher.update(buffer);
        let consumed = buffer.len();
        reader.consume(consumed);
    }
    
    Ok(hasher.finalize().to_hex().to_string())
}

fn format_size(size: u64) -> String {
    const UNITS: &[&str] = &["B", "KB", "MB", "GB", "TB"];
    
    if size == 0 {
        return "0 B".to_string();
    }
    
    let mut size_f = size as f64;
    let mut unit_idx = 0;
    
    while size_f >= 1024.0 && unit_idx < UNITS.len() - 1 {
        size_f /= 1024.0;
        unit_idx += 1;
    }
    
    if unit_idx == 0 {
        format!("{} B", size)
    } else {
        format!("{:.2} {}", size_f, UNITS[unit_idx])
    }
}

fn check_for_duplicates(files: &[FastqFile]) -> HashSet<PathBuf> {
    let mut seen = HashSet::new();
    let mut duplicates = HashSet::new();
    
    for file in files {
        let canonical = file.path.canonicalize().unwrap_or(file.path.clone());
        if !seen.insert(canonical.clone()) {
            duplicates.insert(file.path.clone());
        }
    }
    
    duplicates
}

fn generate_output_filename(input_path: &Path, output_dir: &str) -> Result<PathBuf> {
    let files = find_fastq_files(input_path)?;
    if files.is_empty() {
        return Err(anyhow!("No FASTQ files found to generate output filename"));
    }
    
    let groups = group_files_by_base(&files);
    let first_group = groups.iter().next()
        .ok_or_else(|| anyhow!("No file groups found"))?;
    
    let base_name = &first_group.0;
    let first_file = &first_group.1[0];
    
    let extension = match first_file.compression {
        CompressionType::Gzip => ".all.fastq.gz",
        CompressionType::Bzip2 => ".all.fastq.bz2",
        CompressionType::Xz => ".all.fastq.xz",
        CompressionType::Zip => ".all.fastq.zip",
        CompressionType::Zstd => ".all.fastq.zst",
        CompressionType::None => ".all.fastq",
    };
    
    let output_file = format!("{}{}", base_name, extension);
    let output_path = Path::new(output_dir).join(output_file);
    
    println!("📝 Auto-generated output filename: {}", output_path.display());
    Ok(output_path)
}

fn generate_file_hashes(files: &[FastqFile]) -> Result<BTreeMap<PathBuf, String>> {
    let mut hashes = BTreeMap::new();
    
    println!("\n🔐 Generating Blake3 hashes for input files...");
    let pb = ProgressBar::new(files.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{bar:40.cyan/blue}] {pos}/{len} {msg}")
            .unwrap()
            .progress_chars("#>-")
    );
    
    for (idx, file) in files.iter().enumerate() {
        pb.set_position(idx as u64);
        pb.set_message(format!("Hashing {}", file.path.file_name().unwrap().to_string_lossy()));
        
        let hash = generate_blake3_hash(&file.path)?;
        hashes.insert(file.path.clone(), hash);
    }
    
    pb.finish_with_message("Hashes generated");
    Ok(hashes)
}

fn parse_nanopore_metadata(files: &[FastqFile]) -> Result<BTreeMap<PathBuf, Option<NanoporeMetadata>>> {
    let mut metadata_map = BTreeMap::new();
    
    for file in files {
        let metadata = parse_single_file_metadata(&file.path)?;
        metadata_map.insert(file.path.clone(), metadata);
    }
    
    Ok(metadata_map)
}

fn parse_single_file_metadata(path: &Path) -> Result<Option<NanoporeMetadata>> {
    let file = File::open(path)?;
    let reader: Box<dyn BufRead> = if path.to_string_lossy().ends_with(".gz") {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    
    let mut lines = reader.lines();
    if let Some(Ok(header)) = lines.next() {
        if header.starts_with('@') {
            return parse_nanopore_header(&header);
        }
    }
    
    Ok(None)
}

fn parse_nanopore_header(header: &str) -> Result<Option<NanoporeMetadata>> {
    let re = Regex::new(r"runid=(\S+).*flow_cell_id=(\S+).*start_time=(\S+).*protocol_group_id=(\S+).*sample_id=(\S+).*basecall_model_version_id=(\S+)")?;
    
    if let Some(caps) = re.captures(header) {
        let start_time = DateTime::parse_from_rfc3339(&caps[3])
            .or_else(|_| DateTime::parse_from_str(&caps[3], "%Y-%m-%dT%H:%M:%S%.f%:z"))?;
        
        Ok(Some(NanoporeMetadata {
            runid: caps[1].to_string(),
            flow_cell_id: caps[2].to_string(),
            start_time,
            protocol_group_id: caps[4].to_string(),
            sample_id: caps[5].to_string(),
            basecall_model: caps[6].to_string(),
        }))
    } else {
        Ok(None)
    }
}

fn display_nanopore_summary(metadata_map: &BTreeMap<PathBuf, Option<NanoporeMetadata>>) -> Result<()> {
    let valid_metadata: Vec<_> = metadata_map.values()
        .filter_map(|m| m.as_ref())
        .collect();
    
    if valid_metadata.is_empty() {
        return Ok(());
    }
    
    println!("\n🧬 Oxford Nanopore Metadata Summary:");
    println!("{}", "=".repeat(50));
    
    let runids: HashSet<_> = valid_metadata.iter().map(|m| &m.runid).collect();
    println!("  Run IDs found: {}", runids.len());
    for runid in &runids {
        println!("    - {}", runid);
    }
    
    let flow_cells: HashSet<_> = valid_metadata.iter().map(|m| &m.flow_cell_id).collect();
    println!("  Flow cells: {:?}", flow_cells);
    
    let start_times: Vec<_> = valid_metadata.iter().map(|m| m.start_time).collect();
    if let (Some(first), Some(last)) = (start_times.iter().min(), start_times.iter().max()) {
        println!("  First start time: {}", first);
        println!("  Last start time: {}", last);
        let duration = last.signed_duration_since(*first);
        println!("  Duration span: {} hours {} minutes", 
            duration.num_hours(), 
            duration.num_minutes() % 60);
    }
    
    let protocol_groups: HashSet<_> = valid_metadata.iter().map(|m| &m.protocol_group_id).collect();
    println!("  Protocol groups: {:?}", protocol_groups);
    
    let sample_ids: HashSet<_> = valid_metadata.iter().map(|m| &m.sample_id).collect();
    println!("  Sample IDs: {:?}", sample_ids);
    
    let models: HashSet<_> = valid_metadata.iter().map(|m| &m.basecall_model).collect();
    println!("  Basecall models: {:?}", models);
    
    Ok(())
}