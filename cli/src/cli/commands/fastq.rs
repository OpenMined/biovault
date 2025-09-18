use anyhow::{anyhow, Result};
use blake3::Hasher;
use chrono::{DateTime, FixedOffset};
use indicatif::{ProgressBar, ProgressStyle};
use regex::Regex;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

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
    stats_format: Option<String>,
) -> Result<()> {
    let input_path = Path::new(&input_folder);

    // Temporarily moved - will be after metadata parsing

    if !input_path.exists() {
        return Err(anyhow!("Input folder does not exist: {}", input_folder));
    }

    if !input_path.is_dir() {
        return Err(anyhow!("Input path is not a directory: {}", input_folder));
    }

    println!("Scanning for FASTQ files in: {}", input_folder);

    // Find all FASTQ files
    let fastq_files = find_fastq_files(input_path)?;

    if fastq_files.is_empty() {
        return Err(anyhow!("No FASTQ files found in {}", input_folder));
    }

    // Early duplicate detection based on file size (quick check)
    let size_duplicates = check_size_duplicates(&fastq_files);
    if !size_duplicates.is_empty() {
        println!("\n⚠️  Potential duplicate files detected (same size):");
        for (size, files) in &size_duplicates {
            if files.len() > 1 {
                println!("  Size {}: {} files", format_size(*size), files.len());
                for file in files {
                    println!("    - {}", file.path.file_name().unwrap().to_string_lossy());
                }
            }
        }
        print!("\n❓ Do you want to continue despite potential duplicates? (y/n): ");
        io::stdout().flush()?;
        let mut input = String::new();
        io::stdin().read_line(&mut input)?;
        if input.trim().to_lowercase() != "y" {
            return Err(anyhow!("Aborted: Potential duplicate files detected"));
        }
    }

    // Group files by base name
    let grouped_files = group_files_by_base(&fastq_files);

    // Display file groups
    for (base_name, files) in &grouped_files {
        println!("\n📁 File group: {}", base_name);

        // Check for missing sequence numbers
        let sequence_numbers: Vec<u32> = files.iter().filter_map(|f| f.sequence_number).collect();

        if !sequence_numbers.is_empty() {
            let min_seq = *sequence_numbers.iter().min().unwrap();
            let max_seq = *sequence_numbers.iter().max().unwrap();

            let mut missing = Vec::new();
            for i in min_seq..=max_seq {
                if !sequence_numbers.contains(&i) {
                    missing.push(i);
                }
            }

            println!(
                "  Files: {} (sequences {}-{})",
                files.len(),
                min_seq,
                max_seq
            );

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
        let compression_types: HashSet<CompressionType> =
            files.iter().map(|f| f.compression.clone()).collect();

        if compression_types.len() > 1 {
            println!("  ⚠️  Multiple compression types detected:");
            let mut compression_files: HashMap<CompressionType, Vec<String>> = HashMap::new();
            for file in files {
                let file_name = file.path.file_name().unwrap().to_string_lossy().to_string();
                compression_files
                    .entry(file.compression.clone())
                    .or_default()
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
            return Err(anyhow!(
                "Cannot combine files with different compression types"
            ));
        }
    }

    // Select which group to combine (for now, we'll combine all files)
    let all_files: Vec<FastqFile> = grouped_files
        .into_iter()
        .flat_map(|(_, files)| files)
        .collect();

    // Parse Nanopore metadata from first file of each group
    println!("\n🔬 Checking for Oxford Nanopore metadata...");
    let metadata_map = parse_nanopore_metadata(&all_files)?;

    if !metadata_map.is_empty() {
        display_nanopore_summary(&metadata_map)?;

        // Check for multiple sample IDs
        let sample_ids: HashSet<String> = metadata_map
            .values()
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

    // Determine if output is a directory or file
    let output_path = {
        let out_path = Path::new(&output_file);

        // Check if it's explicitly a directory (ends with /)
        if output_file.ends_with('/') {
            // Create directory if it doesn't exist
            if !out_path.exists() {
                fs::create_dir_all(&output_file)?;
                println!("📁 Created output directory: {}", output_file);
            }
            generate_output_filename_with_metadata(&all_files, &output_file, &metadata_map)?
        }
        // Check if it exists and is a directory
        else if out_path.exists() && out_path.is_dir() {
            generate_output_filename_with_metadata(&all_files, &output_file, &metadata_map)?
        }
        // If it doesn't exist and has no extension, ask if it should be a directory
        else if !out_path.exists() && !output_file.contains('.') {
            print!(
                "\n❓ '{}' doesn't exist and has no extension. Create as directory? (y/n): ",
                output_file
            );
            io::stdout().flush()?;

            let mut input = String::new();
            io::stdin().read_line(&mut input)?;

            if input.trim().to_lowercase() == "y" {
                fs::create_dir_all(&output_file)?;
                println!("📁 Created output directory: {}", output_file);
                generate_output_filename_with_metadata(&all_files, &output_file, &metadata_map)?
            } else {
                // Treat as a file without extension
                PathBuf::from(&output_file)
            }
        }
        // Otherwise treat as a file
        else {
            PathBuf::from(&output_file)
        }
    };

    // Create output directory if it doesn't exist
    if let Some(parent) = output_path.parent() {
        fs::create_dir_all(parent)?;
    }

    // Check if output file already exists
    if output_path.exists() {
        println!(
            "\n❌ Error: Output file already exists: {}",
            output_path.display()
        );
        println!("   Please remove the existing file or choose a different output location.");
        return Err(anyhow!("Output file already exists"));
    }

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
                .progress_chars("#>-"),
        );

        let mut validation_errors = Vec::new();

        for (idx, file) in all_files.iter().enumerate() {
            pb.set_position(idx as u64);
            pb.set_message(format!(
                "Validating {}",
                file.path.file_name().unwrap().to_string_lossy()
            ));

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
            return Err(anyhow!(
                "Validation failed for {} files",
                validation_errors.len()
            ));
        }

        println!("✅ All files validated successfully!");

        // Generate hashes and check for content duplicates
        let file_hashes = generate_file_hashes(&all_files)?;
        let content_duplicates = check_content_duplicates(&file_hashes);
        if !content_duplicates.is_empty() {
            println!("\n⚠️  DUPLICATE FILES DETECTED (identical content):");
            for (hash, files) in &content_duplicates {
                if files.len() > 1 {
                    println!("  Blake3: {}...", &hash[..16]);
                    for file in files {
                        println!("    - {}", file.display());
                    }
                }
            }
            print!("❓ Continue anyway? (y/n): ");
            io::stdout().flush()?;
            let mut input = String::new();
            io::stdin().read_line(&mut input)?;
            if input.trim().to_lowercase() != "y" {
                return Err(anyhow!("Aborted: Duplicate files detected"));
            }
        }

        // Save validation stats to file
        let format = stats_format.as_deref().unwrap_or("tsv");
        let ext = match format {
            "yaml" => "stats.yaml",
            "json" => "stats.json",
            _ => "stats.tsv",
        };
        let stats_file = output_path.with_extension(format!("pre_combine_{}", ext));
        save_stats_to_file(&stats_map, &stats_file, Some(&file_hashes), format)?;
        println!("📄 Validation stats saved to: {}", stats_file.display());
    }

    // Combine files
    println!("\n🔗 Combining FASTQ files...");
    combine_fastq_files(&all_files, &output_path)?;

    // Generate Blake3 hash (only once)
    println!("\n🔐 Generating Blake3 hash...");
    let final_hash = generate_blake3_hash(&output_path)?;
    let hash_file = output_path.with_extension("blake3");
    fs::write(&hash_file, &final_hash)?;
    println!("📄 Blake3 hash saved to: {}", hash_file.display());

    // Validate combined file if validation was requested
    if should_validate {
        println!("\n🔍 Validating combined file...");

        let combined_stats = validate_fastq_file(&output_path)?;

        // Calculate input stats
        let total_input_seqs: u64 = stats_map.values().map(|s| s.num_seqs).sum();
        let input_file_count = stats_map.len();

        // Get actual file sizes (not sequence lengths)
        let total_input_size: u64 = all_files.iter().map(|f| f.size).sum();

        // Get output file size
        let output_metadata = fs::metadata(&output_path)?;
        let output_size = output_metadata.len();

        // Display comparison
        println!("\n📊 Validation Summary:");
        println!("───────────────────────────────────────────");
        println!("INPUT:");
        println!("  Files: {} files", input_file_count);
        println!("  Total size: {}", format_size(total_input_size));
        println!("  Total sequences: {}", total_input_seqs);
        println!();
        println!("OUTPUT:");
        println!("  Files: 1 file");
        println!("  Total size: {}", format_size(output_size));
        println!("  Total sequences: {}", combined_stats.num_seqs);
        println!("───────────────────────────────────────────");

        // Check if sequences match
        if combined_stats.num_seqs != total_input_seqs {
            println!("\n❌ ERROR: Sequence count mismatch!");
            println!(
                "   Expected {} sequences but found {}",
                total_input_seqs, combined_stats.num_seqs
            );
            println!(
                "   Missing {} sequences",
                total_input_seqs.abs_diff(combined_stats.num_seqs)
            );
            return Err(anyhow!(
                "Combined file validation failed: sequence count mismatch"
            ));
        } else {
            println!("\n✅ Validation successful: Input and output sequence counts match!");
            println!("   Average length: {:.2} bp", combined_stats.avg_len);
            println!("   GC content: {:.2}%", combined_stats.gc_content);
        }

        // Save final stats (reuse the hash already calculated)
        let format = stats_format.as_deref().unwrap_or("tsv");
        let ext = match format {
            "yaml" => "stats.yaml",
            "json" => "stats.json",
            _ => "stats.tsv",
        };
        let final_stats_file = output_path.with_extension(ext);
        let mut final_stats = BTreeMap::new();
        final_stats.insert(output_path.to_path_buf(), combined_stats);
        let mut final_hashes = BTreeMap::new();
        final_hashes.insert(output_path.to_path_buf(), final_hash.clone());
        save_stats_to_file(&final_stats, &final_stats_file, Some(&final_hashes), format)?;
        println!("\n📄 Final stats saved to: {}", final_stats_file.display());
    }

    println!(
        "\n✨ FASTQ files successfully combined to: {}",
        output_path.display()
    );

    Ok(())
}

fn find_fastq_files(dir: &Path) -> Result<Vec<FastqFile>> {
    let mut files = Vec::new();

    // Define patterns for FASTQ files
    let patterns = [
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
    files.sort_by(|a, b| match a.base_name.cmp(&b.base_name) {
        std::cmp::Ordering::Equal => a.sequence_number.cmp(&b.sequence_number),
        other => other,
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
    let mut groups: BTreeMap<String, Vec<FastqFile>> = BTreeMap::new();

    for file in files {
        groups
            .entry(file.base_name.clone())
            .or_default()
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
        .arg("-T") // Tab-separated output
        .arg("-a") // All statistics
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
        gc_content: fields[14].parse().unwrap_or(0.0),
    })
}

fn save_stats_to_file(
    stats: &BTreeMap<PathBuf, FastqStats>,
    path: &Path,
    hashes: Option<&BTreeMap<PathBuf, String>>,
    format: &str,
) -> Result<()> {
    match format {
        "yaml" => save_stats_yaml(stats, path, hashes),
        "json" => save_stats_json(stats, path, hashes),
        _ => save_stats_tsv(stats, path, hashes),
    }
}

fn save_stats_tsv(
    stats: &BTreeMap<PathBuf, FastqStats>,
    path: &Path,
    hashes: Option<&BTreeMap<PathBuf, String>>,
) -> Result<()> {
    let mut file = File::create(path)?;

    // Write header
    writeln!(file, "file\tsequences\ttotal_length\tmin_length\tavg_length\tmax_length\tgc_content\tblake3_hash")?;

    // Write data rows
    for (file_path, stats) in stats {
        let hash = hashes
            .and_then(|h| h.get(file_path))
            .map(|h| h.as_str())
            .unwrap_or("");

        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{:.2}\t{}\t{:.2}\t{}",
            file_path.display(),
            stats.num_seqs,
            stats.total_len,
            stats.min_len,
            stats.avg_len,
            stats.max_len,
            stats.gc_content,
            hash
        )?;
    }

    Ok(())
}

fn save_stats_yaml(
    stats: &BTreeMap<PathBuf, FastqStats>,
    path: &Path,
    hashes: Option<&BTreeMap<PathBuf, String>>,
) -> Result<()> {
    use serde::Serialize;

    #[derive(Serialize)]
    struct StatsEntry {
        file: String,
        sequences: u64,
        total_length: u64,
        min_length: u64,
        avg_length: f64,
        max_length: u64,
        gc_content: f64,
        blake3_hash: Option<String>,
    }

    let entries: Vec<StatsEntry> = stats
        .iter()
        .map(|(file_path, stats)| StatsEntry {
            file: file_path.display().to_string(),
            sequences: stats.num_seqs,
            total_length: stats.total_len,
            min_length: stats.min_len,
            avg_length: stats.avg_len,
            max_length: stats.max_len,
            gc_content: stats.gc_content,
            blake3_hash: hashes.and_then(|h| h.get(file_path).cloned()),
        })
        .collect();

    let yaml = serde_yaml::to_string(&entries)?;
    fs::write(path, yaml)?;

    Ok(())
}

fn save_stats_json(
    stats: &BTreeMap<PathBuf, FastqStats>,
    path: &Path,
    hashes: Option<&BTreeMap<PathBuf, String>>,
) -> Result<()> {
    use serde::Serialize;

    #[derive(Serialize)]
    struct StatsEntry {
        file: String,
        sequences: u64,
        total_length: u64,
        min_length: u64,
        avg_length: f64,
        max_length: u64,
        gc_content: f64,
        blake3_hash: Option<String>,
    }

    let entries: Vec<StatsEntry> = stats
        .iter()
        .map(|(file_path, stats)| StatsEntry {
            file: file_path.display().to_string(),
            sequences: stats.num_seqs,
            total_length: stats.total_len,
            min_length: stats.min_len,
            avg_length: stats.avg_len,
            max_length: stats.max_len,
            gc_content: stats.gc_content,
            blake3_hash: hashes.and_then(|h| h.get(file_path).cloned()),
        })
        .collect();

    let json = serde_json::to_string_pretty(&entries)?;
    fs::write(path, json)?;

    Ok(())
}

fn combine_fastq_files(files: &[FastqFile], output_path: &Path) -> Result<()> {
    // Check that all files have the same compression type
    let compression_types: HashSet<CompressionType> =
        files.iter().map(|f| f.compression.clone()).collect();

    if compression_types.len() > 1 {
        return Err(anyhow!(
            "Cannot combine files with different compression types"
        ));
    }

    let compression = files.first().unwrap().compression.clone();

    // Just use direct concatenation - it's faster and preserves original compression
    combine_fastq_files_direct(files, output_path, compression)
}

fn combine_fastq_files_direct(
    files: &[FastqFile],
    output_path: &Path,
    _compression: CompressionType,
) -> Result<()> {
    // Calculate total bytes to process
    #[cfg(not(test))]
    let total_bytes: u64 = files.iter().map(|f| f.size).sum();
    #[cfg(test)]
    let _total_bytes: u64 = files.iter().map(|f| f.size).sum();

    let mut output_file = File::create(output_path)?;

    let pb = {
        #[cfg(test)]
        {
            // Hide progress overhead during tests to keep slow suite leaner
            ProgressBar::hidden()
        }
        #[cfg(not(test))]
        {
            let pb = ProgressBar::new(total_bytes);
            pb.set_style(
                ProgressStyle::default_bar()
                    .template("{spinner:.green} [{bar:40.cyan/blue}] {bytes}/{total_bytes} {bytes_per_sec} | ETA: {eta} | {msg}")
                    .unwrap()
                    .progress_chars("#>-")
            );
            pb
        }
    };

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

    pb.finish_with_message(format!(
        "Files combined - {} in {:.1}s ({}/s)",
        format_size(bytes_processed),
        elapsed.as_secs_f32(),
        throughput
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

fn check_size_duplicates(files: &[FastqFile]) -> HashMap<u64, Vec<&FastqFile>> {
    let mut size_map: HashMap<u64, Vec<&FastqFile>> = HashMap::new();

    for file in files {
        size_map.entry(file.size).or_default().push(file);
    }

    // Only keep sizes with multiple files
    size_map.retain(|_, files| files.len() > 1);
    size_map
}

fn check_content_duplicates(hashes: &BTreeMap<PathBuf, String>) -> HashMap<String, Vec<PathBuf>> {
    let mut hash_map: HashMap<String, Vec<PathBuf>> = HashMap::new();

    for (path, hash) in hashes {
        hash_map.entry(hash.clone()).or_default().push(path.clone());
    }

    // Only keep hashes with multiple files
    hash_map.retain(|_, files| files.len() > 1);
    hash_map
}

fn generate_output_filename_with_metadata(
    files: &[FastqFile],
    output_dir: &str,
    metadata_map: &BTreeMap<PathBuf, Option<NanoporeMetadata>>,
) -> Result<PathBuf> {
    if files.is_empty() {
        return Err(anyhow!("No FASTQ files found to generate output filename"));
    }

    let valid_metadata: Vec<_> = metadata_map.values().filter_map(|m| m.as_ref()).collect();

    let first_file = &files[0];
    let extension = match first_file.compression {
        CompressionType::Gzip => "fastq.gz",
        CompressionType::Bzip2 => "fastq.bz2",
        CompressionType::Xz => "fastq.xz",
        CompressionType::Zip => "fastq.zip",
        CompressionType::Zstd => "fastq.zst",
        CompressionType::None => "fastq",
    };

    let suggested_name = if !valid_metadata.is_empty() {
        // Extract unique values
        let mut flow_cells: Vec<String> = valid_metadata
            .iter()
            .map(|m| m.flow_cell_id.clone())
            .collect::<HashSet<_>>()
            .into_iter()
            .collect();
        flow_cells.sort();

        let mut sample_ids: Vec<String> = valid_metadata
            .iter()
            .map(|m| m.sample_id.clone())
            .collect::<HashSet<_>>()
            .into_iter()
            .collect();
        sample_ids.sort();

        // Get the latest date
        let latest_time = valid_metadata.iter().map(|m| m.start_time).max().unwrap();
        let date_str = latest_time.format("%Y%m%d").to_string();

        // Build filename: flowcells-samples-ONT-date-all.extension
        format!(
            "{}-{}-ONT-{}-all.{}",
            flow_cells.join("-"),
            sample_ids.join("-"),
            date_str,
            extension
        )
    } else {
        // Fallback to base name if no metadata
        let groups = group_files_by_base(files);
        let first_group = groups
            .iter()
            .next()
            .ok_or_else(|| anyhow!("No file groups found"))?;
        let base_name = &first_group.0;
        format!("{}.all.{}", base_name, extension)
    };

    let mut output_path = Path::new(output_dir).join(&suggested_name);

    // Non-interactive mode for tests or when BIOVAULT_NONINTERACTIVE=1
    let noninteractive = cfg!(test)
        || std::env::var("BIOVAULT_NONINTERACTIVE")
            .map(|v| v == "1")
            .unwrap_or(false);

    if !noninteractive {
        println!("\n📝 Suggested output filename: {}", output_path.display());
        print!("Accept this filename? (y/n/e to edit): ");
        io::stdout().flush()?;

        let mut input = String::new();
        io::stdin().read_line(&mut input)?;
        let choice = input.trim().to_lowercase();

        if choice == "n" {
            return Err(anyhow!("Filename generation cancelled"));
        } else if choice == "e" {
            print!("Enter new filename (without path) [{}]: ", suggested_name);
            io::stdout().flush()?;
            input.clear();
            io::stdin().read_line(&mut input)?;
            let new_name = input.trim();

            // Use suggested name if user just presses enter
            let final_name = if new_name.is_empty() {
                &suggested_name
            } else {
                new_name
            };

            output_path = Path::new(output_dir).join(final_name);
            println!("Using filename: {}", output_path.display());
        } else if choice != "y" {
            println!("Using suggested filename: {}", output_path.display());
        }
    }

    Ok(output_path)
}

fn generate_file_hashes(files: &[FastqFile]) -> Result<BTreeMap<PathBuf, String>> {
    let mut hashes = BTreeMap::new();

    println!("\n🔐 Generating Blake3 hashes for input files...");
    let pb = {
        #[cfg(test)]
        {
            ProgressBar::hidden()
        }
        #[cfg(not(test))]
        {
            let pb = ProgressBar::new(files.len() as u64);
            pb.set_style(
                ProgressStyle::default_bar()
                    .template("{spinner:.green} [{bar:40.cyan/blue}] {pos}/{len} {msg}")
                    .unwrap()
                    .progress_chars("#>-"),
            );
            pb
        }
    };

    for (idx, file) in files.iter().enumerate() {
        pb.set_position(idx as u64);
        pb.set_message(format!(
            "Hashing {}",
            file.path.file_name().unwrap().to_string_lossy()
        ));

        let hash = generate_blake3_hash(&file.path)?;
        hashes.insert(file.path.clone(), hash);
    }

    pb.finish_with_message("Hashes generated");
    Ok(hashes)
}

fn parse_nanopore_metadata(
    files: &[FastqFile],
) -> Result<BTreeMap<PathBuf, Option<NanoporeMetadata>>> {
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
            // Try to parse, but don't fail if it doesn't match
            match parse_nanopore_header(&header) {
                Ok(metadata) => return Ok(metadata),
                Err(_) => return Ok(None), // Not a Nanopore file or different format
            }
        }
    }

    Ok(None)
}

fn parse_nanopore_header(header: &str) -> Result<Option<NanoporeMetadata>> {
    // More flexible regex that captures fields individually
    let runid_re = Regex::new(r"runid=(\S+)")?;
    let flow_cell_re = Regex::new(r"flow_cell_id=(\S+)")?;
    let start_time_re = Regex::new(r"start_time=(\S+)")?;
    let protocol_re = Regex::new(r"protocol_group_id=(\S+)")?;
    let sample_re = Regex::new(r"sample_id=(\S+)")?;
    let model_re = Regex::new(r"basecall_model_version_id=(\S+)")?;

    // Check if all required fields are present
    if let (
        Some(runid_cap),
        Some(flow_cap),
        Some(time_cap),
        Some(proto_cap),
        Some(sample_cap),
        Some(model_cap),
    ) = (
        runid_re.captures(header),
        flow_cell_re.captures(header),
        start_time_re.captures(header),
        protocol_re.captures(header),
        sample_re.captures(header),
        model_re.captures(header),
    ) {
        let time_str = &time_cap[1];
        // Try multiple date formats
        let start_time_result = DateTime::parse_from_rfc3339(time_str)
            .or_else(|_| DateTime::parse_from_str(time_str, "%Y-%m-%dT%H:%M:%S%.f%:z"))
            .or_else(|_| DateTime::parse_from_str(time_str, "%Y-%m-%dT%H:%M:%S%.f%#z"))
            .or_else(|_| {
                // Try with manual timezone parsing for formats like -03:00
                if let Some(idx) = time_str.rfind(['+', '-']) {
                    let (dt_part, tz_part) = time_str.split_at(idx);
                    let full_str = format!("{}{}", dt_part, tz_part);
                    DateTime::parse_from_str(&full_str, "%Y-%m-%dT%H:%M:%S%.f%:z")
                } else {
                    // Return an error that will be handled below
                    DateTime::parse_from_str("invalid", "%Y-%m-%dT%H:%M:%S%.f%:z")
                }
            });

        let start_time = match start_time_result {
            Ok(dt) => dt,
            Err(_) => return Ok(None), // Can't parse date, not a valid Nanopore file
        };

        Ok(Some(NanoporeMetadata {
            runid: runid_cap[1].to_string(),
            flow_cell_id: flow_cap[1].to_string(),
            start_time,
            protocol_group_id: proto_cap[1].to_string(),
            sample_id: sample_cap[1].to_string(),
            basecall_model: model_cap[1].to_string(),
        }))
    } else {
        Ok(None)
    }
}

fn display_nanopore_summary(
    metadata_map: &BTreeMap<PathBuf, Option<NanoporeMetadata>>,
) -> Result<()> {
    let valid_metadata: Vec<_> = metadata_map.values().filter_map(|m| m.as_ref()).collect();

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
        println!(
            "  Duration span: {} hours {} minutes",
            duration.num_hours(),
            duration.num_minutes() % 60
        );
    }

    let protocol_groups: HashSet<_> = valid_metadata
        .iter()
        .map(|m| &m.protocol_group_id)
        .collect();
    println!("  Protocol groups: {:?}", protocol_groups);

    let sample_ids: HashSet<_> = valid_metadata.iter().map(|m| &m.sample_id).collect();
    println!("  Sample IDs: {:?}", sample_ids);

    let models: HashSet<_> = valid_metadata.iter().map(|m| &m.basecall_model).collect();
    println!("  Basecall models: {:?}", models);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::TimeZone;
    use tempfile::TempDir;

    #[test]
    fn compression_type_and_format_size() {
        assert_eq!(
            CompressionType::from_path(Path::new("a.fastq")),
            CompressionType::None
        );
        assert_eq!(
            CompressionType::from_path(Path::new("a.fq.gz")),
            CompressionType::Gzip
        );
        assert_eq!(
            CompressionType::from_path(Path::new("a.fastq.bz2")),
            CompressionType::Bzip2
        );
        assert_eq!(
            CompressionType::from_path(Path::new("a.fastq.xz")),
            CompressionType::Xz
        );
        assert_eq!(
            CompressionType::from_path(Path::new("a.fastq.zip")),
            CompressionType::Zip
        );
        assert_eq!(
            CompressionType::from_path(Path::new("a.fastq.zst")),
            CompressionType::Zstd
        );

        // Keep a single small size check
        assert_eq!(format_size(512), "512 B");
    }

    #[test]
    fn size_and_content_duplicates_detection() {
        let f1 = FastqFile {
            path: PathBuf::from("a.fastq"),
            size: 10,
            compression: CompressionType::None,
            sequence_number: Some(1),
            base_name: "a".into(),
        };
        let f2 = FastqFile {
            path: PathBuf::from("b.fastq"),
            size: 10,
            compression: CompressionType::None,
            sequence_number: Some(2),
            base_name: "a".into(),
        };
        let f3 = FastqFile {
            path: PathBuf::from("c.fastq"),
            size: 5,
            compression: CompressionType::None,
            sequence_number: None,
            base_name: "c".into(),
        };
        let binding = [f1.clone(), f2.clone(), f3.clone()];
        let dups = check_size_duplicates(&binding);
        assert!(dups.contains_key(&10));
        assert!(!dups.contains_key(&5));

        let mut hashes = BTreeMap::new();
        hashes.insert(f1.path.clone(), "h1".into());
        hashes.insert(f2.path.clone(), "h1".into());
        hashes.insert(f3.path.clone(), "h2".into());
        let cdups = check_content_duplicates(&hashes);
        assert!(cdups.contains_key("h1"));
        assert!(!cdups.contains_key("h2"));
    }

    #[test]
    fn group_and_sequence_extraction_and_find_files() {
        let tmp = TempDir::new().unwrap();
        // Create files that match patterns and with sequence numbers
        let files = ["sample_001.fastq.gz", "sample_002.fastq.gz", "other.fastq"];
        for name in &files {
            fs::write(tmp.path().join(name), b"@h\nA\n+\n!\n").unwrap();
        }
        let found = find_fastq_files(tmp.path()).unwrap();
        // Expect 3 files
        assert_eq!(found.len(), 3);
        // Grouping by base name should put first two together
        let groups = group_files_by_base(&found);
        assert!(groups.get("sample").unwrap().len() >= 2);
    }

    #[test]
    fn stats_savers_and_blake3() {
        let tmp = TempDir::new().unwrap();
        // Prepare a fake file and hash it
        let file = tmp.path().join("x.fastq");
        fs::write(&file, b"@h\nAC\n+\n!!\n").unwrap();
        let h = generate_blake3_hash(&file).unwrap();
        assert_eq!(h.len(), 64);

        let mut stats = BTreeMap::new();
        stats.insert(
            file.clone(),
            FastqStats {
                num_seqs: 1,
                total_len: 2,
                min_len: 2,
                avg_len: 2.0,
                max_len: 2,
                gc_content: 50.0,
            },
        );
        let mut hashes = BTreeMap::new();
        hashes.insert(file.clone(), h);
        let tsv = tmp.path().join("stats.tsv");
        let yaml = tmp.path().join("stats.yaml");
        let json = tmp.path().join("stats.json");
        save_stats_to_file(&stats, &tsv, Some(&hashes), "tsv").unwrap();
        save_stats_to_file(&stats, &yaml, Some(&hashes), "yaml").unwrap();
        save_stats_to_file(&stats, &json, Some(&hashes), "json").unwrap();
        assert!(tsv.exists() && yaml.exists() && json.exists());
    }

    #[test]
    fn nanopore_header_parsing_and_single_file_metadata() {
        // Valid header
        let header = "@... runid=RID flow_cell_id=FC start_time=2024-01-02T03:04:05+00:00 protocol_group_id=PG sample_id=S basecall_model_version_id=M v=1";
        let md = parse_nanopore_header(header).unwrap();
        assert!(md.is_some());
        // Invalid header
        let md2 = parse_nanopore_header("@no fields").unwrap();
        assert!(md2.is_none());

        // Test parse_single_file_metadata with gz
        let tmp = TempDir::new().unwrap();
        let gz = tmp.path().join("s.fastq.gz");
        {
            let f = File::create(&gz).unwrap();
            // Use fastest compression in tests to reduce overhead
            let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::fast());
            writeln!(enc, "{}", header).unwrap();
            writeln!(enc, "+").unwrap();
            enc.finish().unwrap();
        }
        let parsed = parse_single_file_metadata(&gz).unwrap();
        assert!(parsed.is_some());
    }

    #[test]
    fn output_filename_with_metadata_and_combine_direct() {
        let tmp = TempDir::new().unwrap();
        let f1 = tmp.path().join("a.fastq");
        let f2 = tmp.path().join("b.fastq");
        fs::write(&f1, b"@a\nA\n+\n!\n").unwrap();
        fs::write(&f2, b"@b\nC\n+\n!\n").unwrap();
        let files = vec![
            FastqFile {
                path: f1.clone(),
                size: fs::metadata(&f1).unwrap().len(),
                compression: CompressionType::None,
                sequence_number: None,
                base_name: "g".into(),
            },
            FastqFile {
                path: f2.clone(),
                size: fs::metadata(&f2).unwrap().len(),
                compression: CompressionType::None,
                sequence_number: None,
                base_name: "g".into(),
            },
        ];

        let mut meta = BTreeMap::new();
        meta.insert(
            f1.clone(),
            Some(NanoporeMetadata {
                runid: "r1".into(),
                flow_cell_id: "FC1".into(),
                start_time: FixedOffset::east_opt(0)
                    .unwrap()
                    .with_ymd_and_hms(2024, 1, 2, 3, 4, 5)
                    .unwrap(),
                protocol_group_id: "PG".into(),
                sample_id: "S1".into(),
                basecall_model: "M".into(),
            }),
        );
        meta.insert(
            f2.clone(),
            Some(NanoporeMetadata {
                runid: "r2".into(),
                flow_cell_id: "FC2".into(),
                start_time: FixedOffset::east_opt(0)
                    .unwrap()
                    .with_ymd_and_hms(2024, 1, 3, 0, 0, 0)
                    .unwrap(),
                protocol_group_id: "PG".into(),
                sample_id: "S2".into(),
                basecall_model: "M".into(),
            }),
        );

        let outdir = tmp.path().join("out");
        fs::create_dir_all(&outdir).unwrap();
        let name = generate_output_filename_with_metadata(&files, outdir.to_str().unwrap(), &meta)
            .unwrap();
        assert!(name
            .file_name()
            .unwrap()
            .to_string_lossy()
            .contains("FC1-FC2"));

        // Combine direct and verify concatenation length
        let combined = tmp.path().join("combined.fastq");
        combine_fastq_files_direct(&files, &combined, CompressionType::None).unwrap();
        let content = fs::read(&combined).unwrap();
        assert!(!content.is_empty());
    }

    #[test]
    fn parse_seqkit_stats_output() {
        let header = "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tQ1\tQ2\tQ3\tsum_gap\tN50\tQ20(%)\tGC(%)";
        let data = "a.fastq\tFASTQ\tDNA\t1\t2\t2\t2.0\t2\t0\t0\t0\t0\t0\t0\t50.0";
        let out = format!("{}\n{}\n", header, data);
        let s = parse_seqkit_stats(&out).unwrap();
        assert_eq!(s.num_seqs, 1);
        assert_eq!(s.total_len, 2);
        assert_eq!(s.gc_content, 50.0);
    }

    #[test]
    fn generate_file_hashes_and_output_name_fallback() {
        let tmp = TempDir::new().unwrap();
        let f1 = tmp.path().join("base_1.fastq");
        let f2 = tmp.path().join("base_2.fastq");
        fs::write(&f1, b"@h\nAA\n+\n!!\n").unwrap();
        fs::write(&f2, b"@h\nCC\n+\n!!\n").unwrap();
        let files = vec![
            FastqFile {
                path: f1.clone(),
                size: fs::metadata(&f1).unwrap().len(),
                compression: CompressionType::None,
                sequence_number: Some(1),
                base_name: "base".into(),
            },
            FastqFile {
                path: f2.clone(),
                size: fs::metadata(&f2).unwrap().len(),
                compression: CompressionType::None,
                sequence_number: Some(2),
                base_name: "base".into(),
            },
        ];
        let hashes = generate_file_hashes(&files).unwrap();
        assert_eq!(hashes.len(), 2);

        let meta = BTreeMap::new();
        let outdir = tmp.path().join("out");
        fs::create_dir_all(&outdir).unwrap();
        let name = generate_output_filename_with_metadata(&files, outdir.to_str().unwrap(), &meta)
            .unwrap();
        let s = name.file_name().unwrap().to_string_lossy();
        assert!(s.contains("base.all.fastq"));
    }

    #[test]
    fn combine_fastq_files_mixed_compression_fails() {
        let tmp = TempDir::new().unwrap();
        let f1 = tmp.path().join("a.fastq");
        let f2 = tmp.path().join("b.fastq.gz");
        fs::write(&f1, b"@a\nA\n+\n!\n").unwrap();
        fs::write(&f2, b"@b\nC\n+\n!\n").unwrap();
        let files = vec![
            FastqFile {
                path: f1.clone(),
                size: fs::metadata(&f1).unwrap().len(),
                compression: CompressionType::None,
                sequence_number: None,
                base_name: "g".into(),
            },
            FastqFile {
                path: f2.clone(),
                size: fs::metadata(&f2).unwrap().len(),
                compression: CompressionType::Gzip,
                sequence_number: None,
                base_name: "g".into(),
            },
        ];
        let err = combine_fastq_files(&files, &tmp.path().join("out.fastq"));
        assert!(err.is_err());
    }
}
