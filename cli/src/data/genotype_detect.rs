use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

// Male-specific Y chromosome markers (confirmed in 100% of males, 0% of females)
const MALE_Y_MARKERS: &[&str] = &[
    "rs11575897",
    "rs2534636",
    "i3000043",
    "i3000045",
    "i4000162",
    "rs13303871",
    "rs35284970",
    "rs3895",
    "i4000120",
    "i4000121",
];

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenotypeMetadata {
    pub data_type: String,
    pub source: Option<String>,
    pub grch_version: Option<String>,
    pub row_count: Option<i64>,
    pub chromosome_count: Option<i64>,
    pub inferred_sex: Option<String>,
}

impl Default for GenotypeMetadata {
    fn default() -> Self {
        Self {
            data_type: "Unknown".to_string(),
            source: None,
            grch_version: None,
            row_count: None,
            chromosome_count: None,
            inferred_sex: None,
        }
    }
}

/// Detect if a file is a genotype file and extract metadata (fast version - header only)
pub fn detect_genotype_metadata(file_path: &str) -> Result<GenotypeMetadata> {
    let path = Path::new(file_path);

    if !path.exists() {
        return Ok(GenotypeMetadata::default());
    }

    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Read first 50 lines for header analysis
    let mut header_text = String::new();
    let mut data_lines = Vec::new();
    let mut line_count = 0;
    let max_header_lines = 50;

    for line_result in &mut lines {
        let line = line_result?;
        line_count += 1;

        // Collect lines starting with # or comment markers as header
        if line.starts_with('#') || line.starts_with("//") {
            header_text.push_str(&line);
            header_text.push('\n');
        } else if !line.trim().is_empty() {
            // Non-comment, non-empty line - potential data
            data_lines.push(line);

            // Collect up to 10 data lines for structure validation
            if data_lines.len() >= 10 {
                break;
            }
        }

        if line_count >= max_header_lines {
            break;
        }
    }

    // Check if this looks like a genotype file
    let is_genotype = validate_genotype_structure(&data_lines);

    if !is_genotype {
        return Ok(GenotypeMetadata::default());
    }

    // File is a genotype - extract metadata from headers
    let header_lower = header_text.to_lowercase();

    let source = detect_source(&header_lower, &data_lines);
    let mut grch_version = detect_grch_version(&header_lower);

    if matches!(source.as_deref(), Some("FamilyTreeDNA"))
        && matches!(grch_version.as_deref(), Some("Unknown"))
    {
        grch_version = Some("37".to_string());
    }

    if matches!(source.as_deref(), Some("deCODEme"))
        && matches!(grch_version.as_deref(), Some("Unknown"))
    {
        grch_version = Some("36".to_string());
    }

    Ok(GenotypeMetadata {
        data_type: "Genotype".to_string(),
        source,
        grch_version,
        row_count: None,
        chromosome_count: None,
        inferred_sex: None,
    })
}

/// Analyze genotype file for row count, chromosome count, and inferred sex (slow - reads entire file)
pub fn analyze_genotype_file(file_path: &str) -> Result<GenotypeMetadata> {
    let path = Path::new(file_path);

    if !path.exists() {
        return Ok(GenotypeMetadata::default());
    }

    // First detect source from header (need this for sex inference)
    let metadata = detect_genotype_metadata(file_path)?;
    let source = metadata.source.clone();

    // Count total rows and chromosomes by reading the entire file
    let (row_count, chromosome_count, inferred_sex) =
        count_rows_and_chromosomes(file_path, source.as_deref())?;

    // Return combined metadata (detected + analyzed)
    Ok(GenotypeMetadata {
        data_type: metadata.data_type,
        source: metadata.source,
        grch_version: metadata.grch_version,
        row_count: Some(row_count as i64),
        chromosome_count: Some(chromosome_count as i64),
        inferred_sex: Some(inferred_sex),
    })
}

/// Count data rows and unique chromosomes in a genotype file
fn count_rows_and_chromosomes(
    file_path: &str,
    source: Option<&str>,
) -> Result<(usize, usize, String)> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let mut row_count = 0;
    let mut chromosomes = HashSet::new();
    let mut male_markers_called = 0;

    let is_23andme = source.map(|s| s.contains("23andMe")).unwrap_or(false);

    for line_result in reader.lines() {
        let line = line_result?;

        // Skip comments and empty lines
        if line.starts_with('#') || line.starts_with("//") || line.trim().is_empty() {
            continue;
        }

        let fields = split_fields(&line);

        let Some((chr_idx, _pos_idx, genotype_idx)) = identify_data_columns(&fields) else {
            continue;
        };

        row_count += 1;

        let rsid = fields[0].trim();
        let chr = fields[chr_idx].trim();
        let genotype = fields[genotype_idx].trim();

        chromosomes.insert(chr.to_string());

        if is_23andme
            && (chr == "Y" || chr == "24")
            && MALE_Y_MARKERS.contains(&rsid)
            && !genotype.is_empty()
            && genotype != "--"
            && genotype != "00"
        {
            male_markers_called += 1;
        }
    }

    // Infer sex based on source
    let inferred_sex = if is_23andme {
        // 23andMe: Use male-specific Y marker logic (validated heuristics)
        if male_markers_called >= 5 {
            "Male".to_string()
        } else if chromosomes.contains("X") || chromosomes.contains("23") {
            "Female".to_string()
        } else {
            "Unknown".to_string()
        }
    } else {
        // Other sources: Conservative inference - we don't have validated heuristics yet
        let has_x = chromosomes.contains("X") || chromosomes.contains("23");
        let has_y = chromosomes.contains("Y") || chromosomes.contains("24");

        if has_x && has_y {
            // Both X and Y present - return Unknown since we don't have validated heuristics
            // for determining if Y chromosome data indicates male sex in non-23andMe files
            "Unknown".to_string()
        } else if has_x && !has_y {
            // Only X, no Y - likely female
            "Female".to_string()
        } else {
            "Unknown".to_string()
        }
    };

    Ok((row_count, chromosomes.len(), inferred_sex))
}

/// Validate that data lines match genotype file structure
fn validate_genotype_structure(data_lines: &[String]) -> bool {
    if data_lines.is_empty() {
        return false;
    }

    let mut valid_lines = 0;

    for line in data_lines {
        let parts = split_fields(line);

        if identify_data_columns(&parts).is_some() {
            valid_lines += 1;
        }
    }

    // Consider it a genotype file if at least 70% of data lines are valid
    let threshold = (data_lines.len() * 7) / 10;
    valid_lines >= threshold
}

fn is_valid_chromosome(chr: &str) -> bool {
    // Check for numeric chromosomes 1-22
    if let Ok(num) = chr.parse::<u8>() {
        return (1..=22).contains(&num);
    }

    // Check for sex chromosomes and mitochondrial
    matches!(
        chr.to_uppercase().as_str(),
        "X" | "Y" | "MT" | "M" | "23" | "24"
    )
}

fn is_valid_genotype(genotype: &str) -> bool {
    // Valid genotypes include:
    // - Two-letter genotypes: AA, AG, TT, etc.
    // - Indels: D, I, DD, DI, II
    // - Missing: --, 00
    // - Single letters for haploid: A, T, G, C

    if genotype.is_empty() {
        return false;
    }

    let g = genotype.to_uppercase();

    // Missing data
    if g == "--" || g == "00" {
        return true;
    }

    // Check if all characters are valid nucleotides or indel markers
    for c in g.chars() {
        if !matches!(c, 'A' | 'T' | 'G' | 'C' | 'D' | 'I' | '-' | '0') {
            return false;
        }
    }

    // Valid length (1-2 characters)
    matches!(g.len(), 1 | 2)
}

fn split_fields(line: &str) -> Vec<String> {
    // Strip inline comments (anything after #)
    let line_without_comment = if let Some(comment_pos) = line.find('#') {
        &line[..comment_pos]
    } else {
        line
    };

    let delimiter = if line_without_comment.contains('\t') {
        '\t'
    } else if line_without_comment.contains(',') {
        ','
    } else {
        '\t'
    };

    line_without_comment
        .split(delimiter)
        .map(|part| {
            part.trim()
                .trim_matches('\"')
                .trim_matches('\r')
                .to_string()
        })
        .collect()
}

fn is_header_fields(fields: &[String]) -> bool {
    fields
        .first()
        .map(|value| {
            let lower = value.trim().to_lowercase();
            matches!(lower.as_str(), "rsid" | "name" | "snp" | "marker" | "id")
        })
        .unwrap_or(false)
}

fn identify_data_columns(fields: &[String]) -> Option<(usize, usize, usize)> {
    if fields.len() < 4 || is_header_fields(fields) {
        return None;
    }

    let rsid = fields[0].trim();
    if !rsid.starts_with("rs") && !rsid.starts_with('i') {
        return None;
    }

    for chr_idx in 1..fields.len() {
        if !is_valid_chromosome(fields[chr_idx].trim()) {
            continue;
        }

        for pos_idx in (chr_idx + 1)..fields.len() {
            if fields[pos_idx].trim().parse::<u64>().is_err() {
                continue;
            }

            for (genotype_idx, genotype_value) in fields.iter().enumerate().skip(pos_idx + 1) {
                if is_valid_genotype(genotype_value.trim()) {
                    return Some((chr_idx, pos_idx, genotype_idx));
                }
            }
        }
    }

    None
}

/// Detect source from header text (case-insensitive)
fn detect_source(header_lower: &str, data_lines: &[String]) -> Option<String> {
    if header_lower.contains("23andme") {
        Some("23andMe".to_string())
    } else if header_lower.contains("ancestrydna") || header_lower.contains("ancestry dna") {
        Some("AncestryDNA".to_string())
    } else if header_lower.contains("genes for good") || header_lower.contains("genesforgood") {
        Some("Genes for Good".to_string())
    } else if header_lower.contains("dynamic dna")
        || header_lower.contains("ddna")
        || header_lower.contains("dynamicdnalabs")
    {
        Some("Dynamic DNA".to_string())
    } else if header_lower.contains("living dna") {
        Some("Living DNA".to_string())
    } else if header_lower.contains("myheritage") {
        Some("MyHeritage".to_string())
    } else if header_lower.contains("decodeme")
        || data_lines.iter().any(|line| is_decodeme_header(line))
    {
        Some("deCODEme".to_string())
    } else if data_lines.iter().any(|line| is_familytreedna_header(line)) {
        Some("FamilyTreeDNA".to_string())
    } else {
        Some("Unknown".to_string())
    }
}

fn is_familytreedna_header(line: &str) -> bool {
    let parts = split_fields(line);
    if parts.len() < 4 {
        return false;
    }

    parts[0].eq_ignore_ascii_case("rsid")
        && parts[1].eq_ignore_ascii_case("chromosome")
        && parts[2].eq_ignore_ascii_case("position")
        && (parts[3].eq_ignore_ascii_case("result") || parts[3].eq_ignore_ascii_case("result1"))
}

fn is_decodeme_header(line: &str) -> bool {
    let parts = split_fields(line);
    if parts.len() < 6 {
        return false;
    }

    parts[0].eq_ignore_ascii_case("name")
        && parts[1].eq_ignore_ascii_case("variation")
        && parts[2].eq_ignore_ascii_case("chromosome")
        && parts[3].eq_ignore_ascii_case("position")
        && parts[5].eq_ignore_ascii_case("yourcode")
}

/// Detect GRCh version from header text (case-insensitive)
fn detect_grch_version(header_lower: &str) -> Option<String> {
    // Build 36 / GRCh36 / hg18
    if header_lower.contains("build 36")
        || header_lower.contains("grch36")
        || header_lower.contains("hg18")
    {
        return Some("36".to_string());
    }

    // Build 37 / GRCh37 / hg19
    if header_lower.contains("build 37")
        || header_lower.contains("grch37")
        || header_lower.contains("hg19")
    {
        return Some("37".to_string());
    }

    // Build 38 / GRCh38 / hg38
    if header_lower.contains("build 38")
        || header_lower.contains("grch38")
        || header_lower.contains("hg38")
    {
        return Some("38".to_string());
    }

    Some("Unknown".to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valid_chromosome() {
        assert!(is_valid_chromosome("1"));
        assert!(is_valid_chromosome("22"));
        assert!(is_valid_chromosome("X"));
        assert!(is_valid_chromosome("Y"));
        assert!(is_valid_chromosome("MT"));
        assert!(!is_valid_chromosome("0"));
        assert!(!is_valid_chromosome("23"));
        assert!(!is_valid_chromosome("ABC"));
    }

    #[test]
    fn test_valid_genotype() {
        assert!(is_valid_genotype("AA"));
        assert!(is_valid_genotype("AG"));
        assert!(is_valid_genotype("TT"));
        assert!(is_valid_genotype("--"));
        assert!(is_valid_genotype("00"));
        assert!(is_valid_genotype("D"));
        assert!(is_valid_genotype("I"));
        assert!(is_valid_genotype("DD"));
        assert!(is_valid_genotype("A"));
        assert!(!is_valid_genotype(""));
        assert!(!is_valid_genotype("AAA"));
        assert!(!is_valid_genotype("XY"));
    }

    #[test]
    fn test_detect_source() {
        assert_eq!(
            detect_source("this is from 23andme", &[]),
            Some("23andMe".to_string())
        );
        assert_eq!(
            detect_source("ancestrydna test", &[]),
            Some("AncestryDNA".to_string())
        );
        assert_eq!(
            detect_source("genes for good data", &[]),
            Some("Genes for Good".to_string())
        );
        assert_eq!(
            detect_source(
                "# this data file generated by dynamic dna (ddna) laboratories",
                &[]
            ),
            Some("Dynamic DNA".to_string())
        );
        assert_eq!(
            detect_source("data from ddna", &[]),
            Some("Dynamic DNA".to_string())
        );
        assert_eq!(
            detect_source("https://dynamicdnalabs.com", &[]),
            Some("Dynamic DNA".to_string())
        );
        assert_eq!(
            detect_source("# myheritage dna raw data", &[]),
            Some("MyHeritage".to_string())
        );
        let decodeme_header = "Name,Variation,Chromosome,Position,Strand,YourCode".to_string();
        assert_eq!(
            detect_source("", &[decodeme_header]),
            Some("deCODEme".to_string())
        );
        assert_eq!(
            detect_source("unknown source", &[]),
            Some("Unknown".to_string())
        );
    }

    #[test]
    fn test_detect_source_familytree_via_header_line() {
        let header = "RSID,CHROMOSOME,POSITION,RESULT".to_string();
        assert_eq!(
            detect_source("", &[header]),
            Some("FamilyTreeDNA".to_string())
        );
    }

    #[test]
    fn test_detect_grch_version() {
        assert_eq!(detect_grch_version("build 36"), Some("36".to_string()));
        assert_eq!(detect_grch_version("grch37"), Some("37".to_string()));
        assert_eq!(detect_grch_version("hg38"), Some("38".to_string()));
        assert_eq!(detect_grch_version("build 38"), Some("38".to_string()));
        assert_eq!(detect_grch_version("grch38"), Some("38".to_string()));
        assert_eq!(
            detect_grch_version("chromosomal location realtive to build 38 of the human reference"),
            Some("38".to_string())
        );
        assert_eq!(
            detect_grch_version("no version info"),
            Some("Unknown".to_string())
        );
    }
}
