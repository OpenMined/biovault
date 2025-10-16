use biovault::data::detect_genotype_metadata;
use std::path::PathBuf;

struct SampleExpectation {
    file: &'static str,
    source: &'static str,
    grch: &'static str,
}

fn sample_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/genotype_files")
        .join(name)
}

#[test]
fn detect_genotype_metadata_from_real_world_headers() {
    let expectations = [
        SampleExpectation {
            file: "23andme_grch36_standard_4col.txt",
            source: "23andMe",
            grch: "36",
        },
        SampleExpectation {
            file: "23andme_grch37_standard_4col.txt",
            source: "23andMe",
            grch: "37",
        },
        SampleExpectation {
            file: "23andme_standard_4col.txt",
            source: "23andMe",
            grch: "Unknown",
        },
        SampleExpectation {
            file: "ancestrydna_grch37_v1_split_alleles_5col.txt",
            source: "AncestryDNA",
            grch: "37",
        },
        SampleExpectation {
            file: "ancestrydna_split_alleles_5col.txt",
            source: "AncestryDNA",
            grch: "Unknown",
        },
        SampleExpectation {
            file: "dynamicdna_grch38_extended_7col.txt",
            source: "Dynamic DNA",
            grch: "38",
        },
        SampleExpectation {
            file: "familytreedna_grch37_standard_4col.csv",
            source: "FamilyTreeDNA",
            grch: "37",
        },
        SampleExpectation {
            file: "genesforgood_standard_4col.txt",
            source: "Genes for Good",
            grch: "37",
        },
        SampleExpectation {
            file: "livingdna_grch37_standard_4col.txt",
            source: "Living DNA",
            grch: "37",
        },
        SampleExpectation {
            file: "myheritage_grch37_standard_4col.csv",
            source: "MyHeritage",
            grch: "37",
        },
        SampleExpectation {
            file: "decodeme_grch36_variation_6col.csv",
            source: "deCODEme",
            grch: "36",
        },
        SampleExpectation {
            file: "unknown_standard_4col.txt",
            source: "Unknown",
            grch: "Unknown",
        },
        SampleExpectation {
            file: "unknown_no_header.txt",
            source: "Unknown",
            grch: "Unknown",
        },
    ];

    for expectation in expectations {
        let path = sample_path(expectation.file);
        let path_str = path.to_str().expect("path to string");
        let metadata = detect_genotype_metadata(path_str).expect("metadata detection");

        assert_eq!(
            metadata.data_type, "Genotype",
            "{} should be detected as genotype",
            expectation.file
        );

        assert_eq!(
            metadata.source.as_deref(),
            Some(expectation.source),
            "{} source",
            expectation.file
        );

        assert_eq!(
            metadata.grch_version.as_deref(),
            Some(expectation.grch),
            "{} GRCh version",
            expectation.file
        );
    }
}
