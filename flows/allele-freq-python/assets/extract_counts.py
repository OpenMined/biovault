#!/usr/bin/env python3
import argparse
import gzip
import sys
import os


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract SNP genotype class counts from a VCF.")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--participant", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


def open_vcf(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def is_snp(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1


def log(msg: str, enabled: bool) -> None:
    if enabled:
        print(msg, file=sys.stderr)


def dosages_for_gt(alts, gt):
    if gt in (".", "./.", ".|."):
        return None
    gt_clean = gt.replace("|", "/")
    alleles = gt_clean.split("/")
    if len(alleles) != 2 or "." in alleles:
        return None
    try:
        a = int(alleles[0])
        b = int(alleles[1])
    except ValueError:
        return None
    if a > len(alts) or b > len(alts):
        return None
    alt_dosages = [0] * len(alts)
    allele_counts = {}
    allele_counts[a] = allele_counts.get(a, 0) + 1
    allele_counts[b] = allele_counts.get(b, 0) + 1
    for idx in range(1, len(alts) + 1):
        alt_dosages[idx - 1] = allele_counts.get(idx, 0)
    return alt_dosages


def main() -> int:
    args = parse_args()
    debug = args.debug or ("ALLELE_FREQ_DEBUG" in os.environ)
    try:
        with open(args.output, "w", encoding="utf-8") as out:
            out.write("locus_key\trsid\tparticipant_id\tdosage\n")
            with open_vcf(args.vcf) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 10:
                        continue

                    chrom, pos, rsid, ref, alt_field = parts[0:5]
                    sample = parts[9]

                # Only SNPs
                    alts = alt_field.split(",")
                    if any(not is_snp(ref, alt) for alt in alts):
                        continue

                    gt = sample.split(":")[0]
                    alt_keys = [f"{chrom}-{pos}-{ref}-{alt}" for alt in alts]
                    rsid_val = "" if rsid in (".", "") else rsid

                    if not alt_keys:
                        continue

                # Missing genotype: mark -1 for all rows
                    alt_dosages = dosages_for_gt(alts, gt)
                    if alt_dosages is None:
                        for locus in alt_keys:
                            log(f"{locus} missing=-1", debug)
                            out.write(f"{locus}\t{rsid_val}\t{args.participant}\t-1\n")
                        continue

                    for alt_key, dosage in zip(alt_keys, alt_dosages):
                        log(f"{alt_key} dosage={dosage}", debug)
                        out.write(f"{alt_key}\t{rsid_val}\t{args.participant}\t{dosage}\n")

    except FileNotFoundError:
        print(f"Missing VCF: {args.vcf}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
