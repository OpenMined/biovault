#!/usr/bin/env python3
import os
import sys


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: validate_outputs.py <results_dir>", file=sys.stderr)
        return 1

    results_dir = sys.argv[1]
    module_dir = os.path.join(results_dir, "allele_freq")

    matrix_tsv = os.path.join(results_dir, "dosage_matrix.tsv")
    npz = os.path.join(results_dir, "dosage_matrix.npz")
    loci = os.path.join(results_dir, "locus_index.txt")
    participants = os.path.join(results_dir, "participants.txt")
    vcf = os.path.join(results_dir, "vcf_conversion_results.tsv")
    af = os.path.join(results_dir, "allele_freq.tsv")

    if not os.path.exists(matrix_tsv) and os.path.isdir(module_dir):
        matrix_tsv = os.path.join(module_dir, "dosage_matrix.tsv")
        npz = os.path.join(module_dir, "dosage_matrix.npz")
        loci = os.path.join(module_dir, "locus_index.txt")
        participants = os.path.join(module_dir, "participants.txt")
        vcf = os.path.join(module_dir, "vcf_conversion_results.tsv")
        af = os.path.join(module_dir, "allele_freq.tsv")

    for path in (matrix_tsv, npz, loci, participants, vcf, af):
        if not os.path.exists(path):
            print(f"Missing output: {path}", file=sys.stderr)
            return 1

    with open(loci, encoding="utf-8") as fh:
        lines = [line.rstrip("\n") for line in fh]

    if not lines or not lines[0].startswith("#format=loci-v1"):
        print("Invalid loci header", file=sys.stderr)
        return 1

    try:
        n_loci = int([line for line in lines if line.startswith("#n_loci=")][0].split("=", 1)[1])
    except Exception:
        print("Missing #n_loci header", file=sys.stderr)
        return 1

    if "locus_key" not in lines:
        print("Missing locus_key header", file=sys.stderr)
        return 1

    first_locus_idx = lines.index("locus_key") + 1
    locus_rows = lines[first_locus_idx:]
    if len(locus_rows) != n_loci:
        print(f"n_loci header {n_loci} does not match rows {len(locus_rows)}", file=sys.stderr)
        return 1

    with open(participants, encoding="utf-8") as fh:
        part_lines = [line.strip() for line in fh if line.strip()]
    if not part_lines or part_lines[0] != "participant_id":
        print("Invalid participants header", file=sys.stderr)
        return 1
    participants_list = part_lines[1:]

    with open(matrix_tsv, encoding="utf-8") as fh:
        header = next(fh, "").strip().split("\t")
        if len(header) < 2 or header[0] != "locus_key":
            print("Invalid dosage matrix header", file=sys.stderr)
            return 1
        header_participants = header[1:]
        if header_participants != participants_list:
            print("Participants header does not match participants.txt", file=sys.stderr)
            return 1
        matrix_rows = [line.strip().split("\t") for line in fh if line.strip()]

    loci_in_matrix = {row[0]: row[1:] for row in matrix_rows}

    # Verify allele_freq.tsv is consistent with dosage matrix
    with open(af, encoding="utf-8") as fh:
        af_lines = [line.strip().split("\t") for line in fh if line.strip()]
    if not af_lines:
        print("allele_freq.tsv is empty", file=sys.stderr)
        return 1
    af_header = af_lines[0]
    if af_header[:5] != ["locus_key", "allele_count", "allele_number", "num_homo", "allele_freq"]:
        print("Invalid allele_freq.tsv header", file=sys.stderr)
        return 1
    for row in af_lines[1:]:
        if len(row) < 5:
            print("Invalid allele_freq.tsv row", file=sys.stderr)
            return 1
        locus = row[0]
        if locus not in loci_in_matrix:
            print(f"Allele freq locus missing from matrix: {locus}", file=sys.stderr)
            return 1
        vals = [int(v) for v in loci_in_matrix[locus]]
        observed = [v for v in vals if v != -1]
        n_obs = len(observed)
        allele_number = 2 * n_obs
        allele_count = sum(observed)
        num_homo = sum(1 for v in observed if v == 2)
        af_calc = allele_count / allele_number if allele_number > 0 else 0.0
        try:
            af_file = float(row[4])
        except ValueError:
            print(f"Invalid allele_freq value for {locus}", file=sys.stderr)
            return 1
        if row[1] != str(allele_count) or row[2] != str(allele_number) or row[3] != str(num_homo):
            print(f"Allele freq counts mismatch for {locus}", file=sys.stderr)
            return 1
        if abs(af_calc - af_file) > 1e-6:
            print(f"Allele freq mismatch for {locus}: {af_calc} != {af_file}", file=sys.stderr)
            return 1

    apol1_loci = [
        "22-36265988-T-C",
        "22-36265988-T-G",
        "22-36265860-A-G",
    ]
    thal_loci = [
        "11-5225485-T-C",
        "11-5225660-G-A",
        "11-5225660-G-C",
    ]

    apol1_present = [locus for locus in apol1_loci if locus in loci_in_matrix]
    thal_present = [locus for locus in thal_loci if locus in loci_in_matrix]
    if not apol1_present:
        print("No APOL1 loci found in matrix", file=sys.stderr)
        return 1
    if not thal_present:
        print("No THAL loci found in matrix", file=sys.stderr)
        return 1

    def any_nonzero(loci):
        for locus in loci:
            vals = [int(x) for x in loci_in_matrix[locus]]
            if any(v > 0 for v in vals):
                return True
        return False

    if not any_nonzero(apol1_present):
        print("APOL1 override loci have no non-zero dosages", file=sys.stderr)
        return 1
    if not any_nonzero(thal_present):
        print("THAL override loci have no non-zero dosages", file=sys.stderr)
        return 1

    no_call_freq = float(os.environ.get("ALLELE_FREQ_NO_CALL_FREQ", "0.0") or "0.0")
    if no_call_freq >= 0.1:
        def any_missing(loci):
            for locus in loci:
                vals = [int(x) for x in loci_in_matrix[locus]]
                if any(v == -1 for v in vals):
                    return True
            return False
        if not any_missing(apol1_present + thal_present):
            print("Expected missing (-1) dosages not found for override loci", file=sys.stderr)
            return 1

    try:
        import numpy as np
        arr = np.load(npz)
        if "dosage_mat" not in arr:
            print("Missing dosage_mat in npz", file=sys.stderr)
            return 1
        if arr["dosage_mat"].shape[0] != n_loci:
            print(f"NPZ dosage_mat length {arr['dosage_mat'].shape[0]} != n_loci {n_loci}", file=sys.stderr)
            return 1
        if arr["dosage_mat"].shape[1] != len(participants_list):
            print(f"NPZ dosage_mat cols {arr['dosage_mat'].shape[1]} != participants {len(participants_list)}", file=sys.stderr)
            return 1
    except Exception as exc:
        print(f"Warning: skipping npz validation ({exc})")

    print("✓ allele-freq outputs validated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
