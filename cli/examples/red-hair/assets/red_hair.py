#!/usr/bin/env python3
import argparse
import sys

DESC = """Interpret MC1R rs1805007 (R151C) for red-hair tendency.
Input line format from bcftools query:
  REF<TAB>ALT<TAB>GT<TAB>AD
Example:
  C\tT\t0/1\t18,12
"""


def interpret(line: str, participant: str) -> str:
    parts = line.strip().split("\t")
    if len(parts) < 4:
        return f"No valid variant information found for {participant}."

    ref, alt, gt, ad = parts[0], parts[1], parts[2], parts[3]

    # map GT (0/0,0/1,1/1 or phased) to allele letters using REF/ALT (assume single-alt site)
    gt_simple = gt.replace("|", "/")
    a = gt_simple.split("/")
    alleles = []
    for x in a:
        if x == "0":
            alleles.append(ref)
        elif x == "1":
            alleles.append(alt.split(",")[0])
        else:
            alleles.append("?")
    allele_pair = "(" + ";".join(alleles) + ")"

    # counts string
    ad_parts = ad.split(",")
    counts_str = (
        f"Counts: ({', '.join(ad_parts)})"
        if len(ad_parts) > 1
        else f"Counts: {ad_parts[0]}"
    )

    # rs1805007 is C>T; T is the red-hair–associated allele
    # simple, single-locus heuristic
    if gt_simple in ("0/0",):
        interp = "No rs1805007 variant detected at this site (CC). Not predictive for red hair on its own."
    elif gt_simple in ("0/1", "1/0"):
        interp = "Carrier of rs1805007 T allele (CT). Tends toward fair skin/freckling; red hair possible with other MC1R variants."
    elif gt_simple in ("1/1",):
        interp = "Homozygous rs1805007 T (TT). Strongly associated with red hair and freckling."
    else:
        interp = f"Genotype {gt} not recognized; check depth/quality."

    return f"""MC1R rs1805007 (R151C) — red-hair–associated variant
Participant: {participant}
Site: chr16 (build-dependent; see workflow notes)
Genotype: {allele_pair}  {counts_str}
Interpretation: {interp}
"""


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=DESC)
    p.add_argument("--participant", required=True)
    args = p.parse_args()
    line = sys.stdin.readline()
    print(interpret(line, args.participant))
