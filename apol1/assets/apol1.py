#!/usr/bin/env python3
import sys


def interpret_eye_color(line: str) -> str:
    """
    Input line format from bcftools query:
      REF<TAB>ALT<TAB>GT<TAB>AD
    Example:
      G	A	0/1	10,8
    """

    parts = line.strip().split("\t")
    if len(parts) < 4:
        return "No valid variant information found."

    ref, alt, gt, ad = parts[0], parts[1], parts[2], parts[3]

    # Determine the alleles based on the genotype
    alleles = ";".join(
        ref if a == "0" else alt if a == "1" else "?"
        for a in gt.replace("|", "/").split("/")
    )

    # Format the counts
    ad_parts = ad.split(",")
    if len(ad_parts) == 1:
        counts_str = f"Counts: {ad_parts[0]}"
    else:
        counts_str = f"Counts: ({', '.join(ad_parts)})"

    # Determine the interpretation
    if gt in ("0/0", "0|0"):  # homozygous reference (e.g., G/G)
        interpretation = "Eye color tends toward brown."
    elif gt in ("0/1", "1/0", "0|1", "1|0"):  # heterozygous (e.g., G/A)
        interpretation = "Eye color tends toward brown."
    elif gt in ("1/1", "1|1"):  # homozygous alternate (e.g., A/A)
        interpretation = "Eye color is likely blue."
    else:
        interpretation = f"Genotype {gt} not recognized"

    result = f"""
The rs12913832 variant has the following alleles:
(A;A) yields brown eye color ~80% of the time.
(A;G) also tends toward brown.
(G;G) gives blue eye color ~99% of the time.

This person has:
Genotype: ({alleles})
{counts_str}
Interpretation: {interpretation}
"""

    return result


if __name__ == "__main__":
    # Read line from stdin
    line = sys.stdin.readline()
    print(line)
    # result = interpret_eye_color(line)
    # print(result)
