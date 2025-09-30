import os
import sys

from bioscript.classifier import (
    Classifier,
    DiploidAll,
    DiploidResult,
    DiploidSite,
    GenotypeEnum,
    PhaseMode,
)
from bioscript.reader import load_variants_tsv
from bioscript.types import Alleles, MatchList, MatchType, VariantCall

rs73885319 = VariantCall(rsid="rs73885319", ref=Alleles.A, alt=Alleles.NOT_A)
rs60910145 = VariantCall(rsid="rs60910145", ref=Alleles.T, alt=Alleles.NOT_T)
rs71785313 = VariantCall(
    rsid=["rs71785313", "rs1317778148", "rs143830837"], ref=Alleles.I, alt=Alleles.D
)


# 1) Severity-ordered genotypes (top = most severe)
class APOL1Genotypes(GenotypeEnum):
    G2 = "G2"
    G1 = "G1"
    G0 = "G0"


MISSING = "G-"

g1 = DiploidAll(
    sites=(rs73885319, rs60910145),
    genotype=APOL1Genotypes.G1,
    phase=PhaseMode.CIS,
    unknown_is_none=True,
)

g2 = DiploidSite(
    site=rs71785313,
    genotype=APOL1Genotypes.G2,
    unknown_is_none=True,
)


class APOL1(Classifier):
    def classify(self, matches) -> DiploidResult:
        # Check if any APOL1-related SNPs were found in any match category
        # (reference_matches, variant_matches, or no_call_matches)
        has_apol1_data = len(matches.all_matches) > 0 and any(
            match.variant_call.rsid.matches(rs73885319.rsid)
            or match.variant_call.rsid.matches(rs60910145.rsid)
            or match.variant_call.rsid.matches(rs71785313.rsid)
            for match in matches.all_matches
        )

        # If no APOL1 SNPs found, return missing
        if not has_apol1_data:
            return DiploidResult(MISSING, MISSING)

        g2_count = g2.count(matches)
        g1_count = g1.count(matches)

        # Special case: If there's exactly 1 deletion (het ID) and G1 variants present,
        # the deletion masks G1 on the same haplotype
        if g2_count == 1 and g1_count >= 1:
            # Check if the G2 deletion is heterozygous
            g2_match = None
            for match in matches.all_matches:
                if match.variant_call.rsid.matches(rs71785313.rsid):
                    if match.match_type == MatchType.VARIANT_CALL:
                        g2_match = match
                        break

            if g2_match and g2_match.snp.is_heterozygous():
                # Heterozygous deletion masks one copy of G1
                if g1_count == 2:
                    # Homozygous G1 (GG/GG) + het deletion (ID)
                    # -> One haplotype has G1+G2 (deletion wins) = G2
                    # -> Other haplotype has G1 only = G1
                    # Result: G2/G1 (already correct, don't modify)
                    pass
                elif g1_count == 1:
                    # Heterozygous G1 (het at both sites) + het deletion (ID)
                    # -> Assume CIS: G1 variants on same haplotype as deletion
                    # -> Deletion masks the G1, leaving G2/G0
                    g1_count = 0

        match (g2_count, g1_count):
            case (2, _):
                return DiploidResult(APOL1Genotypes.G2, APOL1Genotypes.G2)
            case (1, y) if y >= 1:
                return DiploidResult(APOL1Genotypes.G2, APOL1Genotypes.G1)
            case (1, _):
                return DiploidResult(APOL1Genotypes.G2, APOL1Genotypes.G0)
            case (0, 2):
                return DiploidResult(APOL1Genotypes.G1, APOL1Genotypes.G1)
            case (0, 1):
                return DiploidResult(APOL1Genotypes.G1, APOL1Genotypes.G0)
            case (0, 0):
                return DiploidResult(APOL1Genotypes.G0, APOL1Genotypes.G0)
            case _:
                return DiploidResult(MISSING, MISSING)


def main(file_path):
    # Retrieve the participant ID from the environment variable
    participant_id = os.environ.get("APOL1_PARTICIPANT_ID", "unknown_participant")

    # Load the variant data from the provided file path
    rows = load_variants_tsv(file_path)
    calls = MatchList(variant_calls=[rs73885319, rs60910145, rs71785313])
    matches = calls.match_rows(rows)

    # Classify the matches using the APOL1 classifier
    apol1 = APOL1()
    result = apol1.classify(matches)

    # Write the genotype information to a file
    output_file = f"apol1_{participant_id}.genotype.txt"
    with open(output_file, "w") as fh:
        fh.write(f"APOL1 genotype: {result}\n")

    print(f"Genotype information written to {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python bs_classifier.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    main(file_path)
