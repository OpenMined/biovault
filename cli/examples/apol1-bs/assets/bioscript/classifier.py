from enum import Enum

from .types import MatchType, Tuple


class GenotypeEnum(Enum):
    pass


class PhaseMode(Enum):
    CIS = "cis"
    ANY = "any"


class DiploidAll:
    def __init__(self, sites, genotype, phase, unknown_is_none):
        self.sites = sites
        self.genotype = genotype
        self.phase = phase
        self.unknown_is_none = unknown_is_none

    def count(self, matches):
        # Collect matches for each required site (use index as key)
        site_matches = {i: [] for i in range(len(self.sites))}

        for match in matches:
            for i, site in enumerate(self.sites):
                if match.variant_call.rsid.matches(site.rsid):
                    if match.match_type == MatchType.VARIANT_CALL:
                        site_matches[i].append(match)

        # Check if all required sites have at least one variant call
        if not all(site_matches[i] for i in range(len(self.sites))):
            return 0

        # For CIS phase, count how many chromosome copies have ALL variants
        # For simplicity, if all sites are heterozygous, assume 1 copy
        # If any site is homozygous for the variant, that's 2 copies (but rare for G1)
        if self.phase == PhaseMode.CIS:
            # Check if all are heterozygous -> 1 haplotype
            # Check if all are homozygous -> 2 haplotypes (both copies have the variant)
            all_heterozygous = all(
                match.snp.is_heterozygous()
                for matches_list in site_matches.values()
                for match in matches_list
            )
            all_homozygous = all(
                match.snp.is_homozygous()
                for matches_list in site_matches.values()
                for match in matches_list
            )

            if all_homozygous:
                return 2  # Both haplotypes have all variants
            elif all_heterozygous:
                return 1  # One haplotype has all variants
            else:
                return 1  # Mixed: assume at least one haplotype has all
        else:  # PhaseMode.ANY
            return 1  # Just need all sites present as variants

    def as_pair(self, matches, missing):
        # Placeholder for actual logic to project to a pair
        return (self.genotype, missing)


class DiploidSite:
    def __init__(self, site, genotype, unknown_is_none):
        self.site = site
        self.genotype = genotype
        self.unknown_is_none = unknown_is_none

    def count(self, matches):
        # Collect matches for this site
        # Since VariantCall can have multiple rsID aliases, we may get duplicate matches
        # for the same physical variant - deduplicate by taking first match only
        site_matches = []
        seen_genotypes = set()

        for match in matches:
            if match.variant_call.rsid.matches(self.site.rsid):
                if match.match_type == MatchType.VARIANT_CALL:
                    # Deduplicate: if we've seen this exact genotype for this site, skip it
                    genotype_key = tuple(match.snp)
                    if genotype_key not in seen_genotypes:
                        site_matches.append(match)
                        seen_genotypes.add(genotype_key)
                    else:
                        # Multiple rsID aliases for same variant
                        print(f"Duplicate match for rsID aliases (genotype: {''.join(str(n.value) for n in match.snp)}) - using first match only", flush=True)

        # Count alleles across deduplicated matches
        # Should only be one match after deduplication
        if len(site_matches) > 1:
            print(f"WARNING: Multiple distinct genotypes found for site {self.site.rsid} - using first", flush=True)

        count = 0
        for match in site_matches[:1]:  # Take only first match
            # Count how many alleles match the variant (alt allele)
            # For homozygous variant (e.g., DD), count = 2
            # For heterozygous variant (e.g., ID), count = 1
            if match.snp[0] in match.variant_call.alt:
                count += 1
            if match.snp[1] in match.variant_call.alt:
                count += 1

        return count


# 3) Switch/match classifier
class Classifier:
    def classify(self, matches) -> Tuple[GenotypeEnum, GenotypeEnum]:
        raise NotImplementedError("Subclasses should implement this method.")


class DiploidResult:
    def __init__(self, genotype1: GenotypeEnum, genotype2: GenotypeEnum):
        self.genotype1 = genotype1
        self.genotype2 = genotype2

    def __str__(self):
        """
        Provides a pretty-printed string representation of the DiploidResult.
        Displays the two genotypes in the format VAL1/VAL2.
        If both genotypes are the same, only display one.
        """
        g1 = self.genotype1.value if hasattr(self.genotype1, 'value') else self.genotype1
        g2 = self.genotype2.value if hasattr(self.genotype2, 'value') else self.genotype2
        return f"{g1}/{g2}"

    def __repr__(self):
        return self.__str__()
