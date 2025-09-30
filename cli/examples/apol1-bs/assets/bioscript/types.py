from dataclasses import dataclass, field
from enum import Enum
from typing import Iterable, Optional, Tuple, Union


class Nucleotide(str, Enum):
    A = "A"  # Adenine
    T = "T"  # Thymine (DNA only)
    U = "U"  # Uracil (RNA only)
    C = "C"  # Cytosine
    G = "G"  # Guanine
    MISSING = "."


class InDel(str, Enum):
    I = "I"  # Insertion
    D = "D"  # Deletion


class MatchType(Enum):
    REFERENCE_CALL = "Reference call"
    VARIANT_CALL = "Variant call"
    NO_CALL = "No call"


class AllelesMeta(type):
    @property
    def A(cls):
        return Alleles(Nucleotide.A)

    @property
    def T(cls):
        return Alleles(Nucleotide.T)

    @property
    def C(cls):
        return Alleles(Nucleotide.C)

    @property
    def G(cls):
        return Alleles(Nucleotide.G)

    @property
    def U(cls):
        return Alleles(Nucleotide.U)

    @property
    def I(cls):
        return Alleles(InDel.I)

    @property
    def D(cls):
        return Alleles(InDel.D)

    @property
    def NOT_A(cls):
        return Alleles(
            {n for n in Nucleotide if n != Nucleotide.A and n != Nucleotide.MISSING}
        )

    @property
    def NOT_T(cls):
        return Alleles(
            {n for n in Nucleotide if n != Nucleotide.T and n != Nucleotide.MISSING}
        )

    @property
    def NOT_U(cls):
        return Alleles(
            {n for n in Nucleotide if n != Nucleotide.U and n != Nucleotide.MISSING}
        )

    @property
    def NOT_C(cls):
        return Alleles(
            {n for n in Nucleotide if n != Nucleotide.C and n != Nucleotide.MISSING}
        )

    @property
    def NOT_G(cls):
        return Alleles(
            {n for n in Nucleotide if n != Nucleotide.G and n != Nucleotide.MISSING}
        )


class Alleles(metaclass=AllelesMeta):
    def __init__(
        self,
        nucleotides: Union[Nucleotide, InDel, Tuple[Union[Nucleotide, InDel], ...]],
    ):
        if isinstance(nucleotides, (Nucleotide, InDel)):
            self.nucleotides = {nucleotides}
        else:
            self.nucleotides = set(nucleotides)

    def __eq__(self, other):
        if isinstance(other, Alleles):
            return self.nucleotides == other.nucleotides
        return False

    def __contains__(self, nucleotide: Union[Nucleotide, InDel]) -> bool:
        return nucleotide in self.nucleotides

    def __len__(self) -> int:
        return len(self.nucleotides)

    def __iter__(self):
        return iter(self.nucleotides)

    def add(self, nucleotide: Union[Nucleotide, InDel]):
        self.nucleotides.add(nucleotide)

    def remove(self, nucleotide: Union[Nucleotide, InDel]):
        self.nucleotides.remove(nucleotide)


@dataclass
class SNP:
    ploidy: str


class DiploidSNP(tuple):
    def __new__(
        cls,
        nucleotide1: Union[Nucleotide, InDel],
        nucleotide2: Union[Nucleotide, InDel],
    ):
        return super().__new__(cls, (nucleotide1, nucleotide2))

    def count(self, nucleotide: Union[Nucleotide, InDel]) -> int:
        return super().count(nucleotide)

    def is_homozygous(self) -> bool:
        return self[0] == self[1]

    def is_heterozygous(self) -> bool:
        return self[0] != self[1]


@dataclass
class VariantMatch:
    variant_call: "VariantCall"
    snp: DiploidSNP
    match_type: MatchType

    def __str__(self):
        """
        Provides a pretty-printed string representation of the VariantMatch for easy reading.
        Displays relevant details of the match.
        """
        rsid = self.variant_call.rsid
        genotype = "".join(nuc.value for nuc in self.snp)
        match_type = self.match_type.name
        return f"{match_type}: {rsid} ref={self.snp[0].value} genotype: {genotype}"

    def __repr__(self):
        return self.__str__()


@dataclass
class RSID:
    aliases: set

    def __init__(self, *aliases: Union[str, "RSID"]):
        self.aliases = set()
        for alias in aliases:
            if isinstance(alias, RSID):
                self.aliases.update(alias.aliases)
            else:
                self.aliases.add(alias)

    def matches(self, rsid: Union[str, "RSID"]) -> bool:
        if isinstance(rsid, RSID):
            return not self.aliases.isdisjoint(rsid.aliases)
        return rsid in self.aliases


@dataclass
class VariantCall:
    rsid: Union[RSID, str, Iterable]
    ploidy: str = "diploid"
    ref: Alleles = field(default_factory=lambda: Alleles({Nucleotide.MISSING}))
    alt: Alleles = field(default_factory=lambda: Alleles({Nucleotide.MISSING}))

    def __post_init__(self):
        if isinstance(self.rsid, str):
            self.rsid = RSID(self.rsid)
        elif isinstance(self.rsid, Iterable) and not isinstance(self.rsid, RSID):
            self.rsid = RSID(*self.rsid)

    def filter_variant_row(self, variant_row: "VariantRow") -> Optional["VariantMatch"]:
        """
        Filters for matching VariantRow based on rsid.
        If rsids match, returns a VariantMatch object with a parsed DiploidSNP; otherwise, returns None.
        """
        if self.rsid.matches(variant_row.rsid):
            if isinstance(variant_row.genotype, str) and len(variant_row.genotype) == 2:
                allele1 = None
                allele2 = None

                # Try to parse first allele
                try:
                    allele1 = Nucleotide(variant_row.genotype[0])
                except ValueError:
                    try:
                        allele1 = InDel(variant_row.genotype[0])
                    except ValueError:
                        print(f"Invalid allele value '{variant_row.genotype[0]}' cast as MISSING '.'", flush=True)
                        allele1 = Nucleotide.MISSING

                # Try to parse second allele
                try:
                    allele2 = Nucleotide(variant_row.genotype[1])
                except ValueError:
                    try:
                        allele2 = InDel(variant_row.genotype[1])
                    except ValueError:
                        print(f"Invalid allele value '{variant_row.genotype[1]}' cast as MISSING '.'", flush=True)
                        allele2 = Nucleotide.MISSING

                diploid_snp = DiploidSNP(allele1, allele2)

                if diploid_snp.is_homozygous() and diploid_snp[0] in self.ref:
                    match_type = MatchType.REFERENCE_CALL
                elif (
                    diploid_snp.is_heterozygous()
                    or diploid_snp[0] in self.alt
                    or diploid_snp[1] in self.alt
                ):
                    match_type = MatchType.VARIANT_CALL
                else:
                    match_type = MatchType.NO_CALL
                return VariantMatch(
                    variant_call=self, snp=diploid_snp, match_type=match_type
                )
        return None


@dataclass
class VariantRow:
    rsid: str
    chromosome: str  # str to allow "X","Y","MT"
    position: int
    genotype: str  # Keep as str until needed
    ploidy: str = "diploid"
    gs: Optional[float] = None
    baf: Optional[float] = None
    lrr: Optional[float] = None

    def __post_init__(self):
        # No need to parse genotype here; it will be parsed in VariantCall if needed
        pass


@dataclass
class MatchList:
    variant_calls: Iterable["VariantCall"]
    reference_matches: list = field(default_factory=list)
    variant_matches: list = field(default_factory=list)
    no_call_matches: list = field(default_factory=list)
    all_matches: list = field(default_factory=list)

    def match_rows(self, variant_rows: Iterable["VariantRow"]) -> "MatchList":
        """
        Iterates through the given variant rows, applies each variant call's filter_variant_row method,
        and collects the outputs into three buckets: reference calls, variant calls, and no calls.
        Also maintains a list of all matches in the order they were added.
        """
        for variant_row in variant_rows:
            for variant_call in self.variant_calls:
                match = variant_call.filter_variant_row(variant_row)
                if match is not None:
                    self.all_matches.append(match)
                    if match.match_type == MatchType.REFERENCE_CALL:
                        self.reference_matches.append(match)
                    elif match.match_type == MatchType.VARIANT_CALL:
                        self.variant_matches.append(match)
                    elif match.match_type == MatchType.NO_CALL:
                        self.no_call_matches.append(match)

        return self

    def __iter__(self):
        return iter(self.all_matches)

    def __getitem__(self, index):
        return self.all_matches[index]

    def iter_reference(self):
        return iter(self.reference_matches)

    def iter_variant(self):
        return iter(self.variant_matches)

    def iter_no_call(self):
        return iter(self.no_call_matches)

    def __str__(self):
        """
        Provides a pretty-printed string representation of the MatchList for easy reading.
        Each match is displayed on its own line with relevant details.
        """
        return "\n".join(str(match) for match in self.all_matches)

    def __repr__(self):
        return self.__str__()
