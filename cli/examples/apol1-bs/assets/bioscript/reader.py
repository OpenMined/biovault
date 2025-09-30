import csv
from typing import Iterator, List, Optional

from .types import VariantRow

REQUIRED = ("rsid", "chromosome", "position", "genotype")
OPTIONAL = ("gs", "baf", "lrr")


def _float_or_none(x: Optional[str]) -> Optional[float]:
    if x is None:
        return None
    x = x.strip()
    if not x:
        return None
    try:
        return float(x)
    except ValueError:
        return None


def _int_or_none(x: Optional[str]) -> Optional[int]:
    if x is None:
        return None
    x = x.strip()
    if not x:
        return None
    try:
        return int(x)
    except ValueError:
        return None


def _extract_header_and_data(lines: List[str]) -> (List[str], Iterator[str]):
    """
    Returns (header_fields, data_lines_iterator).
    Header is the last commented line starting with '#', or the first line if no commented header.
    """
    # Normalize newlines and drop blank lines
    lines = [ln.rstrip("\n\r") for ln in lines if ln.strip()]

    # Find the last commented line (potential commented header)
    commented_header_idx = None
    for i, ln in enumerate(lines):
        if ln.lstrip().startswith("#"):
            commented_header_idx = i

    header_fields = None
    data_start_idx = None
    if commented_header_idx is not None:
        raw = lines[commented_header_idx].lstrip()
        # Strip leading '#' and whitespace
        if raw.startswith("#"):
            raw = raw[1:].lstrip()
        header_fields = next(csv.reader([raw], delimiter="\t"))
        # Data begins at the first non-comment line after the commented header
        data_start_idx = next(
            (
                i
                for i in range(commented_header_idx + 1, len(lines))
                if not lines[i].lstrip().startswith("#")
            ),
            None,
        )
    else:
        # Assume the first line is the header if no commented header is found
        header_fields = next(csv.reader([lines[0]], delimiter="\t"))
        data_start_idx = 1

    if not header_fields or data_start_idx is None:
        raise ValueError("No header found in TSV (considering commented headers).")

    # Prepare iterator over data lines (skip comments anywhere after header)
    def data_iter():
        for ln in lines[data_start_idx:]:
            if ln.lstrip().startswith("#"):
                continue
            yield ln

    return [h.strip() for h in header_fields], data_iter()


def load_variants_tsv(path: str) -> Iterator[VariantRow]:
    """
    Load a TSV with columns:
      required: rsid, chromosome, position, genotype
      optional: gs, baf, lrr
    Supports either a normal header line, or a commented header like:
      # rsid\tchromosome\tposition\tgenotype\tgs\tbaf\tlrr
    Ignores lines starting with '#'. Yields VariantRow objects.
    """
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        all_lines = f.readlines()

    header, data_lines = _extract_header_and_data(all_lines)

    # Validate required columns
    missing = [col for col in REQUIRED if col not in header]
    if missing:
        raise ValueError(f"Missing required column(s): {', '.join(missing)}")

    has_opt = {col: (col in header) for col in OPTIONAL}

    reader = csv.DictReader(data_lines, fieldnames=header, delimiter="\t")
    for row in reader:
        # Required fields
        rsid = (row.get("rsid") or "").strip()
        chrom = (row.get("chromosome") or "").strip()
        pos = _int_or_none(row.get("position"))
        genotype = (row.get("genotype") or "").strip()

        if not rsid or not chrom or pos is None or not genotype:
            raise ValueError(f"Invalid row (required fields empty/bad): {row}")

        # Optional fields
        gs = _float_or_none(row.get("gs")) if has_opt["gs"] else None
        baf = _float_or_none(row.get("baf")) if has_opt["baf"] else None
        lrr = _float_or_none(row.get("lrr")) if has_opt["lrr"] else None

        yield VariantRow(
            rsid=rsid,
            chromosome=chrom,
            position=pos,
            genotype=genotype,
            gs=gs,
            baf=baf,
            lrr=lrr,
        )
