#!/usr/bin/env python3
import argparse
from itertools import product
from pathlib import Path
from typing import List, Tuple, Optional

# -----------------------------------
# Targets (GRCh38)
# -----------------------------------
G1A = ("rs73885319", "22", "36265860")  # A>G (normally), may appear as A>C
G1B = ("rs60910145", "22", "36265988")  # T>G (normally), may appear as T>C
G2C = ("rs71785313", "22", "36266000")  # TTATAA > -
G2A1 = ("rs1317778148", "22", "36266000")  # merged into rs71785313
G2A2 = ("rs143830837", "22", "36266000")  # merged into rs71785313

IRRELEVANT_ROWS = [
    ("rs9701055", "1", "630053", "CC"),
    ("rs9651229", "1", "632287", "CC"),
    ("rs9701872", "1", "632828", "TT"),
]

INDEL_STATES = ["II", "ID", "DD"]  # 0/0, 0/1, 1/1

HEADER = ["rsid", "chromosome", "position", "genotype"]


# -----------------------------------
# Calling helpers (unphased)
# -----------------------------------
def g1_state_from_letters(ref: str, alt: str, gt: str) -> Optional[str]:
    if not gt or len(gt) != 2 or not gt.isalpha():
        return None
    a, b = gt[0].upper(), gt[1].upper()
    ref, alt = ref.upper(), alt.upper()
    s = {a, b}
    if not s.issubset({ref, alt}):
        return None
    if a == ref and b == ref:
        return "hom_ref"
    if a == alt and b == alt:
        return "hom_alt"
    return "het"


def g1_count_from_two(a: Optional[str], b: Optional[str]) -> Tuple[Optional[int], bool]:
    known = (a is not None) and (b is not None)
    if not known:
        return None, False
    if a == "hom_alt" and b == "hom_alt":
        return 2, True
    if (a in ("het", "hom_alt")) and (b in ("het", "hom_alt")):
        return 1, True
    return 0, True


def g2_count_from_indel(indel: Optional[str]) -> Tuple[Optional[int], bool]:
    if indel not in INDEL_STATES:
        return None, False
    if indel == "DD":
        return 2, True
    if indel == "ID":
        return 1, True
    return 0, True


def final_label(
    g1_count: Optional[int], g1_known: bool, g2_count: Optional[int], g2_known: bool
) -> str:
    if g1_known and g2_known:
        gc1, gc2 = g1_count or 0, g2_count or 0
        if gc2 == 2:
            return "G2G2"
        if gc1 == 2:
            return "G1G1"
        if gc1 >= 1 and gc2 >= 1:
            return "G1G2"
        if gc1 == 1:
            return "G1G0"
        if gc2 == 1:
            return "G2G0"
        return "G0G0"
    if g1_known and not g2_known:
        return "G1G-" if (g1_count or 0) >= 1 else "G0G-"
    if g2_known and not g1_known:
        return "G2G-" if (g2_count or 0) >= 1 else "G0G-"
    return "G-G-"


def derive_label(rows: List[Tuple[str, str, str, str]], a_alt: str, b_alt: str) -> str:
    by_id = {r[0]: r for r in rows}
    g1a = (
        g1_state_from_letters(
            "A", a_alt, by_id.get(G1A[0], (None, None, None, None))[3]
        )
        if G1A[0] in by_id
        else None
    )
    g1b = (
        g1_state_from_letters(
            "T", b_alt, by_id.get(G1B[0], (None, None, None, None))[3]
        )
        if G1B[0] in by_id
        else None
    )
    g2_gt = None
    for rid in (G2C[0], G2A1[0], G2A2[0]):
        if rid in by_id:
            g2_gt = by_id[rid][3]
            break
    g1_count, g1_known = g1_count_from_two(g1a, g1b)
    g2_count, g2_known = g2_count_from_indel(g2_gt)
    return final_label(g1_count, g1_known, g2_count, g2_known)


# -----------------------------------
# Writers
# -----------------------------------
def write_tsv(path: Path, rows: List[Tuple[str, str, str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(HEADER) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


def append_participant(rows_accum, rel_path: Path, label: str, idx: int):
    weight = 50 + (idx % 51)
    height = 150 + (idx % 51)
    age = 16 + (idx % 65)
    pid = f"PID{idx:05d}"
    # normalize label with slashes
    mapping = {
        "G0G0": "G0/G0",
        "G1G0": "G1/G0",
        "G2G0": "G2/G0",
        "G1G1": "G1/G1",
        "G2G2": "G2/G2",
        "G1G2": "G1/G2",
        "G0G-": "G0/G-",
        "G1G-": "G1/G-",
        "G2G-": "G2/G-",
        "G-G-": "G-/G-",
    }
    slash_label = mapping.get(label, label)
    rows_accum.append([pid, str(rel_path), weight, height, age, slash_label])


# -----------------------------------
# Row builders (first 4 cols only)
# -----------------------------------
def canon_rows(g1a_gt: str, g1b_gt: str, g2_indel: Optional[str]):
    rows = [(G1A[0], G1A[1], G1A[2], g1a_gt), (G1B[0], G1B[1], G1B[2], g1b_gt)]
    if g2_indel is not None:
        rows.append((G2C[0], G2C[1], G2C[2], g2_indel))
    return rows


def g2_alias_only_rows(g1a_gt: str, g1b_gt: str, g2_indel: str, which="a1"):
    alias = G2A1 if which == "a1" else G2A2
    return [
        (G1A[0], G1A[1], G1A[2], g1a_gt),
        (G1B[0], G1B[1], G1B[2], g1b_gt),
        (alias[0], alias[1], alias[2], g2_indel),
    ]


def g2_both_ids_rows(g1a_gt: str, g1b_gt: str, g2_indel: str):
    return [
        (G1A[0], G1A[1], G1A[2], g1a_gt),
        (G1B[0], G1B[1], G1B[2], g1b_gt),
        (G2C[0], G2C[1], G2C[2], g2_indel),
        (G2A1[0], G2A1[1], G2A1[2], g2_indel),
    ]


def rows_missing(
    which_missing: str,
    g1a_gt: Optional[str],
    g1b_gt: Optional[str],
    g2_indel: Optional[str],
):
    rows: List[Tuple[str, str, str, str]] = []
    if which_missing != "G1A" and g1a_gt is not None:
        rows.append((G1A[0], G1A[1], G1A[2], g1a_gt))
    if which_missing != "G1B" and g1b_gt is not None:
        rows.append((G1B[0], G1B[1], G1B[2], g1b_gt))
    if which_missing != "G2" and g2_indel is not None:
        rows.append((G2C[0], G2C[1], G2C[2], g2_indel))
    return rows


def rows_invalid_gt(
    break_which: str,
    g1a_gt: Optional[str],
    g1b_gt: Optional[str],
    g2_indel: Optional[str],
):
    bad_snp = "A?"
    bad_indel = "IDK"
    rows: List[Tuple[str, str, str, str]] = []
    rows.append(
        (G1A[0], G1A[1], G1A[2], bad_snp if break_which == "G1A" else (g1a_gt or "AA"))
    )
    rows.append(
        (G1B[0], G1B[1], G1B[2], bad_snp if break_which == "G1B" else (g1b_gt or "TT"))
    )
    if break_which == "G2":
        rows.append((G2C[0], G2C[1], G2C[2], bad_indel))
    elif g2_indel is not None:
        rows.append((G2C[0], G2C[1], G2C[2], g2_indel))
    return rows


def rows_irrelevant_only():
    return list(IRRELEVANT_ROWS)


# -----------------------------------
# Main
# -----------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Generate APOL1 genotype test TSVs + participants.csv"
    )
    ap.add_argument("--out", default="apol1_tests", help="Output directory")
    args = ap.parse_args()

    out_root = Path(args.out)
    participants = []
    idx = 1

    # ---- A) All "normal" cases once (ALT=G for both SNPs) ----
    a_alt = b_alt = "G"

    def gt_a(state: str) -> str:
        return "AA" if state == "hom_ref" else ("AG" if state == "het" else "GG")

    def gt_b(state: str) -> str:
        return "TT" if state == "hom_ref" else ("TG" if state == "het" else "GG")

    for s1, s2, indel in product(
        ["hom_ref", "het", "hom_alt"], ["hom_ref", "het", "hom_alt"], INDEL_STATES
    ):
        rows = canon_rows(gt_a(s1), gt_b(s2), indel)
        label = derive_label(rows, a_alt, b_alt)
        rel = Path(label) / f"case_{idx:04d}_canonical.tsv"
        write_tsv(out_root / rel, rows)
        append_participant(participants, rel, label, idx)
        idx += 1

    # ---- B) Exactly one exception per scenario ----
    # 1) G1A ALT=C (exercise strand/assay flip once)
    rows = canon_rows("AC", "TT", "II")  # should be G1 het only -> G1G0
    label = derive_label(rows, "C", "G")
    rel = Path(label) / f"case_{idx:04d}_g1a_alt_C.tsv"
    write_tsv(out_root / rel, rows)
    append_participant(participants, rel, label, idx)
    idx += 1

    # 2) G1B ALT=C
    rows = canon_rows("AA", "TC", "II")  # G1 het at B only -> G1G0
    label = derive_label(rows, "G", "C")
    rel = Path(label) / f"case_{idx:04d}_g1b_alt_C.tsv"
    write_tsv(out_root / rel, rows)
    append_participant(participants, rel, label, idx)
    idx += 1

    # 3) G2 alias-only (a1)
    rows = g2_alias_only_rows("AA", "TT", "ID", which="a1")  # single G2 -> G2G0
    label = derive_label(rows, "G", "G")
    rel = Path(label) / f"case_{idx:04d}_g2_alias_a1.tsv"
    write_tsv(out_root / rel, rows)
    append_participant(participants, rel, label, idx)
    idx += 1

    # 4) G2 alias-only (a2)
    rows = g2_alias_only_rows("AA", "TT", "DD", which="a2")  # two G2 -> G2G2
    label = derive_label(rows, "G", "G")
    rel = Path(label) / f"case_{idx:04d}_g2_alias_a2.tsv"
    write_tsv(out_root / rel, rows)
    append_participant(participants, rel, label, idx)
    idx += 1

    # 5) Canonical + alias duplicated
    rows = g2_both_ids_rows("AG", "TG", "ID")  # both G1 het + 1 G2 -> G1G2
    label = derive_label(rows, "G", "G")
    rel = Path(label) / f"case_{idx:04d}_g2_canonical_plus_alias.tsv"
    write_tsv(out_root / rel, rows)
    append_participant(participants, rel, label, idx)
    idx += 1

    # 6) Missing G2 (unknown G2)
    rows = rows_missing("G2", "AG", "TT", None)  # G1 het only, G2 unknown -> G1G-
    label = derive_label(rows, "G", "G")
    rel = Path(label) / f"case_{idx:04d}_missing_g2.tsv"
    write_tsv(out_root / rel, rows)
    append_participant(participants, rel, label, idx)
    idx += 1

    # 7) Missing G1A
    rows = rows_missing("G1A", None, "TG", "II")  # G1 unknown (need both), G2=0 -> G0G-
    label = derive_label(rows, "G", "G")
    rel = Path(label) / f"case_{idx:04d}_missing_g1a.tsv"
    write_tsv(out_root / rel, rows)
    append_participant(participants, rel, label, idx)
    idx += 1

    # 8) Missing G1B
    rows = rows_missing("G1B", "AG", None, "DD")  # G1 unknown, G2=2 -> G2G-
    label = derive_label(rows, "G", "G")
    rel = Path(label) / f"case_{idx:04d}_missing_g1b.tsv"
    write_tsv(out_root / rel, rows)
    append_participant(participants, rel, label, idx)
    idx += 1

    # 9) Invalid G1A genotype (A?) + G2 known
    rows = rows_invalid_gt("G1A", "AG", "TT", "ID")  # G1 unknown, G2=1 -> G2G-
    label = derive_label(rows, "G", "G")
    rel = Path(label) / f"case_{idx:04d}_invalid_g1a.tsv"
    write_tsv(out_root / rel, rows)
    append_participant(participants, rel, label, idx)
    idx += 1

    # 10) Invalid G2 genotype (IDK)
    rows = rows_invalid_gt("G2", "AG", "TG", "ID")  # G1=1, G2 unknown -> G1G-
    label = derive_label(rows, "G", "G")
    rel = Path(label) / f"case_{idx:04d}_invalid_g2.tsv"
    write_tsv(out_root / rel, rows)
    append_participant(participants, rel, label, idx)
    idx += 1

    # 11) Exactly one G-G- file with irrelevant rsids only
    rows = rows_irrelevant_only()
    label = derive_label(rows, "G", "G")  # will be G-G-
    rel = Path(label) / f"case_{idx:04d}_irrelevant_only.tsv"
    write_tsv(out_root / rel, rows)
    append_participant(participants, rel, label, idx)
    idx += 1

    # ---- participants.csv (deterministic) ----
    csv_path = out_root / "participants.csv"
    with csv_path.open("w", encoding="utf-8") as f:
        f.write(
            "participant_id,genotype_file_path,weight,height,age,expected_genotype\n"
        )
        for row in participants:
            f.write(",".join(map(str, row)) + "\n")

    print(f"✅ Wrote test TSVs + participants.csv under: {out_root.resolve()}")


if __name__ == "__main__":
    main()
