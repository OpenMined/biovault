#!/usr/bin/env python3
import csv
import os
import re
import sys
from typing import Dict, List, Optional, Set, Tuple

# // https://www.ncbi.nlm.nih.gov/snp/?term=rs73885319
# // rs73885319
# // Variant type:SNVAlleles:A>G [Show Flanks]Chromosome:
# // 22:36265860 (GRCh38)
# // 22:36661906 (GRCh37)
# // https://www.ncbi.nlm.nih.gov/snp/?term=rs60910145
# // rs60910145
# // Alleles:T>C,G [Show Flanks]Chromosome:
# // 22:36265988 (GRCh38)
# // 22:36662034 (GRCh37)

# // https://www.ncbi.nlm.nih.gov/snp/rs71785313
# // rs71785313
# // Variant type:INDELAlleles:TTATAA>- [Show Flanks]Chromosome:

# G1 = actually two missense mutations (S342G and I384M) that occur together on the same haplotype.
# G2 = a 6 base-pair deletion (loss of two amino acids, N388 and Y389).
# G0/G1 or G0/G2 → “carriers,” usually no major increased risk
# G1/G1, G1/G2, G2/G2 → “high-risk genotypes” linked to kidney disease progression.


# ---------- Config ----------
TARGET_RSIDS = [
    "rs73885319",  # G1 part A (A>G, sometimes A>C)
    "rs60910145",  # G1 part B (T>G, sometimes T>C)
    "rs71785313",  # G2 (TTATAA > -) canonical
    "rs1317778148",  # merged into rs71785313 (alias)
    "rs143830837",  # merged into rs71785313 (alias)
]

DEFAULT_FIELDS = ["rsid", "chromosome", "position", "genotype", "gs", "baf", "lrr"]

RS_G1_A = "rs73885319"
RS_G1_B = "rs60910145"
RS_G2_CANON = "rs71785313"
RS_G2_ALIASES = {"rs1317778148", "rs143830837"}
RS_G2_ALL = {RS_G2_CANON} | RS_G2_ALIASES

# Authoritative REF and allowed ALT sets for G1 SNPs
KNOWN_REF_ALLOWED_ALTS: Dict[str, Tuple[str, Set[str]]] = {
    RS_G1_A: ("A", {"G", "C"}),  # A>G normally; some assays show A>C
    RS_G1_B: ("T", {"G", "C"}),  # T>G normally; some assays show T>C
}

# Heuristic thresholds (used only for SNP fallback from counts/BAF; NEVER for G2)
HET_MIN, HET_MAX, HOM_ALT_MIN = 0.20, 0.80, 0.80


# ---------- Robust TSV ingestion ----------
def detect_header_and_stream(lines):
    header_fields, data_lines = None, []
    for raw in lines:
        line = raw.lstrip("\ufeff").rstrip("\r\n")
        if not line:
            continue
        if line.startswith("#"):
            maybe = line.lstrip("#").strip()
            if maybe:
                parts = maybe.split("\t")
                lower = [p.strip().lower() for p in parts]
                if {"rsid", "chromosome", "position", "genotype"} <= set(lower):
                    header_fields = parts
            continue
        if header_fields is None:
            parts = line.split("\t")
            lower = [p.strip().lower() for p in parts]
            if {"rsid", "chromosome", "position", "genotype"} <= set(lower):
                header_fields = parts
                continue
        data_lines.append(line)
    return header_fields or DEFAULT_FIELDS, data_lines


def extract_snps(input_stream):
    header_fields, data_lines = detect_header_and_stream(input_stream)
    reader = csv.DictReader(data_lines, fieldnames=header_fields, delimiter="\t")

    def key_for(needed):
        return next((k for k in header_fields if k.lower() == needed), None)

    key_map = {
        k: key_for(k)
        for k in [
            "rsid",
            "chromosome",
            "position",
            "genotype",
            "ref",
            "alt",
            "type",
            "event_gt",
            "ad",
            "dp",
            "baf",
        ]
    }

    results: Dict[str, Optional[Dict[str, str]]] = {rsid: None for rsid in TARGET_RSIDS}
    for row in reader:
        rsid_key = key_map["rsid"]
        if not rsid_key:
            continue
        rsid_val = (row.get(rsid_key) or "").strip()
        if rsid_val in results:
            results[rsid_val] = row
    return header_fields, results, key_map


# ---------- Parsing helpers ----------
def parse_int(x, default=0):
    try:
        return int(x)
    except Exception:
        return default


def parse_float(x, default=None):
    try:
        return float(x)
    except Exception:
        return default


def parse_ad(ad_str):
    if not ad_str:
        return (0, 0)
    parts = [p.strip() for p in ad_str.split(",")]
    return parse_int(parts[0] if len(parts) > 0 else "0"), parse_int(
        parts[1] if len(parts) > 1 else "0"
    )


def normalize_gt_like(genostr):
    """
    Normalize index-style GT strings like '0/1', '1|1' → '0/1' or '1/1'.
    """
    if not genostr:
        return None
    m = re.fullmatch(r"\s*([01])[\/|]([01])\s*", genostr)
    return f"{m.group(1)}/{m.group(2)}" if m else None


def two_letter_with_allowed(
    genostr: Optional[str], ref: str, allowed_alts: Set[str]
) -> Optional[str]:
    """
    Map two-letter genotype using authoritative ref and a SET of allowed alts.
    - If both letters == ref -> '0/0'
    - If one is ref and the other in allowed_alts -> '0/1'
    - If both in allowed_alts (even if different from each other) -> '1/1'
      (treat mixed-alt like 'GC' as non-ref/non-ref; still show exact 'bases=' in evidence)
    - If any letter not in {ref} ∪ allowed_alts -> unknown (None)
    """
    if not genostr or len(genostr) < 2:
        return None
    g = genostr.strip().upper().replace(" ", "")
    if len(g) != 2 or not g.isalpha():
        return None
    a, b = g[0], g[1]
    valid = {ref} | allowed_alts
    if a not in valid or b not in valid:
        return None  # unexpected letter (e.g., N, or wrong locus)
    if a == ref and b == ref:
        return "0/0"
    if a in allowed_alts and b in allowed_alts:
        return "1/1"
    # one ref, one allowed alt
    return "0/1"


# ---------- Calling core ----------
def g2_call_from_text_indel(genostr: Optional[str]) -> Tuple[str, bool]:
    """
    Map 'II'/'ID'/'DD' → ('0/0'|'0/1'|'1/1', True).
    Anything else or missing → ('./.', False).
    """
    if not genostr:
        return "./.", False
    g = genostr.strip().upper()
    if g == "II":
        return "0/0", True
    if g == "ID":
        return "0/1", True
    if g == "DD":
        return "1/1", True
    return "./.", False


def call_from_row_snp(rsid: str, row, key_map) -> Tuple[str, bool]:
    """
    Return ('0/0'|'0/1'|'1/1'|'./.', known_bool) for SNP rows (G1 sites).
    Preference: AD/DP → GT index → two-letter with (ref + allowed_alts) → BAF.
    (We ignore file REF/ALT if missing/wrong; we rely on authoritative sets.)
    """
    if row is None:
        return "./.", False

    def get(k, d=""):
        kk = key_map.get(k)
        return (row.get(kk) if kk else d) if row else d

    dp = parse_int(get("dp"))
    baf = parse_float(get("baf"))
    ad = get("ad")
    geno = get("genotype")

    # Authoritative ref + allowed alts for this rsid (if known)
    ref, allowed = KNOWN_REF_ALLOWED_ALTS.get(rsid, (None, None))

    # 1) AD/DP counts (if present)
    refc, altc = parse_ad(ad)
    if dp <= 0 and (refc + altc) > 0:
        dp = refc + altc
    if dp > 0:
        baf_est = altc / float(dp)
        if baf_est >= HOM_ALT_MIN:
            return "1/1", True
        if HET_MIN <= baf_est < HET_MAX:
            return "0/1", True
        return "0/0", True

    # 2) Index-style GT like '0/1' or '1|1'
    gt_idx = normalize_gt_like(geno)
    if gt_idx:
        return gt_idx, True

    # 3) Two-letter genotype using authoritative ref + allowed_alts
    if ref and allowed:
        letter_gt = two_letter_with_allowed(geno, ref, allowed)
        if letter_gt:
            return letter_gt, True

    # 4) Bare BAF (last resort)
    if baf is not None:
        if baf >= HOM_ALT_MIN:
            return "1/1", True
        if HET_MIN <= baf < HET_MAX:
            return "0/1", True
        return "0/0", True

    return "./.", False


def call_from_row_g2(row, key_map) -> Tuple[str, bool]:
    """
    STRICT G2: trust ONLY the textual 'genotype' field as II/ID/DD.
    """
    if row is None:
        return "./.", False
    geno = row.get(key_map.get("genotype")) if key_map.get("genotype") else None
    return g2_call_from_text_indel(geno)


def g1_count_from_two_calls(ca: str, cb: str) -> Tuple[Optional[int], bool]:
    known = ca in {"0/0", "0/1", "1/1"} and cb in {"0/0", "0/1", "1/1"}
    if not known:
        return None, False
    if ca == "1/1" and cb == "1/1":
        return 2, True
    if (ca in {"0/1", "1/1"}) and (cb in {"0/1", "1/1"}):
        return 1, True
    return 0, True


def g2_count_from_call(c: str) -> Tuple[Optional[int], bool]:
    if c not in {"0/0", "0/1", "1/1"}:
        return None, False
    return (2 if c == "1/1" else 1 if c == "0/1" else 0), True


def final_label(
    g1_count: Optional[int], g1_known: bool, g2_count: Optional[int], g2_known: bool
) -> str:
    if g1_known and g2_known:
        gc1, gc2 = g1_count or 0, g2_count or 0
        if gc2 == 2:
            return "G2/G2"
        if gc1 == 2:
            return "G1/G1"
        if gc1 >= 1 and gc2 >= 1:
            return "G1/G2"
        if gc1 == 1:
            return "G1/G0"
        if gc2 == 1:
            return "G2/G0"
        return "G0/G0"
    if g1_known and not g2_known:
        return "G1/G-" if (g1_count or 0) >= 1 else "G0/G-"
    if g2_known and not g1_known:
        return "G2/G-" if (g2_count or 0) >= 1 else "G0/G-"
    return "G-/G-"


def pick_best_g2(
    results: Dict[str, Optional[Dict[str, str]]],
) -> Tuple[Optional[str], Optional[Dict[str, str]]]:
    if results.get(RS_G2_CANON):
        return RS_G2_CANON, results[RS_G2_CANON]
    for rid in sorted(RS_G2_ALIASES):
        if results.get(rid):
            return rid, results[rid]
    return None, None


# ---------- Main ----------
if __name__ == "__main__":
    header_fields, results, key_map = extract_snps(sys.stdin)

    # --- G1 site calls (SNP path) ---
    best_g1a = results.get(RS_G1_A)
    best_g1b = results.get(RS_G1_B)
    call_g1a, known_g1a = call_from_row_snp(RS_G1_A, best_g1a, key_map)
    call_g1b, known_g1b = call_from_row_snp(RS_G1_B, best_g1b, key_map)

    # --- G2 call (STRICT indel-from-text path, alias-aware) ---
    g2_id_used, best_g2 = pick_best_g2(results)

    # If multiple G2 IDs are present, require agreement among their text genotypes; else unknown.
    g2_texts: List[str] = []
    for rid in RS_G2_ALL:
        row = results.get(rid)
        if not row:
            continue
        geno = row.get(key_map.get("genotype")) if key_map.get("genotype") else None
        if geno:
            g2_texts.append(geno.strip().upper())
    g2_texts = [t for t in g2_texts if t]

    if len(g2_texts) == 0:
        call_g2, known_g2 = "./.", False
    else:
        agree = all(t == g2_texts[0] for t in g2_texts)
        if not agree:
            call_g2, known_g2 = "./.", False
        else:
            call_g2, known_g2 = g2_call_from_text_indel(g2_texts[0])
            if best_g2 is None:
                best_g2 = next(
                    (results[rid] for rid in RS_G2_ALL if results.get(rid)), None
                )
            if g2_id_used is None:
                g2_id_used = (
                    RS_G2_CANON
                    if results.get(RS_G2_CANON)
                    else next(iter(RS_G2_ALIASES), RS_G2_CANON)
                )

    # --- Aggregate → final genotype with unknown propagation ---
    g1_count, g1_known = g1_count_from_two_calls(call_g1a, call_g1b)
    g2_count, g2_known = g2_count_from_call(call_g2)
    apol1 = final_label(g1_count, g1_known, g2_count, g2_known)

    # --- Write outputs with prefix ---
    base = os.environ.get("APOL1_OUT_PREFIX", "apol1_result")

    # 1) rsid/genotype table (include ref/alt columns if present)
    with open(f"{base}.rsid.tsv", "w") as fh:
        fh.write("rsid\tchromosome\tposition\tgenotype\tref\talt\n")
        for rsid in TARGET_RSIDS:
            row = results.get(rsid)
            geno = row.get(key_map.get("genotype"), "--") if row else "--"
            chrom = row.get(key_map.get("chromosome"), "--") if row else "--"
            pos = row.get(key_map.get("position"), "--") if row else "--"
            ref = row.get(key_map.get("ref"), "--") if row else "--"
            alt = row.get(key_map.get("alt"), "--") if row else "--"
            fh.write(f"{rsid}\t{chrom}\t{pos}\t{geno}\t{ref}\t{alt}\n")

    # 2) genotype line only (final APOL1 call)
    with open(f"{base}.genotype.txt", "w") as fh:
        fh.write(f"APOL1 genotype: {apol1}\n")

    # 3) evidence + heuristic (now includes allowed_alts for G1 clarity)
    def brief_row_g1(rsid, row, call):
        if row is None:
            return f"{rsid}\tno row"
        chrom = row.get(key_map.get("chromosome"), "--")
        pos = row.get(key_map.get("position"), "--")
        typ = (row.get(key_map.get("type"), "") or "").upper()
        ref_file = row.get(key_map.get("ref"), "--")
        alt_file = row.get(key_map.get("alt"), "--")
        ad = row.get(key_map.get("ad"), "--")
        dp = row.get(key_map.get("dp"), "--")
        baf = row.get(key_map.get("baf"), "--")
        geno = row.get(key_map.get("genotype"), "--")
        ref_auth, alts_auth = KNOWN_REF_ALLOWED_ALTS.get(rsid, ("?", set()))
        alts_auth_str = (
            ",".join(sorted(alts_auth))
            if isinstance(alts_auth, set)
            else str(alts_auth)
        )
        return (
            f"{rsid}\tpos={chrom}:{pos}\ttype={typ}\t"
            f"GT≈{call}\tbases={geno}\tref={ref_file} alt={alt_file}\t"
            f"AD={ad}\tDP={dp}\tBAF={baf}\tallowed_alts={ref_auth}->{alts_auth_str}"
        )

    def brief_row_g2(rsid, row, call):
        if row is None:
            return f"{rsid}\tno row"
        chrom = row.get(key_map.get("chromosome"), "--")
        pos = row.get(key_map.get("position"), "--")
        typ = (row.get(key_map.get("type"), "") or "").upper()
        ref_file = row.get(key_map.get("ref"), "--")
        alt_file = row.get(key_map.get("alt"), "--")
        ad = row.get(key_map.get("ad"), "--")
        dp = row.get(key_map.get("dp"), "--")
        baf = row.get(key_map.get("baf"), "--")
        geno = row.get(key_map.get("genotype"), "--")
        return (
            f"{rsid}\tpos={chrom}:{pos}\ttype={typ}\t"
            f"GT≈{call}\tbases={geno}\tref={ref_file} alt={alt_file}\t"
            f"AD={ad}\tDP={dp}\tBAF={baf}"
        )

    evidence_lines = [
        brief_row_g1(RS_G1_A, best_g1a, call_g1a),
        brief_row_g1(RS_G1_B, best_g1b, call_g1b),
        brief_row_g2(g2_id_used or RS_G2_CANON, best_g2, call_g2),
        "(Note: rs1317778148 and rs143830837 are merged into rs71785313; all must agree or G2=unknown.)",
    ]

    with open(f"{base}.evidence.txt", "w") as fh:
        fh.write("Calls (evidence):\n")
        for line in evidence_lines:
            fh.write(f"- {line}\n")
        fh.write("\nHeuristic (unphased):\n")
        fh.write(
            "- G1 requires ALT at BOTH rs73885319 (A>G/C) AND rs60910145 (T>G/C).\n"
        )
        fh.write(
            "  * both sites hom-alt -> 2×G1; both non-ref (het/alt) -> 1×G1; else 0×G1.\n"
        )
        fh.write("- G2 is the 6bp deletion at rs71785313 (TTATAA>-).\n")
        fh.write("  * DD -> 2×G2; ID -> 1×G2; II -> 0×G2. (Aliases must agree.)\n")
        fh.write(
            "- Unknown stays unknown (printed as '-'). Final genotype derives from (G1 count, G2 count).\n"
        )
