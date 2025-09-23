#!/usr/bin/env python3
import csv
import os
import re
import sys

# ---------- Config ----------
TARGET_RSIDS = [
    "rs73885319",  # G1 part A (A>G)
    "rs60910145",  # G1 part B (T>G)
    "rs71785313",  # G2 (TTATAA > -)
    "rs1317778148",  # merged into rs71785313
    "rs143830837",  # merged into rs71785313
]

DEFAULT_FIELDS = ["rsid", "chromosome", "position", "genotype", "gs", "baf", "lrr"]

RS_G1_A = "rs73885319"
RS_G1_B = "rs60910145"
RS_G2 = "rs71785313"
RS_MERGED = {"rs1317778148", "rs143830837"}

# Heuristic thresholds for inferring 0/0 vs 0/1 vs 1/1 from BAF/AD
HET_MIN, HET_MAX, HOM_ALT_MIN = 0.20, 0.80, 0.80


# ---------- Helpers ----------
def detect_header_and_stream(lines):
    """
    Accepts a TSV-ish stream that may contain comment lines starting with '#'
    and may or may not contain a header line. Returns (header_fields, data_lines).
    """
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

    results = {rsid: None for rsid in TARGET_RSIDS}
    for row in reader:
        rsid_key = key_map["rsid"]
        if not rsid_key:
            continue
        rsid_val = (row.get(rsid_key) or "").strip()
        if rsid_val in results:
            results[rsid_val] = row
    return header_fields, results, key_map


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


def two_letter_zygosity(genostr, ref=None, alt=None):
    """
    Try to infer 0/0, 0/1, 1/1 from two-letter genotypes (e.g., 'AG', 'TT').
    If REF/ALT are provided, use them; otherwise, fallback to equality check.
    """
    if not genostr or len(genostr) < 2:
        return None
    g = genostr.strip().upper().replace(" ", "")
    if len(g) == 2 and g.isalpha():
        a, b = g[0], g[1]
        if ref and alt:
            ref = ref.upper()
            alt_first = alt.split(",")[0].upper()
            if a == ref and b == ref:
                return "0/0"
            if {a, b} == {ref, alt_first}:
                return "0/1"
            if a == alt_first and b == alt_first:
                return "1/1"
        # no REF/ALT → assume hom letters => 0/0; het letters => 0/1
        return "0/0" if a == b else "0/1"
    return None


def call_from_row(row, key_map):
    """
    Return site-level genotype in index form: '0/0', '0/1', '1/1', or './.' if unknown.
    Uses, in order of preference: TYPE/EVENT_GT → AD/DP → GT (index) → two-letter → BAF.
    """
    if row is None:
        return "./."

    def get(k, d=""):
        kk = key_map.get(k)
        return (row.get(kk) if kk else d) if row else d

    typ, event_gt, dp = (
        (get("type") or "").upper(),
        get("event_gt"),
        parse_int(get("dp")),
    )
    baf, ad, ref, alt, geno = (
        parse_float(get("baf")),
        get("ad"),
        get("ref"),
        get("alt"),
        get("genotype"),
    )

    # 1) Structural-style typed events (DEL/INS) with explicit EVENT_GT
    if typ in ("DEL", "INS") and event_gt in (
        "REF/REF",
        "REF/DEL",
        "DEL/DEL",
        "REF/INS",
        "INS/INS",
    ):
        return {
            "REF/REF": "0/0",
            "REF/DEL": "0/1",
            "DEL/DEL": "1/1",
            "REF/INS": "0/1",
            "INS/INS": "1/1",
        }[event_gt]

    # 2) AD/DP counts
    refc, altc = parse_ad(ad)
    if dp <= 0 and (refc + altc) > 0:
        dp = refc + altc
    if dp > 0:
        baf_est = altc / float(dp)
        if baf_est >= HOM_ALT_MIN:
            return "1/1"
        if HET_MIN <= baf_est < HET_MAX:
            return "0/1"
        return "0/0"

    # 3) Index-style GT like '0/1' or '1|1'
    gt_idx = normalize_gt_like(geno)
    if gt_idx:
        return gt_idx

    # 4) Two-letter genotype (e.g., 'AG'), ideally with REF/ALT
    letter_gt = two_letter_zygosity(geno, ref, alt)
    if letter_gt:
        return letter_gt

    # 5) Bare BAF
    if baf is not None:
        if baf >= HOM_ALT_MIN:
            return "1/1"
        if HET_MIN <= baf < HET_MAX:
            return "0/1"
        return "0/0"

    return "./."


def g1_count_from_two_calls(ca, cb):
    if ca == "1/1" and cb == "1/1":
        return 2
    if (ca in ("0/1", "1/1")) and (cb in ("0/1", "1/1")):
        return 1
    return 0


def g2_count_from_call(c):
    return 2 if c == "1/1" else 1 if c == "0/1" else 0


# ---------- Main ----------
if __name__ == "__main__":
    header_fields, results, key_map = extract_snps(sys.stdin)

    # Pick best records for each signal
    best_g1a, best_g1b = results.get(RS_G1_A), results.get(RS_G1_B)
    G2_CANDIDATES = [RS_G2] + sorted(RS_MERGED)
    g2_id_used = next((rid for rid in G2_CANDIDATES if results.get(rid)), None)
    best_g2 = results.get(g2_id_used) if g2_id_used else None

    # Site-level calls
    call_g1a = call_from_row(best_g1a, key_map)
    call_g1b = call_from_row(best_g1b, key_map)
    call_g2 = call_from_row(best_g2, key_map)

    # Counts → final APOL1 genotype
    g1_count = g1_count_from_two_calls(call_g1a, call_g1b)
    g2_count = g2_count_from_call(call_g2)

    if g2_count == 2:
        apol1 = "G2/G2"
    elif g1_count == 2:
        apol1 = "G1/G1"
    elif g1_count >= 1 and g2_count >= 1:
        apol1 = "G1/G2"
    elif g1_count == 1:
        apol1 = "G1/G0"
    elif g2_count == 1:
        apol1 = "G2/G0"
    else:
        apol1 = "G0/G0"

    # --- Write outputs with prefix ---
    base = os.environ.get("APOL1_OUT_PREFIX", "apol1_result")

    # 1) rsid/genotype table (now includes ref/alt columns, if present)
    with open(f"{base}.rsid.tsv", "w") as fh:
        fh.write("rsid\tchromosome\tposition\tgenotype\tref\talt\n")
        for rsid in TARGET_RSIDS:
            row = results[rsid]
            geno = row.get(key_map.get("genotype"), "--") if row else "--"
            chrom = row.get(key_map.get("chromosome"), "--") if row else "--"
            pos = row.get(key_map.get("position"), "--") if row else "--"
            ref = row.get(key_map.get("ref"), "--") if row else "--"
            alt = row.get(key_map.get("alt"), "--") if row else "--"
            fh.write(f"{rsid}\t{chrom}\t{pos}\t{geno}\t{ref}\t{alt}\n")

    # 2) genotype line only (final APOL1 call)
    with open(f"{base}.genotype.txt", "w") as fh:
        fh.write(f"APOL1 genotype: {apol1}\n")

    # 3) evidence + heuristic (now includes bases=<raw genotype> in each line)
    def brief_row(rsid, row, call, key_map):
        if row is None:
            return f"{rsid}: no row"
        chrom = row.get(key_map.get("chromosome"), "--")
        pos = row.get(key_map.get("position"), "--")
        typ = (row.get(key_map.get("type"), "") or "").upper()
        ref = row.get(key_map.get("ref"), "--")
        alt = row.get(key_map.get("alt"), "--")
        ad = row.get(key_map.get("ad"), "--")
        dp = row.get(key_map.get("dp"), "--")
        baf = row.get(key_map.get("baf"), "--")
        geno = row.get(key_map.get("genotype"), "--")
        return (
            f"{rsid}\tpos={chrom}:{pos}\ttype={typ}\t"
            f"GT≈{call}\tbases={geno}\tref={ref} alt={alt}\tAD={ad}\tDP={dp}\tBAF={baf}"
        )

    explain = [
        brief_row(RS_G1_A, best_g1a, call_g1a, key_map),
        brief_row(RS_G1_B, best_g1b, call_g1b, key_map),
        brief_row(g2_id_used or RS_G2, best_g2, call_g2, key_map),
        "(Note: rs1317778148 and rs143830837 are merged into rs71785313)",
    ]

    with open(f"{base}.evidence.txt", "w") as fh:
        fh.write("Calls (evidence):\n")
        for line in explain:
            fh.write(f"- {line}\n")
        fh.write("\nHeuristic (unphased):\n")
        fh.write("- G1 requires ALT at BOTH rs73885319 (A>G) AND rs60910145 (T>G).\n")
        fh.write(
            "  * both sites hom-alt -> 2×G1; both non-ref (het/alt) -> 1×G1; else 0×G1.\n"
        )
        fh.write("- G2 is the 6bp deletion at rs71785313 (TTATAA>-).\n")
        fh.write("  * 1/1 -> 2×G2; 0/1 -> 1×G2; 0/0 -> 0×G2.\n")
        fh.write("- Final genotype is derived from (G1 count, G2 count).\n")
