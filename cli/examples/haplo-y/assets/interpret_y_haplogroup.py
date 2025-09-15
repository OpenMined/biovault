#!/usr/bin/env python3
import argparse
import sys
from collections import defaultdict

DESC = """Interpret a tiny Y-haplogroup SNP panel from targeted calls.

Inputs:
  (1) Panel TSV (with header):
      locus_id  region  haplogroup  derived_allele  note
      - region is 1-bp like: chrY:20577481-20577481
  (2) Variants TSV from bcftools query:
      CHROM  POS  REF  ALT  GT  AD  DP

Outputs:
  - Human report (.txt)
  - Clean calls table (.tsv) with per-locus allele status

Calling logic (toy, robust):
  - DP < min_dp -> unknown
  - If ALT equals derived_allele:
       GT contains '1'  -> derived
       GT == 0/0 or 0|0 -> reference
    Fallback by AD if GT missing/odd.
"""


def parse_args():
    ap = argparse.ArgumentParser(description=DESC)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--participant", required=True)
    ap.add_argument("--ref-version", required=True, choices=["GRCh37", "GRCh38"])
    ap.add_argument("--out-report", required=True)
    ap.add_argument("--out-variants", required=True)
    ap.add_argument("--min-dp", type=int, default=5)
    return ap.parse_args()


def load_panel(p):
    rows = []
    with open(p, "r", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        idx = {k: i for i, k in enumerate(header)}
        for need in ["locus_id", "region", "haplogroup", "derived_allele"]:
            if need not in idx:
                raise ValueError(f"Panel missing column: {need}")
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n\r").split("\t")
            region = parts[idx["region"]]
            chrom, rng = region.split(":")
            start, end = rng.split("-")
            pos = int(start)
            rows.append(
                {
                    "locus_id": parts[idx["locus_id"]],
                    "chrom": chrom,
                    "pos": pos,
                    "haplogroup": parts[idx["haplogroup"]],
                    "derived_allele": parts[idx["derived_allele"]],
                    "note": parts[idx["note"]] if "note" in idx else "",
                }
            )
    return rows


def load_variants(path):
    sites = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            chrom, pos, ref, alt, gt, ad, dp = line.strip().split("\t")
            pos = int(pos)
            dp = int(dp) if dp.isdigit() else 0
            ad_ref, ad_alt = 0, 0
            if ad and ad != ".":
                bits = ad.split(",")
                if len(bits) >= 2:
                    try:
                        ad_ref = int(bits[0])
                        ad_alt = int(bits[1])
                    except:
                        pass
            sites[(chrom, pos)] = {
                "REF": ref,
                "ALT": "" if alt == "." else alt,
                "GT": gt,
                "AD_REF": ad_ref,
                "AD_ALT": ad_alt,
                "DP": dp,
            }
    return sites


def call_status(site, derived, min_dp):
    if site is None:
        return ("unknown", ".")
    dp, gt, ref, alt = site["DP"], site["GT"], site["REF"], site["ALT"]
    r, a = site["AD_REF"], site["AD_ALT"]
    if dp < min_dp:
        return ("unknown", f"{ref}>{alt or '.'};DP={dp}")
    if alt and derived and alt == derived:
        if "1" in gt.replace("|", "/"):
            return ("derived", f"{ref}>{alt};GT={gt};AD={r},{a};DP={dp}")
        if gt in ("0/0", "0|0"):
            return ("reference", f"{ref}>{alt};GT={gt};AD={r},{a};DP={dp}")
    if a > 0 and a >= max(2, r):
        return ("derived", f"{ref}>{alt or derived or '.'};AD={r},{a};DP={dp}")
    if r > 0 and a == 0:
        return ("reference", f"{ref}>{alt or '.'};AD={r},{a};DP={dp}")
    return ("unknown", f"{ref}>{alt or '.'};GT={gt};AD={r},{a};DP={dp}")


def main():
    args = parse_args()
    panel = load_panel(args.panel)
    sites = load_variants(args.variants)

    per_locus = []
    hap_hits = defaultdict(list)

    for row in panel:
        key = (row["chrom"], row["pos"])
        status, detail = call_status(sites.get(key), row["derived_allele"], args.min_dp)
        per_locus.append(
            {
                "locus_id": row["locus_id"],
                "chrom": row["chrom"],
                "pos": row["pos"],
                "haplogroup": row["haplogroup"],
                "derived_allele": row["derived_allele"],
                "status": status,
                "detail": detail,
                "note": row["note"],
            }
        )
        if status == "derived":
            hap_hits[row["haplogroup"]].append(row["locus_id"])

    best_hap = None
    if hap_hits:
        best_hap = sorted(hap_hits.items(), key=lambda kv: (-len(kv[1]), kv[0]))[0][0]

    with open(args.out_variants, "w", encoding="utf-8") as f:
        f.write(
            "participant\tref_version\tlocus_id\thaplogroup\tchrom\tpos\tderived_allele\tstatus\tdetail\tnote\n"
        )
        for r in per_locus:
            f.write(
                f"{args.participant}\t{args.ref_version}\t{r['locus_id']}\t{r['haplogroup']}\t"
                f"{r['chrom']}\t{r['pos']}\t{r['derived_allele']}\t{r['status']}\t{r['detail']}\t{r['note']}\n"
            )

    with open(args.out_report, "w", encoding="utf-8") as f:
        f.write(f"Participant: {args.participant}\n")
        f.write(f"Reference:   {args.ref_version}\n")
        f.write(f"Panel:       {args.panel}\n\n")
        if best_hap:
            f.write(f"🧬 Y-haplogroup (toy call): {best_hap}\n")
            f.write(
                f"Supporting markers: {', '.join(sorted(set(hap_hits[best_hap])))}\n\n"
            )
        else:
            f.write(
                "🧬 Y-haplogroup (toy call): Unknown (no confident derived markers)\n\n"
            )
        f.write("Per-locus details:\n")
        f.write("locus_id\thaplogroup\tpos\tstatus\tdetail\tnote\n")
        for r in per_locus:
            f.write(
                f"{r['locus_id']}\t{r['haplogroup']}\t{r['chrom']}:{r['pos']}\t{r['status']}\t{r['detail']}\t{r['note']}\n"
            )


if __name__ == "__main__":
    main()
