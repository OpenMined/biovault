"""
Allele Frequency Utilities

Functions for converting between TSV and numpy array formats,
supporting secure aggregation with aligned locus indices.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple, Optional
import json
import numpy as np


@dataclass(frozen=True)
class LocusIndex:
    """Canonical locus ordering for aligned numpy arrays."""
    loci: np.ndarray      # dtype=object (strings), shape (N,)
    rsids: np.ndarray     # dtype=object (strings), shape (N,)
    version: str = "1.0"

    def __len__(self) -> int:
        return len(self.loci)

    def to_pos_map(self) -> dict[str, int]:
        """Return locus -> index mapping."""
        return {str(l): i for i, l in enumerate(self.loci)}

    def save(self, path: str) -> None:
        """Save index to JSON file."""
        data = {
            "version": self.version,
            "n_loci": len(self.loci),
            "loci": [str(l) for l in self.loci],
            "rsids": [str(r) for r in self.rsids]
        }
        with open(path, 'w') as f:
            json.dump(data, f)

    @staticmethod
    def load(path: str) -> "LocusIndex":
        """Load index from JSON file."""
        with open(path) as f:
            data = json.load(f)
        return LocusIndex(
            loci=np.array(data["loci"], dtype=object),
            rsids=np.array(data.get("rsids", [""] * len(data["loci"])), dtype=object),
            version=data.get("version", "1.0")
        )

    @staticmethod
    def from_union(*indices: "LocusIndex", sort: bool = True) -> "LocusIndex":
        """Create union index from multiple indices."""
        locus_to_rsid = {}
        for idx in indices:
            for locus, rsid in zip(idx.loci, idx.rsids):
                if str(locus) not in locus_to_rsid:
                    locus_to_rsid[str(locus)] = str(rsid)
                elif not locus_to_rsid[str(locus)] and rsid:
                    locus_to_rsid[str(locus)] = str(rsid)

        loci = sorted(locus_to_rsid.keys()) if sort else list(locus_to_rsid.keys())
        rsids = [locus_to_rsid[l] for l in loci]
        return LocusIndex(
            loci=np.array(loci, dtype=object),
            rsids=np.array(rsids, dtype=object)
        )


@dataclass
class AlleleFreqData:
    """Allele frequency data aligned to a locus index."""
    ac: np.ndarray       # allele counts, shape (N,)
    an: np.ndarray       # allele numbers, shape (N,)
    n_samples: np.ndarray  # sample counts, shape (N,)
    index: LocusIndex

    def __len__(self) -> int:
        return len(self.index)

    @property
    def af(self) -> np.ndarray:
        """Compute allele frequencies."""
        with np.errstate(divide='ignore', invalid='ignore'):
            af = self.ac / self.an
            af[~np.isfinite(af)] = 0.0
        return af

    def save_npz(self, path: str) -> None:
        """Save arrays to compressed npz file."""
        np.savez_compressed(
            path,
            ac=self.ac,
            an=self.an,
            n_samples=self.n_samples
        )

    @staticmethod
    def load(npz_path: str, index_path: str) -> "AlleleFreqData":
        """Load from npz and index files."""
        index = LocusIndex.load(index_path)
        data = np.load(npz_path)
        return AlleleFreqData(
            ac=data["ac"],
            an=data["an"],
            n_samples=data["n_samples"],
            index=index
        )

    def to_tsv(self, include_zero: bool = True) -> str:
        """Convert to TSV string."""
        lines = ["locus\trsid\tac\tan\taf\tn_samples"]
        af = self.af

        for i in range(len(self.index)):
            if not include_zero and self.ac[i] == 0:
                continue
            locus = self.index.loci[i]
            rsid = self.index.rsids[i]
            lines.append(
                f"{locus}\t{rsid}\t{self.ac[i]}\t{self.an[i]}\t{af[i]:.6f}\t{self.n_samples[i]}"
            )

        return "\n".join(lines)

    def save_tsv(self, path: str, include_zero: bool = True) -> None:
        """Save to TSV file."""
        with open(path, 'w') as f:
            f.write(self.to_tsv(include_zero=include_zero))

    @staticmethod
    def from_tsv(tsv_path: str) -> "AlleleFreqData":
        """Load from TSV file."""
        loci, rsids, ac, an, n_samples = [], [], [], [], []

        with open(tsv_path) as f:
            header = f.readline().strip().split('\t')
            col_idx = {name: i for i, name in enumerate(header)}

            for line in f:
                parts = line.strip().split('\t')
                loci.append(parts[col_idx["locus"]])
                rsids.append(parts[col_idx.get("rsid", 1)] if "rsid" in col_idx else "")
                ac.append(int(parts[col_idx["ac"]]))
                an.append(int(parts[col_idx["an"]]))
                n_samples.append(int(parts[col_idx["n_samples"]]))

        index = LocusIndex(
            loci=np.array(loci, dtype=object),
            rsids=np.array(rsids, dtype=object)
        )
        return AlleleFreqData(
            ac=np.array(ac, dtype=np.int64),
            an=np.array(an, dtype=np.int64),
            n_samples=np.array(n_samples, dtype=np.int64),
            index=index
        )

    def align_to(self, target_index: LocusIndex) -> "AlleleFreqData":
        """
        Realign data to a different (possibly larger) index.
        Missing loci get zeros.
        """
        n = len(target_index)
        new_ac = np.zeros(n, dtype=np.int64)
        new_an = np.zeros(n, dtype=np.int64)
        new_n_samples = np.zeros(n, dtype=np.int64)

        src_map = self.index.to_pos_map()
        tgt_map = target_index.to_pos_map()

        for locus, src_i in src_map.items():
            if locus in tgt_map:
                tgt_i = tgt_map[locus]
                new_ac[tgt_i] = self.ac[src_i]
                new_an[tgt_i] = self.an[src_i]
                new_n_samples[tgt_i] = self.n_samples[src_i]

        return AlleleFreqData(
            ac=new_ac,
            an=new_an,
            n_samples=new_n_samples,
            index=target_index
        )


def aggregate(*datasets: AlleleFreqData) -> AlleleFreqData:
    """
    Aggregate multiple AlleleFreqData objects.
    Creates union index and sums counts.
    """
    if len(datasets) == 0:
        raise ValueError("Need at least one dataset")

    # Build union index
    union_index = LocusIndex.from_union(*[d.index for d in datasets])

    # Align all datasets and sum
    aligned = [d.align_to(union_index) for d in datasets]

    total_ac = sum(d.ac for d in aligned)
    total_an = sum(d.an for d in aligned)
    total_n_samples = sum(d.n_samples for d in aligned)

    return AlleleFreqData(
        ac=total_ac,
        an=total_an,
        n_samples=total_n_samples,
        index=union_index
    )
