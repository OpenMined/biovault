#!/usr/bin/env python3
import atexit
import importlib.util
import os
import subprocess
import tempfile
import unittest
import gzip


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXTRACT_PATH = os.path.join(ROOT, "assets", "extract_counts.py")
MERGE_PATH = os.path.join(ROOT, "assets", "merge_all_batches.py")

extract_counts = load_module(EXTRACT_PATH, "extract_counts")
merge_all_batches = load_module(MERGE_PATH, "merge_all_batches")

LOG_PATH = os.path.join(os.path.dirname(__file__), "output.txt")

with open(LOG_PATH, "w", encoding="utf-8") as fh:
    fh.write("")

EMIT_ENABLED = False


def _print_log_path() -> None:
    print(f"\n[allele-freq tests] log: {LOG_PATH}")


atexit.register(_print_log_path)


def emit(msg: str) -> None:
    if not EMIT_ENABLED:
        return
    print(msg)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(msg + "\n")


class DosageTests(unittest.TestCase):
    def test_biallelic(self):
        alts = ["A"]
        self.assertEqual(extract_counts.dosages_for_gt(alts, "0/0"), [0])
        self.assertEqual(extract_counts.dosages_for_gt(alts, "0/1"), [1])
        self.assertEqual(extract_counts.dosages_for_gt(alts, "1/1"), [2])
        self.assertEqual(extract_counts.dosages_for_gt(alts, "1|0"), [1])

    def test_multiallelic(self):
        alts = ["A", "C"]
        self.assertEqual(extract_counts.dosages_for_gt(alts, "0/2"), [0, 1])
        self.assertEqual(extract_counts.dosages_for_gt(alts, "1/2"), [1, 1])
        self.assertEqual(extract_counts.dosages_for_gt(alts, "2/2"), [0, 2])

    def test_missing_invalid(self):
        alts = ["A", "C"]
        self.assertIsNone(extract_counts.dosages_for_gt(alts, "./."))
        self.assertIsNone(extract_counts.dosages_for_gt(alts, ".|."))
        self.assertIsNone(extract_counts.dosages_for_gt(alts, "A/B"))
        self.assertIsNone(extract_counts.dosages_for_gt(alts, "0/3"))

    def test_extract_counts_vcf_fixture(self):
        vcf_text = "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                "1\t100\trs1\tG\tA\t.\tPASS\t.\tGT\t0/1",
                "1\t101\trs2\tG\tA\t.\tPASS\t.\tGT\t1/1",
                "1\t102\trs3\tG\tA\t.\tPASS\t.\tGT\t0/0",
                "1\t103\trs4\tG\tA,C\t.\tPASS\t.\tGT\t0/2",
                "1\t104\trs5\tG\tA,C\t.\tPASS\t.\tGT\t1/2",
                "1\t105\trs6\tG\tT\t.\tPASS\t.\tGT\t./.",
                "1\t106\trs7\tG\tT\t.\tPASS\t.\tGT\tA/B",
                "1\t107\trs8\tG\tT\t.\tPASS\t.\tGT\t0/3",
                "1\t108\trs9\tG\tGA\t.\tPASS\t.\tGT\t0/1",
                "1\t109\trs10\tG\tA,GA\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        with tempfile.TemporaryDirectory() as tmp:
            vcf_path = os.path.join(tmp, "mini.vcf")
            out_path = os.path.join(tmp, "out.tsv")
            with open(vcf_path, "w", encoding="utf-8") as fh:
                fh.write(vcf_text + "\n")

            result = subprocess.run(
                [
                    "python3",
                    EXTRACT_PATH,
                    "--vcf",
                    vcf_path,
                    "--participant",
                    "p1",
                    "--output",
                    out_path,
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            self.assertEqual(
                result.returncode,
                0,
                msg=f"extract_counts failed: {result.stderr}",
            )

            with open(out_path, encoding="utf-8") as fh:
                out_text = fh.read()
            lines = [line.strip().split("\t") for line in out_text.splitlines() if line.strip()]

            self.assertEqual(lines[0], ["locus_key", "participant_id", "dosage"])
            rows = {row[0]: row[2] for row in lines[1:]}

            expected = {
                "1-100-G-A": "1",
                "1-101-G-A": "2",
                "1-102-G-A": "0",
                "1-103-G-A": "0",
                "1-103-G-C": "1",
                "1-104-G-A": "1",
                "1-104-G-C": "1",
                "1-105-G-T": "-1",
                "1-106-G-T": "-1",
                "1-107-G-T": "-1",
            }
            self.assertEqual(rows, expected)

    def test_genotype_to_vcf_integration(self):
        if subprocess.run(["which", "bvs"], stdout=subprocess.DEVNULL).returncode != 0:
            self.skipTest("bvs not available")
        repo_root = os.path.abspath(os.path.join(ROOT, "../../../.."))
        candidates = [
            os.path.join(repo_root, "data", "genostats.sqlite"),
            os.path.join(repo_root, "tests", "scenarios", "allele-freq", "data", "genostats.sqlite"),
            os.path.join(os.path.expanduser("~/.biovault"), "genostats.sqlite"),
            os.path.join(os.getcwd(), "data", "genostats.sqlite"),
        ]
        sqlite = next((path for path in candidates if os.path.exists(path)), None)
        if not sqlite:
            self.skipTest("genostats.sqlite not available")

        geno_text = "\n".join(
            [
                "# test genotype file",
                "rs33985472\t11\t5225485\tTC\t0.50\t0.10\t0.01",
                "rs33971634\t11\t5225660\tGC\t0.60\t0.20\t0.02",
                "rs73885319\t22\t36265860\tGG\t0.70\t0.30\t0.03",
                "rs60910145\t22\t36265988\t--\t0.80\t0.40\t0.04",
            ]
        )
        with tempfile.TemporaryDirectory() as tmp:
            geno_path = os.path.join(tmp, "mini.txt")
            vcf_path = os.path.join(tmp, "mini.vcf.gz")
            out_path = os.path.join(tmp, "out.tsv")
            with open(geno_path, "w", encoding="utf-8") as fh:
                fh.write(geno_text + "\n")

            result = subprocess.run(
                [
                    "bvs",
                    "genotype-to-vcf",
                    "--input",
                    geno_path,
                    "--output",
                    vcf_path,
                    "--gzip",
                    "--sqlite",
                    sqlite,
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            self.assertEqual(result.returncode, 0, msg=f"bvs genotype-to-vcf failed: {result.stderr}")

            with gzip.open(vcf_path, "rt", encoding="utf-8") as fh:
                vcf_text = "".join(fh.readlines())
            result = subprocess.run(
                [
                    "python3",
                    EXTRACT_PATH,
                    "--vcf",
                    vcf_path,
                    "--participant",
                    "p1",
                    "--output",
                    out_path,
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            self.assertEqual(result.returncode, 0, msg=f"extract_counts failed: {result.stderr}")

            with open(out_path, encoding="utf-8") as fh:
                out_text = fh.read()
            rows = [line.strip().split("\t") for line in out_text.splitlines() if line.strip()]
            data = {row[0]: row[2] for row in rows[1:]}

            expected = {
                "11-5225485-T-C": "1",
                "11-5225485-T-A": "0",
                "11-5225660-N-.": "-1",
                "22-36265860-A-G": "2",
                "22-36265988-T-C": "-1",
                "22-36265988-T-G": "-1",
            }
            for locus, dosage in expected.items():
                self.assertEqual(data.get(locus), dosage, msg=f"Unexpected dosage for {locus}")

    def test_full_pipeline_with_af(self):
        try:
            import numpy  # noqa: F401
        except Exception:
            self.skipTest("numpy not available")
        if subprocess.run(["which", "bvs"], stdout=subprocess.DEVNULL).returncode != 0:
            self.skipTest("bvs not available")
        repo_root = os.path.abspath(os.path.join(ROOT, "../../../.."))
        candidates = [
            os.path.join(repo_root, "data", "genostats.sqlite"),
            os.path.join(repo_root, "tests", "scenarios", "allele-freq", "data", "genostats.sqlite"),
            os.path.join(os.path.expanduser("~/.biovault"), "genostats.sqlite"),
            os.path.join(os.getcwd(), "data", "genostats.sqlite"),
        ]
        sqlite = next((path for path in candidates if os.path.exists(path)), None)
        if not sqlite:
            self.skipTest("genostats.sqlite not available")

        geno_p1 = "\n".join(
            [
                "# p1 genotype file",
                "rs33985472\t11\t5225485\tTC\t0.50\t0.10\t0.01",
                "rs73885319\t22\t36265860\tGG\t0.70\t0.30\t0.03",
                "rs60910145\t22\t36265988\tCG\t0.80\t0.40\t0.04",
            ]
        )
        geno_p2 = "\n".join(
            [
                "# p2 genotype file",
                "rs33985472\t11\t5225485\tTT\t0.40\t0.20\t0.02",
                "rs73885319\t22\t36265860\tAA\t0.60\t0.10\t0.01",
                "rs60910145\t22\t36265988\t--\t0.90\t0.50\t0.05",
            ]
        )
        global EMIT_ENABLED
        EMIT_ENABLED = True
        try:
            emit("\n[Genotype inputs]")
            emit(geno_p1)
            emit(geno_p2)

            with tempfile.TemporaryDirectory() as tmp:
                geno1 = os.path.join(tmp, "p1.txt")
                geno2 = os.path.join(tmp, "p2.txt")
                vcf1 = os.path.join(tmp, "p1.vcf.gz")
                vcf2 = os.path.join(tmp, "p2.vcf.gz")
                c1 = os.path.join(tmp, "p1_counts.tsv")
                c2 = os.path.join(tmp, "p2_counts.tsv")
                batch = os.path.join(tmp, "batch.tsv")
                matrix = os.path.join(tmp, "dosage_matrix.tsv")
                npz = os.path.join(tmp, "dosage_matrix.npz")
                loci = os.path.join(tmp, "locus_index.txt")
                participants = os.path.join(tmp, "participants.txt")
                af_tsv = os.path.join(tmp, "allele_freq.tsv")

                with open(geno1, "w", encoding="utf-8") as fh:
                    fh.write(geno_p1 + "\n")
                with open(geno2, "w", encoding="utf-8") as fh:
                    fh.write(geno_p2 + "\n")

                for geno, vcf in ((geno1, vcf1), (geno2, vcf2)):
                    result = subprocess.run(
                        [
                            "bvs",
                            "genotype-to-vcf",
                            "--input",
                            geno,
                            "--output",
                            vcf,
                            "--gzip",
                            "--sqlite",
                            sqlite,
                        ],
                        check=False,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                    )
                    self.assertEqual(result.returncode, 0, msg=f"bvs genotype-to-vcf failed: {result.stderr}")

                emit("\n[VCF outputs]")
                emit(gzip.open(vcf1, "rt", encoding="utf-8").read())
                emit(gzip.open(vcf2, "rt", encoding="utf-8").read())

                for pid, vcf, out in (("p1", vcf1, c1), ("p2", vcf2, c2)):
                    result = subprocess.run(
                        [
                            "python3",
                            EXTRACT_PATH,
                            "--vcf",
                            vcf,
                            "--participant",
                            pid,
                            "--output",
                            out,
                        ],
                        check=False,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                    )
                    self.assertEqual(result.returncode, 0, msg=f"extract_counts failed: {result.stderr}")

                result = subprocess.run(
                    [
                        "python3",
                        os.path.join(ROOT, "assets", "aggregate_batch.py"),
                        "--batch-id",
                        "batch1",
                        "--counts",
                        f"{c1} {c2}",
                        "--output",
                        batch,
                    ],
                    check=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                self.assertEqual(result.returncode, 0, msg=f"aggregate_batch failed: {result.stderr}")

                result = subprocess.run(
                    [
                        "python3",
                        MERGE_PATH,
                        "--batch-files",
                        batch,
                        "--matrix-tsv",
                        matrix,
                        "--npz",
                        npz,
                        "--loci",
                        loci,
                        "--participants",
                        participants,
                    ],
                    check=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                self.assertEqual(result.returncode, 0, msg=f"merge_all_batches failed: {result.stderr}")

                emit("\n[Counts matrix]")
                with open(matrix, encoding="utf-8") as fh:
                    emit(fh.read())

                with open(matrix, encoding="utf-8") as fh:
                    header = next(fh).strip().split("\t")
                    pids = header[1:]
                    rows = [line.strip().split("\t") for line in fh if line.strip()]

                result = subprocess.run(
                    [
                        "python3",
                        os.path.join(ROOT, "assets", "calc_allele_freq.py"),
                        "--matrix",
                        matrix,
                        "--output",
                        af_tsv,
                        "--emit-equations",
                    ],
                    check=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                self.assertEqual(result.returncode, 0, msg=f"calc_allele_freq failed: {result.stderr}")

                with open(af_tsv, encoding="utf-8") as fh:
                    af_text = fh.read()
                emit("\n[Allele frequencies]")
                emit(af_text)
                af_lines = [line.strip().split("\t") for line in af_text.splitlines() if line.strip()]
                af_rows = {row[0]: row for row in af_lines[1:]}

                expected_af = {
                    "11-5225485-T-C": 0.25,
                    "22-36265860-A-G": 0.50,
                    "22-36265988-T-C": 0.50,
                    "22-36265988-T-G": 0.50,
                }
                for locus, expected in expected_af.items():
                    row = af_rows.get(locus)
                    self.assertIsNotNone(row, msg=f"Missing locus in allele_freq.tsv: {locus}")
                    found = float(row[4])
                    self.assertAlmostEqual(found, expected, places=4, msg=f"AF mismatch for {locus}")
        finally:
            EMIT_ENABLED = False

class MergeTests(unittest.TestCase):
    def test_merge_missing_defaults(self):
        try:
            import numpy  # noqa: F401
        except Exception:
            self.skipTest("numpy not available")
        with tempfile.TemporaryDirectory() as tmp:
            batch1 = os.path.join(tmp, "batch1.tsv")
            batch2 = os.path.join(tmp, "batch2.tsv")
            matrix_tsv = os.path.join(tmp, "matrix.tsv")
            npz = os.path.join(tmp, "matrix.npz")
            loci = os.path.join(tmp, "loci.txt")
            participants = os.path.join(tmp, "participants.txt")

            with open(batch1, "w", encoding="utf-8") as fh:
                fh.write("locus_key\tparticipant_id\tdosage\n")
                fh.write("10:1:G:A\tp1\t2\n")
                fh.write("10:2:G:C\tp1\t0\n")
                fh.write("10:3:G:T\tp1\t-1\n")

            with open(batch2, "w", encoding="utf-8") as fh:
                fh.write("locus_key\tparticipant_id\tdosage\n")
                fh.write("10:2:G:C\tp2\t1\n")
                fh.write("10:3:G:T\tp1\t1\n")

            result = subprocess.run(
                [
                    "python3",
                    MERGE_PATH,
                    "--batch-files",
                    f"{batch1} {batch2}",
                    "--matrix-tsv",
                    matrix_tsv,
                    "--npz",
                    npz,
                    "--loci",
                    loci,
                    "--participants",
                    participants,
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            self.assertEqual(
                result.returncode,
                0,
                msg=f"merge_all_batches failed: {result.stderr}",
            )

            with open(matrix_tsv, encoding="utf-8") as fh:
                out_text = fh.read()
            lines = [line.strip().split("\t") for line in out_text.splitlines() if line.strip()]

            header = lines[0]
            self.assertEqual(header, ["locus_key", "p1", "p2"])
            rows = {row[0]: row[1:] for row in lines[1:]}
            self.assertEqual(rows["10:1:G:A"], ["2", "-1"])
            self.assertEqual(rows["10:2:G:C"], ["0", "1"])
            self.assertEqual(rows["10:3:G:T"], ["1", "-1"])


if __name__ == "__main__":
    unittest.main()
