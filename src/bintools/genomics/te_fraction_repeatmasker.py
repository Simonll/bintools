#!/usr/bin/env python3
"""
Run RepeatMasker on a representative genome and compute the transposable-element
(TE) fraction, for tetrapod groups that have no UCSC RepeatMasker track.

(Formerly run_repeatmasker_te.py.)

Pipeline per genome (accession from manifest):
  1. resolve accession -> assembly name + NCBI FTP path (datasets API)
  2. download *_genomic.fna.gz
  3. compute genome size (sum of sequence lengths)
  4. run RepeatMasker (must be installed; e.g. conda env on the cluster)
  5. parse the .out file -> TE fraction (TE classes LINE/SINE/LTR/DNA/RC,
     intervals merged per (chrom, class))
  6. cleanup: delete the decompressed .fna AND the downloaded .fna.gz
     (both regenerable via the accession) unless --keep-fna

Usage:
  python3 te_fraction_repeatmasker.py --accession GCF_016161935.1 \
      --species Mauremys_reevesii --group Testudines \
      --out outputs --rm "RepeatMasker" --threads 8

RepeatMasker .out columns (space-separated):
  SWscore %div %del %ins query qstart qend qleft strand repeat class/family ...
  The class/family field is like "LINE/L1", "SINE/Alu", "Simple_repeat".
  TE class = part before "/".
"""

import argparse
import gzip
import json
import re
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

if __package__ in (None, ""):  # run as a script: make `bintools` importable
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from bintools.genomics.utils import TE_CLASSES  # noqa: E402,F401  (re-exported)
from bintools.genomics.utils import TE_TSV_HEADER  # noqa: E402
from bintools.genomics.utils import format_te_row  # noqa: E402
from bintools.genomics.utils import parse_rm_out  # noqa: E402

_UA = {"User-Agent": "curl/8.0"}


def datasets_report(accession: str) -> dict:
    url = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/%s/"
        "dataset_report" % accession
    )
    req = urllib.request.Request(url, headers=_UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def ftp_path(accession: str) -> str:
    d = datasets_report(accession)
    rep = d["reports"][0]
    ai = rep["assembly_info"]
    name = ai.get("assembly_name", "")
    # FTP layout: /genomes/all/GCF/016/161/935/GCF_016161935.1_ASM1616193v1/
    a = accession.replace("_", "")
    parts = [a[:3], a[3:6], a[6:9], a[9:12]]
    dirname = f"{accession}_{name}"
    return (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
        + "/".join(parts)
        + "/"
        + dirname
        + "/"
    )


def download_genome(accession: str, workdir: Path) -> Path:
    base = ftp_path(accession)
    req = urllib.request.Request(base, headers=_UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        html = r.read().decode("utf-8", "replace")
    # pick the primary genomic fasta, NOT cds_from_genomic / rna_from_genomic
    cands = re.findall(r'href="([^"]*_genomic\.fna\.gz)"', html)
    cands = [
        c for c in cands if "cds_from_genomic" not in c and "rna_from_genomic" not in c
    ]
    if not cands:
        raise RuntimeError(f"no primary _genomic.fna.gz found at {base}")
    fname = cands[0]
    url = base + fname
    dest = workdir / fname
    if not dest.exists():
        print(f"  downloading {url}", flush=True)
        req = urllib.request.Request(url, headers=_UA)
        with urllib.request.urlopen(req, timeout=3600) as r, open(dest, "wb") as fh:
            shutil.copyfileobj(r, fh)
    return dest


def genome_size(fna_gz: Path) -> int:
    total = 0
    with gzip.open(fna_gz, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def run_repeatmasker(fna_gz: Path, workdir: Path, rm: str, threads: int) -> Path:
    fasta = workdir / (fna_gz.name[:-3])
    if not fasta.exists():
        with gzip.open(fna_gz, "rt") as fin, open(fasta, "w") as fout:
            shutil.copyfileobj(fin, fout)
    outdir = workdir / "rm_out"
    outdir.mkdir(exist_ok=True)
    cmd = [rm, "-pa", str(threads), "-dir", str(outdir), str(fasta)]
    print("  running: " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    return outdir / (fasta.name + ".out")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Download an NCBI genome, run RepeatMasker, compute TE fraction."
    )
    ap.add_argument("--accession", required=True)
    ap.add_argument("--species", default="")
    ap.add_argument("--group", default="")
    ap.add_argument("--out", default="outputs")
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--rm", default="RepeatMasker")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument(
        "--keep-fna",
        action="store_true",
        help="garder le .fna décompressé et le .fna.gz après le parse "
        "(défaut: les supprimer, re-téléchargeables via l'accession)",
    )
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    workdir = Path(args.workdir) if args.workdir else out / "work" / args.accession
    workdir.mkdir(parents=True, exist_ok=True)

    fna = download_genome(args.accession, workdir)
    gsize = genome_size(fna)
    rm_out = run_repeatmasker(fna, workdir, args.rm, args.threads)
    res = parse_rm_out(rm_out)
    frac = res["te_bp"] / gsize if gsize else float("nan")

    # nettoyage: supprimer le .fna décompressé ET le .fna.gz téléchargé
    # (tous deux régénérables via l'accession) sauf --keep-fna
    fasta = workdir / (fna.name[:-3])
    if not args.keep_fna:
        for p in (fasta, fna):
            if p.exists():
                p.unlink()
                print(
                    f"  [cleanup] supprimé {p.name} (re-téléchargeable via {args.accession})",
                    flush=True,
                )

    dest = out / "te_fraction_repeatmasker.tsv"
    if not dest.exists():
        dest.write_text(TE_TSV_HEADER)
    row = format_te_row(
        args.group,
        args.accession,
        args.species,
        gsize,
        res["te_bp"],
        f"{frac:.6f}",
        res["per_class"],
    )
    with dest.open("a") as fh:
        fh.write(row)
    print(
        f"[OK] {args.group}/{args.species} TE={frac:.4f} (genome {gsize})", flush=True
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
