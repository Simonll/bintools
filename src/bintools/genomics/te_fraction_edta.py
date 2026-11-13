#!/usr/bin/env python3
"""
Run EDTA (Extensive de-novo TE Annotator) on a genome and compute the
transposable-element (TE) fraction, for tetrapod groups with no UCSC
RepeatMasker track. De novo, so it captures TEs missed by Dfam-based
RepeatMasker (e.g. large salamander genomes).

(EDTA analogue of te_fraction_repeatmasker.py.)

Pipeline per genome:
  1. (optional) download a genome by NCBI accession (datasets API) into the
     workdir; or pass a local FASTA with --genome.
  2. run EDTA.pl --step <all|anno> --sensitive <0|1> --threads N
     - all  = de novo: build a TE library (RepeatModeler + LTR_retriever)
              then annotate. Slow (weeks on multi-Gb genomes -> benchmark first).
     - anno = fast: annotate with an existing TE library (a *.EDTA.TElib.fa
              in the workdir). Requires that library to already exist.
  3. parse the RepeatMasker-format *.mod.out EDTA produced -> TE fraction
     (classes LINE/SINE/LTR/DNA/RC, intervals merged per (chrom, class)).
  4. append the TE fraction to te_fraction_edta.tsv.

Usage:
  python3 -m bintools.genomics.te_fraction_edta --help
  python3 -m bintools.genomics.te_fraction_edta --genome chr1.fasta \
      --group Lissotriton --species Lissotriton_helveticus \
      --step all --threads 12 --out outputs
  python3 -m bintools.genomics.te_fraction_edta --accession GCF_964261635.1 \
      --group Salamandrininae --species Lissotriton_helveticus --threads 12

Output `*.end` note: EDTA writes `*.<genome>.mod.out` (RepeatMasker-format) and
a `*.EDTA.TEanno.gff3`; we parse the `.out` with the same helper as
RepeatMasker. Non-TE classes (Simple_repeat, Satellite, ...) are excluded.
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
from pathlib import Path

if __package__ in (None, ""):  # run as a script: make `bintools` importable
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from bintools.genomics.te_fraction_repeatmasker import (  # noqa: E402
    datasets_report,
    download_genome,
    genome_size,
)
from bintools.genomics.utils import TE_TSV_HEADER  # noqa: E402
from bintools.genomics.utils import format_te_row  # noqa: E402
from bintools.genomics.utils import parse_rm_out  # noqa: E402


def find_edta_out(genome_fasta: Path, extra_dir: str = "anno") -> Path:
    """Locate the RepeatMasker-format `.out` EDTA produced for `genome_fasta`.

    EDTA 2.x writes the annotation under `<stem>.mod.EDTA.anno/` (plus a
    symlink `<stem>.mod` in the workdir). We glob across likely names so the
    exact layout of a given EDTA release does not matter.
    """
    workdir = genome_fasta.parent
    stem = genome_fasta.stem  # e.g. "chr1.fasta" -> "chr1"
    candidates = [
        workdir / f"{genome_fasta.name}.mod.out",
        workdir / f"{stem}.mod.out",
        workdir / f"{genome_fasta.name}.mod.EDTA.{extra_dir}" / f"{stem}.mod.out",
        workdir / f"{stem}.mod.EDTA.{extra_dir}" / f"{stem}.mod.out",
    ]
    globbed = glob.glob(
        str(workdir / f"*{stem}*.mod{extra_dir}*/{stem}.mod.out")
    ) + glob.glob(str(workdir / f"*{stem}.mod.out"))
    for g in globbed:
        candidates.append(Path(g))
    seen = set()
    for c in candidates:
        if c.exists() and c not in seen:
            return c
    # last resort: anything named *.out under the .EDTA.anno dir
    for d in (workdir / f"{genome_fasta.name}.mod.EDTA.{extra_dir}",):
        for c in sorted(d.glob("*.out")):
            if c.exists() and c not in seen:
                return c
    return None


def run_edta(
    fasta: Path,
    edta: str,
    step: str,
    sensitive: int,
    threads: int,
    anno: int,
) -> Path:
    """Run EDTA.pl and return the `.out` it produced (post-annotation)."""
    cmd = [
        edta,
        "--genome", str(fasta),
        "--species", "others",
        "--step", step,
        "--sensitive", str(sensitive),
        "--threads", str(threads),
        "--anno", str(anno),
    ]
    print("  running: " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    out = find_edta_out(fasta)
    if out is None:
        raise RuntimeError(
            "EDTA finished but no RepeatMasker-format *.out was found. "
            "Look for *<genome>.mod.out or under the .EDTA.anno/ directory."
        )
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Run EDTA and compute TE fraction per genome "
            "(de novo TE annotation), appending to te_fraction_edta.tsv."
        )
    )
    ap.add_argument("--accession", default="", help="NCBI GCF_/GCA_ accession")
    ap.add_argument("--genome", default="", help="local FASTA genome (overrides --accession)")
    ap.add_argument("--species", default="")
    ap.add_argument("--group", default="")
    ap.add_argument("--out", default="outputs")
    ap.add_argument("--workdir", default=None)
    ap.add_argument(
        "--edta",
        default=os.environ.get("EDTA_PL")
        or (
            str(Path(os.environ["EDTA_ENV"]) / "bin" / "EDTA.pl")
            if os.environ.get("EDTA_ENV")
            else "EDTA.pl"
        ),
        help="path to EDTA.pl (default: $EDTA_PL, else $EDTA_ENV/bin/EDTA.pl, "
        "else EDTA.pl from PATH)",
    )
    ap.add_argument("--step", choices=["all", "anno"], default="all",
                    help="all=de novo (construit la librairie, lent) ; "
                         "anno=rapide (exige une *.EDTA.TElib.fa existante)")
    ap.add_argument("--sensitive", type=int, choices=[0, 1], default=0)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--anno", type=int, default=1)
    ap.add_argument(
        "--keep-fna",
        action="store_true",
        help="garder le .fna décompressé et le .fna.gz après le parse "
        "(défaut: les supprimer, re-téléchargeables via l'accession)",
    )
    args = ap.parse_args()

    if not args.genome and not args.accession:
        ap.error("need either --genome <fasta> or --accession <GCF_...>")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    workdir = Path(args.workdir) if args.workdir else out / "work" / (
        args.accession or args.genome
    )
    workdir.mkdir(parents=True, exist_ok=True)

    if args.genome:
        fasta = Path(args.genome)
    else:
        fna = download_genome(args.accession, workdir)
        gsize = genome_size(fna)
        fasta = workdir / (fna.name[:-3])
        if not fasta.exists():
            with __import__("gzip").open(fna, "rt") as fin, open(fasta, "w") as fout:
                shutil.copyfileobj(fin, fout)

    if args.step == "anno":
        # --step anno needs a pre-existing TE library (a *.EDTA.TElib.fa in the
        # workdir). Refuse early if none is found.
        if not glob.glob(str(fasta.parent / "*EDTA*.TElib.fa")):
            ap.error(
                "--step anno requires an existing TE library "
                "(*.EDTA.TElib.fa in the workdir); run --step all first "
                "(de novo) to build it."
            )

    rm_out = run_edta(fasta, args.edta, args.step, args.sensitive, args.threads, args.anno)
    res = parse_rm_out(rm_out)

    if args.genome:
        gsize = sum(len(line.strip()) for line in open(fasta) if not line.startswith(">"))
    frac = res["te_bp"] / gsize if gsize else float("nan")

    # cleanup: same policy as te_fraction_repeatmasker (remove the decompressed
    # .fna and the downloaded .fna.gz, both regenerable, unless --keep-fna)
    if not args.genome and not args.keep_fna:
        fna = workdir / (args.accession and download_genome(args.accession, workdir).name)
        for p in (fasta, fna):
            if p.exists():
                p.unlink()
                print(f"  [cleanup] supprimé {p.name}", flush=True)

    dest = out / "te_fraction_edta.tsv"
    if not dest.exists():
        dest.write_text(TE_TSV_HEADER)
    row = format_te_row(
        args.group,
        args.accession or fasta.name,
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
