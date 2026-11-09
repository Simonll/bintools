#!/usr/bin/env python3
"""
Consolider la fraction TE d'un groupe RepeatMasker déjà exécuté (sans relancer).

Utile pour le 14e groupe (Salamandrininae / Lissotriton helveticus) dont le run
RepeatMasker tourne sur le cluster : dès que le fichier .out est complet, ce
script parse le .out existant et ajoute la ligne au tsv consolidé
te_fraction_repeatmasker.tsv — sans re-télécharger ni re-lancer RepeatMasker.

Usage:
  python3 consolidate_te_fraction.py \
      --out /mnt/disk4/sll/CpG_tetrapods/bioinfo/workflow/outputs/genotypic/te_fraction_repeatmasker.tsv \
      --rm-out /mnt/disk4/sll/CpG_tetrapods/bioinfo/workflow/RM_3026915.TueOct200656282026/GCF_964261635.1_aLisHel1.1_genomic.fna.out \
      --genome-size 23633485539 \
      --group Salamandrininae --accession GCF_964261635.1 --species Lissotriton_helveticus

Le .out de RepeatMasker est le fichier <genome>.out produit dans le répertoire
de sortie du run (pas les .cat/.masked par lots). S'il n'existe pas encore, le
script le signale et sort en code 2 (à relancer plus tard).
"""
import argparse
import sys
from pathlib import Path

if __package__ in (None, ""):  # run as a script: make `bintools` importable
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from bintools.genomics.utils import TE_CLASSES  # noqa: E402,F401  (re-exported)
from bintools.genomics.utils import TE_TSV_HEADER  # noqa: E402
from bintools.genomics.utils import format_te_row  # noqa: E402
from bintools.genomics.utils import parse_rm_out  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Ajoute la fraction TE d'un .out RepeatMasker existant au tsv consolidé."
    )
    ap.add_argument(
        "--out", required=True, help="tsv consolidé te_fraction_repeatmasker.tsv"
    )
    ap.add_argument(
        "--rm-out", required=True, help="chemin du fichier .out de RepeatMasker"
    )
    ap.add_argument("--genome-size", type=int, required=True)
    ap.add_argument("--group", required=True)
    ap.add_argument("--accession", required=True)
    ap.add_argument("--species", default="")
    args = ap.parse_args()

    rm_out = Path(args.rm_out)
    if not rm_out.exists():
        print(
            f"[WAIT] {rm_out.name} n'existe pas encore — run RepeatMasker "
            f"pas terminé. Relance ce script quand le .out est là.",
            file=sys.stderr,
        )
        return 2

    res = parse_rm_out(rm_out)
    gsize = args.genome_size
    frac = res["te_bp"] / gsize if gsize else float("nan")

    dest = Path(args.out)
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
    print(f"[OK] {args.group}/{args.species} TE={frac:.4f} (genome {gsize})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
