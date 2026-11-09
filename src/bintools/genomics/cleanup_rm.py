#!/usr/bin/env python3
"""
Nettoyage post-run RepeatMasker : libère l'espace disque SANS perdre le
travail ni le protocole (science reproductible).

Ce qui est SUPPRIMÉ (régénérable) :
  - <genome>.fna                 (génome brut décompressé, re-téléchargeable via l'accession)
  - <genome>_batch-*.cat/.masked (batches intermédiaires)
  - *.tmp.*, *.err, trfResults-* (fichiers temporaires/vides)

Ce qui est GARDÉ (le travail / le protocole) :
  - <genome>.fna.cat.all.gz      (résultat brut consolidé du run RepeatMasker)
  - <genome>.out et .tbl         (résultat final parsé par consolidate_te_fraction.py)
  - les tsv consolidés (te_fraction*.tsv, manifest, portrait)

SÉCURITÉ : refuse de tourner si le .out final n'existe pas (run pas terminé),
à moins de --force. Sans --yes, ne fait qu'un dry-run.

Usage:
  python3 cleanup_rm.py --rm-dir <dossier_RM_...> \
      --accession GCF_964261635.1
  python3 cleanup_rm.py --rm-dir <dossier_RM_...> \
      --accession GCF_964261635.1 --yes
"""
import argparse
import re
from pathlib import Path

KEEP_GLOBS = ["*.cat.all.gz", "*.out", "*.tbl"]
DELETE_RE = re.compile(
    r"(_batch-.*\.(cat|masked|div|align)|\.fna$|\.tmp\..*$|\.err$|^trfResults|^ncResults|^masks?/)"
)


def plan(rm_dir: Path) -> tuple:
    keep, delete = [], []
    for p in rm_dir.iterdir():
        if not p.is_file():
            continue
        # les trfResults/ncResults sont des temporaires de RepeatMasker (vides) -> DELETE
        if p.name.startswith(("trfResults", "ncResults")):
            delete.append(p)
        elif any(p.match(g) for g in KEEP_GLOBS):
            keep.append(p)
        elif DELETE_RE.search(p.name):
            delete.append(p)
        else:
            keep.append(p)  # fichier inconnu -> on garde par prudence
    return keep, delete


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Nettoyage post-run RepeatMasker (dry-run par défaut)."
    )
    ap.add_argument("--rm-dir", required=True, help="dossier RM_* du run")
    ap.add_argument("--accession", required=True, help="ex GCF_964261635.1")
    ap.add_argument("--yes", action="store_true", help="appliquer (sinon dry-run)")
    ap.add_argument(
        "--force",
        action="store_true",
        help="nettoyer meme si le .out final n'existe pas (dangereux)",
    )
    args = ap.parse_args()

    rd = Path(args.rm_dir)
    if not rd.is_dir():
        print(f"[ERR] {rd} n'est pas un dossier")
        return 1

    out_files = list(rd.glob(f"{args.accession}*.out"))
    out_final = [o for o in out_files if not o.name.startswith("trf")]
    if not out_final and not args.force:
        print(
            f"[BLOQUE] Pas de .out final pour {args.accession} — "
            f"le run n'est probablement pas terminé. (--force pour passer outre)"
        )
        return 2

    keep, delete = plan(rd)
    dels = [
        p
        for p in delete
        if args.accession in p.name or p.name.startswith(("trf", "nc"))
    ]
    del_sz = sum(p.stat().st_size for p in dels)

    print(f"=== {rd} ===")
    print(
        f"REÇU par tamisage : {len(delete)} fichiers (garde tous les autres, {len(keep)})"
    )
    print(f"Supprimerait {len(dels)} fichiers = {del_sz/1e9:.1f} Go")
    print("--- GARDÉS (essentiels) ---")
    for p in sorted([*keep, *[p for p in delete if p not in dels]]):
        print(f"  {p.name}  ({p.stat().st_size/1e6:.1f} Mo)")
    print("--- À SUPPRIMER ---")
    for p in sorted(dels):
        print(f"  {p.name}  ({p.stat().st_size/1e6:.1f} Mo)")

    if not args.yes:
        print("\n[dry-run] Aucun changement. Relance avec --yes pour appliquer.")
        return 0

    freed = 0
    for p in dels:
        sz = p.stat().st_size
        p.unlink()
        freed += sz
    print(f"\n[OK] Supprimé {freed/1e9:.1f} Go. Résultat+tsv intacts.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
