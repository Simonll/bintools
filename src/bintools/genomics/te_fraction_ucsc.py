#!/usr/bin/env python3
"""
Compute transposable-element (TE) fraction per genome from UCSC RepeatMasker
tracks.

(Formerly genotypic_te.py; renamed te_fraction_ucsc.py to make the data
source explicit, by contrast with te_fraction_repeatmasker.py.)

For each UCSC assembly in the manifest:
  - downloads chromInfo.txt.gz (chromosome sizes -> denominator)
  - downloads rmsk.txt.gz      (RepeatMasker annotation, per-repeat rows)

TE fraction = merged bases covered by TE classes / total genome length.

TE classes counted (repClass): LINE, SINE, LTR, DNA, RC.
Non-TE repeat classes (Simple_repeat, Satellite, Low_complexity, RNA, Other,
etc.) are EXCLUDED - the trait is the transposable-element fraction.

Usage:
  python3 te_fraction_ucsc.py --manifest manifest.tsv --out outputs
  python3 -m bintools.genomics.te_fraction_ucsc --manifest manifest.tsv

manifest.tsv (tab-separated, header optional):  genome<TAB>species<TAB>group
Example:
  hg38<TAB>Homo_sapiens<TAB>Primates
  loxAfr3<TAB>Loxodonta_africana<TAB>Afrotheria
"""

import argparse
import gzip
import io
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

if __package__ in (None, ""):  # run as a script: make `bintools` importable
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from bintools.genomics.utils import TE_CLASSES  # noqa: E402
from bintools.genomics.utils import TE_CLASS_ORDER  # noqa: E402
from bintools.genomics.utils import merged_bp  # noqa: E402,F401  (re-exported)

UCSC_BASE = "https://hgdownload.soe.ucsc.edu/goldenPath/{genome}"

_UA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36"}


def chrom_sizes(genome: str) -> dict[str, int]:
    url = f"{UCSC_BASE.format(genome=genome)}/database/chromInfo.txt.gz"
    req = urllib.request.Request(url, headers=_UA)
    with urllib.request.urlopen(req, timeout=60) as resp:
        data = resp.read()
    sizes = {}
    with gzip.open(io.BytesIO(data), "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split()
            if len(f) >= 2 and f[1].isdigit():
                sizes[f[0]] = int(f[1])
    return sizes


def iter_rmsk(genome: str):
    url = f"{UCSC_BASE.format(genome=genome)}/database/rmsk.txt.gz"
    req = urllib.request.Request(url, headers=_UA)
    with urllib.request.urlopen(req, timeout=1200) as resp:
        with gzip.GzipFile(fileobj=resp, mode="rb") as gz:
            for raw in gz:
                line = raw.decode("ascii", "replace").rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                f = line.split()
                if len(f) < 13:
                    continue
                # rmsk: ... index5=genoName 6=genoStart 7=genoEnd 11=repClass
                try:
                    yield f[5], int(f[6]), int(f[7]), f[11]
                except ValueError:
                    continue


def te_fraction(genome: str) -> dict:
    sizes = chrom_sizes(genome)
    genome_size = sum(sizes.values())
    spans = defaultdict(list)
    for chrom, start, end, rep_class in iter_rmsk(genome):
        if rep_class in TE_CLASSES:
            spans[(chrom, rep_class)].append((start, end))
    te_bp = sum(merged_bp(iv) for iv in spans.values())
    frac = te_bp / genome_size if genome_size else float("nan")
    per_class = {
        c: sum(merged_bp(iv) for k, iv in spans.items() if k[1] == c)
        for c in TE_CLASSES
    }
    return {
        "assembly": genome,
        "genome_size": genome_size,
        "te_bp": te_bp,
        "te_fraction": round(frac, 6),
        "per_class": per_class,
    }


def main() -> int:
    ap = argparse.ArgumentParser(
        description="TE fraction per genome from UCSC RepeatMasker (rmsk) tracks."
    )
    ap.add_argument("--manifest", required=True, help="TSV: genome species group")
    ap.add_argument("--out", default="outputs")
    ap.add_argument(
        "--cache",
        default=None,
        help="dir to cache downloaded rmsk.gz (reuse across runs)",
    )
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cache = Path(args.cache) if args.cache else out
    cache.mkdir(parents=True, exist_ok=True)

    rows = []
    with open(args.manifest) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            genome, species = parts[0], parts[1]
            group = parts[2] if len(parts) > 2 else ""
            try:
                r = te_fraction(genome)
                r["species"] = species
                r["group"] = group
                rows.append(r)
                print(
                    f"[OK] {group or genome}/{species} " f"TE={r['te_fraction']:.4f}",
                    flush=True,
                )
            except Exception as exc:  # coverage gaps
                print(
                    f"[SKIP] {genome} ({species}) : {type(exc).__name__} {exc}",
                    file=sys.stderr,
                    flush=True,
                )

    dest = out / "te_fraction.tsv"
    cls = TE_CLASS_ORDER
    with dest.open("w") as fh:
        fh.write(
            "group\tassembly\tspecies\tgenome_size\tte_bp\tte_fraction\t"
            + "\t".join(c + "_pct" for c in cls)
            + "\n"
        )
        for r in rows:
            line = f"{r['group']}\t{r['assembly']}\t{r['species']}\t{r['genome_size']}\t{r['te_bp']}\t{r['te_fraction']}"
            for c in cls:
                line += f"\t{100*r['per_class'].get(c,0)/r['genome_size']:.4f}"
            fh.write(line + "\n")
    print(f"[DONE] wrote {dest} ({len(rows)} genomes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
