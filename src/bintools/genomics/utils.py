"""Shared helpers for bintools.genomics.

- TE (transposable element) helpers shared by te_fraction_ucsc,
  te_fraction_repeatmasker and consolidate_te_fraction (previously duplicated
  verbatim in each script).
- get_dict_of_map: copied from CpG_tetrapods/bioinfo/workflow/scripts/utils.py
  (it does not exist in bintools.utils). Only this function was copied, to
  avoid pulling that module's heavy dependencies (ete3, seaborn, sklearn...).
"""
import glob
import io
from collections import defaultdict
from pathlib import Path
from typing import Dict

# repClass values counted as transposable elements. Non-TE repeat classes
# (Simple_repeat, Satellite, Low_complexity, RNA, Other, ...) are excluded.
TE_CLASSES = {"LINE", "SINE", "LTR", "DNA", "RC"}

# fixed column order for the *_pct columns of the te_fraction TSV outputs
TE_CLASS_ORDER = ["LINE", "SINE", "LTR", "DNA", "RC"]

TE_TSV_HEADER = (
    "group\tassembly\tspecies\tgenome_size\tte_bp\tte_fraction\t"
    + "\t".join(c + "_pct" for c in TE_CLASS_ORDER)
    + "\n"
)


def merged_bp(intervals):
    """Total bp covered by a list of (start, end) intervals, overlaps merged.

    Sorts `intervals` in place.
    """
    intervals.sort()
    total = 0
    cs = ce = None
    for s, e in intervals:
        if cs is None:
            cs, ce = s, e
        elif s <= ce:
            ce = max(ce, e)
        else:
            total += ce - cs
            cs, ce = s, e
    if cs is not None:
        total += ce - cs
    return total


def parse_rm_out(out_file: Path) -> dict:
    """Return {te_bp, per_class} from a RepeatMasker .out file.

    RepeatMasker .out columns (space-separated):
      0=SWscore 1=%div 2=%del 3=%ins 4=query 5=qstart 6=qend
      7=qleft 8=strand 9=repeat 10=class/family ...
    TE class = part of class/family before "/". Intervals are merged per
    (chrom, class).
    """
    spans = defaultdict(list)
    with open(out_file) as fh:
        for line in fh:
            s = line.strip()
            if (
                not s
                or s.startswith("SW")
                or s.startswith("score")
                or s.startswith("There were")
            ):
                continue
            f = line.split()
            if len(f) < 12:
                continue
            chrom, start, end = f[4], int(f[5]), int(f[6])
            class_family = f[10]
            rep_class = class_family.split("/")[0]
            if rep_class in TE_CLASSES:
                spans[(chrom, rep_class)].append((start, end))
    te_bp = sum(merged_bp(iv) for iv in spans.values())
    per_class = {
        c: sum(merged_bp(iv) for k, iv in spans.items() if k[1] == c)
        for c in TE_CLASSES
    }
    return {"te_bp": te_bp, "per_class": per_class}


def format_te_row(
    group: str,
    assembly: str,
    species: str,
    gsize: int,
    te_bp: int,
    frac_str: str,
    per_class: dict,
) -> str:
    """One TSV line (with trailing newline) matching TE_TSV_HEADER."""
    row = f"{group}\t{assembly}\t{species}\t{gsize}\t{te_bp}\t{frac_str}"
    for c in TE_CLASS_ORDER:
        row += f"\t{100 * per_class.get(c, 0) / gsize:.4f}"
    return row + "\n"


def get_dict_of_map(input: str) -> Dict[str, str]:
    """Map geneID -> HGNC symbol from BLAST tabular outputs matching glob `input`.

    geneID is parsed from the file name (last "_" field, before "-M2018");
    the symbol is taken from the best hit (highest bit score) subject id,
    field [-2] when split on ":".
    """
    import pandas as pd

    dict_of_map: Dict[str, str] = {}
    files = glob.glob(input)
    names = [
        "query id",
        "subject id",
        "% identity",
        "alignment length",
        "mismatches",
        "gap opens",
        "q. start",
        "q. end",
        "s. start",
        "s. end",
        "evalue",
        "bit score",
    ]
    for f in files:
        geneID = f.split("/")[-1].split("_")[-1].split("-M2018")[0]
        with open(f, "r") as fh:
            lines = fh.readlines()
            k = 0
            l = lines[k]
            while l.startswith("#"):
                k += 1
                l = lines[k]
            s = ""
            for i in lines[k:-1]:
                s += i.__str__()
            df = pd.read_csv(io.StringIO(s), sep="\t", names=names)
            hgcn = (
                df.sort_values("bit score", ascending=False).iloc[0, 1].split(":")[-2]
            )
            dict_of_map[geneID] = hgcn
    return dict_of_map
