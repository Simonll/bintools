#!/usr/bin/env python3
"""
GO enrichment of tetrapod CpG genes by stability of CpG hypermutability.

Genes are ranked by the variance of their CpG transition rate across the
14 tetrapod groups x 1000 posterior samples (gene_variance.tsv). Low variance
= stable/homogeneous hypermutability.

Two complementary tests:
  A) ORA  (over-representation) via g:Profiler REST, on the LOW- and HIGH-
     variance gene sets (bottom/top percentile), GO BP/MF/CC + KEGG.
  B) GSEA (preranked) via gseapy, on the full ranking by -log10(variance),
     using a GO Biological Process GMT (Baderlab) as the gene-set collection.

Usage:
  python3 go_enrichment_stable_cpg.py --root <project_root> [--cutoff 0.10] [--out DIR]

Inputs (relative to --root unless overridden):
  outputs/blast_out/*.blast           (--blast-glob)  geneID -> symbol
  outputs/genotypic/gene_variance.tsv (--variance)
  data/gmt/<GMT>                      (--gmt)

Outputs (TSV):
  ora_low_gprofiler.tsv    enriched terms in the low-variance (stable) set
  ora_high_gprofiler.tsv   enriched terms in the high-variance set
  gsea_GOBP.tsv            preranked GSEA over the full ranking
  summary.tsv              heads of each result

Dependencies: pandas, requests, gseapy (imported lazily so the module can be
imported / --help'd without them).
"""
import argparse
import sys
from pathlib import Path

if __package__ in (None, ""):  # run as a script: make `bintools` importable
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from bintools.genomics.utils import get_dict_of_map  # noqa: E402

GPROFILER = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
GMT_NAME = "Human_GOBP_AllPathways_noPFOCR_with_GO_iea_October_01_2026_symbol.gmt"


def gost(query, sources=("GO:BP", "GO:MF", "GO:CC", "KEGG"), threshold=0.05):
    import pandas as pd
    import requests

    payload = {
        "organism": "hsapiens",
        "query": list(query),
        "sources": list(sources),
        "user_threshold": threshold,
    }
    r = requests.post(GPROFILER, json=payload, timeout=120)
    r.raise_for_status()
    res = r.json().get("result", [])
    if not res:
        return pd.DataFrame()
    df = pd.DataFrame(res)
    cols = [
        "source",
        "native",
        "name",
        "p_value",
        "significant",
        "intersection_size",
        "term_size",
        "description",
        "intersections",
    ]
    cols = [c for c in cols if c in df.columns]
    return df[cols].sort_values("p_value")


def main():
    ap = argparse.ArgumentParser(
        description="GO enrichment (ORA g:Profiler + preranked GSEA) of stable-CpG genes."
    )
    ap.add_argument(
        "--root",
        default=".",
        help="project root holding outputs/ and data/ (default: current dir). "
        "Replaces the former ROOT = <script dir>/.. which no longer applies "
        "inside the bintools package.",
    )
    ap.add_argument(
        "--cutoff",
        type=float,
        default=0.10,
        help="fraction of extreme genes per oracle (low/high)",
    )
    ap.add_argument(
        "--out", default=None, help="default: <root>/outputs/genotypic/enrichment"
    )
    ap.add_argument(
        "--blast-glob", default=None, help="default: <root>/outputs/blast_out/*.blast"
    )
    ap.add_argument(
        "--variance",
        default=None,
        help="default: <root>/outputs/genotypic/gene_variance.tsv",
    )
    ap.add_argument("--gmt", default=None, help=f"default: <root>/data/gmt/{GMT_NAME}")
    args = ap.parse_args()

    import gseapy
    import pandas as pd

    root = Path(args.root).resolve()
    out = Path(args.out) if args.out else root / "outputs/genotypic/enrichment"
    out.mkdir(parents=True, exist_ok=True)
    blast_glob = args.blast_glob or str(root / "outputs/blast_out/*.blast")
    variance = (
        Path(args.variance)
        if args.variance
        else root / "outputs/genotypic/gene_variance.tsv"
    )
    gmt = Path(args.gmt) if args.gmt else root / "data/gmt" / GMT_NAME

    # --- 1) data: geneID -> symbol + variance ---
    gene2sym = get_dict_of_map(blast_glob)
    var = pd.read_csv(variance, sep="\t")
    var["gene_symbol"] = var["geneID"].map(gene2sym)
    var["score"] = -var["variance"].apply(lambda v: max(v, 1e-9))
    var = var.sort_values("variance")  # ascending: stable first
    symbols = var["gene_symbol"].dropna().tolist()
    print(f"[1] {len(symbols)} genes with symbol")

    # --- 2) ORA: low- vs high-variance gene sets ---
    k = max(1, int(len(symbols) * args.cutoff))
    low, high = symbols[:k], symbols[-k:]
    for name, genes in [("low", low), ("high", high)]:
        d = gost(genes)
        d.to_csv(out / f"ora_{name}_gprofiler.tsv", sep="\t", index=False)
        print(f"[ORA-{name}] {len(genes)} genes, {len(d)} enriched terms")

    # --- 3) GSEA (preranked) over the full ranking by -log10(variance) ---
    rnk = pd.DataFrame({"gene": var["gene_symbol"], "score": var["score"]}).dropna()
    rnk = rnk[~rnk["score"].isin([-float("inf")])]
    prerank = gseapy.prerank(
        rnk=rnk,
        gene_sets=str(gmt),
        outdir=None,
        min_size=5,
        max_size=500,
        permutation_num=1000,
        seed=7,
        threads=4,
        no_plot=True,
    )
    res2d = prerank.res2d.copy()
    res2d.to_csv(out / "gsea_GOBP.tsv", sep="\t", index=False)
    print(f"[GSEA] {len(res2d)} gene-sets reported")

    # --- 4) summary ---
    with open(out / "summary.tsv", "w") as fh:
        fh.write("## ORA low (stable CpG), top5\n")
        fh.write(
            pd.read_csv(out / "ora_low_gprofiler.tsv", sep="\t")
            .head(5)
            .to_csv(sep="\t", index=False)
        )
        fh.write("## GSEA GOBP, top10 (by NES)\n")
        s = res2d.sort_values("NES", key=abs, ascending=False).head(10)
        cols = [c for c in ["Term", "NES", "pval", "fdr", "genes"] if c in s.columns]
        fh.write(s[cols].to_csv(sep="\t", index=False))
    print(f"[DONE] outputs in {out}")


if __name__ == "__main__":
    main()
