# bintools.genomics

Genomics / transposable-element (TE) tools. Each module is importable
(`from bintools.genomics.<module> import ...`) and runnable as a script
(`python3 src/bintools/genomics/<module>.py --help` or
`python3 -m bintools.genomics.<module> --help`).

| Module | Former script | Role | Dependencies |
|---|---|---|---|
| `te_fraction_ucsc.py` | `genotypic_te.py` | TE fraction per genome from UCSC RepeatMasker tracks (`chromInfo` + `rmsk`), for a manifest of UCSC assemblies → `te_fraction.tsv` | stdlib |
| `te_fraction_repeatmasker.py` | `run_repeatmasker_te.py` | Download an NCBI (GCF/GCA) genome, run RepeatMasker, compute TE fraction → appends to `te_fraction_repeatmasker.tsv` | stdlib (+ `RepeatMasker` binary) |
| `consolidate_te_fraction.py` | same | Parse an already-finished RepeatMasker `.out` and append its TE fraction to the consolidated TSV (no re-run) | stdlib |
| `cleanup_rm.py` | same | Post-RepeatMasker cleanup of regenerable intermediates (dry-run unless `--yes`; refuses without a final `.out` unless `--force`) | stdlib |
| `go_enrichment_stable_cpg.py` | same | GO enrichment of genes ranked by CpG-rate variance: ORA via g:Profiler REST + preranked GSEA | pandas, requests, gseapy |
| `find_reps.py` | same | Pick a representative species per tetrapod group (from `<group>.tre`) with a chromosome-level RefSeq assembly (NCBI datasets API) | stdlib |
| `utils.py` | — | Shared helpers: `TE_CLASSES`, `merged_bp`, `parse_rm_out`, `format_te_row`, `get_dict_of_map` (BLAST geneID → symbol; needs pandas) | stdlib (+ pandas for `get_dict_of_map`) |

TE classes counted: LINE, SINE, LTR, DNA, RC (repeat classes such as
Simple_repeat, Satellite, Low_complexity are excluded).
