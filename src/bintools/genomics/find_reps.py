#!/usr/bin/env python3
"""Resolve representative species per tetrapod group using the SPECIES NAMES
in the group .tre files (format Species_name@internal_id), then find a
chromosome-level RefSeq (GCF) assembly via the NCBI datasets API (by name).

The original script ran everything at import time with hard-coded DATA/GROUPS;
it is now wrapped in main() with argparse (same defaults) so it is importable.

Usage:
  python3 find_reps.py [--data DIR] [--groups Testudines Serpentes ...]
"""
import argparse
import json
import re
import time
import urllib.parse
import urllib.request

DATA = "/mnt/disk4/sll/CpG_tetrapods/bioinfo/workflow/data"
GROUPS = [
    "Testudines",
    "Serpentes",
    "Gekkota",
    "Lacertibaenia",
    "Scinciformata",
    "Hyloidea",
    "Microhylidae",
    "Salamandrininae",
]


def species_names(group, data=DATA):
    with open(f"{data}/{group}.tre") as fh:
        s = fh.read()
    # tokens like  Species_name@12345
    return re.findall(r"([A-Za-z0-9_]+)@\d+", s)


def resolve(name):
    q = urllib.parse.quote(name.replace("_", " "))
    url = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/%s/dataset_report"
        "?page_size=10&assembly_level=chromosome" % q
    )
    req = urllib.request.Request(url, headers={"User-Agent": "curl/8.0"})
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.load(r)


def find_representative(group, data=DATA, sleep=0.3):
    """Print and return (name, accession, assembly_name, organism) of the first
    species of `group` with a chromosome-level GCF assembly, or None."""
    names = species_names(group, data)
    print("=== %s (%d species) ===" % (group, len(names)))
    for name in names:
        try:
            d = resolve(name)
        except Exception:
            continue
        for rep in d.get("reports", []):
            acc = rep.get("accession", "")
            ai = rep.get("assembly_info", {})
            org = rep.get("organism", {}).get("organism_name", "")
            if acc.startswith("GCF_"):
                print(
                    "  %-28s %s  %s  (%s)"
                    % (name, acc, ai.get("assembly_name", ""), org)
                )
                return name, acc, ai.get("assembly_name", ""), org
        time.sleep(sleep)
    print("  (no chromosome-level GCF for any species)")
    return None


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Find a chromosome-level RefSeq (GCF) representative per group."
    )
    ap.add_argument("--data", default=DATA, help="dir containing <group>.tre files")
    ap.add_argument("--groups", nargs="+", default=GROUPS)
    args = ap.parse_args()
    for group in args.groups:
        find_representative(group, args.data)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
