import os
from typing import IO
from typing import Any
from typing import Dict
from typing import Iterator
from typing import List
from typing import TextIO

import pandas as pd
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord


class ali:
    def __init__(
        self,
        n_taxa: int,
        n_site: int,
        dict_of_seq: Dict[str, Dict[int, str]],
        genetic_code: str,
    ) -> None:
        self.genetic_code: str = genetic_code
        self.n_taxa: int = n_taxa
        self.n_site: int = n_site
        self.dict_of_seq: Dict[str, Dict[int, str]] = dict_of_seq

    def set_genetic_code(self, table: str):
        self.genetic_code = table

    def get_genetic_code(self) -> str:
        return self.genetic_code

    def get_n_taxa(self) -> int:
        return self.n_taxa

    def get_n_site(self) -> int:
        return self.n_site

    def get_dict_of_seq(self) -> Dict[str, Dict[int, str]]:
        return self.dict_of_seq

    def __iter__(self) -> Iterator[Any]:
        yield self.dict_of_seq.items()

    def get_biopython_align(self) -> MultipleSeqAlignment:
        records: List[SeqRecord] = []
        for k, v in self.dict_of_seq.items():
            records += [SeqRecord(MutableSeq("".join([j for i, j in v.items()])), id=k)]
        return MultipleSeqAlignment(records=records)

    def get_biopython_align_codon2aa(self) -> MultipleSeqAlignment:
        records: List[SeqRecord] = []
        for k, v in self.dict_of_seq.items():
            records += [SeqRecord(MutableSeq("".join([j for i, j in v.items()])), id=k)]
        return MultipleSeqAlignment(records=records)


def sub_sample(seq: str, start: int, stop: int, step: int) -> str:
    return "".join([seq[l] for l in range(start, stop, step)])


def write_fasta(fh: IO[Any], defline: str, content: str) -> bool:
    try:
        fh.write(">" + defline + "\n")
        fh.write(content + "\n")
        return True
    except Exception as e:
        print("Something wrong writing fasta %s" % str(e))
        return False


def write_phylip_header(fh: IO[Any], Ntaxa: int, Nsite: int) -> bool:
    try:
        fh.write(str(Ntaxa) + "\t" + str(Nsite) + "\n")
        return True
    except Exception as e:
        print("Something wrong writing phylip header %s" % str(e))
        return False


def write_phylip_line(fh: IO[Any], defline: str, content: str) -> bool:
    try:
        fh.write(defline + "  " + content + "\n")
        return True
    except Exception as e:
        print("Something wrong writing phylip line %s" % str(e))
        return False


def write_fasta_from_align(filename: str, align: MultipleSeqAlignment) -> bool:
    try:
        for record in align:
            if os.path.exists(filename):
                append_write = "a"  # append if already exists
            else:
                append_write = "w"  # make a new file if not

            with open(filename, append_write) as fh:
                write_fasta(fh=fh, defline=record.id, content=str(record.seq))
        return True
    except Exception as e:
        print("Something wrong writing fasta from MultipleSeqAlignment %s" % str(e))
        return False


def write_phylip(filename: str, align: MultipleSeqAlignment):
    try:
        with open(filename, "w") as fh:
            write_phylip_header(
                fh=fh, Ntaxa=len(align), Nsite=align.get_alignment_length()
            )
            for record in align:
                write_phylip_line(fh=fh, defline=record.id, content=str(record.seq))
        return True
    except Exception as e:
        print("Something wrong writing phylip from MultipleSeqAlignment %s" % str(e))
        return False


def read_ali(fh: IO[Any]):
    lines: List[str] = fh.readlines()
    lines: Iterator[str] = iter(lines)
    dict_of_seq: Dict[str, Dict[int, str]] = {}
    for l in lines:
        if l.startswith("#"):
            pass
        elif l.startswith(">"):
            defline = l.strip()[1:]
            seq = next(lines)
            dict_of_seq[defline.strip()] = {
                i: j for i, j in enumerate(list(seq.strip()))
            }
    df: pd.DataFrame = pd.DataFrame.from_dict(data=dict_of_seq, orient="index")

    return ali(
        n_site=df.shape[1],
        n_taxa=df.shape[0],
        dict_of_seq=dict_of_seq,
        genetic_code="standard",
    )


def read_fasta(fh: IO[Any]) -> ali:
    lines: List[str] = fh.readlines()
    lines: Iterator[str] = iter(lines)
    dict_of_seq: Dict[str, Dict[int, str]] = {}
    for l in lines:
        if l.startswith(">"):
            defline = l.strip()[1:]
            seq = next(lines)
            dict_of_seq[defline.strip()] = {
                i: j for i, j in enumerate(list(seq.strip()))
            }
    df: pd.DataFrame = pd.DataFrame.from_dict(data=dict_of_seq, orient="index")

    return ali(
        n_site=df.shape[1],
        n_taxa=df.shape[0],
        dict_of_seq=dict_of_seq,
        genetic_code="Standard",
    )


def read_phylip(fh: TextIO) -> ali:
    lines: List[str] = fh.readlines()
    n_taxa, n_site = lines[0].split("\t")
    dict_of_seq: Dict[str, Dict[int, str]] = {}
    for l in range(1, len(lines)):
        sp_name, seq = lines[l].split("  ")
        dict_of_seq[sp_name.strip()] = {i: j for i, j in enumerate(list(seq.strip()))}

    return ali(
        n_site=int(n_site),
        n_taxa=int(n_taxa),
        dict_of_seq=dict_of_seq,
        genetic_code="standard",
    )


def find_missing_col(
    ali: ali, tripplet_patterns: List[str], ratio: float = 1
) -> List[int]:
    for i in tripplet_patterns:
        if len(i) != 3:
            raise RuntimeError

    codon_missing_count: Dict[int, int] = {i: 0 for i in range(ali.get_n_site())}

    for sp, dict_of_seq in ali:
        seq: str = "".join(list(dict_of_seq.values()))
        for codon, i, j, k in zip(
            range(len(sub_sample(seq, start=0, stop=len(seq), step=3))),
            sub_sample(seq, start=0, stop=len(seq), step=3),
            sub_sample(seq, start=1, stop=len(seq), step=3),
            sub_sample(seq, start=2, stop=len(seq), step=3),
        ):
            if "".join([i, j, k]) in tripplet_patterns:
                codon_missing_count[codon] += 1
    site_to_keep: List[int] = []
    for l, m in codon_missing_count.items():
        if m <= (1 - ratio) * ali.get_n_taxa():
            site_to_keep += [l]
    return site_to_keep


def codon_to_nuc(list_of_codon: List[int]) -> List[int]:
    list_of_nuc: List[int] = []
    for i in list_of_codon:
        list_of_nuc += [i * 3]
        list_of_nuc += [i * 3 + 1]
        list_of_nuc += [i * 3 + 2]

    return list_of_nuc


def select_positions(dict_ali: Dict[str, Any], positions: List[int]) -> Dict[str, Any]:
    df: pd.DataFrame = pd.DataFrame.from_dict(
        data=dict_ali["alignment"], orient="index"
    )
    positions_nuc: List[int] = codon_to_nuc(list_of_codon=positions)
    df = df.iloc[:, positions_nuc]
    return {
        "genetic_code": dict_ali["genetic_code"],
        "n_taxa": df.shape[0],
        "n_site": df.shape[1],
        "alignment": df.to_dict(orient="index"),
    }
