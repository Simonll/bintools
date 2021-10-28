from typing import Dict
from typing import List

precision: int = 10000
amino_acids_upper: List[str] = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
    "-",
]
amino_acids_lower: List[str] = [
    "a",
    "c",
    "d",
    "e",
    "f",
    "g",
    "h",
    "i",
    "k",
    "l",
    "m",
    "n",
    "p",
    "q",
    "r",
    "s",
    "t",
    "v",
    "w",
    "y",
    "-",
]
dna_upper: List[str] = ["A", "C", "G", "T"]
dna_lower: List[str] = ["a", "c", "g", "t"]
rna_upper: List[str] = ["A", "C", "G", "U"]
rna_lower: List[str] = ["a", "c", "g", "u"]
codons: List[str] = [
    "TTT",
    "TTC",
    "TTA",
    "TTG",
    "TCT",
    "TCC",
    "TCA",
    "TCG",
    "TAT",
    "TAC",
    "TAA",
    "TAG",
    "TGT",
    "TGC",
    "TGA",
    "TGG",
    "CTT",
    "CTC",
    "CTA",
    "CTG",
    "CCT",
    "CCC",
    "CCA",
    "CCG",
    "CAT",
    "CAC",
    "CAA",
    "CAG",
    "CGT",
    "CGC",
    "CGA",
    "CGG",
    "ATT",
    "ATC",
    "ATA",
    "ATG",
    "ACT",
    "ACC",
    "ACA",
    "ACG",
    "AAT",
    "AAC",
    "AAA",
    "AAG",
    "AGT",
    "AGC",
    "AGA",
    "AGG",
    "GTT",
    "GTC",
    "GTA",
    "GTG",
    "GCT",
    "GCC",
    "GCA",
    "GCG",
    "GAT",
    "GAC",
    "GAA",
    "GAG",
    "GGT",
    "GGC",
    "GGA",
    "GGG",
]
codons_position: List[List[int]] = [
    [
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
    ],
    [
        3,
        3,
        3,
        3,
        1,
        1,
        1,
        1,
        0,
        0,
        0,
        0,
        2,
        2,
        2,
        2,
        3,
        3,
        3,
        3,
        1,
        1,
        1,
        1,
        0,
        0,
        0,
        0,
        2,
        2,
        2,
        2,
        3,
        3,
        3,
        3,
        1,
        1,
        1,
        1,
        0,
        0,
        0,
        0,
        2,
        2,
        2,
        2,
        3,
        3,
        3,
        3,
        1,
        1,
        1,
        1,
        0,
        0,
        0,
        0,
        2,
        2,
        2,
        2,
    ],
    [
        3,
        1,
        0,
        2,
        3,
        1,
        0,
        2,
        3,
        1,
        0,
        2,
        3,
        1,
        0,
        2,
        3,
        1,
        0,
        2,
        3,
        1,
        0,
        2,
        3,
        1,
        0,
        2,
        3,
        1,
        0,
        2,
        3,
        1,
        0,
        2,
        3,
        1,
        0,
        2,
        3,
        1,
        0,
        2,
        3,
        1,
        0,
        2,
        3,
        1,
        0,
        2,
        3,
        1,
        0,
        2,
        3,
        1,
        0,
        2,
        3,
        1,
        0,
        2,
    ],
]

nucleotide_rr: List[str] = ["AC", "AG", "AT", "CG", "CT", "GT"]
number_of_amino_acids: int = 20
number_of_nucleotides: int = 4
number_of_codons: int = 64
gen_code_standard: List[int] = [
    4,
    4,
    9,
    9,
    15,
    15,
    15,
    15,
    19,
    19,
    -1,
    -1,
    1,
    1,
    -1,
    18,
    9,
    9,
    9,
    9,
    12,
    12,
    12,
    12,
    6,
    6,
    13,
    13,
    14,
    14,
    14,
    14,
    7,
    7,
    7,
    10,
    16,
    16,
    16,
    16,
    11,
    11,
    8,
    8,
    15,
    15,
    14,
    14,
    17,
    17,
    17,
    17,
    0,
    0,
    0,
    0,
    2,
    2,
    3,
    3,
    5,
    5,
    5,
    5,
]
gen_code_standard_deg_lvl: List[int] = [
    4,
    2,
    2,
    2,
    2,
    4,
    2,
    3,
    2,
    6,
    1,
    2,
    4,
    2,
    6,
    6,
    4,
    4,
    1,
    2,
]
gen_code_mito_invert: List[int] = [
    4,
    4,
    9,
    9,
    15,
    15,
    15,
    15,
    19,
    19,
    -1,
    -1,
    1,
    1,
    18,
    18,
    9,
    9,
    9,
    9,
    12,
    12,
    12,
    12,
    6,
    6,
    13,
    13,
    14,
    14,
    14,
    14,
    7,
    7,
    10,
    10,
    16,
    16,
    16,
    16,
    11,
    11,
    8,
    8,
    15,
    15,
    15,
    15,
    17,
    17,
    17,
    17,
    0,
    0,
    0,
    0,
    2,
    2,
    3,
    3,
    5,
    5,
    5,
    5,
]
gen_code_mito_invert_deg_lvl: List[int] = [
    4,
    2,
    2,
    2,
    2,
    4,
    2,
    2,
    2,
    6,
    2,
    2,
    4,
    2,
    4,
    8,
    4,
    4,
    2,
    2,
]
gen_code_mito_mam: List[int] = [
    4,
    4,
    9,
    9,
    15,
    15,
    15,
    15,
    19,
    19,
    -1,
    -1,
    1,
    1,
    18,
    18,
    9,
    9,
    9,
    9,
    12,
    12,
    12,
    12,
    6,
    6,
    13,
    13,
    14,
    14,
    14,
    14,
    7,
    7,
    10,
    10,
    16,
    16,
    16,
    16,
    11,
    11,
    8,
    8,
    15,
    15,
    -1,
    -1,
    17,
    17,
    17,
    17,
    0,
    0,
    0,
    0,
    2,
    2,
    3,
    3,
    5,
    5,
    5,
    5,
]

dayhoff6_code: List[int] = [3, 5, 2, 2, 4, 3, 1, 0, 1, 0, 0, 2, 3, 2, 1, 3, 3, 0, 4, 4]
dayhoff4_code: List[int] = [3, -1, 2, 2, 0, 3, 1, 0, 1, 0, 0, 2, 3, 2, 1, 3, 3, 0, 2, 2]

TOOSMALL: float = 1e-30
TOOLARGE: float = 500
TOOLARGENEGATIVE: float = -500


dict_of_codon_standard_code_str_int: Dict[str, int] = {
    codon_i: gen_code_standard[i] for i, codon_i in enumerate(codons)
}

dict_of_codon_standard_code_int_str: Dict[int, str] = {
    gen_code_standard[i]: codon_i
    for i, codon_i in enumerate(codons)
    if gen_code_standard[i] != -1
}


def get_diff_codon_position(codon_1: str, codon_2: str):
    Npos, pos = 0, 0

    if codon_1 == codon_2:
        return -1

    for i in range(0, 3):
        if codon_1[i] != codon_2[i]:
            pos = i
            Npos += 1
        if Npos > 1:
            return -1
    return pos


def is_conservative(amino_acid_1: int, amino_acid_2: int) -> bool:

    if dayhoff6_code[amino_acid_1] == dayhoff6_code[amino_acid_2]:
        return True

    return False


def is_synonymous(codon_1: int, codon_2: int, gen_code: List[int]) -> bool:
    amino_acid_1: int = gen_code[codon_1]
    amino_acid_2: int = gen_code[codon_2]
    if amino_acid_1 != -1 and amino_acid_2 != -1 and amino_acid_1 == amino_acid_2:
        return True
    return False


dict_of_nuleotides_str_int: Dict[str, int] = {nuc: i for i, nuc in enumerate(dna_upper)}
