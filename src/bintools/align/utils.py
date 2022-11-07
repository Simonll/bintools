from typing import Dict


def sub_sample(seq: str, start: int, stop: int, step: int) -> str:
    return "".join([seq[l] for l in range(start, stop, step)])


def compute_CpG(seq: str) -> Dict[str, float]:
    dict_of_stats: Dict[str, float] = {}

    for i, j in zip(
        sub_sample(seq=seq, start=0, stop=len(seq), step=1),
        sub_sample(seq=seq, start=1, stop=len(seq), step=1),
    ):
        if (i == "C" or i == "c") and (j == "G" or j == "g"):
            if "CpG" in dict_of_stats:
                dict_of_stats["CpG"] += 1
            else:
                dict_of_stats["CpG"] = 1

        if (i == "T" or i == "t") and (j == "A" or j == "a"):
            if "TpA" in dict_of_stats:
                dict_of_stats["TpA"] += 1
            else:
                dict_of_stats["TpA"] = 1

    for i, j in zip(
        sub_sample(seq=seq, start=0, stop=len(seq), step=3),
        sub_sample(seq=seq, start=1, stop=len(seq), step=3),
    ):
        if (i == "C" or i == "c") and (j == "G" or j == "g"):
            if "CpG_12" in dict_of_stats:
                dict_of_stats["CpG_12"] += 1
            else:
                dict_of_stats["CpG_12"] = 1
        if (i == "T" or i == "t") and (j == "A" or j == "a"):
            if "TpA_12" in dict_of_stats:
                dict_of_stats["TpA_12"] += 1
            else:
                dict_of_stats["TpA_12"] = 1

    for i, j in zip(
        sub_sample(seq=seq, start=1, stop=len(seq), step=3),
        sub_sample(seq=seq, start=2, stop=len(seq), step=3),
    ):
        if (i == "C" or i == "c") and (j == "G" or j == "g"):
            if "CpG_23" in dict_of_stats:
                dict_of_stats["CpG_23"] += 1
            else:
                dict_of_stats["CpG_23"] = 1
        if (i == "T" or i == "t") and (j == "A" or j == "a"):
            if "TpA_23" in dict_of_stats:
                dict_of_stats["TpA_23"] += 1
            else:
                dict_of_stats["TpA_23"] = 1
    for i, j in zip(
        sub_sample(seq=seq, start=2, stop=len(seq), step=3),
        sub_sample(seq=seq, start=3, stop=len(seq), step=3),
    ):
        if (i == "C" or i == "c") and (j == "G" or j == "g"):
            if "CpG_31" in dict_of_stats:
                dict_of_stats["CpG_31"] += 1
            else:
                dict_of_stats["CpG_31"] = 1

        if (i == "T" or i == "t") and (j == "A" or j == "a"):
            if "TpA_31" in dict_of_stats:
                dict_of_stats["TpA_31"] += 1
            else:
                dict_of_stats["TpA_31"] = 1

    return {k: v / (len(seq) - 1) for k, v in dict_of_stats.items()}
