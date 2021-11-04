import datetime
import os
import shlex
from typing import Any
from typing import Dict
from typing import Optional

from bintools.run.utils import joint_kwargs

DOCKER_RUN: str = "docker run --user $(id -u):$(id -g) --rm -v "


def generate_phylobayes_cmd(
    method: str,
    mapping: str,
    logger: Optional[str] = None,
    image: str = "ubuntu20.04/phylobayes",
    **kwargs
) -> Optional[str]:
    cmd: Optional[str] = None
    if method == "pb":
        chainname: str = ""
        if "-chainname" in kwargs:
            chainname = kwargs["-chainname"]
            del kwargs["-chainname"]
        cmd = " ".join(
            [
                DOCKER_RUN,
                mapping,
                image,
                "pb",
                joint_kwargs(**kwargs),
                chainname,
                logger if logger is not None else "",
            ]
        )
    if method == "ppred":
        if "-chainname" in kwargs:
            chainname = kwargs["-chainname"]
            del kwargs["-chainname"]
        cmd = " ".join(
            [
                DOCKER_RUN,
                mapping,
                image,
                "ppred",
                joint_kwargs(**kwargs),
                chainname,
                logger if logger is not None else "",
            ]
        )
    return cmd


def generate_abc_cmd(
    method: str, mapping: str, image: str = "r-base3.6.3/abc", **kwargs
) -> Optional[str]:
    cmd: Optional[str] = None
    if method == "abc":
        cmd = " ".join(
            [
                DOCKER_RUN,
                mapping,
                image,
                joint_kwargs(**kwargs),
            ]
        )
    return cmd


def generate_pb_mpi_cmd(
    method: str,
    mapping: str,
    logger: Optional[str] = None,
    image: str = "ubuntu20.04/pbmpi",
    **kwargs
) -> Optional[str]:
    cmd: Optional[str] = None
    if method == "pb_mpi":
        np: Optional[int] = os.cpu_count()
        if np is None:
            np = 1

        if "-np" in kwargs:
            np = kwargs["-np"]
            del kwargs["-np"]

        chainname: str = ""
        if "-chainname" in kwargs:
            chainname = kwargs["-chainname"]
            del kwargs["-chainname"]

        cmd = " ".join(
            [
                DOCKER_RUN,
                mapping,
                image,
                "mpirun --allow-run-as-root  -np",
                str(np),
                "pb_mpi",
                joint_kwargs(**kwargs),
                chainname,
                logger if logger is not None else "",
            ]
        )

    elif method == "readpb_mpi":
        if "-np" in kwargs:
            np = kwargs["-np"]
            del kwargs["-np"]

        if "-chainname" in kwargs:
            chainname = kwargs["-chainname"]
            del kwargs["-chainname"]

        cmd = " ".join(
            [
                DOCKER_RUN,
                mapping,
                image,
                "mpirun --allow-run-as-root  -np",
                str(np),
                "readpb_mpi",
                joint_kwargs(**kwargs),
                chainname,
                logger if logger is not None else "",
            ]
        )

    else:
        raise NotImplementedError(
            "ERROR: readpb_mpi method %s not implemented yet" % method
        )

    return cmd


def generate_coevol_cmd(
    method: str, mapping: str, image: str = "ubuntu20.04/coevol", **kwargs
) -> Optional[str]:
    cmd: Optional[str] = None
    if method == "coevol":
        cmd = " ".join([DOCKER_RUN, mapping, image, "coevol", joint_kwargs(**kwargs)])

    elif method == "readcoevol":
        cmd = " ".join(
            [DOCKER_RUN, mapping, image, "readcoevol", joint_kwargs(**kwargs)]
        )

    else:
        raise NotImplementedError(
            "ERROR: coevol method %s not implemented yet" % method
        )

    return cmd


def generate_datasets_cmd(
    method: str, mapping: str, image: str = "/opt/datasets", **kwargs
) -> Optional[str]:
    cmd: Optional[str] = None
    if method == "download genome":
        cmd = " ".join(
            [
                image,
                "download genome",
                joint_kwargs(**kwargs),
                "--exclude-protein",
                "--exclude-rna",
            ]
        )

    if method == "summary":
        cmd = " ".join(
            [
                image,
                "summary genome",
                joint_kwargs(**kwargs),
            ]
        )

    else:
        raise NotImplementedError(
            "ERROR: datasets method %s not implemented yet" % method
        )

    return cmd


def generate_blat_cmd(
    method,
    mapping,
    image: str = "quay.io/biocontainers/ucsc-blat:377--h0b8a92a_3",
    **kwargs
):
    # information about all command line details -->http://genome.ucsc.edu/goldenPath/help/blatSpec.html
    cmd: Optional[str] = None
    if method == "blat":
        cmd = " ".join([DOCKER_RUN, mapping, image, "blat", joint_kwargs(**kwargs)])
    else:
        raise NotImplementedError("ERROR: blast method %s not implemented" % method)

    return cmd


def generate_blast_cmd(
    method: str, mapping: str, image: str = "biocontainers/blast:2.2.31", **kwargs
):
    # information about all command line details --> https://www.ncbi.nlm.nih.gov/books/NBK279684/
    cmd: Optional[str] = None
    if method == "makeblastdb":
        cmd = " ".join(
            [DOCKER_RUN, mapping, image, "makeblastdb", joint_kwargs(**kwargs)]
        )

    elif method == "blastn":
        cmd = " ".join([DOCKER_RUN, mapping, image, "blastn", joint_kwargs(**kwargs)])

    elif method == "tblastx":
        cmd = " ".join([DOCKER_RUN, mapping, image, "blastx", joint_kwargs(**kwargs)])

    else:
        raise NotImplementedError("ERROR: blast method %s not implemented" % method)

    return cmd


def generate_fasttree_cmd(
    method,
    aln_fname,
    tree_fname,
    log_fname,
    image: str = "/opt/anaconda3/envs/vert/bin/fasttree",
):
    cmd: Optional[str] = None

    if method == "fasttree":
        cmd = "fasttree -nosupport -lg %s 1> %s 2> %s" % (
            shlex.quote(aln_fname),
            shlex.quote(tree_fname),
            shlex.quote(log_fname),
        )

    return cmd


def generate_alignment_mafft_cmd(
    method: str,
    nthreads: int,
    existing_aln_fname: str,
    seqs_to_align_fname: str,
    aln_fname: str,
    log_fname: str,
):
    # taken from https://github.com/nextstrain/augur/blob/867c5e368e5f621528039fe5d417c85e9e66c0f0/augur/align.py
    cmd: Optional[str] = None
    if method == "mafft":
        if existing_aln_fname:
            cmd = "mafft --add %s --keeplength --reorder --anysymbol --nomemsave --adjustdirection --thread %d %s 1> %s 2> %s" % (
                shlex.quote(seqs_to_align_fname),
                nthreads,
                shlex.quote(existing_aln_fname),
                shlex.quote(aln_fname),
                shlex.quote(log_fname),
            )
        else:
            cmd = "mafft --reorder --anysymbol --nomemsave --adjustdirection --thread %d %s 1> %s 2> %s" % (
                nthreads,
                shlex.quote(seqs_to_align_fname),
                shlex.quote(aln_fname),
                shlex.quote(log_fname),
            )
    #        print("\nusing mafft to align via:\n\t" + cmd +
    #            " \n\n\tKatoh et al, Nucleic Acid Research, vol 30, issue 14"
    #            "\n\thttps://doi.org/10.1093%2Fnar%2Fgkf436\n")

    return cmd


def generate_alignment_prank_cmd(
    method: str, codon: str, seqs_to_align_fname: str, aln_fname: str, log_fname: str
):
    cmd: Optional[str] = None

    if codon:
        cmd = "prank -codon -d=%s -o=%s 1>%s" % (
            shlex.quote(seqs_to_align_fname),
            shlex.quote(aln_fname),
            shlex.quote(log_fname),
        )

    else:
        cmd = "prank -codon -d=%s -o=%s 1>%s" % (
            shlex.quote(seqs_to_align_fname),
            shlex.quote(aln_fname),
            shlex.quote(log_fname),
        )
    return cmd


def generate_simu_cmd(
    method: str,
    mapping: str,
    config_filename: str,
    logger: Optional[str] = None,
    image: str = "ubuntu20.04/lfp",
    **kwargs
) -> Optional[str]:
    cmd: Optional[str] = None
    if method == "M7":
        cmd = " ".join(
            [
                DOCKER_RUN,
                mapping,
                image,
                "codemlM7M8",
                joint_kwargs(**kwargs),
                config_filename,
                logger if logger is not None else "",
            ]
        )
    else:
        raise NotImplementedError(
            "ERROR: simulation method %s not implemented" % method
        )
    return cmd


def generate_codeml_cmd(
    method: str, dict_conf: Dict[str, Any], output_dir: str, codeml_path: str
):

    cmd = ""
    keys = [
        # aaRatefile, fix_rho, rho
        "seqfile",
        "treefile",
        "outfile",
        "NSsites",
        "ncatG",
        "noisy",
        "verbose",
        "runmode",
        "seqtype",
        "CodonFreq",
        "clock",
        "aaDist",
        "model",
        "icode",
        "Mgene",
        "fix_kappa",
        "kappa",
        "fix_omega",
        "omega",
        "fix_alpha",
        "alpha",
        "Malpha",
        "getSE",
        "RateAncestor",
        "fix_blength",
        "method",
        "Small_Diff",
        "cleandata",
        "hkyREV",
    ]

    dict_conf.update(
        {
            "noisy": 0,
            "verbose": 1,
            "runmode": 0,
            "seqtype": 1,
            "clock": 0,
            "aaDist": 0,
            "model": 0,
            "icode": 0,
            "Mgene": 0,
            "fix_kappa": 0,
            "kappa": 2,
            "fix_alpha": 1,
            "alpha": 0,
            "Malpha": 0,
            "getSE": 0,
            "RateAncestor": 0,
            "fix_blength": 0,
            "method": 0,
            "Small_Diff": 0.45e-6,
            "cleandata": 0,
        }
    )

    filename_conf: str = dict_conf["mark"] + ".conf"
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        print("Something wrong when creating dir %s" % output_dir)
        raise RuntimeError

    if not os.path.exists(dict_conf["seqfile"]):
        print("path to alignment doesnt exist %s" % dict_conf["seqfile"])
        raise RuntimeError

    if not os.path.exists(dict_conf["treefile"]):
        print("path to tree doesnt exist %s" % dict_conf["treefile"])
        raise RuntimeError

    if method == "codeml":
        with open(output_dir + filename_conf, "w") as hl:
            for k in keys:
                if k in dict_conf:
                    hl.write(w_codeml(dict_conf=dict_conf, key=k))
                else:
                    print(k)
        cmd = " ".join([codeml_path, filename_conf])

    if cmd is None:
        print("Something wrong with cmd")
        raise RuntimeError

    return cmd


def w_codeml(dict_conf: Dict[str, Any], key):
    if key in dict_conf:
        v = dict_conf[key]
        if type(v) in (float, int):
            v = str(v)
        return key + " = " + v + "\n"
    return None


def get_run_stamp(output_dir: str) -> str:
    return (
        output_dir + "run-" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M") + "/"
    )
