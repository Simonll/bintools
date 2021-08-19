import datetime
import glob
import os
import shlex

shquote = shlex.quote

def print_kwargs(**kwargs):
    print(" ".join([k + " " + v for k, v in kwargs.items()]))


def generate_cabc_cmd(method, **kwargs):
    cmd = ""
    if method == "knn":
        cmd = " ".join([param + " " + value for param, value in kwargs.items()])
    return cmd


def generate_coevol_cmd(method, mapping, **kwargs):
    cmd = ""
    image = " coevol_ubuntu20"
    if method == "coevol":
        cmd = (
            "docker run --user $(id -u):$(id -g) --rm -v "
            + mapping
            + image
            + " coevol "
            + " ".join([param + " " + value for param, value in kwargs.items()])
        )

    elif method == "readcoevol":
        cmd = (
            "docker run --rm -v "
            + mapping
            + image
            + " readcoevol "
            + " ".join([param + " " + value for param, value in kwargs.items()])
        )

    else:
        raise NotImplementedError("ERROR: coevol method %s not implemented" % method)

    return cmd


def generate_datasets_cmd(method, mapping, **kwargs):
    cmd = ""
    image = " /opt/datasets "
    if method == "download":
        cmd = (
            image
            + " download genome "
            + " ".join([param + " " + value for param, value in kwargs.items()])
            + " --exclude-protein --exclude-rna "
        )

    return cmd


def generate_blast_cmd(method, mapping, **kwargs):
    # information about all command line details --> https://www.ncbi.nlm.nih.gov/books/NBK279684/
    cmd = ""
    image = " biocontainers/blast:2.2.31 "
    if method == "makeblastdb":
        cmd = (
            "docker run --rm -v "
            + mapping
            + image
            + " makeblastdb "
            + " ".join([param + " " + value for param, value in kwargs.items()])
        )

    elif method == "blastn":
        cmd = (
            "docker run --rm -v "
            + mapping
            + image
            + " blastn "
            + " ".join([param + " " + value for param, value in kwargs.items()])
        )

    elif method == "tblastx":
        cmd = (
            "docker run --rm -v "
            + mapping
            + image
            + " tblastx "
            + " ".join([param + " " + value for param, value in kwargs.items()])
        )

    else:
        raise NotImplementedError("ERROR: blast method %s not implemented" % method)

    return cmd


def generate_fasttree_cmd(method, aln_fname, tree_fname, log_fname):
    cmd = ""
    image = "/opt/anaconda3/envs/vert/bin/fasttree"

    if method == "fasttree":
        cmd = "fasttree -nosupport -lg %s 1> %s 2> %s" % (
            shquote(aln_fname),
            shquote(tree_fname),
            shquote(log_fname),
        )

    return cmd


def generate_alignment_mafft_cmd(
    method, nthreads, existing_aln_fname, seqs_to_align_fname, aln_fname, log_fname
):
    # taken from https://github.com/nextstrain/augur/blob/867c5e368e5f621528039fe5d417c85e9e66c0f0/augur/align.py
    cmd = ""
    if method == "mafft":
        if existing_aln_fname:
            cmd = "mafft --add %s --keeplength --reorder --anysymbol --nomemsave --adjustdirection --thread %d %s 1> %s 2> %s" % (
                shquote(seqs_to_align_fname),
                nthreads,
                shquote(existing_aln_fname),
                shquote(aln_fname),
                shquote(log_fname),
            )
        else:
            cmd = "mafft --reorder --anysymbol --nomemsave --adjustdirection --thread %d %s 1> %s 2> %s" % (
                nthreads,
                shquote(seqs_to_align_fname),
                shquote(aln_fname),
                shquote(log_fname),
            )
    #        print("\nusing mafft to align via:\n\t" + cmd +
    #            " \n\n\tKatoh et al, Nucleic Acid Research, vol 30, issue 14"
    #            "\n\thttps://doi.org/10.1093%2Fnar%2Fgkf436\n")

    return cmd


def generate_alignment_prank_cmd(
    method, codon, seqs_to_align_fname, aln_fname, log_fname
):
    if codon:
        cmd = "prank -codon -d=%s -o=%s 1>%s" % (
            shquote(seqs_to_align_fname),
            shquote(aln_fname),
            shquote(log_fname),
        )

    else:
        cmd = "prank -codon -d=%s -o=%s 1>%s" % (
            shquote(seqs_to_align_fname),
            shquote(aln_fname),
            shquote(log_fname),
        )
    return cmd


def generate_simu_cmd(method, dict_conf, output_dir, log_fname):
    lfp_path = "/opt/LikelihoodFreePhylogenetics/data/LFP"
    codemlM7M8_path = "/opt/LikelihoodFreePhylogenetics/data/codemlM7M8"

    filename = (
        dict_conf["geneID"].split(".")[0]
        + "-"
        + dict_conf["model"]
        + "-"
        + dict_conf["mark"]
    )
    cmd = ""
    os.makedirs(output_dir, exist_ok=True)
    if not os.path.exists(dict_conf["aln_path"]):
        print("path to alignment doesnt exist %s" % dict_conf["aln_path"])

    if not os.path.exists(dict_conf["chain_path"] + ".chain"):
        print("path to pb_mpi doesnt exist %s" % dict_conf["chain_path"])

    if method == "simulation":
        with open(output_dir + filename + ".conf", "w") as hl:
            hl.write("\t".join(["#SUMMARIES"] + dict_conf["ss"]))
            hl.write("\n")
            hl.write("\t".join(["#PARAM"] + dict_conf["param"]))
            hl.write("\n")
            hl.write("\t".join(["#MAP"] + dict_conf["map"]))
            hl.write("\n")
            hl.write("\t".join(["#SAMPLING" + " " + dict_conf["sampling"]]))
            hl.write("\n")
            # hl.write("\t".join(["#CHAIN"+" "+dict_conf["chain_path"]]))
            # hl.write("\n")
            hl.write("\t".join(["#RUN" + " " + dict_conf["nrun"]]))
            hl.write("\n")
            hl.write("\t".join(["#NTHREADS" + " " + str(dict_conf["nthreads"])]))
            hl.write("\n")
            hl.write("#TRANS no")
            hl.write("\n")
            hl.write("\t".join(["#OUTPUT"] + [output_dir + filename]))
            hl.write("\n")
            hl.write(
                "\t".join(
                    ["#LOCALPARAM"]
                    + dict_conf["localparam"]
                    + ["-d " + dict_conf["aln_path"]]
                    + ["-chain " + dict_conf["chain_path"]]
                )
            )
        cmd = " ".join(
            [lfp_path] + [dict_conf["args"]] + [output_dir + filename + ".conf"]
        )
    elif method == "M7" or method == "M8" or method == "M8a":
        with open(output_dir + filename + ".conf", "w") as hl:
            hl.write("\t".join(["#SUMMARIES"] + dict_conf["ss"]))
            hl.write("\n")
            hl.write("\t".join(["#PARAM"] + dict_conf["param"]))
            hl.write("\n")
            hl.write("\t".join(["#MAP"] + dict_conf["map"]))
            hl.write("\n")
            hl.write("\t".join(["#SAMPLING" + " " + dict_conf["sampling"]]))
            hl.write("\n")
            hl.write("\t".join(["#RUN" + " " + dict_conf["nrun"]]))
            hl.write("\n")
            hl.write("\t".join(["#NTHREADS" + " " + str(1)]))
            hl.write("\n")
            hl.write("#TRANS no")
            hl.write("\n")
            hl.write("\t".join(["#OUTPUT"] + [output_dir + filename]))
            hl.write("\n")
            hl.write(
                "\t".join(
                    ["#LOCALPARAM"]
                    + dict_conf["localparam"]
                    + ["-d " + dict_conf["aln_path"]]
                    + ["-chain " + dict_conf["chain_path"]]
                )
            )
        cmd = " ".join(
            [codemlM7M8_path] + [dict_conf["args"]] + [output_dir + filename + ".conf"]
        )

    cmd = cmd + " 2>" + output_dir + log_fname

    return cmd


def generate_codeml_cmd(method, dict_conf, output_dir, codeml_path):

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

    filename_conf = dict_conf["mark"] + ".conf"
    os.makedirs(output_dir, exist_ok=True)

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

    return cmd


def w_codeml(dict_conf, key):
    if key in dict_conf:

        v = dict_conf[key]
        if type(v) in (float, int):
            v = str(v)
        return key + " = " + v + "\n"
    return None


def get_run_stamp(output_dir):
    return (
        output_dir + "run-" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M") + "/"
    )


def generate_experiment_MG_freeF1x4_GTR_beta_l(list_of_cores):
    """
    simulated under MG - GTR mix of omega values {(0.5,0.5, "beta_u"), (10,10, "beta_norm"),(1,1, "beta_unif"),(0.9,0.1, "beta_neutral"),(0.1,0.9, "beta_l")}
    analyzed M0-MG-HKY
    """
    codeml_path = "/home/sll/paml4.9j/bin/codeml"
    sel = ["M7", "M8", "M8a"]
    mut = ["HKY", "GTR"]
    input_dir = "/home/sll/forward_sim/pseudo_data/"
    output_dir = get_run_stamp("/home/sll/forward_sim/")
    os.makedirs(output_dir, exist_ok=True)
    tree = input_dir + "Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    files = glob.glob(input_dir + "*GTR*.fasta")

    for f in files:
        if "beta_l" in f and ("CpG-1-" in f or "CpG-2-" in f or "CpG-8-" in f):
            for s in sel:
                for m in mut:
                    prefix = f.split("/")[-1].split(".fasta")[0]
                    mark = s + "-" + m
                    if not os.path.exists(input_dir):
                        print("input_dir missing" % input_dir)
                        raise RuntimeError

                    if not os.path.exists(output_dir):
                        print("output_dir missing" % output_dir)
                        raise RuntimeError

                    dict_conf = {
                        "prefix": prefix,
                        "mark": mark,
                        "seqfile": f,
                        "CodonFreq": 4,
                        "treefile": tree,
                        "outfile": output_dir + prefix + "-" + mark + "/" + mark,
                        "NSsites": 7
                        if s in ["M7"]
                        else 8
                        if s in ["M8", "M8a"]
                        else "",  # 7:M7, M8
                        "ncatG": 10,
                        "hkyREV": 1 if m in ["GTR"] else 0 if m in ["HKY"] else "",
                        "fix_omega": 0
                        if s in ["M7", "M8"]
                        else 1
                        if s in ["M8a"]
                        else "",
                        "omega": 1,
                    }

                    cmd = generate_codeml_cmd(
                        method="codeml",
                        dict_conf=dict_conf,
                        output_dir=output_dir + prefix + "-" + mark + "/",
                        codeml_path=codeml_path,
                    )
                    append_write = ""
                    cmdfile = output_dir + mark + ".conf"
                    if os.path.exists(cmdfile):
                        append_write = "a"
                    else:
                        append_write = "w"

                    with open(cmdfile, append_write) as fh:
                        fh.write("cd " + output_dir + prefix + "-" + mark + "/")
                        fh.write("\n")
                        fh.write(
                            "qsub "
                            + " ".join(["-q " + i for i in list_of_cores])
                            + " -cwd launch_codeml.sh"
                        )
                        fh.write("\n")

                    with open(
                        output_dir + prefix + "-" + mark + "/" + "launch_codeml.sh", "w"
                    ) as hl:
                        hl.write(cmd)


def generate_experiment_MG_freeF1x4_HKY_beta_neutral(list_of_cores):
    """
    simulated under MG - HKY mix of omega values {(0.5,0.5, "beta_u"), (10,10, "beta_norm"),(1,1, "beta_unif"),(0.9,0.1, "beta_neutral")}
    analyzed M0-MG-HKY
    """
    codeml_path = "/home/sll/paml4.9j/bin/codeml"
    sel = ["M7", "M8", "M8a"]
    mut = ["HKY", "GTR"]
    input_dir = "/home/sll/forward_sim/pseudo_data/"
    output_dir = get_run_stamp("/home/sll/forward_sim/")
    os.makedirs(output_dir, exist_ok=True)
    tree = input_dir + "Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    files = glob.glob(input_dir + "*HKY*.fasta")

    for f in files:
        if "beta_neutral" in f and ("CpG-1-" in f or "CpG-2-" in f or "CpG-8-" in f):
            for s in sel:
                for m in mut:
                    prefix = f.split("/")[-1].split(".fasta")[0]
                    mark = s + "-" + m
                    if not os.path.exists(input_dir):
                        print("input_dir missing" % input_dir)
                        raise RuntimeError

                    if not os.path.exists(output_dir):
                        print("output_dir missing" % output_dir)
                        raise RuntimeError

                    dict_conf = {
                        "prefix": prefix,
                        "mark": mark,
                        "seqfile": f,
                        "CodonFreq": 4,
                        "treefile": tree,
                        "outfile": output_dir + prefix + "-" + mark + "/" + mark,
                        "NSsites": 7
                        if s in ["M7"]
                        else 8
                        if s in ["M8", "M8a"]
                        else "",  # 7:M7, M8
                        "ncatG": 10,
                        "hkyREV": 1 if m in ["GTR"] else 0 if m in ["HKY"] else "",
                        "fix_omega": 0
                        if s in ["M7", "M8"]
                        else 1
                        if s in ["M8a"]
                        else "",
                        "omega": 1,
                    }

                    cmd = generate_codeml_cmd(
                        method="codeml",
                        dict_conf=dict_conf,
                        output_dir=output_dir + prefix + "-" + mark + "/",
                        codeml_path=codeml_path,
                    )
                    append_write = ""
                    cmdfile = output_dir + mark + ".conf"
                    if os.path.exists(cmdfile):
                        append_write = "a"
                    else:
                        append_write = "w"

                    with open(cmdfile, append_write) as fh:
                        fh.write("cd " + output_dir + prefix + "-" + mark + "/")
                        fh.write("\n")
                        fh.write(
                            "qsub "
                            + " ".join(["-q " + i for i in list_of_cores])
                            + " -cwd launch_codeml.sh"
                        )
                        fh.write("\n")

                    with open(
                        output_dir + prefix + "-" + mark + "/" + "launch_codeml.sh", "w"
                    ) as hl:
                        hl.write(cmd)


def generate_experiment_MG_F1x4_HKY():
    """
    simulated under MG - HKY mix of omega values {(0.5,0.5, "beta_u"), (10,10, "beta_norm"),(1,1, "beta_unif"),(0.9,0.1, "beta_neutral")}
    analyzed M0-MG-HKY
    """
    codeml_path = "/home/sll/paml4.9j/bin/codeml"
    sel = ["M7", "M8", "M8a"]
    mut = ["HKY", "GTR"]
    input_dir = "/home/sll/forward_sim/pseudo_data/"
    output_dir = get_run_stamp("/home/sll/forward_sim/")
    os.makedirs(output_dir, exist_ok=True)
    tree = input_dir + "Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    files = glob.glob(input_dir + "*HKY*.fasta")

    for f in files:
        if "beta_norm" in f and ("CpG-1-" in f or "CpG-2-" in f or "CpG-8-" in f):
            for s in sel:
                for m in mut:
                    prefix = f.split("/")[-1].split(".fasta")[0]
                    mark = s + "-" + m
                    if not os.path.exists(input_dir):
                        print("input_dir missing" % input_dir)
                        raise RuntimeError

                    if not os.path.exists(output_dir):
                        print("output_dir missing" % output_dir)
                        raise RuntimeError

                    dict_conf = {
                        "prefix": prefix,
                        "mark": mark,
                        "seqfile": f,
                        "CodonFreq": 1,
                        "treefile": tree,
                        "outfile": output_dir + prefix + "-" + mark + "/" + mark,
                        "NSsites": 7
                        if s in ["M7"]
                        else 8
                        if s in ["M8", "M8a"]
                        else "",  # 7:M7, M8
                        "ncatG": 10,
                        "hkyREV": 1 if m in ["GTR"] else 0 if m in ["HKY"] else "",
                        "fix_omega": 0
                        if s in ["M7", "M8"]
                        else 1
                        if s in ["M8a"]
                        else "",
                        "omega": 1,
                    }

                    cmd = generate_codeml_cmd(
                        method="codeml",
                        dict_conf=dict_conf,
                        output_dir=output_dir + prefix + "-" + mark + "/",
                        codeml_path=codeml_path,
                    )
                    append_write = ""
                    cmdfile = output_dir + mark + ".conf"
                    if os.path.exists(cmdfile):
                        append_write = "a"
                    else:
                        append_write = "w"

                    with open(cmdfile, append_write) as fh:
                        fh.write("cd " + output_dir + prefix + "-" + mark + "/")
                        fh.write("\n")
                        fh.write("qsub -cwd launch_codeml.sh")
                        fh.write("\n")

                    with open(
                        output_dir + prefix + "-" + mark + "/" + "launch_codeml.sh", "w"
                    ) as hl:
                        hl.write(cmd)


def generate_experiment_MG_F1x4_HKY_with_fixed():
    """
    simulated under MG - HKY mix of omega values {(0.2, "fixed_0.2"), (0.5, "fixed_0.5"),(0.8, "fixed_0.8")}
    analyzed M0-MG-HKY
    """
    codeml_path = "/home/sll/paml4.9j/bin/codeml"
    sel = ["M0"]
    mut = ["HKY", "GTR"]
    input_dir = "/home/sll/forward_sim/pseudo_data/M0/"
    output_dir = get_run_stamp("/home/sll/forward_sim/")
    os.makedirs(output_dir, exist_ok=True)
    tree = input_dir + "Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    files = glob.glob(input_dir + "*.fasta")

    for f in files:
        if (
            "fixed" in f
            and "-HKY-" in f
            and ("CpG-1-" in f or "CpG-2-" in f or "CpG-8-" in f)
        ):
            for s in sel:
                for m in mut:
                    prefix = f.split("/")[-1].split(".fasta")[0]
                    mark = s + "-" + m
                    if not os.path.exists(input_dir):
                        print("input_dir missing" % input_dir)
                        raise RuntimeError

                    if not os.path.exists(output_dir):
                        print("output_dir missing" % output_dir)
                        raise RuntimeError

                    dict_conf = {
                        "prefix": prefix,
                        "mark": mark,
                        "seqfile": f,
                        "CodonFreq": 1,
                        "treefile": tree,
                        "outfile": output_dir + prefix + "-" + mark + "/" + mark,
                        "NSsites": 7
                        if s in ["M7"]
                        else 8
                        if s in ["M8", "M8a"]
                        else 0
                        if s in ["M0"]
                        else "",  # 7:M7, M8
                        "ncatG": 10,
                        "hkyREV": 1 if m in ["GTR"] else 0 if m in ["HKY"] else "",
                        "fix_omega": 0
                        if s in ["M7", "M8", "M0"]
                        else 1
                        if s in ["M8a"]
                        else "",
                        "omega": 1,
                    }

                    cmd = generate_codeml_cmd(
                        method="codeml",
                        dict_conf=dict_conf,
                        output_dir=output_dir + prefix + "-" + mark + "/",
                        codeml_path=codeml_path,
                    )
                    append_write = ""
                    cmdfile = output_dir + mark + ".conf"
                    if os.path.exists(cmdfile):
                        append_write = "a"
                    else:
                        append_write = "w"

                    with open(cmdfile, append_write) as fh:
                        fh.write("cd " + output_dir + prefix + "-" + mark + "/")
                        fh.write("\n")
                        fh.write("qsub -cwd launch_codeml.sh")
                        fh.write("\n")

                    with open(
                        output_dir + prefix + "-" + mark + "/" + "launch_codeml.sh", "w"
                    ) as hl:
                        hl.write(cmd)


def generate_experiment_MG_freeF1x4_GTR_with_fixed(list_of_cores):
    """
    simulated under MG - GTR choice of omega values {(0.2, "fixed_0.2"), (0.5, "fixed_0.5"),(0.8, "fixed_0.8")}
    analyzed M0-freeMG-GTR
    """
    codeml_path = "/home/sll/paml4.9j/bin/codeml"
    sel = ["M0"]
    mut = ["HKY", "GTR"]
    input_dir = "/home/sll/forward_sim/pseudo_data/M0/"
    output_dir = get_run_stamp("/home/sll/forward_sim/")
    os.makedirs(output_dir, exist_ok=True)
    tree = input_dir + "Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    files = glob.glob(input_dir + "*.fasta")

    for f in files:
        if (
            "fixed" in f
            and "-GTR-" in f
            and ("CpG-1-" in f or "CpG-2-" in f or "CpG-8-" in f)
        ):
            for s in sel:
                for m in mut:
                    prefix = f.split("/")[-1].split(".fasta")[0]
                    mark = s + "-" + m
                    if not os.path.exists(input_dir):
                        print("input_dir missing" % input_dir)
                        raise RuntimeError

                    if not os.path.exists(output_dir):
                        print("output_dir missing" % output_dir)
                        raise RuntimeError

                    dict_conf = {
                        "prefix": prefix,
                        "mark": mark,
                        "seqfile": f,
                        "CodonFreq": 4,
                        "treefile": tree,
                        "outfile": output_dir + prefix + "-" + mark + "/" + mark,
                        "NSsites": 7
                        if s in ["M7"]
                        else 8
                        if s in ["M8", "M8a"]
                        else 0
                        if s in ["M0"]
                        else "",  # 7:M7, M8
                        "ncatG": 10,
                        "hkyREV": 1 if m in ["GTR"] else 0 if m in ["HKY"] else "",
                        "fix_omega": 0
                        if s in ["M7", "M8", "M0"]
                        else 1
                        if s in ["M8a"]
                        else "",
                        "omega": 1,
                    }

                    cmd = generate_codeml_cmd(
                        method="codeml",
                        dict_conf=dict_conf,
                        output_dir=output_dir + prefix + "-" + mark + "/",
                        codeml_path=codeml_path,
                    )
                    append_write = ""
                    cmdfile = output_dir + mark + ".conf"
                    if os.path.exists(cmdfile):
                        append_write = "a"
                    else:
                        append_write = "w"

                    with open(cmdfile, append_write) as fh:
                        fh.write("cd " + output_dir + prefix + "-" + mark + "/")
                        fh.write("\n")
                        fh.write(
                            "qsub "
                            + " ".join(["-q " + i for i in list_of_cores])
                            + " -cwd launch_codeml.sh"
                        )
                        fh.write("\n")

                    with open(
                        output_dir + prefix + "-" + mark + "/" + "launch_codeml.sh", "w"
                    ) as hl:
                        hl.write(cmd)


def generate_experiment_MG_freeF1x4_HKY_with_fixed(list_of_cores):
    """
    simulated under MG - HKY mix of omega values {(0.2, "fixed_0.2"), (0.5, "fixed_0.5"),(0.8, "fixed_0.8")}
    analyzed M0-freeMG-HKY
    """
    codeml_path = "/home/sll/paml4.9j/bin/codeml"
    sel = ["M0"]
    mut = ["HKY", "GTR"]
    input_dir = "/home/sll/forward_sim/pseudo_data/M0_fix_0/"
    output_dir = get_run_stamp("/home/sll/forward_sim/")
    os.makedirs(output_dir, exist_ok=True)
    tree = input_dir + "Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    files = glob.glob(input_dir + "*.fasta")

    for f in files:
        if (
            "fixed" in f
            and "-HKY-" in f
            and ("CpG-1-" in f or "CpG-2-" in f or "CpG-8-" in f)
        ):
            for s in sel:
                for m in mut:
                    prefix = f.split("/")[-1].split(".fasta")[0]
                    mark = s + "-" + m
                    if not os.path.exists(input_dir):
                        print("input_dir missing" % input_dir)
                        raise RuntimeError

                    if not os.path.exists(output_dir):
                        print("output_dir missing" % output_dir)
                        raise RuntimeError

                    dict_conf = {
                        "prefix": prefix,
                        "mark": mark,
                        "seqfile": f,
                        "CodonFreq": 4,
                        "treefile": tree,
                        "outfile": output_dir + prefix + "-" + mark + "/" + mark,
                        "NSsites": 7
                        if s in ["M7"]
                        else 8
                        if s in ["M8", "M8a"]
                        else 0
                        if s in ["M0"]
                        else "",  # 7:M7, M8
                        "ncatG": 10,
                        "hkyREV": 1 if m in ["GTR"] else 0 if m in ["HKY"] else "",
                        "fix_omega": 0
                        if s in ["M7", "M8", "M0"]
                        else 1
                        if s in ["M8a"]
                        else "",
                        "omega": 1,
                    }

                    cmd = generate_codeml_cmd(
                        method="codeml",
                        dict_conf=dict_conf,
                        output_dir=output_dir + prefix + "-" + mark + "/",
                        codeml_path=codeml_path,
                    )
                    append_write = ""
                    cmdfile = output_dir + mark + ".conf"
                    if os.path.exists(cmdfile):
                        append_write = "a"
                    else:
                        append_write = "w"

                    with open(cmdfile, append_write) as fh:
                        fh.write("cd " + output_dir + prefix + "-" + mark + "/")
                        fh.write("\n")
                        fh.write(
                            "qsub "
                            + " ".join(["-q " + i for i in list_of_cores])
                            + " -cwd launch_codeml.sh"
                        )
                        fh.write("\n")

                    with open(
                        output_dir + prefix + "-" + mark + "/" + "launch_codeml.sh", "w"
                    ) as hl:
                        hl.write(cmd)


def generate_experiment_MG_freeF1x4_HKY_with_fixed_CpG_4(list_of_cores):
    """
    simulated under MG - HKY mix of omega values {(0.2, "fixed_0.2"), (0.5, "fixed_0.5"),(0.8, "fixed_0.8")}
    analyzed M0-freeMG-HKY
    """
    codeml_path = "/home/sll/paml4.9j/bin/codeml"
    sel = ["M0"]
    mut = ["HKY", "GTR"]
    input_dir = "/home/sll/forward_sim/pseudo_data/M0_fix_0/"
    output_dir = get_run_stamp("/home/sll/forward_sim/")
    os.makedirs(output_dir, exist_ok=True)
    tree = input_dir + "Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    files = glob.glob(input_dir + "*.fasta")

    for f in files:
        if "fixed" in f and "-HKY-" in f and ("CpG-4-" in f):
            for s in sel:
                for m in mut:
                    prefix = f.split("/")[-1].split(".fasta")[0]
                    mark = s + "-" + m
                    if not os.path.exists(input_dir):
                        print("input_dir missing" % input_dir)
                        raise RuntimeError

                    if not os.path.exists(output_dir):
                        print("output_dir missing" % output_dir)
                        raise RuntimeError

                    dict_conf = {
                        "prefix": prefix,
                        "mark": mark,
                        "seqfile": f,
                        "CodonFreq": 4,
                        "treefile": tree,
                        "outfile": output_dir + prefix + "-" + mark + "/" + mark,
                        "NSsites": 7
                        if s in ["M7"]
                        else 8
                        if s in ["M8", "M8a"]
                        else 0
                        if s in ["M0"]
                        else "",  # 7:M7, M8
                        "ncatG": 10,
                        "hkyREV": 1 if m in ["GTR"] else 0 if m in ["HKY"] else "",
                        "fix_omega": 0
                        if s in ["M7", "M8", "M0"]
                        else 1
                        if s in ["M8a"]
                        else "",
                        "omega": 1,
                    }

                    cmd = generate_codeml_cmd(
                        method="codeml",
                        dict_conf=dict_conf,
                        output_dir=output_dir + prefix + "-" + mark + "/",
                        codeml_path=codeml_path,
                    )
                    append_write = ""
                    cmdfile = output_dir + mark + ".conf"
                    if os.path.exists(cmdfile):
                        append_write = "a"
                    else:
                        append_write = "w"

                    with open(cmdfile, append_write) as fh:
                        fh.write("cd " + output_dir + prefix + "-" + mark + "/")
                        fh.write("\n")
                        fh.write(
                            "qsub "
                            + " ".join(["-q " + i for i in list_of_cores])
                            + " -cwd launch_codeml.sh"
                        )
                        fh.write("\n")

                    with open(
                        output_dir + prefix + "-" + mark + "/" + "launch_codeml.sh", "w"
                    ) as hl:
                        hl.write(cmd)


def generate_experiment_GY_F3x4_HKY_with_fixed(list_of_cores):
    """
    simulated under GY F3 X 4 - HKY mix of omega values {(0.2, "fixed_0.2"), (0.5, "fixed_0.5"),(0.8, "fixed_0.8")}
    analyzed M0-MG-HKY
    """
    codeml_path = "/home/sll/paml4.9j/bin/codeml"
    sel = ["M0"]
    mut = ["HKY", "GTR"]
    input_dir = "/home/sll/forward_sim/pseudo_data/M0/"
    output_dir = get_run_stamp("/home/sll/forward_sim/")
    os.makedirs(output_dir, exist_ok=True)
    tree = input_dir + "Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    files = glob.glob(input_dir + "*.fasta")

    for f in files:
        if (
            "fixed" in f
            and "-HKY-" in f
            and ("CpG-1-" in f or "CpG-2-" in f or "CpG-8-" in f)
        ):
            for s in sel:
                for m in mut:
                    prefix = f.split("/")[-1].split(".fasta")[0]
                    mark = s + "-" + m
                    if not os.path.exists(input_dir):
                        print("input_dir missing" % input_dir)
                        raise RuntimeError

                    if not os.path.exists(output_dir):
                        print("output_dir missing" % output_dir)
                        raise RuntimeError

                    dict_conf = {
                        "prefix": prefix,
                        "mark": mark,
                        "seqfile": f,
                        "CodonFreq": 2,
                        "treefile": tree,
                        "outfile": output_dir + prefix + "-" + mark + "/" + mark,
                        "NSsites": 7
                        if s in ["M7"]
                        else 8
                        if s in ["M8", "M8a"]
                        else 0
                        if s in ["M0"]
                        else "",  # 7:M7, M8
                        "ncatG": 10,
                        "hkyREV": 1 if m in ["GTR"] else 0 if m in ["HKY"] else "",
                        "fix_omega": 0
                        if s in ["M7", "M8", "M0"]
                        else 1
                        if s in ["M8a"]
                        else "",
                        "omega": 1,
                    }

                    cmd = generate_codeml_cmd(
                        method="codeml",
                        dict_conf=dict_conf,
                        output_dir=output_dir + prefix + "-" + mark + "/",
                        codeml_path=codeml_path,
                    )
                    append_write = ""
                    cmdfile = output_dir + mark + ".conf"
                    if os.path.exists(cmdfile):
                        append_write = "a"
                    else:
                        append_write = "w"

                    with open(cmdfile, append_write) as fh:
                        fh.write("cd " + output_dir + prefix + "-" + mark + "/")
                        fh.write("\n")
                        fh.write(
                            "qsub "
                            + " ".join(["-q " + i for i in list_of_cores])
                            + " -cwd launch_codeml.sh"
                        )
                        fh.write("\n")

                    with open(
                        output_dir + prefix + "-" + mark + "/" + "launch_codeml.sh", "w"
                    ) as hl:
                        hl.write(cmd)


def generate_experiment_GY_F61_HKY_with_fixed(list_of_cores):
    """
    simulated under GY F61 - HKY mix of omega values {(0.2, "fixed_0.2"), (0.5, "fixed_0.5"),(0.8, "fixed_0.8")}
    analyzed M0-MG-HKY
    """
    codeml_path = "/home/sll/paml4.9j/bin/codeml"
    sel = ["M0"]
    mut = ["HKY", "GTR"]
    input_dir = "/home/sll/forward_sim/pseudo_data/M0/"
    output_dir = get_run_stamp("/home/sll/forward_sim/")
    os.makedirs(output_dir, exist_ok=True)
    tree = input_dir + "Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    files = glob.glob(input_dir + "*.fasta")

    for f in files:
        if (
            "fixed" in f
            and "-HKY-" in f
            and ("CpG-1-" in f or "CpG-2-" in f or "CpG-8-" in f)
        ):
            for s in sel:
                for m in mut:
                    prefix = f.split("/")[-1].split(".fasta")[0]
                    mark = s + "-" + m
                    if not os.path.exists(input_dir):
                        print("input_dir missing" % input_dir)
                        raise RuntimeError

                    if not os.path.exists(output_dir):
                        print("output_dir missing" % output_dir)
                        raise RuntimeError

                    dict_conf = {
                        "prefix": prefix,
                        "mark": mark,
                        "seqfile": f,
                        "CodonFreq": 3,
                        "treefile": tree,
                        "outfile": output_dir + prefix + "-" + mark + "/" + mark,
                        "NSsites": 7
                        if s in ["M7"]
                        else 8
                        if s in ["M8", "M8a"]
                        else 0
                        if s in ["M0"]
                        else "",  # 7:M7, M8
                        "ncatG": 10,
                        "hkyREV": 1 if m in ["GTR"] else 0 if m in ["HKY"] else "",
                        "fix_omega": 0
                        if s in ["M7", "M8", "M0"]
                        else 1
                        if s in ["M8a"]
                        else "",
                        "omega": 1,
                    }

                    cmd = generate_codeml_cmd(
                        method="codeml",
                        dict_conf=dict_conf,
                        output_dir=output_dir + prefix + "-" + mark + "/",
                        codeml_path=codeml_path,
                    )
                    append_write = ""
                    cmdfile = output_dir + mark + ".conf"
                    if os.path.exists(cmdfile):
                        append_write = "a"
                    else:
                        append_write = "w"

                    with open(cmdfile, append_write) as fh:
                        fh.write("cd " + output_dir + prefix + "-" + mark + "/")
                        fh.write("\n")
                        fh.write(
                            "qsub "
                            + " ".join(["-q " + i for i in list_of_cores])
                            + " -cwd launch_codeml.sh"
                        )
                        fh.write("\n")

                    with open(
                        output_dir + prefix + "-" + mark + "/" + "launch_codeml.sh", "w"
                    ) as hl:
                        hl.write(cmd)


def get_list_of_cores():
    list_of_cores = ["all.q@compute-0-" + str(i) + ".local" for i in range(16)]
    return list_of_cores


def get_sub_list_of_cores(sub_list_of_cores):
    list_of_cores = get_list_of_cores()
    list_of_cores_ = []
    for i in sub_list_of_cores:
        if "all.q@compute-0-" + str(i) + ".local" in list_of_cores:
            list_of_cores_ += ["all.q@compute-0-" + str(i) + ".local"]

    return list_of_cores_


def generate_experiment_MG_freeF1x4_HKY_with_mix(list_of_cores):
    """
    simulated under MG F1x4 - HKY mix of omega values {([0.1, 0.2, 0.3], "mix_1"), ([0.4, 0.5, 0.6], "mix_2"),([0.7, 0.8, 0.9], "mix_3"), ([0.2, 0.5, 0.7], "mix_4"), ([0.5, 0.7, 0.9], "mix_5")}
    analyzed M0-MG-HKY
    """
    codeml_path = "/home/sll/paml4.9j/bin/codeml"
    sel = ["M0"]
    mut = ["HKY", "GTR"]
    input_dir = "/home/sll/forward_sim/pseudo_data/M0_mix_0/"
    output_dir = get_run_stamp("/home/sll/forward_sim/")
    os.makedirs(output_dir, exist_ok=True)
    tree = input_dir + "Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    files = glob.glob(input_dir + "*.fasta")

    for f in files:
        if (
            "mix" in f
            and "-HKY-" in f
            and ("CpG-1-" in f or "CpG-2-" in f or "CpG-8-" in f)
        ):
            for s in sel:
                for m in mut:
                    prefix = f.split("/")[-1].split(".fasta")[0]
                    mark = s + "-" + m
                    if not os.path.exists(input_dir):
                        print("input_dir missing" % input_dir)
                        raise RuntimeError

                    if not os.path.exists(output_dir):
                        print("output_dir missing" % output_dir)
                        raise RuntimeError

                    dict_conf = {
                        "prefix": prefix,
                        "mark": mark,
                        "seqfile": f,
                        "CodonFreq": 4,
                        "treefile": tree,
                        "outfile": output_dir + prefix + "-" + mark + "/" + mark,
                        "NSsites": 7
                        if s in ["M7"]
                        else 8
                        if s in ["M8", "M8a"]
                        else 0
                        if s in ["M0"]
                        else "",  # 7:M7, M8
                        "ncatG": 10,
                        "hkyREV": 1 if m in ["GTR"] else 0 if m in ["HKY"] else "",
                        "fix_omega": 0
                        if s in ["M7", "M8", "M0"]
                        else 1
                        if s in ["M8a"]
                        else "",
                        "omega": 1,
                    }

                    cmd = generate_codeml_cmd(
                        method="codeml",
                        dict_conf=dict_conf,
                        output_dir=output_dir + prefix + "-" + mark + "/",
                        codeml_path=codeml_path,
                    )
                    append_write = ""
                    cmdfile = output_dir + mark + ".conf"
                    if os.path.exists(cmdfile):
                        append_write = "a"
                    else:
                        append_write = "w"

                    with open(cmdfile, append_write) as fh:
                        fh.write("cd " + output_dir + prefix + "-" + mark + "/")
                        fh.write("\n")
                        fh.write(
                            "qsub "
                            + " ".join(["-q " + i for i in list_of_cores])
                            + " -cwd launch_codeml.sh"
                        )
                        fh.write("\n")

                    with open(
                        output_dir + prefix + "-" + mark + "/" + "launch_codeml.sh", "w"
                    ) as hl:
                        hl.write(cmd)


def generate_experiment_MG_freeF1x4_HKY_with_mix_M7M8M8a(list_of_cores):
    """
    simulated under MG F1x4 - HKY mix of omega values {([0.1, 0.2, 0.3], "mix_1"), ([0.4, 0.5, 0.6], "mix_2"),([0.7, 0.8, 0.9], "mix_3"), ([0.2, 0.5, 0.7], "mix_4"), ([0.5, 0.7, 0.9], "mix_5")}
    analyzed M0-MG-HKY
    """
    codeml_path = "/home/sll/paml4.9j/bin/codeml"
    sel = ["M7", "M8", "M8a"]
    mut = ["HKY", "GTR"]
    input_dir = "/home/sll/forward_sim/pseudo_data/M0_mix_0/"
    output_dir = get_run_stamp("/home/sll/forward_sim/")
    os.makedirs(output_dir, exist_ok=True)
    tree = input_dir + "Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    files = glob.glob(input_dir + "*.fasta")

    for f in files:
        if (
            "mix" in f
            and "-HKY-" in f
            and ("CpG-1-" in f or "CpG-2-" in f or "CpG-8-" in f)
        ):
            for s in sel:
                for m in mut:
                    prefix = f.split("/")[-1].split(".fasta")[0]
                    mark = s + "-" + m
                    if not os.path.exists(input_dir):
                        print("input_dir missing" % input_dir)
                        raise RuntimeError

                    if not os.path.exists(output_dir):
                        print("output_dir missing" % output_dir)
                        raise RuntimeError

                    dict_conf = {
                        "prefix": prefix,
                        "mark": mark,
                        "seqfile": f,
                        "CodonFreq": 4,
                        "treefile": tree,
                        "outfile": output_dir + prefix + "-" + mark + "/" + mark,
                        "NSsites": 7
                        if s in ["M7"]
                        else 8
                        if s in ["M8", "M8a"]
                        else 0
                        if s in ["M0"]
                        else "",  # 7:M7, M8
                        "ncatG": 10,
                        "hkyREV": 1 if m in ["GTR"] else 0 if m in ["HKY"] else "",
                        "fix_omega": 0
                        if s in ["M7", "M8", "M0"]
                        else 1
                        if s in ["M8a"]
                        else "",
                        "omega": 1,
                    }

                    cmd = generate_codeml_cmd(
                        method="codeml",
                        dict_conf=dict_conf,
                        output_dir=output_dir + prefix + "-" + mark + "/",
                        codeml_path=codeml_path,
                    )
                    append_write = ""
                    cmdfile = output_dir + mark + ".conf"
                    if os.path.exists(cmdfile):
                        append_write = "a"
                    else:
                        append_write = "w"

                    with open(cmdfile, append_write) as fh:
                        fh.write("cd " + output_dir + prefix + "-" + mark + "/")
                        fh.write("\n")
                        fh.write(
                            "qsub "
                            + " ".join(["-q " + i for i in list_of_cores])
                            + " -cwd launch_codeml.sh"
                        )
                        fh.write("\n")

                    with open(
                        output_dir + prefix + "-" + mark + "/" + "launch_codeml.sh", "w"
                    ) as hl:
                        hl.write(cmd)


def generate_experiment_MG_freeF1x4_HKY_with_mix_M7M8M8a_CpG_4(list_of_cores):
    """
    simulated under MG F1x4 - HKY mix of omega values {([0.1, 0.2, 0.3], "mix_1"), ([0.4, 0.5, 0.6], "mix_2"),([0.7, 0.8, 0.9], "mix_3"), ([0.2, 0.5, 0.7], "mix_4"), ([0.5, 0.7, 0.9], "mix_5")}
    analyzed M0-MG-HKY
    """
    codeml_path = "/home/sll/paml4.9j/bin/codeml"
    sel = ["M7", "M8", "M8a"]
    mut = ["HKY", "GTR"]
    input_dir = "/home/sll/forward_sim/pseudo_data/M0_mix_0/"
    output_dir = get_run_stamp("/home/sll/forward_sim/")
    os.makedirs(output_dir, exist_ok=True)
    tree = input_dir + "Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    files = glob.glob(input_dir + "*.fasta")

    for f in files:
        if "mix" in f and "-HKY-" in f and ("CpG-4-" in f):
            for s in sel:
                for m in mut:
                    prefix = f.split("/")[-1].split(".fasta")[0]
                    mark = s + "-" + m
                    if not os.path.exists(input_dir):
                        print("input_dir missing" % input_dir)
                        raise RuntimeError

                    if not os.path.exists(output_dir):
                        print("output_dir missing" % output_dir)
                        raise RuntimeError

                    dict_conf = {
                        "prefix": prefix,
                        "mark": mark,
                        "seqfile": f,
                        "CodonFreq": 4,
                        "treefile": tree,
                        "outfile": output_dir + prefix + "-" + mark + "/" + mark,
                        "NSsites": 7
                        if s in ["M7"]
                        else 8
                        if s in ["M8", "M8a"]
                        else 0
                        if s in ["M0"]
                        else "",  # 7:M7, M8
                        "ncatG": 10,
                        "hkyREV": 1 if m in ["GTR"] else 0 if m in ["HKY"] else "",
                        "fix_omega": 0
                        if s in ["M7", "M8", "M0"]
                        else 1
                        if s in ["M8a"]
                        else "",
                        "omega": 1,
                    }

                    cmd = generate_codeml_cmd(
                        method="codeml",
                        dict_conf=dict_conf,
                        output_dir=output_dir + prefix + "-" + mark + "/",
                        codeml_path=codeml_path,
                    )
                    append_write = ""
                    cmdfile = output_dir + mark + ".conf"
                    if os.path.exists(cmdfile):
                        append_write = "a"
                    else:
                        append_write = "w"

                    with open(cmdfile, append_write) as fh:
                        fh.write("cd " + output_dir + prefix + "-" + mark + "/")
                        fh.write("\n")
                        fh.write(
                            "qsub "
                            + " ".join(["-q " + i for i in list_of_cores])
                            + " -cwd launch_codeml.sh"
                        )
                        fh.write("\n")

                    with open(
                        output_dir + prefix + "-" + mark + "/" + "launch_codeml.sh", "w"
                    ) as hl:
                        hl.write(cmd)


def generate_experiment_MG_F1x4_HKY_with_mix_root_1(list_of_cores):
    """
    simulated under MG F1x4 - HKY mix of omega values {([0.1, 0.2, 0.3], "mix_1"), ([0.4, 0.5, 0.6], "mix_2"),([0.7, 0.8, 0.9], "mix_3")} and root length = 1 (not stationary)
    analyzed M0-MG-HKY
    """
    codeml_path = "/home/sll/paml4.9j/bin/codeml"
    sel = ["M0"]
    mut = ["HKY", "GTR"]
    input_dir = "/home/sll/forward_sim/pseudo_data/M0_mix_root_1/"
    output_dir = get_run_stamp("/home/sll/forward_sim/")
    os.makedirs(output_dir, exist_ok=True)
    tree = input_dir + "Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    files = glob.glob(input_dir + "*.fasta")

    for f in files:
        if (
            "mix" in f
            and "-HKY-" in f
            and ("CpG-1-" in f or "CpG-2-" in f or "CpG-8-" in f)
        ):
            for s in sel:
                for m in mut:
                    prefix = f.split("/")[-1].split(".fasta")[0]
                    mark = s + "-" + m
                    if not os.path.exists(input_dir):
                        print("input_dir missing" % input_dir)
                        raise RuntimeError

                    if not os.path.exists(output_dir):
                        print("output_dir missing" % output_dir)
                        raise RuntimeError

                    dict_conf = {
                        "prefix": prefix,
                        "mark": mark,
                        "seqfile": f,
                        "CodonFreq": 1,
                        "treefile": tree,
                        "outfile": output_dir + prefix + "-" + mark + "/" + mark,
                        "NSsites": 7
                        if s in ["M7"]
                        else 8
                        if s in ["M8", "M8a"]
                        else 0
                        if s in ["M0"]
                        else "",  # 7:M7, M8
                        "ncatG": 10,
                        "hkyREV": 1 if m in ["GTR"] else 0 if m in ["HKY"] else "",
                        "fix_omega": 0
                        if s in ["M7", "M8", "M0"]
                        else 1
                        if s in ["M8a"]
                        else "",
                        "omega": 1,
                    }

                    cmd = generate_codeml_cmd(
                        method="codeml",
                        dict_conf=dict_conf,
                        output_dir=output_dir + prefix + "-" + mark + "/",
                        codeml_path=codeml_path,
                    )
                    append_write = ""
                    cmdfile = output_dir + mark + ".conf"
                    if os.path.exists(cmdfile):
                        append_write = "a"
                    else:
                        append_write = "w"

                    with open(cmdfile, append_write) as fh:
                        fh.write("cd " + output_dir + prefix + "-" + mark + "/")
                        fh.write("\n")
                        fh.write(
                            "qsub "
                            + " ".join(["-q " + i for i in list_of_cores])
                            + " -cwd launch_codeml.sh"
                        )
                        fh.write("\n")

                    with open(
                        output_dir + prefix + "-" + mark + "/" + "launch_codeml.sh", "w"
                    ) as hl:
                        hl.write(cmd)


def generate_experiment_MG_F1x4_GTR_with_mix(list_of_cores):
    """
    simulated under MG F1x4 - GTR mix of omega values {([0.1, 0.2, 0.3], "mix_1"), ([0.4, 0.5, 0.6], "mix_2"),([0.7, 0.8, 0.9], "mix_3")}
    analyzed M0-MG-HKY
    """
    codeml_path = "/home/sll/paml4.9j/bin/codeml"
    sel = ["M0"]
    mut = ["HKY", "GTR"]
    input_dir = "/home/sll/forward_sim/pseudo_data/M0_mix/"
    output_dir = get_run_stamp("/home/sll/forward_sim/")
    os.makedirs(output_dir, exist_ok=True)
    tree = input_dir + "Mammalia-39sp-CAT_unrooted_withoutsupport.tre"
    files = glob.glob(input_dir + "*.fasta")

    for f in files:
        if (
            "mix" in f
            and "-GTR-" in f
            and ("CpG-1-" in f or "CpG-2-" in f or "CpG-8-" in f)
        ):
            for s in sel:
                for m in mut:
                    prefix = f.split("/")[-1].split(".fasta")[0]
                    mark = s + "-" + m
                    if not os.path.exists(input_dir):
                        print("input_dir missing" % input_dir)
                        raise RuntimeError

                    if not os.path.exists(output_dir):
                        print("output_dir missing" % output_dir)
                        raise RuntimeError

                    dict_conf = {
                        "prefix": prefix,
                        "mark": mark,
                        "seqfile": f,
                        "CodonFreq": 1,
                        "treefile": tree,
                        "outfile": output_dir + prefix + "-" + mark + "/" + mark,
                        "NSsites": 7
                        if s in ["M7"]
                        else 8
                        if s in ["M8", "M8a"]
                        else 0
                        if s in ["M0"]
                        else "",  # 7:M7, M8
                        "ncatG": 10,
                        "hkyREV": 1 if m in ["GTR"] else 0 if m in ["HKY"] else "",
                        "fix_omega": 0
                        if s in ["M7", "M8", "M0"]
                        else 1
                        if s in ["M8a"]
                        else "",
                        "omega": 1,
                    }

                    cmd = generate_codeml_cmd(
                        method="codeml",
                        dict_conf=dict_conf,
                        output_dir=output_dir + prefix + "-" + mark + "/",
                        codeml_path=codeml_path,
                    )
                    append_write = ""
                    cmdfile = output_dir + mark + ".conf"
                    if os.path.exists(cmdfile):
                        append_write = "a"
                    else:
                        append_write = "w"

                    with open(cmdfile, append_write) as fh:
                        fh.write("cd " + output_dir + prefix + "-" + mark + "/")
                        fh.write("\n")
                        fh.write(
                            "qsub "
                            + " ".join(["-q " + i for i in list_of_cores])
                            + " -cwd launch_codeml.sh"
                        )
                        fh.write("\n")

                    with open(
                        output_dir + prefix + "-" + mark + "/" + "launch_codeml.sh", "w"
                    ) as hl:
                        hl.write(cmd)
