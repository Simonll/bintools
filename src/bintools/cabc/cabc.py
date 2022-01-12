import argparse
import pickle
from typing import List
from typing import Optional

import numpy as np
import pandas as pd
from scipy.spatial import distance
from sklearn.base import BaseEstimator

from bintools.cabc.utils import transform


def generate_cabc_conf(method: str, **kwargs) -> List[str]:
    conf: Optional[List[str]] = None
    if method in ["M0", "M7", "M8"]:
        conf = []
        conf += ["#SUMMARIES" + "\t" + kwargs["ss"]]
        conf += ["#PARAM" + "\t" + kwargs["param"]]
        conf += ["#MAP" + "\t" + kwargs["map"]]
        conf += ["#SAMPLING" + "\t" + kwargs["sampling"]]
        conf += ["#RUN" + "\t" + kwargs["nrun"]]
        conf += ["#NTHREADS" + "\t" + kwargs["nthreads"]]
        conf += ["#TRANS no"]
        conf += ["#OUTPUT" + "\t" + kwargs["output"]]
        conf += [
            "\t".join(
                ["#LOCALPARAM"]
                + [kwargs["localparam"]]
                + ["-d" + "\t" + kwargs["align"]]
                + ["-chain" + "\t" + kwargs["chainname"]]
            )
        ]
    else:
        raise NotImplementedError(
            "ERROR: simulation method %s not implemented" % method
        )
    return conf


def read(input_file: str) -> Optional[List[str]]:
    try:
        with open(input_file, "r") as fh:
            lines: List[str] = fh.readlines()
            lines = [l.strip() for l in lines if not l.startswith("#")]
        return lines
    except:
        print("something wrong: %s" % input_file)
        return None


def recover_r_abc_results(
    knn_file: str, model_preprocessing_params_file: str, output_dir: str, prefix: str
) -> bool:
    try:
        model_preprocessing_params: Optional[BaseEstimator] = None
        with open(model_preprocessing_params_file, "rb") as fh:
            model_preprocessing_params = pickle.load(fh)
    except Exception as e:
        print("something wrong with %s" % model_preprocessing_params_file)
        return False
    try:
        df_knn: pd.DataFrame = pd.read_feather(knn_file)
    except Exception as e:
        print("something wrong with %s" % knn_file)
        return False
    try:
        pd.DataFrame(
            data=model_preprocessing_params.inverse_transform(df_knn),
        ).to_csv(output_dir + prefix + "-df_knn_params.tsv", sep="\t")
    except Exception as e:
        print("something wrong with %s" % output_dir + prefix + "-df_knn_params.tsv")
        return False
    return True


def generate_files_r_abc(
    simu_space_file: str,
    knn: int,
    params: List[str],
    ss: List[str],
    preprocessing_ss: str,
    preprocessing_params: str,
    true_ss_file: str,
    output_dir: str,
    prefix: str,
) -> bool:

    list_of_params: List[str] = params
    list_of_ss: List[str] = ss

    print("Here the list of params and summary statistics")
    print(list_of_params)
    print(list_of_ss)

    df_simu_space: pd.DataFrame = pd.read_csv(
        simu_space_file, sep="\t", index_col=False
    )

    print("simulation space dimensions", np.shape(df_simu_space))

    col: List[str] = list(df_simu_space.columns)
    df_simu_space.drop(
        [i for i in col if i not in list_of_params + list_of_ss + ["chainID"]],
        axis=1,
        inplace=True,
    )
    print("shape of reference table", np.shape(df_simu_space))
    df_simu_space.dropna(inplace=True)
    print("shape of reference table after droping na", np.shape(df_simu_space))

    # reading realdata
    df_true_ss: pd.DataFrame = pd.read_csv(true_ss_file, sep="\t", index_col=False)
    df_true_ss.drop(
        [i for i in list(df_true_ss.columns) if (i not in list_of_ss)],
        axis=1,
        inplace=True,
    )
    print("realdata table ", np.shape(df_true_ss))
    random_row: int = np.random.randint(0, df_true_ss.shape[0] - 1)
    df_true_ss = df_true_ss.iloc[random_row : random_row + 1, :]
    print("realdata choosing randomly one row ", df_true_ss.shape)

    model_preprocessing_ss: Optional[BaseEstimator] = transform(
        trans_name=preprocessing_ss, df=df_simu_space[list_of_ss]
    )

    if model_preprocessing_ss is None:
        print("Something wrong with pre-processing of summary statistics")
        raise RuntimeError

    model_preprocessing_params: Optional[BaseEstimator] = transform(
        trans_name=preprocessing_params, df=df_simu_space[list_of_params]
    )

    if model_preprocessing_params is None:
        print("Something wrong with pre-processing of params")
        raise RuntimeError

    try:
        with open(output_dir + prefix + "-model_preprocessing_ss.pickle", "wb") as fh:
            pickle.dump(model_preprocessing_ss, fh, pickle.HIGHEST_PROTOCOL)
    except Exception as e:
        print(
            "something wrong while dumping %s on %s"
            % (prefix + "model_preprocessing_ss", output_dir)
        )

    try:
        with open(
            output_dir + prefix + "-model_preprocessing_params.pickle", "wb"
        ) as fh:
            pickle.dump(model_preprocessing_params, fh, pickle.HIGHEST_PROTOCOL)
    except Exception as e:
        print(
            "something wrong while dumping %s on %s"
            % (prefix + "model_preprocessing_params", output_dir)
        )

    assert df_true_ss.shape[1] == df_simu_space[list_of_ss].shape[1]

    df_simu_space["metric"] = distance.cdist(
        model_preprocessing_ss.transform(
            df_simu_space[list_of_ss],
        ),
        model_preprocessing_ss.transform(df_true_ss),
        lambda u, v: ((u - v) ** 2).sum(),
    )

    df_simu_space.sort_values(by="metric", inplace=True)
    df_simu_space.reset_index(inplace=True)
    df_simu_space = df_simu_space.iloc[:knn, :]
    assert df_simu_space.shape[0] == knn
    try:
        pd.DataFrame(
            data=model_preprocessing_params.transform(df_simu_space[list_of_params]),
            columns=list_of_params,
        ).to_feather(output_dir + prefix + "-df_simu_space_knn_params.feather")
        pd.DataFrame(
            data=model_preprocessing_ss.transform(df_simu_space[list_of_ss]),
            columns=list_of_ss,
        ).to_feather(output_dir + prefix + "-df_simu_space_knn_ss.feather")
        pd.DataFrame(
            data=model_preprocessing_ss.transform(df_true_ss[list_of_ss]),
            columns=list_of_ss,
        ).reset_index()[list_of_ss].to_feather(
            output_dir + prefix + "-df_true_ss.feather"
        )
        # http://onnx.ai/sklearn-onnx/ for better portability and durability

    except Exception as e:
        print("Something wrong saving files for abc r package %s" % str(e))

        return False

    return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="argument", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--simu_space",
        type=str,
        required=True,
        description="simulation space defined by model parameters' values along the summary statistics' values",
    )

    parser.add_argument(
        "--knn",
        type=int,
        required=True,
        description="k-nearest neighbors to be kept according to a metric",
    )

    parser.add_argument(
        "--params",
        type=str,
        required=True,
        description="list of parameters to be considered",
    )

    parser.add_argument(
        "--ss",
        type=str,
        required=True,
        description="list of summary statistics to be considered",
    )

    parser.add_argument(
        "--output_dir", type=str, required=True, description="output directory"
    )

    parser.add_argument(
        "--trans_fct_ss",
        type=str,
        required=True,
        description="function for preprocessing summary statistics",
    )

    parser.add_argument(
        "--trans_fct_params",
        type=str,
        required=True,
        description="function for preprocessing parameter",
    )

    parser.add_argument(
        "--true_ss_file",
        type=str,
        required=True,
        description="summary statistics values computed from true data",
    )

    parser.add_argument(
        "--prefix",
        type=str,
        required=True,
        description="prefix used to id analyses",
    )

    args = parser.parse_args()

    generate_files_r_abc(
        simu_space_file=args.simulation_space_file,
        knn=args.knn,
        params=args.params,
        ss=args.ss,
        preprocessing_ss=args.trans_fct_ss,
        preprocessing_params=args.trans_fct_params,
        true_ss_file=args.true_ss_file,
        prefix=args.prefix,
        output_dir=args.output_dir,
    )
