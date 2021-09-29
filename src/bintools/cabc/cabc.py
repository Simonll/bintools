import argparse
from typing import List
from typing import Optional

import numpy as np
import pandas as pd
from scipy.spatial import distance
from sklearn.base import BaseEstimator

from bintools.cabc.utils import transform


def generate_cabc_conf(method: str, **kwargs) -> List[str]:
    conf: Optional[List[str]] = None
    if method == "M7":
        conf = []
        conf += ["#SUMMARIES" + "\t" + kwargs["ss"]]
        conf += ["#PARAM" + "\t" + kwargs["param"]]
        conf += ["#MAP" + "\t" + kwargs["map"]]
        conf += ["#SAMPLING" + "\t" + kwargs["sampling"]]
        conf += ["#RUN" + kwargs["nrun"]]
        conf += ["#NTHREADS" + kwargs["nthreads"]]
        conf += ["#TRANS no"]
        conf += ["#OUTPUT" + kwargs["output"]]
        conf += [
            "\t".join(
                ["#LOCALPARAM"]
                + [kwargs["localparam"]]
                + ["-d " + kwargs["align"]]
                + ["-chain " + kwargs["chainname"]]
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


def generate_files_r_abc(
    simu_space_file: str,
    knn: int,
    params_file: str,
    ss_file: str,
    preprocessing_ss: str,
    preprocessing_params: str,
    reg_model: str,
    true_ss_file: str,
    input_dir: str,
    output_dir: str,
) -> bool:

    list_of_params: Optional[List[str]] = read(input_file=input_dir + params_file)
    list_of_ss: Optional[List[str]] = read(input_file=input_dir + ss_file)

    if list_of_params is None or list_of_ss is None:
        print(
            "Something wrong with params %s or summary statistics %ss"
            % (params_file, ss_file)
        )
        raise RuntimeError

    print("Here the list of params and summary statistics")
    print(list_of_params)
    print(list_of_ss)

    df_simu_space: pd.DataFrame = pd.read_csv(
        input_dir + simu_space_file, sep="\t", index_col=False
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
    df_true_ss: pd.DataFrame = pd.read_csv(
        input_dir + true_ss_file, sep="\t", index_col=False
    )
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
        ).to_feather(output_dir + "df_simu_space_knn_params.feather")
        pd.DataFrame(
            data=model_preprocessing_ss.transform(df_simu_space[list_of_ss]),
            columns=list_of_ss,
        ).to_feather(output_dir + "df_simu_space_knn_ss.feather")
        pd.DataFrame(
            data=model_preprocessing_ss.transform(df_true_ss[list_of_ss]),
            columns=list_of_ss,
        ).reset_index()[list_of_ss].to_feather(output_dir + "df_true_ss.feather")

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
        "--params_file",
        type=str,
        required=True,
        description="file with list of parameters to be considered",
    )

    parser.add_argument(
        "--ss_file",
        type=str,
        required=True,
        description="file with list of summary statistics to be considered",
    )

    parser.add_argument(
        "--input_dir", type=str, required=True, description="input directory"
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
        "--reg_model",
        type=str,
        required=True,
        description="regression model to be use to get closer to the true posterior distribution",
    )

    parser.add_argument(
        "--true_ss_file",
        type=str,
        required=True,
        description="summary statistics values computed from true data",
    )

    args = parser.parse_args()

    generate_files_r_abc(
        simu_space_file=args.simulation_space_file,
        knn=args.knn,
        params_file=args.params_file,
        ss_file=args.ss_file,
        preprocessing_ss=args.trans_fct_ss,
        preprocessing_params=args.trans_fct_params,
        reg_model=args.reg_model,
        true_ss_file=args.true_ss_file,
        input_dir=args.input_dir,
        output_dir=args.output_dir,
    )
