import argparse
from typing import List
from typing import Optional

import numpy as np
import pandas as pd
from scipy.spatial import distance
from sklearn.base import BaseEstimator

from bintools.cabc.utils import transform


def read(input_file: str) -> Optional[List[str]]:
    list_: List[str] = []
    try:
        with open(input_file, "r") as fh:
            lines = fh.readlines()
            list_ = [l.strip() for l in lines if not l.startswith("#")]
        return list_
    except:
        print("something wrong: %s" % input_file)
        return None


def prepare_files_for_abc_r_package(
    df_simu_space: pd.DataFrame,
    df_true_ss: pd.DataFrame,
    knn: int,
    list_of_ss: List[str],
    list_of_params: List[str],
    model_preprocessing_ss: BaseEstimator,
    model_preprocessing_params: BaseEstimator,
    output: str,
) -> bool:

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
        ).to_feather(output + "/df_simu_space_knn_params.feather")
        pd.DataFrame(
            data=model_preprocessing_ss.transform(df_simu_space[list_of_ss]),
            columns=list_of_ss,
        ).to_feather(output + "/df_simu_space_knn_ss.feather")
        pd.DataFrame(
            data=model_preprocessing_ss.transform(df_true_ss[list_of_ss]),
            columns=list_of_ss,
        ).reset_index()[list_of_ss].to_feather(output + "/df_true_ss.feather")

    except Exception as e:
        print("Something wrong saving files for abc r package %s" % str(e))

        return False

    return True


def generate_out_dir_cavc_step_0() -> Optional[str]:
    """
    step_0 consists of recovering gene for analysis
    """


def generate_out_dir_cavc_step_1() -> Optional[str]:
    """
    step_1 consists of generating posterior distribution using phylobayes-mpi
    """


def generate_output_dir_cabc_step_2() -> Optional[str]:
    """
    step_2 consists of simulating reference table using likelihoodfreephylogenetics
    """


def generate_output_dir_cabc_step_3(
    geneID: str,
    reg_model: str,
    preprocessing_ss: str,
    preprocessing_params,
    hcorr: str,
    kernel: str,
) -> Optional[str]:
    """
    step_3 consists of preparing files for ABC r package analysis by applying rejection sampling
    """
    output_dir: str = "step_3-"
    if reg_model == "loclinear":
        output_dir += (
            geneID
            + "-LRM"  # linear regression model
            + "-SS_"
            + preprocessing_ss
            + "-PARAMS_"
            + preprocessing_params
            + "-"
            + hcorr
            + "-"
            + kernel
        )
    elif reg_model == "neuralnet":
        output_dir += (
            geneID
            + "-NNRM"  # neural network regression model
            + "-SS_"
            + preprocessing_ss
            + "-PARAMS_"
            + preprocessing_params
            + "-"
            + hcorr
            + "-"
            + kernel
        )
    elif reg_model == "randomforest":
        output_dir += (
            geneID
            + "RFRM"  # random forest regression model
            + "-SS_"
            + preprocessing_ss
            + "-PARAMS_"
            + preprocessing_params
        )

    else:
        print("Something wrong %s" % reg_model)
        return None

    return output_dir


def get_closer_to_true_post(
    simu_space_file: str,
    knn: int,
    params_file: str,
    ss_file: str,
    output: str,
    preprocessing_ss: str,
    preprocessing_params: str,
    reg_model: str,
    true_ss_file: str,
) -> bool:

    list_of_params: Optional[List[str]] = read(params_file)
    list_of_ss: Optional[List[str]] = read(ss_file)

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

    assert df_true_ss.shape[1] == df_simu_space[list_of_ss].shape[1]

    prepare_files_for_abc_r_package(
        df_simu_space=df_simu_space,
        df_true_ss=df_true_ss,
        knn=knn,
        list_of_ss=list_of_ss,
        list_of_params=list_of_params,
        model_preprocessing_ss=model_preprocessing_ss,
        model_preprocessing_params=model_preprocessing_params,
        output="/data/",
    )

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

    get_closer_to_true_post(
        simu_space_file=args.simulation_space_file,
        knn=args.knn,
        params_file=args.params_file,
        ss_file=args.ss_file,
        output=args.output,
        preprocessing_ss=args.trans_fct_ss,
        preprocessing_params=args.trans_fct_params,
        reg_model=args.reg_model,
        true_ss_file=args.true_ss_file,
    )
