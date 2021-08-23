import argparse
from typing import List
from typing import Optional

import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects import r
from rpy2.robjects.packages import importr

importr("abc")

from scipy.spatial import distance
from sklearn.base import BaseEstimator

from bintools.cabc.utils import transform
from utils.utils import check_path


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


def wrapper_abc_r_package(
    df_simu_space: pd.DataFrame,
    df_true_ss: pd.DataFrame,
    method: str,
    knn: int,
    transf: str,
    hcorr: str,
    kernel: str,
    list_of_ss: List[str],
    list_of_params: List[str],
    model_preprocessing_ss: BaseEstimator,
    model_preprocessing_params: BaseEstimator,
) -> Optional[pd.DataFrame]:

    df_knn: Optional[pd.DataFrame] = None

    if method == "neuralnet" or method == "loclinear":
        df_simu_space["metric"] = distance.cdist(
            model_preprocessing_ss.transform(df_simu_space[list_of_ss]),
            model_preprocessing_ss.transform(df_true_ss[list_of_ss]),
            lambda u, v: ((u - v) ** 2).sum(),
        )
        pandas2ri.activate()
        df_simu_space.sort_values(by="metric", inplace=True)
        df_simu_space.reset_index(inplace=True)
        df_simu_space = df_simu_space.loc[:knn, :]
        r.assign("method", method)
        r.assign("transf", transf)
        r.assign("hcorr", hcorr)
        r.assign("kernel", kernel)
        r.assign(
            "simu_space_param",
            ro.conversion.py2rpy(
                model_preprocessing_params.transform(df_simu_space[list_of_params])
            ),
        )
        r.assign(
            "simu_space_ss",
            ro.conversion.py2rpy(
                model_preprocessing_ss.transform(df_simu_space[list_of_ss])
            ),
        )
        r.assign(
            "true_ss",
            ro.conversion.py2rpy(
                model_preprocessing_ss.transform(df_true_ss.loc[:, list_of_ss])
            ),
        )
        r("print(true_ss_r[1:10])")

        if method == "neuralnet":
            try:
                r(
                    "df_knn <- abc(target = true_ss, param = simu_space_param, sumstat = simu_pace_ss, tol = 1, method = method, transf = transf, hcorr = hcorr, kernel = kernel, sizenet = 1)$adj.values"
                )
            except Exception as e:
                print("Something wrong %s" % str(e))
        elif method == "loclinear":
            try:
                r(
                    "df_knn <- abc(target = true_ss, param = simu_space_param, sumstat = simu_pace_ss, tol = 1, method = method, transf = transf, hcorr = hcorr, kernel = kernel)$adj.values"
                )
            except Exception as e:
                print("Something wrong %s" % str(e))

        df_knn = pd.DataFrame(
            model_preprocessing_params.inverse_transform(r.df_knn),
            columns=list_of_params,
        )
    return df_knn


def generate_output_filename(
    method: str,
    output: str,
    preprocessing_ss: str,
    preprocessing_params,
    hcorr: str,
    kernel: str,
) -> Optional[str]:

    if not check_path(path=output):
        raise RuntimeError

    output_: Optional[str] = None
    if method == "loclinear":
        output_ = (
            output
            + "-"
            + "LRM"  # linear regression model
            + "-SS_"
            + preprocessing_ss
            + "-PARAMS_"
            + preprocessing_params
            + "-"
            + hcorr
            + "-"
            + kernel
        )
    elif method == "neuralnet":
        output_ = (
            output
            + "-"
            + "NNRM"  # neural network regression model
            + "-SS_"
            + preprocessing_ss
            + "-PARAMS_"
            + preprocessing_params
            + "-"
            + hcorr
            + "-"
            + kernel
        )
    elif method == "randomforest":
        output_ = (
            output
            + "-"
            + "RFRM"  # random forest regression model
            + "-SS_"
            + preprocessing_ss
            + "-PARAMS_"
            + preprocessing_params
        )

    else:
        print("Something wrong %s" % method)
        return None

    return output_


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

    method: str = reg_model
    transf: str = "none"
    hcorr: str = "TRUE"
    kernel: str = "epanechnikov"

    output_: Optional[str] = generate_output_filename(
        method=method,
        output=output,
        preprocessing_ss=preprocessing_ss,
        preprocessing_params=preprocessing_params,
        hcorr=hcorr,
        kernel=kernel,
    )

    if output_ is None:
        print("Something wrong with %s" % output_)
        raise RuntimeError

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
    df_realdata = df_true_ss.iloc[random_row : random_row + 1, :]
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

    df_knn: Optional[pd.DataFrame] = wrapper_abc_r_package(
        df_simu_space=df_simu_space,
        df_true_ss=df_true_ss,
        method=method,
        knn=knn,
        transf=transf,
        hcorr=hcorr,
        kernel=kernel,
        list_of_ss=list_of_ss,
        list_of_params=list_of_params,
        model_preprocessing_ss=model_preprocessing_ss,
        model_preprocessing_params=model_preprocessing_params,
    )

    if df_knn is None:
        print("Something wrong with regression model results")
        raise RuntimeError

    df_knn["chainID"] = df_simu_space["chainID"]
    print(output_ + "." + "knn")
    df_knn[list_of_params + ["chainID"]].to_csv(
        output_ + "." + "knn", sep="\t", index=False, float_format="%.8f"
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
