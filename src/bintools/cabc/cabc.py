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
from sklearn.preprocessing import FunctionTransformer
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import PowerTransformer
from sklearn.preprocessing import QuantileTransformer
from sklearn.preprocessing import StandardScaler

from src.bintools.run.utils import check_path


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
    df_refTab: pd.DataFrame,
    df_realdata: pd.DataFrame,
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
        df_refTab["distance_1"] = distance.cdist(
            model_preprocessing_ss.transform(df_refTab[list_of_ss]),
            model_preprocessing_ss.transform(df_realdata[list_of_ss]),
            lambda u, v: ((u - v) ** 2).sum(),
        )
        pandas2ri.activate()
        df_refTab.sort_values(by="distance_1", inplace=True)
        df_refTab.reset_index(inplace=True)
        df_refTab = df_refTab.loc[:knn, :]
        refTab_param_r = ro.conversion.py2rpy(
            model_preprocessing_params.transform(df_refTab[list_of_params])
        )
        refTab_ss_r = ro.conversion.py2rpy(
            model_preprocessing_ss.transform(df_refTab[list_of_ss])
        )
        true_ss_r = ro.conversion.py2rpy(
            model_preprocessing_ss.transform(df_realdata.loc[:, list_of_ss])
        )

        r.assign("method", method)
        r.assign("transf", transf)
        r.assign("hcorr", hcorr)
        r.assign("kernel", kernel)
        r.assign("simu_space_param", refTab_param_r)
        r.assign("simu_space_ss", refTab_ss_r)
        r.assign("true_ss", true_ss_r)
        r("print(true_ss_r[1:10])")

        if method == "neuralnet":
            try:
                r(
                    "abc_r_knn <- abc(target = true_ss, param = simu_space_param, sumstat = simu_pace_ss, tol = 1, method = method, transf = transf, hcorr = hcorr, kernel = kernel, sizenet = 1)$adj.values"
                )
            except Exception as e:
                print("Something wrong %s" % str(e))
        elif method == "loclinear":
            try:
                r(
                    "abc_r_knn <- abc(target = true_ss, param = simu_space_param, sumstat = simu_pace_ss, tol = 1, method = method, transf = transf, hcorr = hcorr, kernel = kernel)$adj.values"
                )
            except Exception as e:
                print("Something wrong %s" % str(e))

        df_knn = pd.DataFrame(
            model_preprocessing_params.inverse_transform(r.abc_r_knn),
            columns=list_of_params,
        )
    return df_knn


def transform(trans_name: str, df: pd.DataFrame) -> Optional[BaseEstimator]:
    pre_process: Optional[BaseEstimator] = None
    if trans_name == "qt":
        pre_process = QuantileTransformer()
        pre_process.fit(df)
    if trans_name == "boxcox":
        pre_process = PowerTransformer(method="box-cox", standardize=False)
        pre_process.fit(df)
    if trans_name == "minmax":
        pre_process = MinMaxScaler()
        pre_process.fit(df)
    elif trans_name == "standard":
        pre_process = StandardScaler()
        pre_process.fit(df)
    elif trans_name == "log":
        pre_process = FunctionTransformer(func=np.log, inverse_func=np.exp)
        pre_process.fit(df)
    elif trans_name == "log2":
        pre_process = FunctionTransformer(func=np.log2, inverse_func=np.exp2)
        pre_process.fit(df)
    elif trans_name == "log10":
        pre_process = FunctionTransformer(
            func=np.log10, inverse_func=lambda x: np.power(x, 10)
        )
        pre_process.fit(df)
    elif trans_name == "none":
        pre_process = FunctionTransformer(func=lambda x: x, inverse_func=lambda x: x)
        pre_process.fit(df)
    else:
        print("something wrong with transformation name: %s" % trans_name)

    return pre_process


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
    sim_space: str,
    knn: int,
    params_file: str,
    ss_file: str,
    output: str,
    preprocessing_ss: str,
    preprocessing_params: str,
    reg_model: str,
    realdata_ss: str,
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

    df_refTab: pd.DataFrame = pd.read_csv(sim_space, sep="\t", index_col=False)

    print("ref table shape", np.shape(df_refTab))

    col: List[str] = list(df_refTab.columns)
    df_refTab.drop(
        [i for i in col if i not in list_of_params + list_of_ss + ["chainID"]],
        axis=1,
        inplace=True,
    )
    print("shape of reference table", np.shape(df_refTab))
    df_refTab.dropna(inplace=True)
    print("shape of reference table after droping na", np.shape(df_refTab))

    # reading realdata
    df_realdata: pd.DataFrame = pd.read_csv(realdata_ss, sep="\t", index_col=False)
    df_realdata.drop(
        [i for i in list(df_realdata.columns) if (i not in list_of_ss)],
        axis=1,
        inplace=True,
    )
    print("realdata table ", np.shape(df_realdata))
    random_row: int = np.random.randint(0, df_realdata.shape[0] - 1)
    df_realdata = df_realdata.iloc[random_row : random_row + 1, :]
    print("realdata choosing randomly one row ", df_realdata.shape)

    model_preprocessing_ss: Optional[BaseEstimator] = transform(
        trans_name=preprocessing_ss, df=df_refTab[list_of_ss]
    )

    if model_preprocessing_ss is None:
        print("Something wrong with pre-processing of summary statistics")
        raise RuntimeError

    model_preprocessing_params: Optional[BaseEstimator] = transform(
        trans_name=preprocessing_params, df=df_refTab[list_of_params]
    )

    if model_preprocessing_params is None:
        print("Something wrong with pre-processing of params")
        raise RuntimeError

    assert df_realdata.shape[1] == df_refTab[list_of_ss].shape[1]

    df_knn: Optional[pd.DataFrame] = wrapper_abc_r_package(
        df_refTab=df_refTab,
        df_realdata=df_realdata,
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

    df_knn["chainID"] = df_refTab["chainID"]
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
        "--simulation_space",
        type=str,
        required=True,
        description="simulation space defined by model parameters' values along the summary statistics' values",
    )

    parser.add_argument(
        "--knn",
        type=int,
        required=True,
        description="k-nearest neighbors to be kept",
    )

    parser.add_argument(
        "--parameters",
        type=str,
        required=True,
        description="file with list of parameters to be considered",
    )

    parser.add_argument(
        "--summary_statistics",
        type=str,
        required=True,
        description="file with list of summary statistics to be considered",
    )

    parser.add_argument("--output", type=str, required=True, description="output")

    parser.add_argument(
        "--preprocessing_of_summary_statistics",
        type=str,
        required=True,
        description="model for preprocessing summary statistics",
    )

    parser.add_argument(
        "--preprocessing_of_parameters",
        type=str,
        required=True,
        description="model for preprocessing parameters",
    )

    parser.add_argument(
        "--regression_model",
        type=str,
        required=True,
        description="regression model to be use to get closer to the true posterior distribution",
    )

    parser.add_argument(
        "--realdata_summary_statistics_file",
        type=str,
        required=True,
        description="summary statistics values computed from realdata",
    )

    args = parser.parse_args()

    get_closer_to_true_post(
        sim_space=args.simulation_space,
        knn=args.knn,
        params_file=args.parameters,
        ss_file=args.summary_statistics,
        output=args.output,
        preprocessing_ss=args.preprocessing_of_summary_statistics,
        preprocessing_params=args.preprocessing_of_parameters,
        reg_model=args.regression_model,
        realdata_ss=args.realdata_summary_statistics,
    )
