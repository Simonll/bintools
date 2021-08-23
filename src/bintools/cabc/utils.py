from typing import Optional

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator
from sklearn.preprocessing import FunctionTransformer
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import PowerTransformer
from sklearn.preprocessing import QuantileTransformer
from sklearn.preprocessing import StandardScaler


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
