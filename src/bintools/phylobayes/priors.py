from typing import Dict
from typing import List

import numpy as np
import scipy as sc
import scipy.stats


def prior_M7M8_mix(N: int, mixture: List[float]) -> Dict[int, float]:
    """
    generates site-heterogeneous omega values from mixture of omega values
    """
    dict_of_omega: Dict[int, float] = {
        i + 1: np.random.choice(a=mixture, replace=False) for i in range(N)
    }
    return dict_of_omega


def prior_M7M8(N: int, p_0: float, p: float, q: float, l: float) -> Dict[int, float]:
    """
    generates site-heterogeneous omega values from beta distribution and exponential
    distribution when p_0 == 1, this corresponds to M7 model from CodeML
    """
    dict_of_omega: Dict[int, float] = {}
    for i in range(N):
        u_float: float = sc.stats.uniform(0, 1).rvs()
        if u_float > p_0:
            # l is the rate (lambda) of the exponential: scale = 1 / l.
            # (expon(l) would set loc=l, i.e. l + Exp(1), which is not intended)
            dict_of_omega[i + 1] = sc.stats.expon(scale=1.0 / l).rvs()
        else:
            dict_of_omega[i + 1] = sc.stats.beta(p, q).rvs()
    return dict_of_omega


def prior_M7M8_fix(N: int, value: float) -> Dict[int, float]:
    """
    generates site specific omega values with one value
    """
    dict_of_omega: Dict[int, float] = {i + 1: value for i in range(N)}
    return dict_of_omega
