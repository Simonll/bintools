from typing import Dict
from typing import List

import numpy as np
import scipy as sc

# def get_value_from_list(list_of_values: List[float]) -> float:
#     """
#     samples values from list
#     """
#     u_float: float = sc.stats.uniform(0, 1).rvs()
#     k: int = 1
#     gibbs: float = k / len(list_of_values)
#     while gibbs < u_float:
#         k += 1
#         gibbs = k / len(list_of_values)
#     k -= 1
#     return list_of_values[k]


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
    dict_site: Dict[int, float] = {}
    for i in range(N):
        u_float: float = sc.stats.uniform(0, 1).rvs()
        if u_float > p_0:
            dict_site[i + 1] = sc.stats.expon(l).rvs()
        else:
            dict_site[i + 1] = sc.stats.beta(p, q).rvs()
    return dict_site


def prior_M7M8_fix(N: int, value: float) -> Dict[int, float]:
    """
    generates site specific omega values with one value
    """
    dict_site: Dict[int, float] = {i + 1: value for i in range(N)}
    return dict_site
