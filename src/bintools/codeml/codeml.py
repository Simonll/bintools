from typing import Any
from typing import Dict
from typing import List


def generate_codeml_conf(**kwargs) -> List[str]:
    conf: List[str] = []
    for k, v in kwargs.items():
        if w_codeml(kwargs, k) is not None:
            conf += [w_codeml(kwargs, k)]
    else:
        print("something wrong with %s " % (k))

    return conf


def w_codeml(dict_conf: Dict[str, Any], key):
    if key in dict_conf:
        v = dict_conf[key]
        if type(v) in (float, int):
            v = str(v)
        return key + " = " + v + "\n"
    return None
