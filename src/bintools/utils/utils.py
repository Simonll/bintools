import datetime
import os
import pathlib
from typing import Any
from typing import Dict
from typing import Optional
from typing import Union

import yaml


def check_path(path: Union[str, pathlib.Path]) -> bool:
    path_: Optional[str] = path_to_str(path=path)
    if path_ is not None:
        if not os.path.exists(path_):
            print("something wrong %s" % path)
            return False
        return True
    return False


def path_to_str(path: Union[str, pathlib.Path]) -> Optional[str]:
    if isinstance(path, (pathlib.WindowsPath, pathlib.PosixPath, pathlib.PurePath)):
        return str(path)
    elif type(path) == str:
        return path
    return None


def compute_delta_time(t0: datetime.datetime, msg: str = "") -> str:
    t0_date_time = datetime.datetime.combine(t0.date(), t0.time())
    t1_date_time = datetime.datetime.combine(
        datetime.date.today(), datetime.datetime.now().time()
    )
    delta_time: float = (t1_date_time - t0_date_time).total_seconds()
    unit_of_time: str = "s"
    if delta_time > 120:
        delta_time = delta_time / 60
        unit_of_time = "min"
    elif delta_time < 1e-2:
        delta_time = delta_time * 1000
        unit_of_time = "ms"
    elif delta_time < 1e-4:
        delta_time = delta_time * 1e6
        unit_of_time = "us"

    print("%.2f %s spent in %s" % (delta_time, unit_of_time, msg))
    return "%.2f%s" % (delta_time, unit_of_time)


def get_yaml_config(yaml_file: pathlib.Path) -> Optional[Dict[Any, Any]]:
    with open(yaml_file.__str__(), "r") as stream:
        try:
            dict_of_config: Dict[Any, Any] = yaml.safe_load(stream)
            return dict_of_config
        except yaml.YAMLError as e:
            print(str(e))
            return None
