import os
import pathlib
from typing import Optional
from typing import Union


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
