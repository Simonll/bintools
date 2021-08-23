import os
from pathlib import Path
from typing import Optional
from typing import Union


def check_path(path: Union[str, Path]) -> bool:
    path_: Optional[str] = path_to_str(path=path)
    if path_ is not None:
        if not os.path.exists(path_):
            print("something wrong %s" % path)
            return False
        return True
    return False


def path_to_str(path: Union[str, Path]) -> Optional[str]:
    if type(path) == Path:
        return str(path)
    elif type(path) == str:
        return path
    return None
