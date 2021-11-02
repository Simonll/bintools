import os
import pathlib
from typing import Any
from typing import Dict
from typing import Optional

from bintools.utils.utils import check_path
from bintools.utils.utils import get_yaml_config


class workflow:
    def __init__(self, config_file: pathlib.Path, root_local: pathlib.Path) -> None:
        """
        generic construct for phylogenetic workflow
        """

        if not check_path(config_file):
            print("something wrong with config file %s" % config_file)
            raise RuntimeError
        dict_of_dirs: Optional[Dict[Any, Any]] = get_yaml_config(yaml_file=config_file)
        if dict_of_dirs is not None:
            self.dict_of_dirs: Dict[Any, Any] = dict_of_dirs

        self.dict_of_dirs["root_local"] = root_local.__str__()
        self.root_local_root_dir: str = self.dict_of_dirs["root_local"]
        self.mapping: str = (
            self.get_dir(type="local", dir="exp")[:-1]
            + ":"
            + self.get_dir(type="mapped", dir="exp")
        )
        if self.make_dirs():
            print("workflow dirs generated")

    def make_dirs(self) -> bool:
        for k in self.dict_of_dirs.keys():
            try:
                os.makedirs(self.get_dir(type="local", dir=k), exist_ok=True)
            except Exception as e:
                print("something wrong with %s %s " % (k, e))
                return False
        return True

    def get_dir(self, type: str, dir: str):

        if type not in ["local", "mapped"]:
            print("Something wrong with type %s" % type)
            raise RuntimeError

        if dir not in list(self.dict_of_dirs.keys()):
            print("Something wrong with type %s" % dir)
            raise RuntimeError

        if type == "local":
            return self.root_local_root_dir[:-1] + self.dict_of_dirs[dir]

        return self.dict_of_dirs[dir]

    def generate_log_file(self, output: str, dir: str) -> str:
        return "2> " + self.get_dir(type="local", dir=dir) + output + ".log"
