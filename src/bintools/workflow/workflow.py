import os
import pathlib
from typing import Any
from typing import Dict
from typing import List
from typing import Optional

from bintools.utils.utils import check_path
from bintools.utils.utils import get_yaml_config


class workflow:
    def __init__(
        self, config_file: pathlib.Path, root_local: pathlib.Path, test: bool = False
    ) -> None:
        """
        generic construct for phylogenetic workflow
        """

        if not check_path(config_file):
            print("something wrong with config file %s" % config_file)
            raise RuntimeError

        self.list_of_preexisting_dirs: List[str] = ["exp", "data", "scripts", "test"]
        self.dict_of_dirs: Dict[str, str] = {}
        workflow_config: Optional[Dict[Any, Any]] = get_yaml_config(
            yaml_file=config_file
        )
        self.root_local_dir: str = root_local.__str__() + "/"
        if workflow_config is not None:
            try:
                self.exp_dir: str = workflow_config["exp"]
            except Exception as e:
                print("something wrong with %s dir" % "exp")
                raise RuntimeError

            try:
                self.test_dir: str = workflow_config["test"]
            except Exception as e:
                print("something wrong with %s dir" % "test")

            try:
                self.scripts_dir: str = workflow_config["scripts"]
            except Exception as e:
                print("something wrong with %s dir" % "scripts")

            try:
                self.data_dir: str = workflow_config["data"]
            except Exception as e:
                print("something wrong with %s dir" % "data")

            for k, v in workflow_config.items():
                if k not in self.list_of_preexisting_dirs:
                    self.dict_of_dirs[k] = v

        self.test: bool = test
        if self.make_dirs():
            print("workflow dirs generated")

    def make_dirs(self) -> bool:
        for k in self.dict_of_dirs.keys():
            if k not in ["root_local", "scripts", "data"]:
                try:
                    os.makedirs(self.get_dir(type="local", dir=k), exist_ok=True)
                    print("dir: %s created" % self.get_dir(type="local", dir=k))
                except Exception as e:
                    print("something wrong with %s %s " % (k, str(e)))
                    return False
        return True

    def get_dir(self, type: str, dir: str)-> str:

        if type not in ["local", "mapped"]:
            print("Something wrong with type %s" % type)
            raise RuntimeError

        if dir not in list(self.dict_of_dirs.keys()) + self.list_of_preexisting_dirs:
            print("Something wrong with type %s" % dir)
            raise RuntimeError

        if type == "local":
            if dir == "scripts":
                return self.root_local_dir[:-1] + self.scripts_dir
            if dir == "data":
                return self.root_local_dir[:-1] + self.data_dir
            if dir == "exp":
                return self.root_local_dir[:-1] + self.exp_dir
            if self.test:
                return (
                    self.root_local_dir[:-1]
                    + self.test_dir[:-1]
                    + self.dict_of_dirs[dir]
                )
            else:
                return (
                    self.root_local_dir[:-1]
                    + self.exp_dir[:-1]
                    + self.dict_of_dirs[dir]
                )
        else:
            if dir in ["scripts"]:
                print("something wrong with %s dir" % dir)
                raise RuntimeError
            if dir == "scripts":
                return self.scripts_dir
            if dir == "data":
                return self.data_dir
            if dir == "exp":
                return self.exp_dir
            if self.test:
                return self.test_dir[:-1] + self.dict_of_dirs[dir]
            else:
                return self.exp_dir[:-1] + self.dict_of_dirs[dir]

    def generate_log_file(self, output: str, dir: str) -> str:
        return "2> " + self.get_dir(type="local", dir=dir) + output + ".log"
