import abc
from io import TextIOWrapper
from itertools import takewhile
from pathlib import Path
from typing import Dict
from typing import List
from typing import Union

import numpy as np

from utils.utils import check_path


class posterior:
    """
    abstract class of posterior distribution
    """

    def __init__(self, mcmc_path: Path, burnin: int):
        self.mcmc_path: Path = mcmc_path
        self.burnin: int = burnin
        assert check_path(self.mcmc_path)
        assert burnin >= 0

    @abc.abstractmethod
    def parse_mcmc(self):
        """
        parses mcmc
        """
        raise RuntimeError()

    @abc.abstractmethod
    def sample(self):
        """
        draw values
        """
        raise RuntimeError()


class posterior_MGTRtsCpG_SNCatAA(posterior):
    """
    need to read MCMC + ABC parameter values
    """

    def __init__(self, mcmc_path: Path, abc_path: Path, burnin: int):
        super().__init__(mcmc_path=mcmc_path, burnin=burnin)

        self.list_of_trees: List[str] = []
        self.list_of_phi: List[List[float]] = []
        self.list_of_rho: List[List[float]] = []
        self.list_of_omega: List[float] = []
        self.number_of_profiles: List[int] = []
        self.list_of_aa_profiles: List[List[List[float]]] = []
        self.list_of_alloc: List[List[int]] = []
        self.parse_mcmc()

    def parse_mcmc(self) -> bool:
        """
        parses mcmc
        """
        try:
            with open(self.mcmc_path, "r") as hl:
                lines = hl.readlines()
                k = 0
                chainID = 0
                for line in lines:
                    if k == 0:
                        self.list_of_trees += [line.strip()]
                    if k == 3:
                        self.list_of_phi += [
                            np.fromstring(line, dtype=float, sep="\t").tolist()
                        ]
                    if k == 5:
                        self.list_of_rho += [
                            np.fromstring(line, dtype=float, sep="\t").tolist()
                        ]
                    if k == 10:
                        self.list_of_omega += np.fromstring(
                            line, dtype=float, sep="\t"
                        ).tolist()

                    if k == 12:
                        self.number_of_profiles += np.fromstring(
                            line, dtype=int, sep="\t"
                        ).to_list()
                    k += 1
                    if k == 15:
                        k_profile = 0
                        cur_list_of_profiles: List[List[float]] = []
                        while k_profile < self.number_of_profiles[chainID]:
                            cur_list_of_profiles += [
                                np.fromstring(line, dtype=float, sep="\t").tolist()
                            ]
                            k_profile += 1
                        self.list_of_aa_profiles += [cur_list_of_profiles]
                        k += 1
                    if k == self.number_of_profiles[chainID] + 15:
                        self.list_of_alloc += [
                            np.fromstring(line, dtype=int, sep="\t").tolist()
                        ]
                        k = 0
                        chainID += 1

            return True
        except Exception as e:
            print("something wrong %s when parsing %s" % (str(e), self.mcmc_path))
            return False


class posterior_M0_GTR(posterior):
    """
    posterior distribution of M0-GTR F1 x 4 (Muse and Goat 1994) generated under Phylobayes-MPI
    """

    def __init__(self, mcmc_path: Path, burnin: int):
        super().__init__(mcmc_path=mcmc_path, burnin=burnin)

        self.list_of_trees: List[str] = []
        self.list_of_phi: List[List[float]] = []
        self.list_of_rho: List[List[float]] = []
        self.list_of_omega: List[float] = []
        self.parse_mcmc()

    def parse_mcmc(self) -> bool:
        """
        parses mcmc
        """
        try:
            with open(self.mcmc_path, "r") as hl:
                lines = hl.readlines()
                k = 0
                for line in lines:
                    if k == 0:
                        self.list_of_trees += [line.strip()]
                    if k == 3:
                        self.list_of_phi += [
                            np.fromstring(line, dtype=float, sep="\t").tolist()
                        ]
                    if k == 5:
                        self.list_of_rho += [
                            np.fromstring(line, dtype=float, sep="\t").tolist()
                        ]
                    if k == 9:
                        self.list_of_omega += np.fromstring(
                            line, dtype=float, sep="\t"
                        ).tolist()

                    k += 1
                    if k == 15:
                        k = 0

            return True
        except Exception as e:
            print("something wrong %s when parsing %s" % (str(e), self.mcmc_path))
            return False

    def sample(self) -> Dict[str, Union[float, str]]:
        """
        samples paramter values from burnin to end of mcmc and return dictionary
        """
        u_int: int = np.random.randint(self.burnin, len(self.list_of_rho))
        dict_of_params: Dict[str, Union[float, str]] = {
            "phi_A": self.list_of_phi[u_int][0],
            "phi_C": self.list_of_phi[u_int][1],
            "phi_G": self.list_of_phi[u_int][2],
            "phi_T": self.list_of_phi[u_int][3],
            "rho_AC": self.list_of_rho[u_int][0],
            "rho_AG": self.list_of_rho[u_int][1],
            "rho_AT": self.list_of_rho[u_int][2],
            "rho_CG": self.list_of_rho[u_int][3],
            "rho_CT": self.list_of_rho[u_int][4],
            "rho_TG": self.list_of_rho[u_int][5],
            "omega": self.list_of_omega[u_int],
            "tree": self.list_of_trees[u_int],
        }

        return dict_of_params

    def sample_hky(self):
        """
        samples paramter values from burnin to end of mcmc and return dictionary with GTR transformed to HKY
        """
        u_int: int = np.random.randint(self.burnin, len(self.list_of_rho))
        tv: float = np.sum([self.list_of_rho[u_int][i] for i in [0, 2, 3, 5]])
        ts: float = np.sum([self.list_of_rho[u_int][i] for i in [1, 4]])
        ts /= 2.0
        tv /= 4.0
        assert np.isclose(a=ts * 2 + tv * 4, b=1, rtol=0.001)
        dict_of_params: Dict[str, Union[float, str]] = {
            "phi_A": self.list_of_phi[u_int][0],
            "phi_C": self.list_of_phi[u_int][1],
            "phi_G": self.list_of_phi[u_int][2],
            "phi_T": self.list_of_phi[u_int][3],
            "rho_AC": tv,
            "rho_AG": ts,
            "rho_AT": tv,
            "rho_CG": tv,
            "rho_CT": ts,
            "rho_TG": tv,
            "omega": self.list_of_omega[u_int],
            "tree": self.list_of_trees[u_int],
        }

        return dict_of_params

    def write_values(
        self,
        dict_of_params: Dict[str, Union[float, str, Dict[int, float]]],
        file_handler: TextIOWrapper,
    ) -> bool:
        """
        writes parameter values from dict of params
        """
        try:
            file_handler.write(str(dict_of_params["tree"]))
            file_handler.write("\n")
            file_handler.write(
                "\t".join(
                    [
                        str(dict_of_params["phi_A"]),
                        str(dict_of_params["phi_C"]),
                        str(dict_of_params["phi_G"]),
                        str(dict_of_params["phi_T"]),
                    ]
                )
            )
            file_handler.write("\n")
            file_handler.write(
                "\t".join(
                    [
                        str(dict_of_params["rho_AC"]),
                        str(dict_of_params["rho_AG"]),
                        str(dict_of_params["rho_AT"]),
                        str(dict_of_params["rho_CG"]),
                        str(dict_of_params["rho_CT"]),
                        str(dict_of_params["rho_TG"]),
                    ]
                )
            )
            file_handler.write("\n")
            if (
                "omega_site" in dict_of_params
                and type(dict_of_params["omega_site"]) == Dict[int, float]
            ):
                file_handler.write(
                    "\t".join(
                        [str(i) for i in list(dict_of_params["omega_site"].values())]
                    )
                )
                file_handler.write("\n")
            return True
        except Exception as e:
            print("something wrong when writing paramter values %s" % str(e))
            return False

    @classmethod
    def read(cls, input_file: Path):
        if not check_path(input_file):
            raise RuntimeError

        with open(str(input_file), "r") as file_handler:

            lines = iter(file_handler.readlines())
            for line in takewhile(lambda x: x is not None, lines):
                cls.list_of_trees += [line.strip()]

                line = next(lines)
                cls.list_of_phi += [
                    np.fromstring(line.strip(), dtype=float, sep="\t").tolist()
                ]

                line = next(lines)
                cls.list_of_rho += [
                    np.fromstring(line.strip(), dtype=float, sep="\t").tolist()
                ]

                line = next(lines)
                cls.list_of_omega += [
                    np.fromstring(line.strip(), dtype=float, sep="\t").tolist()
                ]
