import abc
from io import TextIOWrapper
from itertools import takewhile
from pathlib import Path
from typing import Dict
from typing import Iterator
from typing import List
from typing import Optional
from typing import Union

import numpy as np

from bintools.utils.utils import check_path


class posterior:
    """
    abstract class of posterior distribution
    """

    def __init__(
        self,
        list_of_sampleID: List[int],
        list_of_trees: List[str],
        list_of_phi: List[List[float]],
        list_of_rho: List[List[float]],
        list_of_omega: List[float],
        burnin: int = 0,
    ) -> None:
        self.list_of_sampleID: List[int] = list_of_sampleID
        self.list_of_trees: List[str] = list_of_trees
        self.list_of_phi: List[List[float]] = list_of_phi
        self.list_of_rho: List[List[float]] = list_of_rho
        self.list_of_omega: List[float] = list_of_omega
        self.burnin: int = burnin

    def __sizeof__(self) -> int:
        return len(self.list_of_sampleID)

    @abc.abstractmethod
    def parse_mcmc(self, mcmc_path: Path, burnin: int):
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


def chunck_chain(lines, chunck_size):
    for i, j in zip(
        range(0, (len(lines) - chunck_size), chunck_size),
        range(chunck_size, len(lines), chunck_size),
    ):
        yield lines[i:j]


class posterior_MGTRtsCpG_SNCatAA(posterior):
    """
    need to read MCMC + ABC parameter values
    """

    def __init__(
        self,
        list_of_sampleID: List[int],
        list_of_trees: List[str],
        list_of_phi: List[List[float]],
        list_of_rho: List[List[float]],
        list_of_omega: List[float],
        number_of_aa_profiles: int,
        list_of_aa_profiles: List[List[List[float]]],
        list_of_alloc: List[List[int]],
        burnin: int = 0,
    ):
        super().__init__(
            list_of_sampleID=list_of_sampleID,
            list_of_trees=list_of_trees,
            list_of_phi=list_of_phi,
            list_of_rho=list_of_rho,
            list_of_omega=list_of_omega,
            burnin=burnin,
        )
        self.number_of_aa_profiles: int = number_of_aa_profiles
        self.list_of_aa_profiles: List[List[List[float]]] = list_of_aa_profiles
        self.list_of_alloc: List[List[int]] = list_of_alloc

    def parse_mcmc(
        self, mcmc_path: Path, burnin: int = 0
    ) -> Optional["posterior_MGTRtsCpG_SNCatAA"]:
        """
        parses mcmc
        """
        list_of_sampleID: List[int] = []
        list_of_trees: List[str] = []
        list_of_phi: List[List[float]] = []
        list_of_rho: List[List[float]] = []
        list_of_omega: List[float] = []
        number_of_aa_profiles: int = 0
        list_of_aa_profiles: List[List[List[float]]] = []
        list_of_alloc: List[List[int]] = []
        try:
            with open(mcmc_path, "r") as hl:
                lines: List[str] = hl.readlines()
                number_of_aa_profiles = int(lines[11])
                for sampleID, chunck in enumerate(
                    chunck_chain(lines=lines, chunck_size=number_of_aa_profiles + 15)
                ):
                    list_of_sampleID += [sampleID]
                    try:
                        list_of_trees += [chunck[0].strip()]
                    except Exception as e:
                        print(
                            "something wrong with %s %s" % ("parsing the tree", str(e))
                        )

                    try:
                        list_of_phi += [
                            np.fromstring(chunck[3], dtype=float, sep="\t").tolist()
                        ]
                    except Exception as e:
                        print(
                            "something wrong with %s %s" % ("parsing the phi", str(e))
                        )

                    try:
                        list_of_rho += [
                            np.fromstring(chunck[5], dtype=float, sep="\t").tolist()
                        ]
                    except Exception as e:
                        print(
                            "something wrong with %s %s" % ("parsing the rho", str(e))
                        )

                    try:
                        list_of_omega += [float(chunck[9].strip())]
                    except Exception as e:
                        print(
                            "something wrong with %s %s" % ("parsing the omega", str(e))
                        )

                    try:
                        cur_list_of_aa_profiles: List[List[float]] = []
                        for i in range(14, (number_of_aa_profiles + 14)):
                            cur_list_of_aa_profiles += [
                                np.fromstring(chunck[i], dtype=float, sep="\t").tolist()
                            ]
                        list_of_aa_profiles += [cur_list_of_aa_profiles]
                    except Exception as e:
                        print(
                            "something wrong with %s %s"
                            % ("parsing the aa profiles", str(e))
                        )

                    try:
                        list_of_alloc += [
                            np.fromstring(
                                chunck[self.number_of_aa_profiles + 14],
                                dtype=int,
                                sep="\t",
                            ).tolist()
                        ]
                    except Exception as e:
                        print(
                            "something wrong with %s %s"
                            % ("parsing the aa profile allocation", str(e))
                        )

            return posterior_MGTRtsCpG_SNCatAA(
                list_of_sampleID=list_of_sampleID,
                list_of_trees=list_of_trees,
                list_of_phi=list_of_phi,
                list_of_rho=list_of_rho,
                list_of_omega=list_of_omega,
                number_of_aa_profiles=number_of_aa_profiles,
                list_of_aa_profiles=list_of_aa_profiles,
                list_of_alloc=list_of_alloc,
                burnin=burnin,
            )
        except Exception as e:
            print("something wrong %s when parsing %s" % (str(e), mcmc_path.__str__()))
            return None

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


class posterior_M0_GTR(posterior):
    """
    posterior distribution of M0-GTR F1 x 4 (Muse and Goat 1994) generated under Phylobayes-MPI
    """

    def __init__(
        self,
        list_of_sampleID: List[int],
        list_of_trees: List[str],
        list_of_phi: List[List[float]],
        list_of_rho: List[List[float]],
        list_of_omega: List[float],
        burnin: int = 0,
    ):
        super().__init__(
            list_of_sampleID=list_of_sampleID,
            list_of_trees=list_of_trees,
            list_of_phi=list_of_phi,
            list_of_rho=list_of_rho,
            list_of_omega=list_of_omega,
            burnin=burnin,
        )

    @classmethod
    def parse_mcmc(
        cls, mcmc_path: Path, burnin: int = 0
    ) -> Optional["posterior_M0_GTR"]:
        """
        parses mcmc
        """
        list_of_sampleID: List[int] = []
        list_of_trees: List[str] = []
        list_of_phi: List[List[float]] = []
        list_of_rho: List[List[float]] = []
        list_of_omega: List[float] = []

        try:
            with open(mcmc_path, "r") as hl:
                lines: List[str] = hl.readlines()
                k: int = 0
                sampleID: int = 0
                for line in lines:
                    if k == 0:
                        list_of_trees += [line.strip()]
                    if k == 3:
                        list_of_phi += [
                            np.fromstring(line, dtype=float, sep="\t").tolist()
                        ]
                    if k == 5:
                        list_of_rho += [
                            np.fromstring(line, dtype=float, sep="\t").tolist()
                        ]
                    if k == 9:
                        list_of_omega += np.fromstring(
                            line, dtype=float, sep="\t"
                        ).tolist()

                    k += 1
                    if k == 15:
                        list_of_sampleID += [sampleID]
                        sampleID += 1
                        k = 0

            return posterior_M0_GTR(
                list_of_sampleID=list_of_sampleID,
                list_of_trees=list_of_trees,
                list_of_phi=list_of_phi,
                list_of_rho=list_of_rho,
                list_of_omega=list_of_omega,
                burnin=burnin,
            )
        except Exception as e:
            print("something wrong %s when parsing %s" % (str(e), mcmc_path.__str__()))
            return None

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
    def read(cls, input_file: Path) -> Optional["posterior_M0_GTR"]:
        list_of_sampleID: List[int] = []
        list_of_trees: List[str] = []
        list_of_phi: List[List[float]] = []
        list_of_rho: List[List[float]] = []
        list_of_omega: List[float] = []
        if not check_path(input_file):
            raise RuntimeError

        try:
            with open(str(input_file), "r") as file_handler:

                lines: Iterator[str] = iter(file_handler.readlines())
                sampleID: int = 0
                for line in takewhile(lambda x: x is not None, lines):
                    list_of_trees += [line.strip()]

                    line = next(lines)
                    list_of_phi += [
                        np.fromstring(line.strip(), dtype=float, sep="\t").tolist()
                    ]

                    line = next(lines)
                    list_of_rho += [
                        np.fromstring(line.strip(), dtype=float, sep="\t").tolist()
                    ]

                    line = next(lines)
                    list_of_omega += [
                        np.fromstring(line.strip(), dtype=float, sep="\t").tolist()
                    ]
                    list_of_sampleID += [sampleID]
                    sampleID += 1

                return posterior_M0_GTR(
                    list_of_sampleID=list_of_sampleID,
                    list_of_trees=list_of_trees,
                    list_of_phi=list_of_phi,
                    list_of_rho=list_of_rho,
                    list_of_omega=list_of_omega,
                )

        except Exception as e:
            print("something wrong %s, %s" % (input_file, str(e)))
            return None
