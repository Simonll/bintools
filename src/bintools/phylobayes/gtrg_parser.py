from io import TextIOWrapper
from itertools import takewhile
from pathlib import Path
from typing import Dict
from typing import Iterator
from typing import List
from typing import Optional
from typing import Union

import numpy as np

from bintools.phylobayes.mcmc_parser import posterior
from bintools.utils.utils import check_path


class posterior_GTRG(posterior):
    """
    posterior distribution of GTR Tavaré 1986
    """

    def __init__(
        self,
        list_of_sampleID: List[int],
        list_of_trees: List[str],
        list_of_phi: List[List[float]],
        list_of_rho: List[List[float]],
        list_of_alpha: List[float],
        burnin: int = 0,
    ):
        super().__init__(
            list_of_sampleID=list_of_sampleID,
            list_of_trees=list_of_trees,
            list_of_phi=list_of_phi,
            list_of_rho=list_of_rho,
            burnin=burnin,
        )
        self.list_of_alpha: List[float] = list_of_alpha

    @classmethod
    def parse_mcmc(cls, mcmc_path: Path, burnin: int = 0) -> Optional["posterior_GTRG"]:
        """
        parses mcmc
        """
        list_of_sampleID: List[int] = []
        list_of_trees: List[str] = []
        list_of_phi: List[List[float]] = []
        list_of_rho: List[List[float]] = []
        list_of_alpha: List[float] = []
        try:
            with open(mcmc_path, "r") as hl:
                lines: List[str] = hl.readlines()
                k: int = 0
                sampleID: int = 0
                for line in lines:
                    if k == 0:
                        list_of_trees += [line.strip()]
                    if k == 3:
                        list_of_alpha += [
                            np.fromstring(line, dtype=float, sep="\t").tolist()[0]
                        ]
                    if k == 7:
                        list_of_rho += [
                            np.fromstring(line, dtype=float, sep="\t").tolist()
                        ]
                    if k == 9:
                        list_of_phi += [
                            np.fromstring(line, dtype=float, sep="\t").tolist()
                        ]
                    k += 1
                    if k == 11:
                        list_of_sampleID += [sampleID]
                        sampleID += 1
                        k = 0

            return posterior_GTRG(
                list_of_sampleID=list_of_sampleID,
                list_of_trees=list_of_trees,
                list_of_phi=list_of_phi,
                list_of_rho=list_of_rho,
                list_of_alpha=list_of_alpha,
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
            "tree": self.list_of_trees[u_int],
            "alpha": self.list_of_alpha[u_int],
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
            "tree": self.list_of_trees[u_int],
            "alpha": self.list_of_alpha[u_int],
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
            file_handler.write(str(dict_of_params["alpha"]))
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
            return True
        except Exception as e:
            print("something wrong when writing paramter values %s" % str(e))
            return False

    @classmethod
    def read(cls, input_file: Path) -> Optional["posterior_GTRG"]:
        list_of_sampleID: List[int] = []
        list_of_trees: List[str] = []
        list_of_phi: List[List[float]] = []
        list_of_rho: List[List[float]] = []
        list_of_alpha: List[float] = []
        if not check_path(input_file):
            raise RuntimeError

        try:
            with open(str(input_file), "r") as file_handler:

                lines: Iterator[str] = iter(file_handler.readlines())
                sampleID: int = 0
                for line in takewhile(lambda x: x is not None, lines):
                    list_of_trees += [line.strip()]
                    line = next(lines)
                    list_of_alpha += [
                        np.fromstring(line.strip(), dtype=float, sep="\t").tolist()[0]
                    ]
                    line = next(lines)
                    list_of_phi += [
                        np.fromstring(line.strip(), dtype=float, sep="\t").tolist()
                    ]

                    line = next(lines)
                    list_of_rho += [
                        np.fromstring(line.strip(), dtype=float, sep="\t").tolist()
                    ]

                    list_of_sampleID += [sampleID]
                    sampleID += 1

                return posterior_GTRG(
                    list_of_sampleID=list_of_sampleID,
                    list_of_trees=list_of_trees,
                    list_of_phi=list_of_phi,
                    list_of_rho=list_of_rho,
                    list_of_alpha=list_of_alpha,
                )

        except Exception as e:
            print("something wrong %s, %s" % (input_file, str(e)))
            return None
