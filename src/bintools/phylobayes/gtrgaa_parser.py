from io import TextIOWrapper
from itertools import takewhile
from pathlib import Path
from typing import Dict
from typing import Iterator
from typing import List
from typing import Optional
from typing import Tuple
from typing import Union

import numpy as np

from bintools.phylobayes.mcmc_parser import posterior
from bintools.utils.utils import check_path


class posterior_GTRGAA(posterior):
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
        self.list_of_phiID: List[str] = [
            "A",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "K",
            "L",
            "M",
            "N",
            "P",
            "Q",
            "R",
            "S",
            "T",
            "V",
            "W",
            "Y",
        ]
        self.list_of_rhoID: List[str] = [
            "VW",
            "AF",
            "RS",
            "EI",
            "DL",
            "HW",
            "CG",
            "AS",
            "IR",
            "AV",
            "FK",
            "LS",
            "CE",
            "GY",
            "EL",
            "CR",
            "EK",
            "TW",
            "HL",
            "VY",
            "AI",
            "IS",
            "QY",
            "PV",
            "AK",
            "GR",
            "MR",
            "GI",
            "CT",
            "TV",
            "FG",
            "AY",
            "HR",
            "HT",
            "DR",
            "IV",
            "CN",
            "CY",
            "CQ",
            "LQ",
            "IK",
            "IQ",
            "NR",
            "RV",
            "LW",
            "CV",
            "NS",
            "LN",
            "FM",
            "GN",
            "EV",
            "FQ",
            "DF",
            "IM",
            "ER",
            "CL",
            "TY",
            "QW",
            "FY",
            "DT",
            "SV",
            "PR",
            "KQ",
            "FN",
            "HP",
            "GW",
            "RY",
            "LR",
            "SW",
            "AD",
            "PY",
            "CS",
            "FL",
            "FV",
            "NV",
            "KW",
            "NP",
            "CK",
            "CP",
            "NT",
            "MT",
            "NW",
            "EH",
            "LP",
            "PT",
            "DK",
            "EF",
            "GS",
            "NY",
            "DV",
            "FW",
            "KS",
            "CD",
            "CI",
            "DQ",
            "DM",
            "KP",
            "HV",
            "EW",
            "FS",
            "FP",
            "SY",
            "GP",
            "EQ",
            "IT",
            "EG",
            "GH",
            "DN",
            "ET",
            "HK",
            "NQ",
            "AH",
            "LM",
            "MV",
            "AQ",
            "QT",
            "AT",
            "IN",
            "HQ",
            "AL",
            "GM",
            "AW",
            "GK",
            "MY",
            "CH",
            "EM",
            "KT",
            "PW",
            "QS",
            "QR",
            "GQ",
            "LT",
            "MN",
            "AR",
            "MW",
            "FI",
            "DH",
            "AG",
            "MP",
            "RW",
            "CW",
            "AN",
            "HY",
            "DG",
            "AM",
            "KV",
            "KY",
            "KM",
            "PS",
            "AP",
            "RT",
            "WY",
            "LV",
            "QV",
            "CF",
            "HS",
            "IP",
            "EP",
            "PQ",
            "ST",
            "CM",
            "IY",
            "HN",
            "HI",
            "AE",
            "AC",
            "EY",
            "GV",
            "KR",
            "FR",
            "FT",
            "GT",
            "LY",
            "MQ",
            "DI",
            "FH",
            "DY",
            "EN",
            "KN",
            "HM",
            "DP",
            "IW",
            "DW",
            "MS",
            "GL",
            "DS",
            "ES",
            "IL",
            "KL",
            "DE",
        ]

    @classmethod
    def parse_mcmc(
        cls, mcmc_path: Path, burnin: int = 0
    ) -> Optional["posterior_GTRGAA"]:
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

            return posterior_GTRGAA(
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
            "tree": self.list_of_trees[u_int],
            "alpha": self.list_of_alpha[u_int],
        }
        dict_of_params.update(
            {
                "phi_" + phi: self.list_of_phi[u_int][i]
                for i, phi in enumerate(self.list_of_phiID)
            }
        )

        dict_of_params.update(
            {
                "rho_" + rho: self.list_of_phi[u_int][i]
                for i, rho in enumerate(self.list_of_rhoID)
            }
        )
        return dict_of_params

    def get_posterior_mean(self):
        """
        return posterior mean
        """
        dict_of_params: Dict[str, Tuple[float, float]] = {
            "alpha": (
                np.mean([i for i in self.list_of_alpha]),
                np.std([i for i in self.list_of_alpha]),
            ),
        }
        dict_of_params.update(
            {
                "phi_"
                + phi: (
                    np.mean([i[j] / np.sum(i) for i in self.list_of_phi]),
                    np.std([i[j] / np.sum(i) for i in self.list_of_phi]),
                )
                for j, phi in enumerate(self.list_of_phiID)
            }
        )
        dict_of_params.update(
            {
                "rho_"
                + rho: (
                    np.mean([i[j] / np.sum(i) for i in self.list_of_rho]),
                    np.std([i[j] / np.sum(i) for i in self.list_of_rho]),
                )
                for j, rho in enumerate(self.list_of_rhoID)
            }
        )

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
                    [str(dict_of_params["phi_" + phi]) for phi in self.list_of_phiID]
                )
            )
            file_handler.write("\n")
            file_handler.write(
                "\t".join(
                    [str(dict_of_params["rho_" + rho]) for rho in self.list_of_rhoID]
                )
            )
            file_handler.write("\n")
            return True
        except Exception as e:
            print("something wrong when writing paramter values %s" % str(e))
            return False

    @classmethod
    def read(cls, input_file: Path) -> Optional["posterior_GTRGAA"]:
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

                return posterior_GTRGAA(
                    list_of_sampleID=list_of_sampleID,
                    list_of_trees=list_of_trees,
                    list_of_phi=list_of_phi,
                    list_of_rho=list_of_rho,
                    list_of_alpha=list_of_alpha,
                )

        except Exception as e:
            print("something wrong %s, %s" % (input_file, str(e)))
            return None
