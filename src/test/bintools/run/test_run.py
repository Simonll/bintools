import unittest
from typing import Any
from typing import Dict
from typing import Optional

from bintools.run.run import generate_pb_mpi_cmd


class TestRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self) -> None:
        pass

    def tearDown(self) -> None:
        pass

    def test_generate_pb_mpi_cmd(self) -> None:
        kwargs: Dict[str, Any] = {
            "-gtr": "",
            "-T": "",
            "-s": "",
            "-np": "4",
            "-chainname": "test_a",
        }
        cmd: Optional[str] = generate_pb_mpi_cmd(
            method="pb_mpi", mapping="/tmp/data/:/data/", **kwargs
        )


if __name__ == "__main__":
    unittest.main()
