import os
import sys

_BT_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _BT_SRC not in sys.path:
    sys.path.insert(0, _BT_SRC)
