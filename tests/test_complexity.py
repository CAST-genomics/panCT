import os
from pathlib import Path

DATADIR = Path(__file__).parent.joinpath("data")

from panct.complexity import *


def test_basic():
    assert True
