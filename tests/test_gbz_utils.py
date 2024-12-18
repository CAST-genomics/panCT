import os
from panct.gbz_utils import *
from pathlib import Path
import pytest
from panct.logging import getLogger

DATADIR = Path(__file__).parent.joinpath("data")


def test_check_gbzfile():
    log = getLogger(name="complexity", level="INFO")

    # Non-existing file
    assert not check_gbzfile(Path("/path/to/bad/file"), log)

    # Check existing file
    assert check_gbzfile(DATADIR / "basic.gbz", log)

    # Check existing file but not indexed
    # TODO - add rest of checks once add gbz-base to
    # list of test dependencies
    # check_gbzfile(DATADIR / "basic_noseq.gbz", log)
