import os
import shutil
import filecmp
from pathlib import Path

from typer.testing import CliRunner

from panct.__main__ import app

runner = CliRunner()

DATADIR = Path(__file__).parent.joinpath("data")

def test_basic(capfd):
    file = DATADIR / "basic.gfa"
    tmp_file = Path("test.walk")

    # copy the file so that we don't affect anything in the tests/data directory
    shutil.copy(str(file), str(tmp_file))

    cmd = f"walks {tmp_file}"
    result = runner.invoke(app, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == ""
    # check that the output .walk file is the same as the file in tests/data/
    assert filecmp.cmp(tmp_file, DATADIR / "basic.walk")
    # TODO: # check that the output.walk.gz.tbi file exists
    # assert tmp_file.with_suffix(".walk.gz.tbi").is_file()
    # assert result.exit_code == 0

    tmp_file.unlink()
    # tmp_file.with_suffix(".walk.gz").unlink()
    # tmp_file.with_suffix(".walk.gz.tbi").unlink()
