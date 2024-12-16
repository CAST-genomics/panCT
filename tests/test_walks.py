import os
import gzip
import shutil
import filecmp
from pathlib import Path

from typer.testing import CliRunner

from panct.__main__ import app

runner = CliRunner()

DATADIR = Path(__file__).parent.joinpath("data")


def test_basic_wo_gz(capfd):
    in_file = DATADIR / "basic.gfa"
    out_file = Path("test.walk")
    exp_file = in_file.with_suffix(".walk")

    # create a simple test.walk file
    cmd = f"walks -o {out_file} {in_file}"
    result = runner.invoke(app, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == ""
    # check that the output .walk file is the same as the file in tests/data/
    assert filecmp.cmp(out_file, exp_file)
    assert result.exit_code == 0

    out_file.unlink()


def test_basic_gz(capfd):
    in_file = DATADIR / "basic.gfa"
    out_file = Path("test.walk.gz")
    exp_file = in_file.with_suffix(".walk")

    # copy the file so that we don't affect anything in the tests/data directory
    tmp_in_file = out_file.with_suffix("").with_suffix(".gfa")
    shutil.copy(str(in_file), tmp_in_file)

    # by default: we also create a gz file and its index
    cmd = f"walks {tmp_in_file}"
    result = runner.invoke(app, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == ""
    with gzip.open(out_file, "rb") as f:
        out_file_content = f.read().decode('utf-8')
    with open(exp_file, "r") as f:
        exp_file_content = f.read()
    # check that the output .walk file is the same as the file in tests/data/
    assert out_file_content == exp_file_content
    # check that an index was also automatically generated
    assert out_file.with_suffix(".gz.tbi").is_file()
    assert result.exit_code == 0

    out_file.unlink()
    tmp_in_file.unlink()
    out_file.with_suffix(".gz.tbi").unlink()
