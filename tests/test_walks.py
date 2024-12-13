import os
from pathlib import Path

from typer.testing import CliRunner

from panct.__main__ import main

runner = CliRunner()

DATADIR = Path(__file__).parent.joinpath("data")

def test_basic(capfd):
    file = DATADIR / "basic.gfa"
    tmp_file = Path("test.walk")

    # copy the file so that we don't affect anything in the tests/data directory
    shutil.copy(str(file), str(tmp_file))

    cmd = f"walks {tmp_file}"
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == ""
    # check that the output .walk file is the same as the file in tests/data/
    with Data.hook_compressed(tmp_file.with_suffix(".walk.gz"), mode="rt") as haps:
        haps = filter(lambda l: not l.startswith("#"), haps.read().splitlines())
        with Data.hook_compressed(file.with_suffix(".walk.gz"), mode="rt") as expected:
            exp = filter(lambda l: not l.startswith("#"), expected.read().splitlines())
            assert list(haps) == list(exp)
    # check that the output .walk.gz.tbi file exists
    assert tmp_file.with_suffix(".walk.gz.tbi").is_file()
    assert result.exit_code == 0

    tmp_file.unlink()
    tmp_file.with_suffix(".walk.gz").unlink()
    tmp_file.with_suffix(".walk.gz.tbi").unlink()
