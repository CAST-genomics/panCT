from typer.testing import CliRunner
from panct.__main__ import app
from pathlib import Path

DATADIR = Path(__file__).parent.joinpath("data")
runner = CliRunner()

def test_complexity_cli():
    short_gfa_test_file = DATADIR / "short_test.gfa"
    short_tsv_test_file = DATADIR / "short_test_node_map.tsv"
    output_file = DATADIR / "test_output.txt"
    result = runner.invoke(app, ["complexity", str(short_gfa_test_file), str(short_tsv_test_file), str(output_file)])

    assert result.exit_code == 0
    assert "Complexity" in result.output
