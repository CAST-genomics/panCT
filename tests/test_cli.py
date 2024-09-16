from typer.testing import CliRunner
from panct.__main__ import app

runner = CliRunner()

def test_complexity_cli():
    short_gfa_test_file = "tests/data/short_test.gfa"
    short_tsv_test_file = "tests/data/short_test_node_map.tsv"
    output_file = "tests/data/test_output.txt"

    result = runner.invoke(app, ["complexity", short_gfa_test_file, short_tsv_test_file])

    assert result.exit_code == 0
    assert "Complexity" in result.output
