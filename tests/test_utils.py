import os
from panct.utils import *
from pathlib import Path
import pytest

DATADIR = Path(__file__).parent.joinpath("data")


def test_parse_region_string():
    # Valid region
    region = parse_region_string("chr1:123-456")
    assert region.chrom == "chr1"
    assert region.start == 123
    assert region.end == 456

    # Invalid region
    with pytest.raises(ValueError):
        parse_region_string("invalid region string")
    with pytest.raises(ValueError):
        parse_region_string(1234)

    # start>end
    with pytest.raises(ValueError):
        parse_region_string("1:8-3")


def test_parse_regions_file():
    # Valid regions file
    regions_list = parse_regions_file(os.path.join(DATADIR, "valid_regions.bed"))
    assert len(regions_list) == 2
    assert regions_list[0].chrom == "chr1"
    assert regions_list[0].start == 123
    assert regions_list[0].end == 456
    assert regions_list[1].chrom == "chrX"
    assert regions_list[1].start == 45678
    assert regions_list[1].end == 89098

    # Malformatted regions files
    with pytest.raises(ValueError):
        parse_regions_file(os.path.join(DATADIR, "invalid_regions1.bed"))
    with pytest.raises(ValueError):
        parse_regions_file(os.path.join(DATADIR, "invalid_regions2.bed"))
    with pytest.raises(ValueError):
        parse_regions_file(os.path.join(DATADIR, "invalid_regions3.bed"))
