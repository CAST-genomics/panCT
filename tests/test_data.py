import os
from pathlib import Path

import pytest

from panct.data import Region, Regions

DATADIR = Path(__file__).parent.joinpath("data")


class TestRegions:
    def test_parse_region_string(self):
        # Valid region
        region = Region.load("chr1:123-456")
        assert region.chrom == "chr1"
        assert region.start == 123
        assert region.end == 456

        # Invalid region
        with pytest.raises(ValueError):
            Region.load("invalid region string")
        with pytest.raises(ValueError):
            Region.load(1234)

        # start>end
        with pytest.raises(ValueError):
            Region.load("1:8-3")

    def test_parse_regions_file(self):
        # Valid regions file
        regions_list = Regions.load(DATADIR / "valid_regions.bed")
        assert len(regions_list) == 2
        assert regions_list[0].chrom == "chr1"
        assert regions_list[0].start == 123
        assert regions_list[0].end == 456
        assert regions_list[1].chrom == "chrX"
        assert regions_list[1].start == 45678
        assert regions_list[1].end == 89098

        # Malformatted regions files
        with pytest.raises(ValueError):
            Regions.load(DATADIR / "invalid_regions1.bed")
        with pytest.raises(ValueError):
            Regions.load(DATADIR / "invalid_regions2.bed")
        with pytest.raises(ValueError):
            Regions.load(DATADIR / "invalid_regions3.bed")
