import os
from pathlib import Path

import pytest

from panct.data import Region, Regions, Walks

DATADIR = Path(__file__).parent.joinpath("data")


class TestRegions:
    def test_parse_region_string(self):
        # Valid region
        region = Region.read("chr1:123-456")
        assert region.chrom == "chr1"
        assert region.start == 123
        assert region.end == 456

        # Invalid region
        with pytest.raises(ValueError):
            Region.read("invalid region string")
        with pytest.raises(ValueError):
            Region.read(1234)

        # start>end
        with pytest.raises(ValueError):
            Region.read("1:8-3")

    def test_parse_regions_file(self):
        # Valid regions file
        regions_list = Regions.read(DATADIR / "valid_regions.bed")
        assert len(regions_list) == 2
        assert regions_list[0].chrom == "chr1"
        assert regions_list[0].start == 123
        assert regions_list[0].end == 456
        assert regions_list[1].chrom == "chrX"
        assert regions_list[1].start == 45678
        assert regions_list[1].end == 89098

        # Malformatted regions files
        with pytest.raises(ValueError):
            Regions.read(DATADIR / "invalid_regions1.bed")
        with pytest.raises(ValueError):
            Regions.read(DATADIR / "invalid_regions2.bed")
        with pytest.raises(ValueError):
            Regions.read(DATADIR / "invalid_regions3.bed")


class TestWalks:
    def _get_dummy_walks(self):
        data = {}
        data[1] = set(("GRCh38", "samp1", "samp2"))
        data[2] = set(("GRCh38", "samp1"))
        return Walks(data=data)

    def test_parse_walks_file(self):
        expected = self._get_dummy_walks()

        nodes = Walks.read(DATADIR / "basic.walk")
        assert nodes.data == expected.data

        nodes = Walks.read(DATADIR / "basic.walk", region="1-2")
        assert nodes.data == expected.data

        nodes = Walks.read(DATADIR / "basic.walk", region="1-")
        assert nodes.data == expected.data

        nodes = Walks.read(DATADIR / "basic.walk", region="-2")
        assert nodes.data == expected.data

        nodes = Walks.read(DATADIR / "basic.walk.gz")
        assert nodes.data == expected.data

        nodes = Walks.read(DATADIR / "basic.walk.gz", region="1-2")
        assert nodes.data == expected.data

        nodes = Walks.read(DATADIR / "basic.walk.gz", region="1-")
        assert nodes.data == expected.data

        nodes = Walks.read(DATADIR / "basic.walk.gz", region="-2")
        assert nodes.data == expected.data

        del expected.data[2]

        nodes = Walks.read(DATADIR / "basic.walk", region="1-1")
        assert nodes.data == expected.data

        nodes = Walks.read(DATADIR / "basic.walk.gz", region="1-1")
        assert nodes.data == expected.data
