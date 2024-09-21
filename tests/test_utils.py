from panct.utils import *
import pytest

def test_parse_region():
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