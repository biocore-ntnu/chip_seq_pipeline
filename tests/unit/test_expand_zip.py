
import pytest

import pandas as pd

from io import StringIO


from utils.helpers import expand_zip

@pytest.fixture
def zip_dict():
    return {"g1": [1, 2, 3], "g2": [4, 5, 6]}


@pytest.fixture
def regular_dict():
    return {"prefix": "hoo", "regions": "tss tes".split()}

@pytest.fixture
def template():
    return "{prefix}/data/{regions}/{g1}vs{g2}"


@pytest.fixture
def expected_result():
    return ["hoo/data/tes/1vs4",
            "hoo/data/tes/2vs5",
            "hoo/data/tes/3vs6",
            "hoo/data/tss/1vs4",
            "hoo/data/tss/2vs5",
            "hoo/data/tss/3vs6"]



@pytest.mark.unit
def test_expand_zip(zip_dict, regular_dict, template, expected_result):

    r = expand_zip(template, zip_dict, regular_dict)

    assert set(r) == set(expected_result)
