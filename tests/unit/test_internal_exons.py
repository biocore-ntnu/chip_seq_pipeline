import pytest

import pandas as pd

from io import StringIO

from scripts.compute_internal_exons import compute_internal_exons



@pytest.fixture
def df():

    c = """Chromosome Start End Name Type Strand Transcript ExonNumber
chrX 38149334 38220924 SRPX gene - NM_001170750 -1
chrX 38149334 38149894 SRPX exon - NM_001170750 1
chrX 38154461 38154583 SRPX exon - NM_001170750 2
chrX 38156895 38157029 SRPX exon - NM_001170750 3
chrX 38160016 38160196 SRPX exon - NM_001170750 4
chrX 38160932 38161054 SRPX exon - NM_001170750 5
chrX 38164768 38164895 SRPX exon - NM_001170750 6
chrX 38171880 38172057 SRPX exon - NM_001170750 7
chrX 38174159 38174351 SRPX exon - NM_001170750 8
chrX 38220695 38220924 SRPX exon - NM_001170750 9"""

    return pd.read_table(StringIO(c), sep="\s+")


@pytest.fixture
def expected_result():

    result_df = """Chromosome     Start       End  Name  Type Strand    Transcript  ExonNumber
chrX  38154461  38154583  SRPX  exon      -  NM_001170750           2
chrX  38156895  38157029  SRPX  exon      -  NM_001170750           3
chrX  38160016  38160196  SRPX  exon      -  NM_001170750           4
chrX  38160932  38161054  SRPX  exon      -  NM_001170750           5
chrX  38164768  38164895  SRPX  exon      -  NM_001170750           6
chrX  38171880  38172057  SRPX  exon      -  NM_001170750           7
chrX  38174159  38174351  SRPX  exon      -  NM_001170750           8"""

    df = pd.read_table(StringIO(result_df), sep="\s+")

    return df


@pytest.mark.unit
def test_internal_exons(df, expected_result):
    result = compute_internal_exons(df)

    # print("indata", exons.to_csv(sep=" "))

    print("result", result.to_csv(sep=" "))
    print("expected_result", expected_result.to_csv(sep=" "))

    print("result", result.index)
    print("expected_result", expected_result.index)

    print("result", result.columns)
    print("expected_result", expected_result.columns)


    print("result", result.dtypes)
    print("expected_result", expected_result.dtypes)

    assert result.equals(expected_result)
