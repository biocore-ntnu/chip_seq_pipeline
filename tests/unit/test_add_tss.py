import pytest

import pandas as pd

from io import StringIO

from scripts.add_tss import add_tss



@pytest.fixture
def df():

    c = """Chromosome      Start   End     Name    Type    Strand  Transcript      ExonNumber
chrX    0        200        SRPX    gene    -       NM_001170751    -1
chrX    10        20        SRPX    exon    -       NM_001170751    1
chrX    200        400        SRPX    gene    +       NM_001170750    -1
chrX    210        230        SRPX    exon    +       NM_001170750    1"""

    return pd.read_table(StringIO(c), sep="\s+")


@pytest.fixture
def expected_result():

    result_df = """Chromosome Start End Name Type Strand Transcript ExonNumber
chrX 0 200 SRPX gene - NM_001170751 -1
chrX 10 20 SRPX exon - NM_001170751 1
chrX 200 400 SRPX gene + NM_001170750 -1
chrX 210 230 SRPX exon + NM_001170750 1
chrX -4800 5200 SRPX tss + NM_001170750 -1
chrX -4800 5200 SRPX tss - NM_001170751 -1
chrX -4600 5400 SRPX tes + NM_001170750 -1
chrX -5000 5000 SRPX tes - NM_001170751 -1"""

    df = pd.read_table(StringIO(result_df), sep="\s+")

    return df


@pytest.mark.unit
def test_add_tss(df, expected_result):
    result = add_tss(df, 5000)

    print("result", result.to_csv(sep=" "))
    print("expected_result", expected_result.to_csv(sep=" "))

    print("result", result.index)
    print("expected_result", expected_result.index)

    print("result", result.columns)
    print("expected_result", expected_result.columns)

    print("result", result.dtypes)
    print("expected_result", expected_result.dtypes)

    assert result.equals(expected_result)
