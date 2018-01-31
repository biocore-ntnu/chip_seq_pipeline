import pytest

import pandas as pd

from io import StringIO

from scripts.remove_duplicate_tss import remove_duplicate_tss

@pytest.fixture
def df():

    c = """Chromosome      Start   End     Name    Type    Strand  Transcript      ExonNumber
chrX    38149334        38220924        SRPX    gene    -       NM_001170751    -1
chrX    38149334        38149894        SRPX    exon    -       NM_001170751    1
chrX    38154461        38154583        SRPX    exon    -       NM_001170751    2
chrX    38156895        38157029        SRPX    exon    -       NM_001170751    3
chrX    38160016        38160196        SRPX    exon    -       NM_001170751    4
chrX    38160932        38161054        SRPX    exon    -       NM_001170751    5
chrX    38164768        38164895        SRPX    exon    -       NM_001170751    6
chrX    38174159        38174351        SRPX    exon    -       NM_001170751    7
chrX    38178284        38178344        SRPX    exon    -       NM_001170751    8
chrX    38149334        38220924        SRPX    gene    -       NM_001170750    -1
chrX    38149334        38149894        SRPX    exon    -       NM_001170750    1
chrX    38154461        38154583        SRPX    exon    -       NM_001170750    2
chrX    38156895        38157029        SRPX    exon    -       NM_001170750    3
chrX    38160016        38160196        SRPX    exon    -       NM_001170750    4
chrX    38160932        38161054        SRPX    exon    -       NM_001170750    5
chrX    38164768        38164895        SRPX    exon    -       NM_001170750    6
chrX    38171880        38172057        SRPX    exon    -       NM_001170750    7
chrX    38174159        38174351        SRPX    exon    -       NM_001170750    8
chrX    38220695        38220924        SRPX    exon    -       NM_001170750    9"""

    return pd.read_table(StringIO(c), sep="\s+")


@pytest.fixture
def expected_result():

    c = """Chromosome      Start   End     Name    Type    Strand  Transcript      ExonNumber
chrX    38149334        38220924        SRPX    gene    -       NM_001170751    -1
chrX    38149334        38149894        SRPX    exon    -       NM_001170751    1
chrX    38154461        38154583        SRPX    exon    -       NM_001170751    2
chrX    38156895        38157029        SRPX    exon    -       NM_001170751    3
chrX    38160016        38160196        SRPX    exon    -       NM_001170751    4
chrX    38160932        38161054        SRPX    exon    -       NM_001170751    5
chrX    38164768        38164895        SRPX    exon    -       NM_001170751    6
chrX    38174159        38174351        SRPX    exon    -       NM_001170751    7
chrX    38178284        38178344        SRPX    exon    -       NM_001170751    8"""

    return pd.read_table(StringIO(c), sep="\s+")

@pytest.fixture
def df_simple():

    c = """Chromosome      Start   End     Name    Type    Strand  Transcript      ExonNumber
chrX    0        200        SRPX    gene    -       NM_001170751    -1
chrX    10        20        SRPX    exon    -       NM_001170751    1
chrX    200        400        SRPX    gene    +       NM_001170750    -1
chrX    210        230        SRPX    exon    +       NM_001170750    1"""

    return pd.read_table(StringIO(c), sep="\s+")


@pytest.fixture
def expected_result_simple():

    c = """Chromosome      Start   End     Name    Type    Strand  Transcript      ExonNumber
chrX    200        400        SRPX    gene    +       NM_001170750    -1
chrX    210        230        SRPX    exon    +       NM_001170750    1"""

    return pd.read_table(StringIO(c), sep="\s+")


@pytest.mark.unit
def test_remove_overlapping_tss(df, expected_result):

    result = remove_duplicate_tss(df)

    print(result)
    print(expected_result)

    assert result.equals(expected_result)

@pytest.mark.unit
def test_remove_overlapping_tss_simple(df_simple, expected_result_simple):

    result = remove_duplicate_tss(df_simple)

    print(result)
    print(expected_result_simple)
    assert expected_result_simple.equals(result)
