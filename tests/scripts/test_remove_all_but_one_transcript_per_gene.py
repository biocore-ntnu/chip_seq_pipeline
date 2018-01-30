import pytest

import sys
from io import StringIO

import pandas as pd

from scripts.remove_all_but_one_transcript_per_gene import remove_all_but_one_transcript_per_gene

@pytest.mark.unit
def test_remove_all_but_one_transcript_per_gene_ucsc(ucsc_df, ucsc_result):

    df = ucsc_df

    result = remove_all_but_one_transcript_per_gene(df)

    print(result.to_csv(sep=" ",index=False))

    assert result.equals(ucsc_result)



@pytest.fixture
def ucsc_df():

    c = """Chromosome      Start   End     Name    Type    Strand  Transcript      ExonNumber
chrX    38149334        38220924        SRPX    transcript    -       NM_001170751    -1
chrX    38149334        38149894        SRPX    exon    -       NM_001170751    1
chrX    38154461        38154583        SRPX    exon    -       NM_001170751    2
chrX    38156895        38157029        SRPX    exon    -       NM_001170751    3
chrX    38160016        38160196        SRPX    exon    -       NM_001170751    4
chrX    38160932        38161054        SRPX    exon    -       NM_001170751    5
chrX    38164768        38164895        SRPX    exon    -       NM_001170751    6
chrX    38174159        38174351        SRPX    exon    -       NM_001170751    7
chrX    38178284        38178344        SRPX    exon    -       NM_001170751    8
chrX    38149334        38220924        SRPX    transcript    -       NM_001170750    -1
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
def ucsc_result():

    c = """Chromosome Start End Name Type Strand Transcript ExonNumber
chrX 38149334 38220924 SRPX transcript - NM_001170750 -1
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
def gencode_df():
    c = """Chromosome      Start   End     Name    Type    Strand  Transcript      ExonNumber
chr1    11868   12227   ENSG00000223972.5       exon    +       ENST00000456328.2       1
chr1    11868   14409   ENSG00000223972.5       transcript      +       ENST00000456328.2       -1
chr1    12009   12057   ENSG00000223972.5       exon    +       ENST00000450305.2       1
chr1    12009   13670   ENSG00000223972.5       transcript      +       ENST00000450305.2       -1
chr1    12178   12227   ENSG00000223972.5       exon    +       ENST00000450305.2       2
chr1    12612   12697   ENSG00000223972.5       exon    +       ENST00000450305.2       3
chr1    12612   12721   ENSG00000223972.5       exon    +       ENST00000456328.2       2
chr1    12974   13052   ENSG00000223972.5       exon    +       ENST00000450305.2       4
chr1    13220   13374   ENSG00000223972.5       exon    +       ENST00000450305.2       5
chr1    13220   14409   ENSG00000223972.5       exon    +       ENST00000456328.2       3
chr1    13452   13670   ENSG00000223972.5       exon    +       ENST00000450305.2       6"""

    return pd.read_table(StringIO(c), sep="\s+")


@pytest.fixture
def gencode_result():
    c = """Chromosome      Start   End     Name    Type    Strand  Transcript      ExonNumber
chr1 12009 12057 ENSG00000223972.5 exon + ENST00000450305.2 1
chr1 12009 13670 ENSG00000223972.5 transcript + ENST00000450305.2 -1
chr1 12178 12227 ENSG00000223972.5 exon + ENST00000450305.2 2
chr1 12612 12697 ENSG00000223972.5 exon + ENST00000450305.2 3
chr1 12974 13052 ENSG00000223972.5 exon + ENST00000450305.2 4
chr1 13220 13374 ENSG00000223972.5 exon + ENST00000450305.2 5
chr1 13452 13670 ENSG00000223972.5 exon + ENST00000450305.2 6"""

    return pd.read_table(StringIO(c), sep="\s+")

def test_remove_all_but_one_transcript_per_gene_gencode(gencode_df, gencode_result):

    df = gencode_df

    result = remove_all_but_one_transcript_per_gene(df)

    print(result.to_csv(sep=" ",index=False))

    assert result.equals(gencode_result)
