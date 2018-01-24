import pytest

import pandas as pd

from io import StringIO

from scripts.compute_internal_exons import compute_internal_exons

@pytest.fixture
def exons():

    exons_df = """Chromosome Start End Name Score Strand Gene
3 chr21 9907189 9909277 exon:ENST00000400754.3:3 . - ENSG00000188681.7
4 chr21 9907193 9908432 exon:ENST00000416067.1:4 . - ENSG00000188681.7
5 chr21 9907246 9907462 exon:ENST00000559466.1:6 . - ENSG00000188681.7
6 chr21 9908277 9908432 exon:ENST00000559466.1:5 . - ENSG00000188681.7
7 chr21 9909046 9909277 exon:ENST00000416067.1:3 . - ENSG00000188681.7
8 chr21 9909046 9909277 exon:ENST00000559466.1:4 . - ENSG00000188681.7
9 chr21 9909960 9910081 exon:ENST00000559466.1:3 . - ENSG00000188681.7
10 chr21 9910441 9910515 exon:ENST00000559466.1:2 . - ENSG00000188681.7
11 chr21 9911927 9912442 exon:ENST00000559466.1:1 . - ENSG00000188681.7
14 chr21 9966321 9966380 exon:ENST00000400754.3:2 . - ENSG00000188681.7
15 chr21 9966321 9966380 exon:ENST00000416067.1:2 . - ENSG00000188681.7
16 chr21 9968515 9968548 exon:ENST00000400754.3:1 . - ENSG00000188681.7
17 chr21 9968515 9968585 exon:ENST00000416067.1:1 . - ENSG00000188681.7"""

    df = pd.read_table(StringIO(exons_df), sep=" ")

    return df

@pytest.fixture
def expected_result():
    result_df = """Chromosome Start End Name Score Strand Gene
0 chr21 9910441 9910515 exon:ENST00000559466.1:2 . - ENSG00000188681.7
1 chr21 9966321 9966380 exon:ENST00000416067.1:2 . - ENSG00000188681.7
2 chr21 9909046 9909277 exon:ENST00000416067.1:3 . - ENSG00000188681.7
3 chr21 9909960 9910081 exon:ENST00000559466.1:3 . - ENSG00000188681.7
4 chr21 9907193 9908432 exon:ENST00000416067.1:4 . - ENSG00000188681.7
5 chr21 9909046 9909277 exon:ENST00000559466.1:4 . - ENSG00000188681.7
6 chr21 9908277 9908432 exon:ENST00000559466.1:5 . - ENSG00000188681.7
7 chr21 9966321 9966380 exon:ENST00000400754.3:2 . - ENSG00000188681.7"""


    df = pd.read_table(StringIO(result_df), sep=" ")

    return df

# old with transcript info
# result_df = """Chromosome	Start	End	Name	Score	Strand	Gene	Transcript
# chr1	12178	12227	exon:ENST00000450305.2:2	.	+	ENST00000450305	2
# chr1	12612	12697	exon:ENST00000450305.2:3	.	+	ENST00000450305	2
# chr1	12974	13052	exon:ENST00000450305.2:4	.	+	ENST00000450305	2
# chr1	13220	13374	exon:ENST00000450305.2:5	.	+	ENST00000450305	2
# chr1	12612	12721	exon:ENST00000456328.2:2	.	+	ENST00000456328	2"""


@pytest.mark.unit
def test_internal_exons(exons, expected_result):
    result = compute_internal_exons(exons)

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
