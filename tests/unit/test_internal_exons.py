import pytest

import pandas as pd

from io import StringIO

from scripts.compute_internal_exons import compute_internal_exons

@pytest.fixture
def exons():
    exons_df = """chr1	11868	12227	exon:ENST00000456328.2:1	.	+
chr1	12009	12057	exon:ENST00000450305.2:1	.	+
chr1	12178	12227	exon:ENST00000450305.2:2	.	+
chr1	12612	12697	exon:ENST00000450305.2:3	.	+
chr1	12612	12721	exon:ENST00000456328.2:2	.	+
chr1	12974	13052	exon:ENST00000450305.2:4	.	+
chr1	13220	13374	exon:ENST00000450305.2:5	.	+
chr1	13220	14409	exon:ENST00000456328.2:3	.	+
chr1	13452	13670	exon:ENST00000450305.2:6	.	+
chr1	14403	14501	exon:ENST00000488147.1:11	.	-"""

    df = pd.read_table(StringIO(exons_df), header=None, names="Chromosome Start End Name Score Strand".split())

    return df

@pytest.fixture
def expected_result():
    result_df = """Chromosome	Start	End	Name	Score	Strand	Gene	Transcript
chr1	12178	12227	exon:ENST00000450305.2:2	.	+	ENST00000450305	2
chr1	12612	12697	exon:ENST00000450305.2:3	.	+	ENST00000450305	2
chr1	12974	13052	exon:ENST00000450305.2:4	.	+	ENST00000450305	2
chr1	13220	13374	exon:ENST00000450305.2:5	.	+	ENST00000450305	2
chr1	12612	12721	exon:ENST00000456328.2:2	.	+	ENST00000456328	2"""

    df = pd.read_table(StringIO(result_df))

    return df


@pytest.mark.unit
def test_internal_exons(exons, expected_result):
    result = compute_internal_exons(exons)

    print("result", result.to_csv(sep="\t"))
    print("expected_result", expected_result.to_csv(sep="\t"))

    print("result", result.index)
    print("expected_result", expected_result.index)

    print("result", result.columns)
    print("expected_result", expected_result.columns)


    print("result", result.dtypes)
    print("expected_result", expected_result.dtypes)

    assert result.equals(expected_result)
