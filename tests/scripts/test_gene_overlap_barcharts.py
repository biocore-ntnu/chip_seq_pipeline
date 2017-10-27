import pytest

import sys
from io import StringIO

import pandas as pd

if sys.version_info[0] < 3:
    from scripts.find_peak_gene_overlaps import (create_intervaltrees,
                                                 find_peak_gene_overlaps,
                                                 parse_overlap_dataframe,
                                                 create_barchart_data)


@pytest.fixture
def genes(tmpdir):
    contents = """chr1    20   40   exon:ENST00000456328.2:1        exon    +
chr1    220   240   exon:ENST00000450305.2:1        exon    +
chr1    260   280   exon:ENST00000450305.2:2        exon    +
chr1    0   100   ENSG00000223972.5       gene    +
chr1    200   300   ENSG00000243485.5       gene    +
chr1    0    10   ENSG00000223972.5       tss     +
chr1    200   210   ENSG00000243485.5       tss     +
chr1    90   100   ENSG00000223972.5       tes     +
chr1    290   300   ENSG00000243485.5       tes     +"""

    f = tmpdir.join("genes.txt")
    f.write(contents)

    return str(f)

@pytest.fixture
def peaks_string():

    contents = u"""chr1    3   5
chr1    12   14
chr1    200   300
chr1    240   297"""

    return contents


@pytest.fixture
def peaks(peaks_string):

    return pd.read_table(StringIO(peaks_string), header=None, sep="\s+")


@pytest.fixture
def intervaltrees(genes):

    genome = create_intervaltrees(genes)

    return genome

@pytest.fixture
def expected_result_find_overlaps():

    contents = u"""Chromosome  Start  End  Peak Region
0        chr1      3    5     0   gene
1        chr1      3    5     0    tss
2        chr1     12   14     1   gene
3        chr1    200  300     2   gene
4        chr1    200  300     2    tss
5        chr1    200  300     2   exon
6        chr1    200  300     2   exon
7        chr1    200  300     2    tes
8        chr1    240  297     3   gene
9        chr1    240  297     3   exon
10       chr1    240  297     3    tes"""

    return pd.read_table(StringIO(contents), header=0, sep="\s+")


@pytest.mark.py27
def test_find_peak_gene_overlaps(intervaltrees, peaks, expected_result_find_overlaps):

    result = find_peak_gene_overlaps(intervaltrees, peaks)

    assert result.equals(expected_result_find_overlaps)


@pytest.fixture
def expected_result_parse_overlap_dataframe():

    contents = u"""Region Counts
intergenic  1
tes         1
tss         2"""

    return pd.read_table(StringIO(contents), header=0, sep="\s+")


@pytest.mark.py27
def test_parse_overlap_dataframe(expected_result_find_overlaps, expected_result_parse_overlap_dataframe):

    counts = parse_overlap_dataframe(expected_result_find_overlaps)

    print(counts)
    print(expected_result_parse_overlap_dataframe)

    assert counts.equals(expected_result_parse_overlap_dataframe)


@pytest.fixture
def peak_file(tmpdir, peaks_string):

    f = tmpdir.join("peaks.txt")
    f.write(peaks_string)

    return str(f)

@pytest.fixture
def expected_result_create_barchart_data():

    contents = u"""Region  Counts   Label
0  intergenic       1  Sample1
1         tes       1  Sample1
2         tss       1  Sample1"""

    return pd.read_table(StringIO(contents), header=0, sep="\s+")

@pytest.mark.py27
def test_create_barchart_data(genes, peak_file, expected_result_create_barchart_data):

    result = create_barchart_data(genes, peak_file, "Sample1")
    print(result)
    print(expected_result_create_barchart_data)

    assert result.equals(expected_result_create_barchart_data)
