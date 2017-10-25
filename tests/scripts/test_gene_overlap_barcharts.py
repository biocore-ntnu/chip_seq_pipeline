import pytest

import sys
from io import StringIO

import pandas as pd

if sys.version_info[0] > 3:
    from scripts.find_peak_gene_overlaps import create_intervaltree


@pytest.fixture
def genes(tmpdir):
    contents = """chr1    11868   12000   exon:ENST00000456328.2:1        exon    +
chr1    12009   12057   exon:ENST00000450305.2:1        exon    +
chr1    12178   12227   exon:ENST00000450305.2:2        exon    +
chr1    11868   14409   ENSG00000223972.5       gene    +
chr1    29553   31109   ENSG00000243485.5       gene    +
chr1    8868    14868   ENSG00000223972.5       tss     +
chr1    26553   32553   ENSG00000243485.5       tss     +
chr1    11409   17409   ENSG00000223972.5       tes     +
chr1    28109   34109   ENSG00000243485.5       tes     +"""

    f = tmpdir.join("genes.txt")
    f.write(contents)

    return str(f)


@pytest.fixture
def peaks():

    contents = u"""chr1    12009   12057
chr1    12178   12227
chr1    11868   14409"""

    return pd.read_table(StringIO(contents), header=None, sep="\s+")


@pytest.fixture
def intervaltree(genes):

    genome = create_intervaltree(genes)

    return genome

@pytest.mark.py27
def test_gene_overlap_barcharts(intervaltree, peaks):

    find_gene_peak_overlaps()

    assert 0
