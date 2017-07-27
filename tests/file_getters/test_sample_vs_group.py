import pytest

import pandas as pd
from io import StringIO

from utils.file_getters import sample_vs_group

@pytest.mark.current
@pytest.mark.unit
def test_sample_vs_group(sample_sheet, expected_result):

    result = sample_vs_group(sample_sheet, "ChIP")

    print()
    print(result)
    print(expected_result)

    assert result.equals(expected_result)


@pytest.fixture
def expected_result():
    contents = """Sample OtherGroup
0    GENE1_KO_ChIP_1    GENE2_KO
1    GENE1_KO_ChIP_1         WT
2    GENE1_KO_ChIP_2    GENE2_KO
3    GENE1_KO_ChIP_2         WT
4    GENE1_KO_ChIP_3    GENE2_KO
5    GENE1_KO_ChIP_3         WT
6   GENE2_KO_ChIP_1     GENE1_KO
7   GENE2_KO_ChIP_1         WT
8   GENE2_KO_ChIP_2     GENE1_KO
9   GENE2_KO_ChIP_2         WT
10  GENE2_KO_ChIP_3     GENE1_KO
11  GENE2_KO_ChIP_3         WT
12       WT_ChIP_1     GENE1_KO
13       WT_ChIP_1    GENE2_KO
14       WT_ChIP_2     GENE1_KO
15       WT_ChIP_2    GENE2_KO
16       WT_ChIP_3     GENE1_KO
17       WT_ChIP_3    GENE2_KO"""

    return pd.read_table(StringIO(contents), sep="\s+", index_col=0)


@pytest.fixture
def sample_sheet():

    contents = """File    Name    Group   ChIP    Mate
/home/endrebas/genomes/barbara/fastq/1_S1_L001_R1_001.fastq.gz  WT_Input_1      WT      Input   1
/home/endrebas/genomes/barbara/fastq/1_S1_L002_R1_001.fastq.gz  WT_Input_1      WT      Input   1
/home/endrebas/genomes/barbara/fastq/1_S1_L003_R1_001.fastq.gz  WT_Input_1      WT      Input   1
/home/endrebas/genomes/barbara/fastq/1_S1_L004_R1_001.fastq.gz  WT_Input_1      WT      Input   1
/home/endrebas/genomes/barbara/fastq/2_S2_L001_R1_001.fastq.gz  GENE1_KO_Input_1  GENE1_KO  Input   1
/home/endrebas/genomes/barbara/fastq/2_S2_L002_R1_001.fastq.gz  GENE1_KO_Input_1  GENE1_KO  Input   1
/home/endrebas/genomes/barbara/fastq/2_S2_L003_R1_001.fastq.gz  GENE1_KO_Input_1  GENE1_KO  Input   1
/home/endrebas/genomes/barbara/fastq/2_S2_L004_R1_001.fastq.gz  GENE1_KO_Input_1  GENE1_KO  Input   1
/home/endrebas/genomes/barbara/fastq/3_S3_L001_R1_001.fastq.gz  GENE2_KO_Input_1 GENE2_KO Input   1
/home/endrebas/genomes/barbara/fastq/3_S3_L002_R1_001.fastq.gz  GENE2_KO_Input_1 GENE2_KO Input   1
/home/endrebas/genomes/barbara/fastq/3_S3_L003_R1_001.fastq.gz  GENE2_KO_Input_1 GENE2_KO Input   1
/home/endrebas/genomes/barbara/fastq/3_S3_L004_R1_001.fastq.gz  GENE2_KO_Input_1 GENE2_KO Input   1
/home/endrebas/genomes/barbara/fastq/4_S4_L001_R1_001.fastq.gz  WT_ChIP_1       WT      ChIP    1
/home/endrebas/genomes/barbara/fastq/4_S4_L002_R1_001.fastq.gz  WT_ChIP_1       WT      ChIP    1
/home/endrebas/genomes/barbara/fastq/4_S4_L003_R1_001.fastq.gz  WT_ChIP_1       WT      ChIP    1
/home/endrebas/genomes/barbara/fastq/4_S4_L004_R1_001.fastq.gz  WT_ChIP_1       WT      ChIP    1
/home/endrebas/genomes/barbara/fastq/5_S5_L001_R1_001.fastq.gz  GENE1_KO_ChIP_1   GENE1_KO  ChIP    1
/home/endrebas/genomes/barbara/fastq/5_S5_L002_R1_001.fastq.gz  GENE1_KO_ChIP_1   GENE1_KO  ChIP    1
/home/endrebas/genomes/barbara/fastq/5_S5_L003_R1_001.fastq.gz  GENE1_KO_ChIP_1   GENE1_KO  ChIP    1
/home/endrebas/genomes/barbara/fastq/5_S5_L004_R1_001.fastq.gz  GENE1_KO_ChIP_1   GENE1_KO  ChIP    1
/home/endrebas/genomes/barbara/fastq/6_S6_L001_R1_001.fastq.gz  GENE2_KO_ChIP_1  GENE2_KO ChIP    1
/home/endrebas/genomes/barbara/fastq/6_S6_L002_R1_001.fastq.gz  GENE2_KO_ChIP_1  GENE2_KO ChIP    1
/home/endrebas/genomes/barbara/fastq/6_S6_L003_R1_001.fastq.gz  GENE2_KO_ChIP_1  GENE2_KO ChIP    1
/home/endrebas/genomes/barbara/fastq/6_S6_L004_R1_001.fastq.gz  GENE2_KO_ChIP_1  GENE2_KO ChIP    1
/home/endrebas/genomes/barbara/fastq/7_S7_L001_R1_001.fastq.gz  WT_Input_2      WT      Input   1
/home/endrebas/genomes/barbara/fastq/7_S7_L002_R1_001.fastq.gz  WT_Input_2      WT      Input   1
/home/endrebas/genomes/barbara/fastq/7_S7_L003_R1_001.fastq.gz  WT_Input_2      WT      Input   1
/home/endrebas/genomes/barbara/fastq/7_S7_L004_R1_001.fastq.gz  WT_Input_2      WT      Input   1
/home/endrebas/genomes/barbara/fastq/8_S8_L001_R1_001.fastq.gz  GENE1_KO_Input_2  GENE1_KO  Input   1
/home/endrebas/genomes/barbara/fastq/8_S8_L002_R1_001.fastq.gz  GENE1_KO_Input_2  GENE1_KO  Input   1
/home/endrebas/genomes/barbara/fastq/8_S8_L003_R1_001.fastq.gz  GENE1_KO_Input_2  GENE1_KO  Input   1
/home/endrebas/genomes/barbara/fastq/8_S8_L004_R1_001.fastq.gz  GENE1_KO_Input_2  GENE1_KO  Input   1
/home/endrebas/genomes/barbara/fastq/9_S9_L001_R1_001.fastq.gz  GENE2_KO_Input_2 GENE2_KO Input   1
/home/endrebas/genomes/barbara/fastq/9_S9_L002_R1_001.fastq.gz  GENE2_KO_Input_2 GENE2_KO Input   1
/home/endrebas/genomes/barbara/fastq/9_S9_L003_R1_001.fastq.gz  GENE2_KO_Input_2 GENE2_KO Input   1
/home/endrebas/genomes/barbara/fastq/9_S9_L004_R1_001.fastq.gz  GENE2_KO_Input_2 GENE2_KO Input   1
/home/endrebas/genomes/barbara/fastq/10_S10_L001_R1_001.fastq.gz        WT_ChIP_2       WT      ChIP    1
/home/endrebas/genomes/barbara/fastq/10_S10_L002_R1_001.fastq.gz        WT_ChIP_2       WT      ChIP    1
/home/endrebas/genomes/barbara/fastq/10_S10_L003_R1_001.fastq.gz        WT_ChIP_2       WT      ChIP    1
/home/endrebas/genomes/barbara/fastq/10_S10_L004_R1_001.fastq.gz        WT_ChIP_2       WT      ChIP    1
/home/endrebas/genomes/barbara/fastq/11_S11_L001_R1_001.fastq.gz        GENE1_KO_ChIP_2   GENE1_KO  ChIP    1
/home/endrebas/genomes/barbara/fastq/11_S11_L002_R1_001.fastq.gz        GENE1_KO_ChIP_2   GENE1_KO  ChIP    1
/home/endrebas/genomes/barbara/fastq/11_S11_L003_R1_001.fastq.gz        GENE1_KO_ChIP_2   GENE1_KO  ChIP    1
/home/endrebas/genomes/barbara/fastq/11_S11_L004_R1_001.fastq.gz        GENE1_KO_ChIP_2   GENE1_KO  ChIP    1
/home/endrebas/genomes/barbara/fastq/12_S12_L001_R1_001.fastq.gz        GENE2_KO_ChIP_2  GENE2_KO ChIP    1
/home/endrebas/genomes/barbara/fastq/12_S12_L002_R1_001.fastq.gz        GENE2_KO_ChIP_2  GENE2_KO ChIP    1
/home/endrebas/genomes/barbara/fastq/12_S12_L003_R1_001.fastq.gz        GENE2_KO_ChIP_2  GENE2_KO ChIP    1
/home/endrebas/genomes/barbara/fastq/12_S12_L004_R1_001.fastq.gz        GENE2_KO_ChIP_2  GENE2_KO ChIP    1
/home/endrebas/genomes/barbara/fastq/13_S13_L001_R1_001.fastq.gz        WT_Input_3      WT      Input   1
/home/endrebas/genomes/barbara/fastq/13_S13_L002_R1_001.fastq.gz        WT_Input_3      WT      Input   1
/home/endrebas/genomes/barbara/fastq/13_S13_L003_R1_001.fastq.gz        WT_Input_3      WT      Input   1
/home/endrebas/genomes/barbara/fastq/13_S13_L004_R1_001.fastq.gz        WT_Input_3      WT      Input   1
/home/endrebas/genomes/barbara/fastq/14_S14_L001_R1_001.fastq.gz        GENE1_KO_Input_3  GENE1_KO  Input   1
/home/endrebas/genomes/barbara/fastq/14_S14_L002_R1_001.fastq.gz        GENE1_KO_Input_3  GENE1_KO  Input   1
/home/endrebas/genomes/barbara/fastq/14_S14_L003_R1_001.fastq.gz        GENE1_KO_Input_3  GENE1_KO  Input   1
/home/endrebas/genomes/barbara/fastq/14_S14_L004_R1_001.fastq.gz        GENE1_KO_Input_3  GENE1_KO  Input   1
/home/endrebas/genomes/barbara/fastq/15_S15_L001_R1_001.fastq.gz        GENE2_KO_Input_3 GENE2_KO Input   1
/home/endrebas/genomes/barbara/fastq/15_S15_L002_R1_001.fastq.gz        GENE2_KO_Input_3 GENE2_KO Input   1
/home/endrebas/genomes/barbara/fastq/15_S15_L003_R1_001.fastq.gz        GENE2_KO_Input_3 GENE2_KO Input   1
/home/endrebas/genomes/barbara/fastq/15_S15_L004_R1_001.fastq.gz        GENE2_KO_Input_3 GENE2_KO Input   1
/home/endrebas/genomes/barbara/fastq/16_S16_L001_R1_001.fastq.gz        WT_ChIP_3       WT      ChIP    1
/home/endrebas/genomes/barbara/fastq/16_S16_L002_R1_001.fastq.gz        WT_ChIP_3       WT      ChIP    1
/home/endrebas/genomes/barbara/fastq/16_S16_L003_R1_001.fastq.gz        WT_ChIP_3       WT      ChIP    1
/home/endrebas/genomes/barbara/fastq/16_S16_L004_R1_001.fastq.gz        WT_ChIP_3       WT      ChIP    1
/home/endrebas/genomes/barbara/fastq/17_S17_L001_R1_001.fastq.gz        GENE1_KO_ChIP_3   GENE1_KO  ChIP    1
/home/endrebas/genomes/barbara/fastq/17_S17_L002_R1_001.fastq.gz        GENE1_KO_ChIP_3   GENE1_KO  ChIP    1
/home/endrebas/genomes/barbara/fastq/17_S17_L003_R1_001.fastq.gz        GENE1_KO_ChIP_3   GENE1_KO  ChIP    1
/home/endrebas/genomes/barbara/fastq/17_S17_L004_R1_001.fastq.gz        GENE1_KO_ChIP_3   GENE1_KO  ChIP    1
/home/endrebas/genomes/barbara/fastq/18_S18_L001_R1_001.fastq.gz        GENE2_KO_ChIP_3  GENE2_KO ChIP    1
/home/endrebas/genomes/barbara/fastq/18_S18_L002_R1_001.fastq.gz        GENE2_KO_ChIP_3  GENE2_KO ChIP    1
/home/endrebas/genomes/barbara/fastq/18_S18_L003_R1_001.fastq.gz        GENE2_KO_ChIP_3  GENE2_KO ChIP    1
/home/endrebas/genomes/barbara/fastq/18_S18_L004_R1_001.fastq.gz        GENE2_KO_ChIP_3  GENE2_KO ChIP    1"""

    df = pd.read_table(StringIO(contents), sep="\s+")

    return df
