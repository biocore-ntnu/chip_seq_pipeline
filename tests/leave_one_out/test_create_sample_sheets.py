import pytest

from io import StringIO
import pandas as pd

from leave_one_out.create_sample_sheets import create_sample_sheet

@pytest.fixture
def sample_sheet():

    contents = """File	Name	Group	ChIP	Mate
/home/endrebas/genomes/barbara/fastq/1_S1_L001_R1_001.fastq.gz	WT_Input_1	WT	Input	1
/home/endrebas/genomes/barbara/fastq/1_S1_L002_R1_001.fastq.gz	WT_Input_1	WT	Input	1
/home/endrebas/genomes/barbara/fastq/1_S1_L003_R1_001.fastq.gz	WT_Input_1	WT	Input	1
/home/endrebas/genomes/barbara/fastq/1_S1_L004_R1_001.fastq.gz	WT_Input_1	WT	Input	1
/home/endrebas/genomes/barbara/fastq/2_S2_L001_R1_001.fastq.gz	GENE1_KO_Input_1	GENE1_KO	Input	1
/home/endrebas/genomes/barbara/fastq/2_S2_L002_R1_001.fastq.gz	GENE1_KO_Input_1	GENE1_KO	Input	1
/home/endrebas/genomes/barbara/fastq/2_S2_L003_R1_001.fastq.gz	GENE1_KO_Input_1	GENE1_KO	Input	1
/home/endrebas/genomes/barbara/fastq/2_S2_L004_R1_001.fastq.gz	GENE1_KO_Input_1	GENE1_KO	Input	1
/home/endrebas/genomes/barbara/fastq/3_S3_L001_R1_001.fastq.gz	GENE2_KO_Input_1	GENE2_KO	Input	1
/home/endrebas/genomes/barbara/fastq/3_S3_L002_R1_001.fastq.gz	GENE2_KO_Input_1	GENE2_KO	Input	1
/home/endrebas/genomes/barbara/fastq/3_S3_L003_R1_001.fastq.gz	GENE2_KO_Input_1	GENE2_KO	Input	1
/home/endrebas/genomes/barbara/fastq/3_S3_L004_R1_001.fastq.gz	GENE2_KO_Input_1	GENE2_KO	Input	1
/home/endrebas/genomes/barbara/fastq/4_S4_L001_R1_001.fastq.gz	WT_ChIP_1	WT	ChIP	1
/home/endrebas/genomes/barbara/fastq/4_S4_L002_R1_001.fastq.gz	WT_ChIP_1	WT	ChIP	1
/home/endrebas/genomes/barbara/fastq/4_S4_L003_R1_001.fastq.gz	WT_ChIP_1	WT	ChIP	1
/home/endrebas/genomes/barbara/fastq/4_S4_L004_R1_001.fastq.gz	WT_ChIP_1	WT	ChIP	1
/home/endrebas/genomes/barbara/fastq/5_S5_L001_R1_001.fastq.gz	GENE1_KO_ChIP_1	GENE1_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/5_S5_L002_R1_001.fastq.gz	GENE1_KO_ChIP_1	GENE1_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/5_S5_L003_R1_001.fastq.gz	GENE1_KO_ChIP_1	GENE1_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/5_S5_L004_R1_001.fastq.gz	GENE1_KO_ChIP_1	GENE1_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/6_S6_L001_R1_001.fastq.gz	GENE2_KO_ChIP_1	GENE2_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/6_S6_L002_R1_001.fastq.gz	GENE2_KO_ChIP_1	GENE2_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/6_S6_L003_R1_001.fastq.gz	GENE2_KO_ChIP_1	GENE2_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/6_S6_L004_R1_001.fastq.gz	GENE2_KO_ChIP_1	GENE2_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/7_S7_L001_R1_001.fastq.gz	WT_Input_2	WT	Input	1
/home/endrebas/genomes/barbara/fastq/7_S7_L002_R1_001.fastq.gz	WT_Input_2	WT	Input	1
/home/endrebas/genomes/barbara/fastq/7_S7_L003_R1_001.fastq.gz	WT_Input_2	WT	Input	1
/home/endrebas/genomes/barbara/fastq/7_S7_L004_R1_001.fastq.gz	WT_Input_2	WT	Input	1
/home/endrebas/genomes/barbara/fastq/8_S8_L001_R1_001.fastq.gz	GENE1_KO_Input_2	GENE1_KO	Input	1
/home/endrebas/genomes/barbara/fastq/8_S8_L002_R1_001.fastq.gz	GENE1_KO_Input_2	GENE1_KO	Input	1
/home/endrebas/genomes/barbara/fastq/8_S8_L003_R1_001.fastq.gz	GENE1_KO_Input_2	GENE1_KO	Input	1
/home/endrebas/genomes/barbara/fastq/8_S8_L004_R1_001.fastq.gz	GENE1_KO_Input_2	GENE1_KO	Input	1
/home/endrebas/genomes/barbara/fastq/9_S9_L001_R1_001.fastq.gz	GENE2_KO_Input_2	GENE2_KO	Input	1
/home/endrebas/genomes/barbara/fastq/9_S9_L002_R1_001.fastq.gz	GENE2_KO_Input_2	GENE2_KO	Input	1
/home/endrebas/genomes/barbara/fastq/9_S9_L003_R1_001.fastq.gz	GENE2_KO_Input_2	GENE2_KO	Input	1
/home/endrebas/genomes/barbara/fastq/9_S9_L004_R1_001.fastq.gz	GENE2_KO_Input_2	GENE2_KO	Input	1
/home/endrebas/genomes/barbara/fastq/10_S10_L001_R1_001.fastq.gz	WT_ChIP_2	WT	ChIP	1
/home/endrebas/genomes/barbara/fastq/10_S10_L002_R1_001.fastq.gz	WT_ChIP_2	WT	ChIP	1
/home/endrebas/genomes/barbara/fastq/10_S10_L003_R1_001.fastq.gz	WT_ChIP_2	WT	ChIP	1
/home/endrebas/genomes/barbara/fastq/10_S10_L004_R1_001.fastq.gz	WT_ChIP_2	WT	ChIP	1
/home/endrebas/genomes/barbara/fastq/11_S11_L001_R1_001.fastq.gz	GENE1_KO_ChIP_2	GENE1_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/11_S11_L002_R1_001.fastq.gz	GENE1_KO_ChIP_2	GENE1_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/11_S11_L003_R1_001.fastq.gz	GENE1_KO_ChIP_2	GENE1_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/11_S11_L004_R1_001.fastq.gz	GENE1_KO_ChIP_2	GENE1_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/12_S12_L001_R1_001.fastq.gz	GENE2_KO_ChIP_2	GENE2_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/12_S12_L002_R1_001.fastq.gz	GENE2_KO_ChIP_2	GENE2_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/12_S12_L003_R1_001.fastq.gz	GENE2_KO_ChIP_2	GENE2_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/12_S12_L004_R1_001.fastq.gz	GENE2_KO_ChIP_2	GENE2_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/13_S13_L001_R1_001.fastq.gz	WT_Input_3	WT	Input	1
/home/endrebas/genomes/barbara/fastq/13_S13_L002_R1_001.fastq.gz	WT_Input_3	WT	Input	1
/home/endrebas/genomes/barbara/fastq/13_S13_L003_R1_001.fastq.gz	WT_Input_3	WT	Input	1
/home/endrebas/genomes/barbara/fastq/13_S13_L004_R1_001.fastq.gz	WT_Input_3	WT	Input	1
/home/endrebas/genomes/barbara/fastq/14_S14_L001_R1_001.fastq.gz	GENE1_KO_Input_3	GENE1_KO	Input	1
/home/endrebas/genomes/barbara/fastq/14_S14_L002_R1_001.fastq.gz	GENE1_KO_Input_3	GENE1_KO	Input	1
/home/endrebas/genomes/barbara/fastq/14_S14_L003_R1_001.fastq.gz	GENE1_KO_Input_3	GENE1_KO	Input	1
/home/endrebas/genomes/barbara/fastq/14_S14_L004_R1_001.fastq.gz	GENE1_KO_Input_3	GENE1_KO	Input	1
/home/endrebas/genomes/barbara/fastq/15_S15_L001_R1_001.fastq.gz	GENE2_KO_Input_3	GENE2_KO	Input	1
/home/endrebas/genomes/barbara/fastq/15_S15_L002_R1_001.fastq.gz	GENE2_KO_Input_3	GENE2_KO	Input	1
/home/endrebas/genomes/barbara/fastq/15_S15_L003_R1_001.fastq.gz	GENE2_KO_Input_3	GENE2_KO	Input	1
/home/endrebas/genomes/barbara/fastq/15_S15_L004_R1_001.fastq.gz	GENE2_KO_Input_3	GENE2_KO	Input	1
/home/endrebas/genomes/barbara/fastq/16_S16_L001_R1_001.fastq.gz	WT_ChIP_3	WT	ChIP	1
/home/endrebas/genomes/barbara/fastq/16_S16_L002_R1_001.fastq.gz	WT_ChIP_3	WT	ChIP	1
/home/endrebas/genomes/barbara/fastq/16_S16_L003_R1_001.fastq.gz	WT_ChIP_3	WT	ChIP	1
/home/endrebas/genomes/barbara/fastq/16_S16_L004_R1_001.fastq.gz	WT_ChIP_3	WT	ChIP	1
/home/endrebas/genomes/barbara/fastq/17_S17_L001_R1_001.fastq.gz	GENE1_KO_ChIP_3	GENE1_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/17_S17_L002_R1_001.fastq.gz	GENE1_KO_ChIP_3	GENE1_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/17_S17_L003_R1_001.fastq.gz	GENE1_KO_ChIP_3	GENE1_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/17_S17_L004_R1_001.fastq.gz	GENE1_KO_ChIP_3	GENE1_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/18_S18_L001_R1_001.fastq.gz	GENE2_KO_ChIP_3	GENE2_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/18_S18_L002_R1_001.fastq.gz	GENE2_KO_ChIP_3	GENE2_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/18_S18_L003_R1_001.fastq.gz	GENE2_KO_ChIP_3	GENE2_KO	ChIP	1
/home/endrebas/genomes/barbara/fastq/18_S18_L004_R1_001.fastq.gz	GENE2_KO_ChIP_3	GENE2_KO	ChIP	1"""

    return pd.read_table(StringIO(contents), sep="\s+", header=0, index_col=None)


@pytest.fixture
def expected_result():
    contents = """File Name Group ChIP Mate OriginalGroup
0 /home/endrebas/genomes/barbara/fastq/2_S2_L001_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_1_lo Input 1 GENE1_KO
1 /home/endrebas/genomes/barbara/fastq/2_S2_L002_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_1_lo Input 1 GENE1_KO
2 /home/endrebas/genomes/barbara/fastq/2_S2_L003_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_1_lo Input 1 GENE1_KO
3 /home/endrebas/genomes/barbara/fastq/2_S2_L004_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_1_lo Input 1 GENE1_KO
4 /home/endrebas/genomes/barbara/fastq/8_S8_L001_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_1_lo Input 1 GENE1_KO
5 /home/endrebas/genomes/barbara/fastq/8_S8_L002_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_1_lo Input 1 GENE1_KO
6 /home/endrebas/genomes/barbara/fastq/8_S8_L003_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_1_lo Input 1 GENE1_KO
7 /home/endrebas/genomes/barbara/fastq/8_S8_L004_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_1_lo Input 1 GENE1_KO
8 /home/endrebas/genomes/barbara/fastq/11_S11_L001_R1_001.fastq.gz GENE1_KO_ChIP_2 GENE1_KO_ChIP_1_lo ChIP 1 GENE1_KO
9 /home/endrebas/genomes/barbara/fastq/11_S11_L002_R1_001.fastq.gz GENE1_KO_ChIP_2 GENE1_KO_ChIP_1_lo ChIP 1 GENE1_KO
10 /home/endrebas/genomes/barbara/fastq/11_S11_L003_R1_001.fastq.gz GENE1_KO_ChIP_2 GENE1_KO_ChIP_1_lo ChIP 1 GENE1_KO
11 /home/endrebas/genomes/barbara/fastq/11_S11_L004_R1_001.fastq.gz GENE1_KO_ChIP_2 GENE1_KO_ChIP_1_lo ChIP 1 GENE1_KO
12 /home/endrebas/genomes/barbara/fastq/14_S14_L001_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_1_lo Input 1 GENE1_KO
13 /home/endrebas/genomes/barbara/fastq/14_S14_L002_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_1_lo Input 1 GENE1_KO
14 /home/endrebas/genomes/barbara/fastq/14_S14_L003_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_1_lo Input 1 GENE1_KO
15 /home/endrebas/genomes/barbara/fastq/14_S14_L004_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_1_lo Input 1 GENE1_KO
16 /home/endrebas/genomes/barbara/fastq/17_S17_L001_R1_001.fastq.gz GENE1_KO_ChIP_3 GENE1_KO_ChIP_1_lo ChIP 1 GENE1_KO
17 /home/endrebas/genomes/barbara/fastq/17_S17_L002_R1_001.fastq.gz GENE1_KO_ChIP_3 GENE1_KO_ChIP_1_lo ChIP 1 GENE1_KO
18 /home/endrebas/genomes/barbara/fastq/17_S17_L003_R1_001.fastq.gz GENE1_KO_ChIP_3 GENE1_KO_ChIP_1_lo ChIP 1 GENE1_KO
19 /home/endrebas/genomes/barbara/fastq/17_S17_L004_R1_001.fastq.gz GENE1_KO_ChIP_3 GENE1_KO_ChIP_1_lo ChIP 1 GENE1_KO
20 /home/endrebas/genomes/barbara/fastq/2_S2_L001_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_2_lo Input 1 GENE1_KO
21 /home/endrebas/genomes/barbara/fastq/2_S2_L002_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_2_lo Input 1 GENE1_KO
22 /home/endrebas/genomes/barbara/fastq/2_S2_L003_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_2_lo Input 1 GENE1_KO
23 /home/endrebas/genomes/barbara/fastq/2_S2_L004_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_2_lo Input 1 GENE1_KO
24 /home/endrebas/genomes/barbara/fastq/5_S5_L001_R1_001.fastq.gz GENE1_KO_ChIP_1 GENE1_KO_ChIP_2_lo ChIP 1 GENE1_KO
25 /home/endrebas/genomes/barbara/fastq/5_S5_L002_R1_001.fastq.gz GENE1_KO_ChIP_1 GENE1_KO_ChIP_2_lo ChIP 1 GENE1_KO
26 /home/endrebas/genomes/barbara/fastq/5_S5_L003_R1_001.fastq.gz GENE1_KO_ChIP_1 GENE1_KO_ChIP_2_lo ChIP 1 GENE1_KO
27 /home/endrebas/genomes/barbara/fastq/5_S5_L004_R1_001.fastq.gz GENE1_KO_ChIP_1 GENE1_KO_ChIP_2_lo ChIP 1 GENE1_KO
28 /home/endrebas/genomes/barbara/fastq/8_S8_L001_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_2_lo Input 1 GENE1_KO
29 /home/endrebas/genomes/barbara/fastq/8_S8_L002_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_2_lo Input 1 GENE1_KO
30 /home/endrebas/genomes/barbara/fastq/8_S8_L003_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_2_lo Input 1 GENE1_KO
31 /home/endrebas/genomes/barbara/fastq/8_S8_L004_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_2_lo Input 1 GENE1_KO
32 /home/endrebas/genomes/barbara/fastq/14_S14_L001_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_2_lo Input 1 GENE1_KO
33 /home/endrebas/genomes/barbara/fastq/14_S14_L002_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_2_lo Input 1 GENE1_KO
34 /home/endrebas/genomes/barbara/fastq/14_S14_L003_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_2_lo Input 1 GENE1_KO
35 /home/endrebas/genomes/barbara/fastq/14_S14_L004_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_2_lo Input 1 GENE1_KO
36 /home/endrebas/genomes/barbara/fastq/17_S17_L001_R1_001.fastq.gz GENE1_KO_ChIP_3 GENE1_KO_ChIP_2_lo ChIP 1 GENE1_KO
37 /home/endrebas/genomes/barbara/fastq/17_S17_L002_R1_001.fastq.gz GENE1_KO_ChIP_3 GENE1_KO_ChIP_2_lo ChIP 1 GENE1_KO
38 /home/endrebas/genomes/barbara/fastq/17_S17_L003_R1_001.fastq.gz GENE1_KO_ChIP_3 GENE1_KO_ChIP_2_lo ChIP 1 GENE1_KO
39 /home/endrebas/genomes/barbara/fastq/17_S17_L004_R1_001.fastq.gz GENE1_KO_ChIP_3 GENE1_KO_ChIP_2_lo ChIP 1 GENE1_KO
40 /home/endrebas/genomes/barbara/fastq/2_S2_L001_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_3_lo Input 1 GENE1_KO
41 /home/endrebas/genomes/barbara/fastq/2_S2_L002_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_3_lo Input 1 GENE1_KO
42 /home/endrebas/genomes/barbara/fastq/2_S2_L003_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_3_lo Input 1 GENE1_KO
43 /home/endrebas/genomes/barbara/fastq/2_S2_L004_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_3_lo Input 1 GENE1_KO
44 /home/endrebas/genomes/barbara/fastq/5_S5_L001_R1_001.fastq.gz GENE1_KO_ChIP_1 GENE1_KO_ChIP_3_lo ChIP 1 GENE1_KO
45 /home/endrebas/genomes/barbara/fastq/5_S5_L002_R1_001.fastq.gz GENE1_KO_ChIP_1 GENE1_KO_ChIP_3_lo ChIP 1 GENE1_KO
46 /home/endrebas/genomes/barbara/fastq/5_S5_L003_R1_001.fastq.gz GENE1_KO_ChIP_1 GENE1_KO_ChIP_3_lo ChIP 1 GENE1_KO
47 /home/endrebas/genomes/barbara/fastq/5_S5_L004_R1_001.fastq.gz GENE1_KO_ChIP_1 GENE1_KO_ChIP_3_lo ChIP 1 GENE1_KO
48 /home/endrebas/genomes/barbara/fastq/8_S8_L001_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_3_lo Input 1 GENE1_KO
49 /home/endrebas/genomes/barbara/fastq/8_S8_L002_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_3_lo Input 1 GENE1_KO
50 /home/endrebas/genomes/barbara/fastq/8_S8_L003_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_3_lo Input 1 GENE1_KO
51 /home/endrebas/genomes/barbara/fastq/8_S8_L004_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_3_lo Input 1 GENE1_KO
52 /home/endrebas/genomes/barbara/fastq/11_S11_L001_R1_001.fastq.gz GENE1_KO_ChIP_2 GENE1_KO_ChIP_3_lo ChIP 1 GENE1_KO
53 /home/endrebas/genomes/barbara/fastq/11_S11_L002_R1_001.fastq.gz GENE1_KO_ChIP_2 GENE1_KO_ChIP_3_lo ChIP 1 GENE1_KO
54 /home/endrebas/genomes/barbara/fastq/11_S11_L003_R1_001.fastq.gz GENE1_KO_ChIP_2 GENE1_KO_ChIP_3_lo ChIP 1 GENE1_KO
55 /home/endrebas/genomes/barbara/fastq/11_S11_L004_R1_001.fastq.gz GENE1_KO_ChIP_2 GENE1_KO_ChIP_3_lo ChIP 1 GENE1_KO
56 /home/endrebas/genomes/barbara/fastq/14_S14_L001_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_3_lo Input 1 GENE1_KO
57 /home/endrebas/genomes/barbara/fastq/14_S14_L002_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_3_lo Input 1 GENE1_KO
58 /home/endrebas/genomes/barbara/fastq/14_S14_L003_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_3_lo Input 1 GENE1_KO
59 /home/endrebas/genomes/barbara/fastq/14_S14_L004_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_3_lo Input 1 GENE1_KO
60 /home/endrebas/genomes/barbara/fastq/3_S3_L001_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_1_lo Input 1 GENE2_KO
61 /home/endrebas/genomes/barbara/fastq/3_S3_L002_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_1_lo Input 1 GENE2_KO
62 /home/endrebas/genomes/barbara/fastq/3_S3_L003_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_1_lo Input 1 GENE2_KO
63 /home/endrebas/genomes/barbara/fastq/3_S3_L004_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_1_lo Input 1 GENE2_KO
64 /home/endrebas/genomes/barbara/fastq/9_S9_L001_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_1_lo Input 1 GENE2_KO
65 /home/endrebas/genomes/barbara/fastq/9_S9_L002_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_1_lo Input 1 GENE2_KO
66 /home/endrebas/genomes/barbara/fastq/9_S9_L003_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_1_lo Input 1 GENE2_KO
67 /home/endrebas/genomes/barbara/fastq/9_S9_L004_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_1_lo Input 1 GENE2_KO
68 /home/endrebas/genomes/barbara/fastq/12_S12_L001_R1_001.fastq.gz GENE2_KO_ChIP_2 GENE2_KO_ChIP_1_lo ChIP 1 GENE2_KO
69 /home/endrebas/genomes/barbara/fastq/12_S12_L002_R1_001.fastq.gz GENE2_KO_ChIP_2 GENE2_KO_ChIP_1_lo ChIP 1 GENE2_KO
70 /home/endrebas/genomes/barbara/fastq/12_S12_L003_R1_001.fastq.gz GENE2_KO_ChIP_2 GENE2_KO_ChIP_1_lo ChIP 1 GENE2_KO
71 /home/endrebas/genomes/barbara/fastq/12_S12_L004_R1_001.fastq.gz GENE2_KO_ChIP_2 GENE2_KO_ChIP_1_lo ChIP 1 GENE2_KO
72 /home/endrebas/genomes/barbara/fastq/15_S15_L001_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_1_lo Input 1 GENE2_KO
73 /home/endrebas/genomes/barbara/fastq/15_S15_L002_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_1_lo Input 1 GENE2_KO
74 /home/endrebas/genomes/barbara/fastq/15_S15_L003_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_1_lo Input 1 GENE2_KO
75 /home/endrebas/genomes/barbara/fastq/15_S15_L004_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_1_lo Input 1 GENE2_KO
76 /home/endrebas/genomes/barbara/fastq/18_S18_L001_R1_001.fastq.gz GENE2_KO_ChIP_3 GENE2_KO_ChIP_1_lo ChIP 1 GENE2_KO
77 /home/endrebas/genomes/barbara/fastq/18_S18_L002_R1_001.fastq.gz GENE2_KO_ChIP_3 GENE2_KO_ChIP_1_lo ChIP 1 GENE2_KO
78 /home/endrebas/genomes/barbara/fastq/18_S18_L003_R1_001.fastq.gz GENE2_KO_ChIP_3 GENE2_KO_ChIP_1_lo ChIP 1 GENE2_KO
79 /home/endrebas/genomes/barbara/fastq/18_S18_L004_R1_001.fastq.gz GENE2_KO_ChIP_3 GENE2_KO_ChIP_1_lo ChIP 1 GENE2_KO
80 /home/endrebas/genomes/barbara/fastq/3_S3_L001_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_2_lo Input 1 GENE2_KO
81 /home/endrebas/genomes/barbara/fastq/3_S3_L002_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_2_lo Input 1 GENE2_KO
82 /home/endrebas/genomes/barbara/fastq/3_S3_L003_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_2_lo Input 1 GENE2_KO
83 /home/endrebas/genomes/barbara/fastq/3_S3_L004_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_2_lo Input 1 GENE2_KO
84 /home/endrebas/genomes/barbara/fastq/6_S6_L001_R1_001.fastq.gz GENE2_KO_ChIP_1 GENE2_KO_ChIP_2_lo ChIP 1 GENE2_KO
85 /home/endrebas/genomes/barbara/fastq/6_S6_L002_R1_001.fastq.gz GENE2_KO_ChIP_1 GENE2_KO_ChIP_2_lo ChIP 1 GENE2_KO
86 /home/endrebas/genomes/barbara/fastq/6_S6_L003_R1_001.fastq.gz GENE2_KO_ChIP_1 GENE2_KO_ChIP_2_lo ChIP 1 GENE2_KO
87 /home/endrebas/genomes/barbara/fastq/6_S6_L004_R1_001.fastq.gz GENE2_KO_ChIP_1 GENE2_KO_ChIP_2_lo ChIP 1 GENE2_KO
88 /home/endrebas/genomes/barbara/fastq/9_S9_L001_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_2_lo Input 1 GENE2_KO
89 /home/endrebas/genomes/barbara/fastq/9_S9_L002_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_2_lo Input 1 GENE2_KO
90 /home/endrebas/genomes/barbara/fastq/9_S9_L003_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_2_lo Input 1 GENE2_KO
91 /home/endrebas/genomes/barbara/fastq/9_S9_L004_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_2_lo Input 1 GENE2_KO
92 /home/endrebas/genomes/barbara/fastq/15_S15_L001_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_2_lo Input 1 GENE2_KO
93 /home/endrebas/genomes/barbara/fastq/15_S15_L002_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_2_lo Input 1 GENE2_KO
94 /home/endrebas/genomes/barbara/fastq/15_S15_L003_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_2_lo Input 1 GENE2_KO
95 /home/endrebas/genomes/barbara/fastq/15_S15_L004_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_2_lo Input 1 GENE2_KO
96 /home/endrebas/genomes/barbara/fastq/18_S18_L001_R1_001.fastq.gz GENE2_KO_ChIP_3 GENE2_KO_ChIP_2_lo ChIP 1 GENE2_KO
97 /home/endrebas/genomes/barbara/fastq/18_S18_L002_R1_001.fastq.gz GENE2_KO_ChIP_3 GENE2_KO_ChIP_2_lo ChIP 1 GENE2_KO
98 /home/endrebas/genomes/barbara/fastq/18_S18_L003_R1_001.fastq.gz GENE2_KO_ChIP_3 GENE2_KO_ChIP_2_lo ChIP 1 GENE2_KO
99 /home/endrebas/genomes/barbara/fastq/18_S18_L004_R1_001.fastq.gz GENE2_KO_ChIP_3 GENE2_KO_ChIP_2_lo ChIP 1 GENE2_KO
100 /home/endrebas/genomes/barbara/fastq/3_S3_L001_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_3_lo Input 1 GENE2_KO
101 /home/endrebas/genomes/barbara/fastq/3_S3_L002_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_3_lo Input 1 GENE2_KO
102 /home/endrebas/genomes/barbara/fastq/3_S3_L003_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_3_lo Input 1 GENE2_KO
103 /home/endrebas/genomes/barbara/fastq/3_S3_L004_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_3_lo Input 1 GENE2_KO
104 /home/endrebas/genomes/barbara/fastq/6_S6_L001_R1_001.fastq.gz GENE2_KO_ChIP_1 GENE2_KO_ChIP_3_lo ChIP 1 GENE2_KO
105 /home/endrebas/genomes/barbara/fastq/6_S6_L002_R1_001.fastq.gz GENE2_KO_ChIP_1 GENE2_KO_ChIP_3_lo ChIP 1 GENE2_KO
106 /home/endrebas/genomes/barbara/fastq/6_S6_L003_R1_001.fastq.gz GENE2_KO_ChIP_1 GENE2_KO_ChIP_3_lo ChIP 1 GENE2_KO
107 /home/endrebas/genomes/barbara/fastq/6_S6_L004_R1_001.fastq.gz GENE2_KO_ChIP_1 GENE2_KO_ChIP_3_lo ChIP 1 GENE2_KO
108 /home/endrebas/genomes/barbara/fastq/9_S9_L001_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_3_lo Input 1 GENE2_KO
109 /home/endrebas/genomes/barbara/fastq/9_S9_L002_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_3_lo Input 1 GENE2_KO
110 /home/endrebas/genomes/barbara/fastq/9_S9_L003_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_3_lo Input 1 GENE2_KO
111 /home/endrebas/genomes/barbara/fastq/9_S9_L004_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_3_lo Input 1 GENE2_KO
112 /home/endrebas/genomes/barbara/fastq/12_S12_L001_R1_001.fastq.gz GENE2_KO_ChIP_2 GENE2_KO_ChIP_3_lo ChIP 1 GENE2_KO
113 /home/endrebas/genomes/barbara/fastq/12_S12_L002_R1_001.fastq.gz GENE2_KO_ChIP_2 GENE2_KO_ChIP_3_lo ChIP 1 GENE2_KO
114 /home/endrebas/genomes/barbara/fastq/12_S12_L003_R1_001.fastq.gz GENE2_KO_ChIP_2 GENE2_KO_ChIP_3_lo ChIP 1 GENE2_KO
115 /home/endrebas/genomes/barbara/fastq/12_S12_L004_R1_001.fastq.gz GENE2_KO_ChIP_2 GENE2_KO_ChIP_3_lo ChIP 1 GENE2_KO
116 /home/endrebas/genomes/barbara/fastq/15_S15_L001_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_3_lo Input 1 GENE2_KO
117 /home/endrebas/genomes/barbara/fastq/15_S15_L002_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_3_lo Input 1 GENE2_KO
118 /home/endrebas/genomes/barbara/fastq/15_S15_L003_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_3_lo Input 1 GENE2_KO
119 /home/endrebas/genomes/barbara/fastq/15_S15_L004_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_3_lo Input 1 GENE2_KO
120 /home/endrebas/genomes/barbara/fastq/1_S1_L001_R1_001.fastq.gz WT_Input_1 WT_ChIP_1_lo Input 1 WT
121 /home/endrebas/genomes/barbara/fastq/1_S1_L002_R1_001.fastq.gz WT_Input_1 WT_ChIP_1_lo Input 1 WT
122 /home/endrebas/genomes/barbara/fastq/1_S1_L003_R1_001.fastq.gz WT_Input_1 WT_ChIP_1_lo Input 1 WT
123 /home/endrebas/genomes/barbara/fastq/1_S1_L004_R1_001.fastq.gz WT_Input_1 WT_ChIP_1_lo Input 1 WT
124 /home/endrebas/genomes/barbara/fastq/7_S7_L001_R1_001.fastq.gz WT_Input_2 WT_ChIP_1_lo Input 1 WT
125 /home/endrebas/genomes/barbara/fastq/7_S7_L002_R1_001.fastq.gz WT_Input_2 WT_ChIP_1_lo Input 1 WT
126 /home/endrebas/genomes/barbara/fastq/7_S7_L003_R1_001.fastq.gz WT_Input_2 WT_ChIP_1_lo Input 1 WT
127 /home/endrebas/genomes/barbara/fastq/7_S7_L004_R1_001.fastq.gz WT_Input_2 WT_ChIP_1_lo Input 1 WT
128 /home/endrebas/genomes/barbara/fastq/10_S10_L001_R1_001.fastq.gz WT_ChIP_2 WT_ChIP_1_lo ChIP 1 WT
129 /home/endrebas/genomes/barbara/fastq/10_S10_L002_R1_001.fastq.gz WT_ChIP_2 WT_ChIP_1_lo ChIP 1 WT
130 /home/endrebas/genomes/barbara/fastq/10_S10_L003_R1_001.fastq.gz WT_ChIP_2 WT_ChIP_1_lo ChIP 1 WT
131 /home/endrebas/genomes/barbara/fastq/10_S10_L004_R1_001.fastq.gz WT_ChIP_2 WT_ChIP_1_lo ChIP 1 WT
132 /home/endrebas/genomes/barbara/fastq/13_S13_L001_R1_001.fastq.gz WT_Input_3 WT_ChIP_1_lo Input 1 WT
133 /home/endrebas/genomes/barbara/fastq/13_S13_L002_R1_001.fastq.gz WT_Input_3 WT_ChIP_1_lo Input 1 WT
134 /home/endrebas/genomes/barbara/fastq/13_S13_L003_R1_001.fastq.gz WT_Input_3 WT_ChIP_1_lo Input 1 WT
135 /home/endrebas/genomes/barbara/fastq/13_S13_L004_R1_001.fastq.gz WT_Input_3 WT_ChIP_1_lo Input 1 WT
136 /home/endrebas/genomes/barbara/fastq/16_S16_L001_R1_001.fastq.gz WT_ChIP_3 WT_ChIP_1_lo ChIP 1 WT
137 /home/endrebas/genomes/barbara/fastq/16_S16_L002_R1_001.fastq.gz WT_ChIP_3 WT_ChIP_1_lo ChIP 1 WT
138 /home/endrebas/genomes/barbara/fastq/16_S16_L003_R1_001.fastq.gz WT_ChIP_3 WT_ChIP_1_lo ChIP 1 WT
139 /home/endrebas/genomes/barbara/fastq/16_S16_L004_R1_001.fastq.gz WT_ChIP_3 WT_ChIP_1_lo ChIP 1 WT
140 /home/endrebas/genomes/barbara/fastq/1_S1_L001_R1_001.fastq.gz WT_Input_1 WT_ChIP_2_lo Input 1 WT
141 /home/endrebas/genomes/barbara/fastq/1_S1_L002_R1_001.fastq.gz WT_Input_1 WT_ChIP_2_lo Input 1 WT
142 /home/endrebas/genomes/barbara/fastq/1_S1_L003_R1_001.fastq.gz WT_Input_1 WT_ChIP_2_lo Input 1 WT
143 /home/endrebas/genomes/barbara/fastq/1_S1_L004_R1_001.fastq.gz WT_Input_1 WT_ChIP_2_lo Input 1 WT
144 /home/endrebas/genomes/barbara/fastq/4_S4_L001_R1_001.fastq.gz WT_ChIP_1 WT_ChIP_2_lo ChIP 1 WT
145 /home/endrebas/genomes/barbara/fastq/4_S4_L002_R1_001.fastq.gz WT_ChIP_1 WT_ChIP_2_lo ChIP 1 WT
146 /home/endrebas/genomes/barbara/fastq/4_S4_L003_R1_001.fastq.gz WT_ChIP_1 WT_ChIP_2_lo ChIP 1 WT
147 /home/endrebas/genomes/barbara/fastq/4_S4_L004_R1_001.fastq.gz WT_ChIP_1 WT_ChIP_2_lo ChIP 1 WT
148 /home/endrebas/genomes/barbara/fastq/7_S7_L001_R1_001.fastq.gz WT_Input_2 WT_ChIP_2_lo Input 1 WT
149 /home/endrebas/genomes/barbara/fastq/7_S7_L002_R1_001.fastq.gz WT_Input_2 WT_ChIP_2_lo Input 1 WT
150 /home/endrebas/genomes/barbara/fastq/7_S7_L003_R1_001.fastq.gz WT_Input_2 WT_ChIP_2_lo Input 1 WT
151 /home/endrebas/genomes/barbara/fastq/7_S7_L004_R1_001.fastq.gz WT_Input_2 WT_ChIP_2_lo Input 1 WT
152 /home/endrebas/genomes/barbara/fastq/13_S13_L001_R1_001.fastq.gz WT_Input_3 WT_ChIP_2_lo Input 1 WT
153 /home/endrebas/genomes/barbara/fastq/13_S13_L002_R1_001.fastq.gz WT_Input_3 WT_ChIP_2_lo Input 1 WT
154 /home/endrebas/genomes/barbara/fastq/13_S13_L003_R1_001.fastq.gz WT_Input_3 WT_ChIP_2_lo Input 1 WT
155 /home/endrebas/genomes/barbara/fastq/13_S13_L004_R1_001.fastq.gz WT_Input_3 WT_ChIP_2_lo Input 1 WT
156 /home/endrebas/genomes/barbara/fastq/16_S16_L001_R1_001.fastq.gz WT_ChIP_3 WT_ChIP_2_lo ChIP 1 WT
157 /home/endrebas/genomes/barbara/fastq/16_S16_L002_R1_001.fastq.gz WT_ChIP_3 WT_ChIP_2_lo ChIP 1 WT
158 /home/endrebas/genomes/barbara/fastq/16_S16_L003_R1_001.fastq.gz WT_ChIP_3 WT_ChIP_2_lo ChIP 1 WT
159 /home/endrebas/genomes/barbara/fastq/16_S16_L004_R1_001.fastq.gz WT_ChIP_3 WT_ChIP_2_lo ChIP 1 WT
160 /home/endrebas/genomes/barbara/fastq/1_S1_L001_R1_001.fastq.gz WT_Input_1 WT_ChIP_3_lo Input 1 WT
161 /home/endrebas/genomes/barbara/fastq/1_S1_L002_R1_001.fastq.gz WT_Input_1 WT_ChIP_3_lo Input 1 WT
162 /home/endrebas/genomes/barbara/fastq/1_S1_L003_R1_001.fastq.gz WT_Input_1 WT_ChIP_3_lo Input 1 WT
163 /home/endrebas/genomes/barbara/fastq/1_S1_L004_R1_001.fastq.gz WT_Input_1 WT_ChIP_3_lo Input 1 WT
164 /home/endrebas/genomes/barbara/fastq/4_S4_L001_R1_001.fastq.gz WT_ChIP_1 WT_ChIP_3_lo ChIP 1 WT
165 /home/endrebas/genomes/barbara/fastq/4_S4_L002_R1_001.fastq.gz WT_ChIP_1 WT_ChIP_3_lo ChIP 1 WT
166 /home/endrebas/genomes/barbara/fastq/4_S4_L003_R1_001.fastq.gz WT_ChIP_1 WT_ChIP_3_lo ChIP 1 WT
167 /home/endrebas/genomes/barbara/fastq/4_S4_L004_R1_001.fastq.gz WT_ChIP_1 WT_ChIP_3_lo ChIP 1 WT
168 /home/endrebas/genomes/barbara/fastq/7_S7_L001_R1_001.fastq.gz WT_Input_2 WT_ChIP_3_lo Input 1 WT
169 /home/endrebas/genomes/barbara/fastq/7_S7_L002_R1_001.fastq.gz WT_Input_2 WT_ChIP_3_lo Input 1 WT
170 /home/endrebas/genomes/barbara/fastq/7_S7_L003_R1_001.fastq.gz WT_Input_2 WT_ChIP_3_lo Input 1 WT
171 /home/endrebas/genomes/barbara/fastq/7_S7_L004_R1_001.fastq.gz WT_Input_2 WT_ChIP_3_lo Input 1 WT
172 /home/endrebas/genomes/barbara/fastq/10_S10_L001_R1_001.fastq.gz WT_ChIP_2 WT_ChIP_3_lo ChIP 1 WT
173 /home/endrebas/genomes/barbara/fastq/10_S10_L002_R1_001.fastq.gz WT_ChIP_2 WT_ChIP_3_lo ChIP 1 WT
174 /home/endrebas/genomes/barbara/fastq/10_S10_L003_R1_001.fastq.gz WT_ChIP_2 WT_ChIP_3_lo ChIP 1 WT
175 /home/endrebas/genomes/barbara/fastq/10_S10_L004_R1_001.fastq.gz WT_ChIP_2 WT_ChIP_3_lo ChIP 1 WT
176 /home/endrebas/genomes/barbara/fastq/13_S13_L001_R1_001.fastq.gz WT_Input_3 WT_ChIP_3_lo Input 1 WT
177 /home/endrebas/genomes/barbara/fastq/13_S13_L002_R1_001.fastq.gz WT_Input_3 WT_ChIP_3_lo Input 1 WT
178 /home/endrebas/genomes/barbara/fastq/13_S13_L003_R1_001.fastq.gz WT_Input_3 WT_ChIP_3_lo Input 1 WT
179 /home/endrebas/genomes/barbara/fastq/13_S13_L004_R1_001.fastq.gz WT_Input_3 WT_ChIP_3_lo Input 1 WT"""

    return pd.read_table(StringIO(contents), sep=" ", index_col=0)


def test_create_sample_sheets(sample_sheet, expected_result):

    df = create_sample_sheet(sample_sheet)

    df.to_csv("corr_res", sep=" ")

    assert df.equals(expected_result)
