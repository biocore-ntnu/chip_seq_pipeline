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
    contents = """File Name Group ChIP Mate
0 /home/endrebas/genomes/barbara/fastq/2_S2_L001_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_1_lo Input 1
1 /home/endrebas/genomes/barbara/fastq/2_S2_L002_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_1_lo Input 1
2 /home/endrebas/genomes/barbara/fastq/2_S2_L003_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_1_lo Input 1
3 /home/endrebas/genomes/barbara/fastq/2_S2_L004_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_1_lo Input 1
4 /home/endrebas/genomes/barbara/fastq/8_S8_L001_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_1_lo Input 1
5 /home/endrebas/genomes/barbara/fastq/8_S8_L002_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_1_lo Input 1
6 /home/endrebas/genomes/barbara/fastq/8_S8_L003_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_1_lo Input 1
7 /home/endrebas/genomes/barbara/fastq/8_S8_L004_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_1_lo Input 1
8 /home/endrebas/genomes/barbara/fastq/11_S11_L001_R1_001.fastq.gz GENE1_KO_ChIP_2 GENE1_KO_ChIP_1_lo ChIP 1
9 /home/endrebas/genomes/barbara/fastq/11_S11_L002_R1_001.fastq.gz GENE1_KO_ChIP_2 GENE1_KO_ChIP_1_lo ChIP 1
10 /home/endrebas/genomes/barbara/fastq/11_S11_L003_R1_001.fastq.gz GENE1_KO_ChIP_2 GENE1_KO_ChIP_1_lo ChIP 1
11 /home/endrebas/genomes/barbara/fastq/11_S11_L004_R1_001.fastq.gz GENE1_KO_ChIP_2 GENE1_KO_ChIP_1_lo ChIP 1
12 /home/endrebas/genomes/barbara/fastq/14_S14_L001_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_1_lo Input 1
13 /home/endrebas/genomes/barbara/fastq/14_S14_L002_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_1_lo Input 1
14 /home/endrebas/genomes/barbara/fastq/14_S14_L003_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_1_lo Input 1
15 /home/endrebas/genomes/barbara/fastq/14_S14_L004_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_1_lo Input 1
16 /home/endrebas/genomes/barbara/fastq/17_S17_L001_R1_001.fastq.gz GENE1_KO_ChIP_3 GENE1_KO_ChIP_1_lo ChIP 1
17 /home/endrebas/genomes/barbara/fastq/17_S17_L002_R1_001.fastq.gz GENE1_KO_ChIP_3 GENE1_KO_ChIP_1_lo ChIP 1
18 /home/endrebas/genomes/barbara/fastq/17_S17_L003_R1_001.fastq.gz GENE1_KO_ChIP_3 GENE1_KO_ChIP_1_lo ChIP 1
19 /home/endrebas/genomes/barbara/fastq/17_S17_L004_R1_001.fastq.gz GENE1_KO_ChIP_3 GENE1_KO_ChIP_1_lo ChIP 1
20 /home/endrebas/genomes/barbara/fastq/2_S2_L001_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_2_lo Input 1
21 /home/endrebas/genomes/barbara/fastq/2_S2_L002_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_2_lo Input 1
22 /home/endrebas/genomes/barbara/fastq/2_S2_L003_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_2_lo Input 1
23 /home/endrebas/genomes/barbara/fastq/2_S2_L004_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_2_lo Input 1
24 /home/endrebas/genomes/barbara/fastq/5_S5_L001_R1_001.fastq.gz GENE1_KO_ChIP_1 GENE1_KO_ChIP_2_lo ChIP 1
25 /home/endrebas/genomes/barbara/fastq/5_S5_L002_R1_001.fastq.gz GENE1_KO_ChIP_1 GENE1_KO_ChIP_2_lo ChIP 1
26 /home/endrebas/genomes/barbara/fastq/5_S5_L003_R1_001.fastq.gz GENE1_KO_ChIP_1 GENE1_KO_ChIP_2_lo ChIP 1
27 /home/endrebas/genomes/barbara/fastq/5_S5_L004_R1_001.fastq.gz GENE1_KO_ChIP_1 GENE1_KO_ChIP_2_lo ChIP 1
28 /home/endrebas/genomes/barbara/fastq/8_S8_L001_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_2_lo Input 1
29 /home/endrebas/genomes/barbara/fastq/8_S8_L002_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_2_lo Input 1
30 /home/endrebas/genomes/barbara/fastq/8_S8_L003_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_2_lo Input 1
31 /home/endrebas/genomes/barbara/fastq/8_S8_L004_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_2_lo Input 1
32 /home/endrebas/genomes/barbara/fastq/14_S14_L001_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_2_lo Input 1
33 /home/endrebas/genomes/barbara/fastq/14_S14_L002_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_2_lo Input 1
34 /home/endrebas/genomes/barbara/fastq/14_S14_L003_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_2_lo Input 1
35 /home/endrebas/genomes/barbara/fastq/14_S14_L004_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_2_lo Input 1
36 /home/endrebas/genomes/barbara/fastq/17_S17_L001_R1_001.fastq.gz GENE1_KO_ChIP_3 GENE1_KO_ChIP_2_lo ChIP 1
37 /home/endrebas/genomes/barbara/fastq/17_S17_L002_R1_001.fastq.gz GENE1_KO_ChIP_3 GENE1_KO_ChIP_2_lo ChIP 1
38 /home/endrebas/genomes/barbara/fastq/17_S17_L003_R1_001.fastq.gz GENE1_KO_ChIP_3 GENE1_KO_ChIP_2_lo ChIP 1
39 /home/endrebas/genomes/barbara/fastq/17_S17_L004_R1_001.fastq.gz GENE1_KO_ChIP_3 GENE1_KO_ChIP_2_lo ChIP 1
40 /home/endrebas/genomes/barbara/fastq/2_S2_L001_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_3_lo Input 1
41 /home/endrebas/genomes/barbara/fastq/2_S2_L002_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_3_lo Input 1
42 /home/endrebas/genomes/barbara/fastq/2_S2_L003_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_3_lo Input 1
43 /home/endrebas/genomes/barbara/fastq/2_S2_L004_R1_001.fastq.gz GENE1_KO_Input_1 GENE1_KO_ChIP_3_lo Input 1
44 /home/endrebas/genomes/barbara/fastq/5_S5_L001_R1_001.fastq.gz GENE1_KO_ChIP_1 GENE1_KO_ChIP_3_lo ChIP 1
45 /home/endrebas/genomes/barbara/fastq/5_S5_L002_R1_001.fastq.gz GENE1_KO_ChIP_1 GENE1_KO_ChIP_3_lo ChIP 1
46 /home/endrebas/genomes/barbara/fastq/5_S5_L003_R1_001.fastq.gz GENE1_KO_ChIP_1 GENE1_KO_ChIP_3_lo ChIP 1
47 /home/endrebas/genomes/barbara/fastq/5_S5_L004_R1_001.fastq.gz GENE1_KO_ChIP_1 GENE1_KO_ChIP_3_lo ChIP 1
48 /home/endrebas/genomes/barbara/fastq/8_S8_L001_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_3_lo Input 1
49 /home/endrebas/genomes/barbara/fastq/8_S8_L002_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_3_lo Input 1
50 /home/endrebas/genomes/barbara/fastq/8_S8_L003_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_3_lo Input 1
51 /home/endrebas/genomes/barbara/fastq/8_S8_L004_R1_001.fastq.gz GENE1_KO_Input_2 GENE1_KO_ChIP_3_lo Input 1
52 /home/endrebas/genomes/barbara/fastq/11_S11_L001_R1_001.fastq.gz GENE1_KO_ChIP_2 GENE1_KO_ChIP_3_lo ChIP 1
53 /home/endrebas/genomes/barbara/fastq/11_S11_L002_R1_001.fastq.gz GENE1_KO_ChIP_2 GENE1_KO_ChIP_3_lo ChIP 1
54 /home/endrebas/genomes/barbara/fastq/11_S11_L003_R1_001.fastq.gz GENE1_KO_ChIP_2 GENE1_KO_ChIP_3_lo ChIP 1
55 /home/endrebas/genomes/barbara/fastq/11_S11_L004_R1_001.fastq.gz GENE1_KO_ChIP_2 GENE1_KO_ChIP_3_lo ChIP 1
56 /home/endrebas/genomes/barbara/fastq/14_S14_L001_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_3_lo Input 1
57 /home/endrebas/genomes/barbara/fastq/14_S14_L002_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_3_lo Input 1
58 /home/endrebas/genomes/barbara/fastq/14_S14_L003_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_3_lo Input 1
59 /home/endrebas/genomes/barbara/fastq/14_S14_L004_R1_001.fastq.gz GENE1_KO_Input_3 GENE1_KO_ChIP_3_lo Input 1
60 /home/endrebas/genomes/barbara/fastq/3_S3_L001_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_1_lo Input 1
61 /home/endrebas/genomes/barbara/fastq/3_S3_L002_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_1_lo Input 1
62 /home/endrebas/genomes/barbara/fastq/3_S3_L003_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_1_lo Input 1
63 /home/endrebas/genomes/barbara/fastq/3_S3_L004_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_1_lo Input 1
64 /home/endrebas/genomes/barbara/fastq/9_S9_L001_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_1_lo Input 1
65 /home/endrebas/genomes/barbara/fastq/9_S9_L002_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_1_lo Input 1
66 /home/endrebas/genomes/barbara/fastq/9_S9_L003_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_1_lo Input 1
67 /home/endrebas/genomes/barbara/fastq/9_S9_L004_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_1_lo Input 1
68 /home/endrebas/genomes/barbara/fastq/12_S12_L001_R1_001.fastq.gz GENE2_KO_ChIP_2 GENE2_KO_ChIP_1_lo ChIP 1
69 /home/endrebas/genomes/barbara/fastq/12_S12_L002_R1_001.fastq.gz GENE2_KO_ChIP_2 GENE2_KO_ChIP_1_lo ChIP 1
70 /home/endrebas/genomes/barbara/fastq/12_S12_L003_R1_001.fastq.gz GENE2_KO_ChIP_2 GENE2_KO_ChIP_1_lo ChIP 1
71 /home/endrebas/genomes/barbara/fastq/12_S12_L004_R1_001.fastq.gz GENE2_KO_ChIP_2 GENE2_KO_ChIP_1_lo ChIP 1
72 /home/endrebas/genomes/barbara/fastq/15_S15_L001_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_1_lo Input 1
73 /home/endrebas/genomes/barbara/fastq/15_S15_L002_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_1_lo Input 1
74 /home/endrebas/genomes/barbara/fastq/15_S15_L003_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_1_lo Input 1
75 /home/endrebas/genomes/barbara/fastq/15_S15_L004_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_1_lo Input 1
76 /home/endrebas/genomes/barbara/fastq/18_S18_L001_R1_001.fastq.gz GENE2_KO_ChIP_3 GENE2_KO_ChIP_1_lo ChIP 1
77 /home/endrebas/genomes/barbara/fastq/18_S18_L002_R1_001.fastq.gz GENE2_KO_ChIP_3 GENE2_KO_ChIP_1_lo ChIP 1
78 /home/endrebas/genomes/barbara/fastq/18_S18_L003_R1_001.fastq.gz GENE2_KO_ChIP_3 GENE2_KO_ChIP_1_lo ChIP 1
79 /home/endrebas/genomes/barbara/fastq/18_S18_L004_R1_001.fastq.gz GENE2_KO_ChIP_3 GENE2_KO_ChIP_1_lo ChIP 1
80 /home/endrebas/genomes/barbara/fastq/3_S3_L001_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_2_lo Input 1
81 /home/endrebas/genomes/barbara/fastq/3_S3_L002_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_2_lo Input 1
82 /home/endrebas/genomes/barbara/fastq/3_S3_L003_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_2_lo Input 1
83 /home/endrebas/genomes/barbara/fastq/3_S3_L004_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_2_lo Input 1
84 /home/endrebas/genomes/barbara/fastq/6_S6_L001_R1_001.fastq.gz GENE2_KO_ChIP_1 GENE2_KO_ChIP_2_lo ChIP 1
85 /home/endrebas/genomes/barbara/fastq/6_S6_L002_R1_001.fastq.gz GENE2_KO_ChIP_1 GENE2_KO_ChIP_2_lo ChIP 1
86 /home/endrebas/genomes/barbara/fastq/6_S6_L003_R1_001.fastq.gz GENE2_KO_ChIP_1 GENE2_KO_ChIP_2_lo ChIP 1
87 /home/endrebas/genomes/barbara/fastq/6_S6_L004_R1_001.fastq.gz GENE2_KO_ChIP_1 GENE2_KO_ChIP_2_lo ChIP 1
88 /home/endrebas/genomes/barbara/fastq/9_S9_L001_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_2_lo Input 1
89 /home/endrebas/genomes/barbara/fastq/9_S9_L002_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_2_lo Input 1
90 /home/endrebas/genomes/barbara/fastq/9_S9_L003_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_2_lo Input 1
91 /home/endrebas/genomes/barbara/fastq/9_S9_L004_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_2_lo Input 1
92 /home/endrebas/genomes/barbara/fastq/15_S15_L001_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_2_lo Input 1
93 /home/endrebas/genomes/barbara/fastq/15_S15_L002_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_2_lo Input 1
94 /home/endrebas/genomes/barbara/fastq/15_S15_L003_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_2_lo Input 1
95 /home/endrebas/genomes/barbara/fastq/15_S15_L004_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_2_lo Input 1
96 /home/endrebas/genomes/barbara/fastq/18_S18_L001_R1_001.fastq.gz GENE2_KO_ChIP_3 GENE2_KO_ChIP_2_lo ChIP 1
97 /home/endrebas/genomes/barbara/fastq/18_S18_L002_R1_001.fastq.gz GENE2_KO_ChIP_3 GENE2_KO_ChIP_2_lo ChIP 1
98 /home/endrebas/genomes/barbara/fastq/18_S18_L003_R1_001.fastq.gz GENE2_KO_ChIP_3 GENE2_KO_ChIP_2_lo ChIP 1
99 /home/endrebas/genomes/barbara/fastq/18_S18_L004_R1_001.fastq.gz GENE2_KO_ChIP_3 GENE2_KO_ChIP_2_lo ChIP 1
100 /home/endrebas/genomes/barbara/fastq/3_S3_L001_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_3_lo Input 1
101 /home/endrebas/genomes/barbara/fastq/3_S3_L002_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_3_lo Input 1
102 /home/endrebas/genomes/barbara/fastq/3_S3_L003_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_3_lo Input 1
103 /home/endrebas/genomes/barbara/fastq/3_S3_L004_R1_001.fastq.gz GENE2_KO_Input_1 GENE2_KO_ChIP_3_lo Input 1
104 /home/endrebas/genomes/barbara/fastq/6_S6_L001_R1_001.fastq.gz GENE2_KO_ChIP_1 GENE2_KO_ChIP_3_lo ChIP 1
105 /home/endrebas/genomes/barbara/fastq/6_S6_L002_R1_001.fastq.gz GENE2_KO_ChIP_1 GENE2_KO_ChIP_3_lo ChIP 1
106 /home/endrebas/genomes/barbara/fastq/6_S6_L003_R1_001.fastq.gz GENE2_KO_ChIP_1 GENE2_KO_ChIP_3_lo ChIP 1
107 /home/endrebas/genomes/barbara/fastq/6_S6_L004_R1_001.fastq.gz GENE2_KO_ChIP_1 GENE2_KO_ChIP_3_lo ChIP 1
108 /home/endrebas/genomes/barbara/fastq/9_S9_L001_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_3_lo Input 1
109 /home/endrebas/genomes/barbara/fastq/9_S9_L002_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_3_lo Input 1
110 /home/endrebas/genomes/barbara/fastq/9_S9_L003_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_3_lo Input 1
111 /home/endrebas/genomes/barbara/fastq/9_S9_L004_R1_001.fastq.gz GENE2_KO_Input_2 GENE2_KO_ChIP_3_lo Input 1
112 /home/endrebas/genomes/barbara/fastq/12_S12_L001_R1_001.fastq.gz GENE2_KO_ChIP_2 GENE2_KO_ChIP_3_lo ChIP 1
113 /home/endrebas/genomes/barbara/fastq/12_S12_L002_R1_001.fastq.gz GENE2_KO_ChIP_2 GENE2_KO_ChIP_3_lo ChIP 1
114 /home/endrebas/genomes/barbara/fastq/12_S12_L003_R1_001.fastq.gz GENE2_KO_ChIP_2 GENE2_KO_ChIP_3_lo ChIP 1
115 /home/endrebas/genomes/barbara/fastq/12_S12_L004_R1_001.fastq.gz GENE2_KO_ChIP_2 GENE2_KO_ChIP_3_lo ChIP 1
116 /home/endrebas/genomes/barbara/fastq/15_S15_L001_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_3_lo Input 1
117 /home/endrebas/genomes/barbara/fastq/15_S15_L002_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_3_lo Input 1
118 /home/endrebas/genomes/barbara/fastq/15_S15_L003_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_3_lo Input 1
119 /home/endrebas/genomes/barbara/fastq/15_S15_L004_R1_001.fastq.gz GENE2_KO_Input_3 GENE2_KO_ChIP_3_lo Input 1
120 /home/endrebas/genomes/barbara/fastq/1_S1_L001_R1_001.fastq.gz WT_Input_1 WT_ChIP_1_lo Input 1
121 /home/endrebas/genomes/barbara/fastq/1_S1_L002_R1_001.fastq.gz WT_Input_1 WT_ChIP_1_lo Input 1
122 /home/endrebas/genomes/barbara/fastq/1_S1_L003_R1_001.fastq.gz WT_Input_1 WT_ChIP_1_lo Input 1
123 /home/endrebas/genomes/barbara/fastq/1_S1_L004_R1_001.fastq.gz WT_Input_1 WT_ChIP_1_lo Input 1
124 /home/endrebas/genomes/barbara/fastq/7_S7_L001_R1_001.fastq.gz WT_Input_2 WT_ChIP_1_lo Input 1
125 /home/endrebas/genomes/barbara/fastq/7_S7_L002_R1_001.fastq.gz WT_Input_2 WT_ChIP_1_lo Input 1
126 /home/endrebas/genomes/barbara/fastq/7_S7_L003_R1_001.fastq.gz WT_Input_2 WT_ChIP_1_lo Input 1
127 /home/endrebas/genomes/barbara/fastq/7_S7_L004_R1_001.fastq.gz WT_Input_2 WT_ChIP_1_lo Input 1
128 /home/endrebas/genomes/barbara/fastq/10_S10_L001_R1_001.fastq.gz WT_ChIP_2 WT_ChIP_1_lo ChIP 1
129 /home/endrebas/genomes/barbara/fastq/10_S10_L002_R1_001.fastq.gz WT_ChIP_2 WT_ChIP_1_lo ChIP 1
130 /home/endrebas/genomes/barbara/fastq/10_S10_L003_R1_001.fastq.gz WT_ChIP_2 WT_ChIP_1_lo ChIP 1
131 /home/endrebas/genomes/barbara/fastq/10_S10_L004_R1_001.fastq.gz WT_ChIP_2 WT_ChIP_1_lo ChIP 1
132 /home/endrebas/genomes/barbara/fastq/13_S13_L001_R1_001.fastq.gz WT_Input_3 WT_ChIP_1_lo Input 1
133 /home/endrebas/genomes/barbara/fastq/13_S13_L002_R1_001.fastq.gz WT_Input_3 WT_ChIP_1_lo Input 1
134 /home/endrebas/genomes/barbara/fastq/13_S13_L003_R1_001.fastq.gz WT_Input_3 WT_ChIP_1_lo Input 1
135 /home/endrebas/genomes/barbara/fastq/13_S13_L004_R1_001.fastq.gz WT_Input_3 WT_ChIP_1_lo Input 1
136 /home/endrebas/genomes/barbara/fastq/16_S16_L001_R1_001.fastq.gz WT_ChIP_3 WT_ChIP_1_lo ChIP 1
137 /home/endrebas/genomes/barbara/fastq/16_S16_L002_R1_001.fastq.gz WT_ChIP_3 WT_ChIP_1_lo ChIP 1
138 /home/endrebas/genomes/barbara/fastq/16_S16_L003_R1_001.fastq.gz WT_ChIP_3 WT_ChIP_1_lo ChIP 1
139 /home/endrebas/genomes/barbara/fastq/16_S16_L004_R1_001.fastq.gz WT_ChIP_3 WT_ChIP_1_lo ChIP 1
140 /home/endrebas/genomes/barbara/fastq/1_S1_L001_R1_001.fastq.gz WT_Input_1 WT_ChIP_2_lo Input 1
141 /home/endrebas/genomes/barbara/fastq/1_S1_L002_R1_001.fastq.gz WT_Input_1 WT_ChIP_2_lo Input 1
142 /home/endrebas/genomes/barbara/fastq/1_S1_L003_R1_001.fastq.gz WT_Input_1 WT_ChIP_2_lo Input 1
143 /home/endrebas/genomes/barbara/fastq/1_S1_L004_R1_001.fastq.gz WT_Input_1 WT_ChIP_2_lo Input 1
144 /home/endrebas/genomes/barbara/fastq/4_S4_L001_R1_001.fastq.gz WT_ChIP_1 WT_ChIP_2_lo ChIP 1
145 /home/endrebas/genomes/barbara/fastq/4_S4_L002_R1_001.fastq.gz WT_ChIP_1 WT_ChIP_2_lo ChIP 1
146 /home/endrebas/genomes/barbara/fastq/4_S4_L003_R1_001.fastq.gz WT_ChIP_1 WT_ChIP_2_lo ChIP 1
147 /home/endrebas/genomes/barbara/fastq/4_S4_L004_R1_001.fastq.gz WT_ChIP_1 WT_ChIP_2_lo ChIP 1
148 /home/endrebas/genomes/barbara/fastq/7_S7_L001_R1_001.fastq.gz WT_Input_2 WT_ChIP_2_lo Input 1
149 /home/endrebas/genomes/barbara/fastq/7_S7_L002_R1_001.fastq.gz WT_Input_2 WT_ChIP_2_lo Input 1
150 /home/endrebas/genomes/barbara/fastq/7_S7_L003_R1_001.fastq.gz WT_Input_2 WT_ChIP_2_lo Input 1
151 /home/endrebas/genomes/barbara/fastq/7_S7_L004_R1_001.fastq.gz WT_Input_2 WT_ChIP_2_lo Input 1
152 /home/endrebas/genomes/barbara/fastq/13_S13_L001_R1_001.fastq.gz WT_Input_3 WT_ChIP_2_lo Input 1
153 /home/endrebas/genomes/barbara/fastq/13_S13_L002_R1_001.fastq.gz WT_Input_3 WT_ChIP_2_lo Input 1
154 /home/endrebas/genomes/barbara/fastq/13_S13_L003_R1_001.fastq.gz WT_Input_3 WT_ChIP_2_lo Input 1
155 /home/endrebas/genomes/barbara/fastq/13_S13_L004_R1_001.fastq.gz WT_Input_3 WT_ChIP_2_lo Input 1
156 /home/endrebas/genomes/barbara/fastq/16_S16_L001_R1_001.fastq.gz WT_ChIP_3 WT_ChIP_2_lo ChIP 1
157 /home/endrebas/genomes/barbara/fastq/16_S16_L002_R1_001.fastq.gz WT_ChIP_3 WT_ChIP_2_lo ChIP 1
158 /home/endrebas/genomes/barbara/fastq/16_S16_L003_R1_001.fastq.gz WT_ChIP_3 WT_ChIP_2_lo ChIP 1
159 /home/endrebas/genomes/barbara/fastq/16_S16_L004_R1_001.fastq.gz WT_ChIP_3 WT_ChIP_2_lo ChIP 1
160 /home/endrebas/genomes/barbara/fastq/1_S1_L001_R1_001.fastq.gz WT_Input_1 WT_ChIP_3_lo Input 1
161 /home/endrebas/genomes/barbara/fastq/1_S1_L002_R1_001.fastq.gz WT_Input_1 WT_ChIP_3_lo Input 1
162 /home/endrebas/genomes/barbara/fastq/1_S1_L003_R1_001.fastq.gz WT_Input_1 WT_ChIP_3_lo Input 1
163 /home/endrebas/genomes/barbara/fastq/1_S1_L004_R1_001.fastq.gz WT_Input_1 WT_ChIP_3_lo Input 1
164 /home/endrebas/genomes/barbara/fastq/4_S4_L001_R1_001.fastq.gz WT_ChIP_1 WT_ChIP_3_lo ChIP 1
165 /home/endrebas/genomes/barbara/fastq/4_S4_L002_R1_001.fastq.gz WT_ChIP_1 WT_ChIP_3_lo ChIP 1
166 /home/endrebas/genomes/barbara/fastq/4_S4_L003_R1_001.fastq.gz WT_ChIP_1 WT_ChIP_3_lo ChIP 1
167 /home/endrebas/genomes/barbara/fastq/4_S4_L004_R1_001.fastq.gz WT_ChIP_1 WT_ChIP_3_lo ChIP 1
168 /home/endrebas/genomes/barbara/fastq/7_S7_L001_R1_001.fastq.gz WT_Input_2 WT_ChIP_3_lo Input 1
169 /home/endrebas/genomes/barbara/fastq/7_S7_L002_R1_001.fastq.gz WT_Input_2 WT_ChIP_3_lo Input 1
170 /home/endrebas/genomes/barbara/fastq/7_S7_L003_R1_001.fastq.gz WT_Input_2 WT_ChIP_3_lo Input 1
171 /home/endrebas/genomes/barbara/fastq/7_S7_L004_R1_001.fastq.gz WT_Input_2 WT_ChIP_3_lo Input 1
172 /home/endrebas/genomes/barbara/fastq/10_S10_L001_R1_001.fastq.gz WT_ChIP_2 WT_ChIP_3_lo ChIP 1
173 /home/endrebas/genomes/barbara/fastq/10_S10_L002_R1_001.fastq.gz WT_ChIP_2 WT_ChIP_3_lo ChIP 1
174 /home/endrebas/genomes/barbara/fastq/10_S10_L003_R1_001.fastq.gz WT_ChIP_2 WT_ChIP_3_lo ChIP 1
175 /home/endrebas/genomes/barbara/fastq/10_S10_L004_R1_001.fastq.gz WT_ChIP_2 WT_ChIP_3_lo ChIP 1
176 /home/endrebas/genomes/barbara/fastq/13_S13_L001_R1_001.fastq.gz WT_Input_3 WT_ChIP_3_lo Input 1
177 /home/endrebas/genomes/barbara/fastq/13_S13_L002_R1_001.fastq.gz WT_Input_3 WT_ChIP_3_lo Input 1
178 /home/endrebas/genomes/barbara/fastq/13_S13_L003_R1_001.fastq.gz WT_Input_3 WT_ChIP_3_lo Input 1
179 /home/endrebas/genomes/barbara/fastq/13_S13_L004_R1_001.fastq.gz WT_Input_3 WT_ChIP_3_lo Input 1"""

    return pd.read_table(StringIO(contents), sep=" ", index_col=0)


def test_create_sample_sheets(sample_sheet, expected_result):

    df = create_sample_sheet(sample_sheet)

    assert df.equals(expected_result)
