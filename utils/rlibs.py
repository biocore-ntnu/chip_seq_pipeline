import pandas as pd

from io import StringIO

def _txdb_df():
    contents = """Genome GeneType LibraryName
ce11 refGene TxDb.Celegans.UCSC.ce11.refGene
ce6 ensGene TxDb.Celegans.UCSC.ce6.ensGene
dm3 ensGene TxDb.Dmelanogaster.UCSC.dm3.ensGene
dm6 ensGene TxDb.Dmelanogaster.UCSC.dm6.ensGene
danRer10 refGene TxDb.Drerio.UCSC.danRer10.refGene
hg18 knownGene TxDb.Hsapiens.UCSC.hg18.knownGene
hg19 knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
hg19 lincRNAsTranscripts TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts
hg38 knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
mm10 ensGene TxDb.Mmusculus.UCSC.mm10.ensGene
mm10 knownGene TxDb.Mmusculus.UCSC.mm10.knownGene
mm9 knownGene TxDb.Mmusculus.UCSC.mm9.knownGene
rn4 ensGene TxDb.Rnorvegicus.UCSC.rn4.ensGene
rn5 refGene TxDb.Rnorvegicus.UCSC.rn5.refGene
rn6 refGene TxDb.Rnorvegicus.UCSC.rn6.refGene
sacCer2 sgdGene TxDb.Scerevisiae.UCSC.sacCer2.sgdGene
sacCer3 sgdGene TxDb.Scerevisiae.UCSC.sacCer3.sgdGene"""

    df = pd.read_table(StringIO(contents), sep=" ", header=0)

    return df

txdb_df = _txdb_df()

# def org_df():

#     contents = """
# org.Bt.eg.db
# org.Ce.eg.db
# org.Cf.eg.db
# org.Dm.eg.db
# org.Dr.eg.db
# org.Gg.eg.db
# org.Hs.eg.db
# org.Mm.eg.db
# org.Rn.eg.db
# org.Sc.sgd.db
# org.Ss.eg.db
#     """


    # df = pd.read_table(StringIO, sep=" ", header=0)

    # return df
