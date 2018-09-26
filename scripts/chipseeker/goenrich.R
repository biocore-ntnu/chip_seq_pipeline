
library(ChIPseeker)

library(clusterProfiler)

library(org.Hs.eg.db)

annos = readRDS(snakemake@input[[1]])
genes = unique(as.data.frame(anno)$geneId)

data(geneList, package="DOSE")

ego <- enrichGO(gene          = genes,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

outdf = as.data.frame(ego)
write.table(snakemake@output[[1]], sep=" ")
## pdf(snakemake@output[[1]])
## plotAnnoBar(annos)
## dev.off()
