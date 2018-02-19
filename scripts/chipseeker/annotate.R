file = snakemake@input[[1]]

tss_length = snakemake@config[["barchart_tss_length"]]

txdb_name = snakemake@params[["txdb"]]

suppressWarnings(suppressMessages(require(txdb_name, character.only=TRUE)))

suppressWarnings(suppressMessages(require(ChIPseeker)))

peakAnno <- annotatePeak(file, tssRegion=c(-tss_length, tss_length),
                         TxDb=eval(parse(text=txdb_name)))

saveRDS(peakAnno, snakemake@output[[1]])
