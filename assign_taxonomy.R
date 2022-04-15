library(dada2); packageVersion("dada2")
seqtab <- readRDS("./seqtab_all_rerun_defaultEE.rds")

taxa <- assignTaxonomy(seqtab, "/media/henk/TheRAID5/GTDB/r202/DADA2/GTDB_bac120_arc122_ssu_r202_Genus_no_euks.fa.gz", multithread=TRUE)

head(taxa)

saveRDS(taxa, "./tax_final_GTDB_defaultEE.rds")
