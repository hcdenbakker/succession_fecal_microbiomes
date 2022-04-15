library(dada2); packageVersion("dada2")
# Merge multiple runs (if necessary)

st1 <- readRDS("./run2/table_run2.rds")
st2 <- readRDS("./run3/table_run3.rds")
st3 <- readRDS("./run4/table_run4.rds")
st4 <- readRDS("./run5/table_run5.rds")
st5 <- readRDS("./run6/table_run6.rds")
st.all <- mergeSequenceTables(st1, st2, st3, st4, st5)
dim(st.all)
# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
dim(seqtab)
sum(seqtab)/sum(st.all)

#remove short amplicons
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 397:428]

# Write to disk
saveRDS(seqtab, "./seqtab_final_all_rerun.rds")
