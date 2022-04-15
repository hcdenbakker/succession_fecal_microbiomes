#Example of a preprocessing file, this was run for each of the 5 runs for the 16S amplicon
#data

library(dada2); packageVersion("dada2")
#library(phyloseq)
#library(ggplot2)

#The following lines of code are necessary for the preprocessing of the raw reads
# I have blocked them out with "'''"

path <- "./fastq"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
#plotQualityProfile(fnFs[4:5])
#plotQualityProfile(fnRs[4:5])
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#trim reads so they overlap, and primers are excluded (trimLeft arguement)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft=20)

out

errF <- learnErrors(filtFs, multithread=TRUE, verbose = TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, verbose = TRUE)

#plotErrors(errF, nominalQ=TRUE)

#consider one of the pooling options here!
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = "pseudo")

dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = "pseudo")

dadaFs[[4]]

#merging
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample

head(mergers[[4]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

table(nchar(getSequences(seqtab.nochim)))

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

saveRDS(seqtab.nochim, "table_run2_defaultEE.rds")
