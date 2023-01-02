


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.13")
library(dada2);packageVersion("dada2")
## [1] '1.20.0'
library(Biostrings); packageVersion("Biostrings")
## [1] '2.60.2'
library(ShortRead); packageVersion("ShortRead")
## [1] '1.50.0'
library(ggplot2); packageVersion("ggplot2")
## [1] '3.3.5'
library(reshape2); packageVersion("reshape2")
## [1] '1.4.4'
library(gridExtra); packageVersion("gridExtra")
## [1] '2.3'
#BiocManager::install("phyloseq")
library(phyloseq); packageVersion("phyloseq")
## [1] '1.36.0'

path1 <- "D:/16S_seqs_PacBio/" # CHANGE ME to location of the First Replicate fastq files
path.out <- "D:/16S_seqs_PacBio/Figures/"
path.rds <- "D:/16S_seqs_PacBio/RDS/"
fns1 <- list.files(path1, pattern="fastq", full.names=TRUE)
x <- readFastq("D:/16S_seqs_PacBio/")

F27 <- "AGRGTTYGATYMTGGCTCAG"
R1492 <- "RGYTACCTTGTTACGACTT"
rc <- dada2:::rc
theme_set(theme_bw())

## https://urldefense.proofpoint.com/v2/url?u=https-3A__benjjneb.github.io_LRASManuscript_LRASms-5Ffecal.html&d=DwIGAg&c=sJ6xIWYx-zLMB3EPkvcnVg&r=1ag44WqDjpGP-PgOfmXDVzD8kytJpCm5jahIwIWsSXE&m=eavk7tyFIigW4s74NueOiRORmtxLaCr0qqI15pm0FuI&s=UVDMCvJJ4lrn7TtPoQdTMCWlY6So7Q0YJ4JV8O1kZrI&e= 

## remove primers
nops1 <- file.path(path1, "noprimers", basename(fns1))
prim1 <- removePrimers(fns1, nops1, primer.fwd=F27, primer.rev=dada2:::rc(R1492), orient=TRUE)

## filter
lens.fn <- lapply(nops1, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)

filts1 <- file.path(path1, "noprimers", "filtered", basename(fns1))
track1 <- filterAndTrim(nops1, filts1, minQ=3, minLen=1000, maxLen=1600, maxN=0, rm.phix=FALSE, maxEE=2)
track1

## dada2
drp1 <- derepFastq(filts1, verbose=TRUE)

## Learn errors
err1 <- learnErrors(drp1, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)
saveRDS(err1, file="soil_err1.RDS") ## change me to your save

## Plot errors
plotErrors(err1)

## Denoise
dd1 <- dada(drp1, err=err1, BAND_SIZE=32, multithread=TRUE)
saveRDS(dd1, file = "soil_dd1.rds")

cbind(ccs=prim1[,1], primers=prim1[,2], filtered=track1[,2], denoised=sapply(dd1, function(x) sum(x$denoised)))

############################################
## Sequence table
st1 <- makeSequenceTable(dd1); dim(st1)

## import the denoise file from the first run and call it dd
dd <- readRDS("soil_dd1.rds")
stcomb <- makeSequenceTable(dd)

## check rownames, should still be 126 samples
rownames(stcomb)

## Assign taxonomy, download appropriate reference database from here - https://urldefense.proofpoint.com/v2/url?u=https-3A__benjjneb.github.io_dada2_training.html&d=DwIGAg&c=sJ6xIWYx-zLMB3EPkvcnVg&r=1ag44WqDjpGP-PgOfmXDVzD8kytJpCm5jahIwIWsSXE&m=eavk7tyFIigW4s74NueOiRORmtxLaCr0qqI15pm0FuI&s=EZBbAJcvYR5Vg0uxPWXHqX-_EuAdcAWJKFHCah3reow&e= 
tax1 <- assignTaxonomy(stcomb, "silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE) # Slowest part
#tax1[,"Genus"] <- gsub("Escherichia/Shigella", "Escherichia", tax1[,"Genus"]) # Reformat to be compatible with other data sources
head(unname(tax1))

## Check Chimeras
bim1 <- isBimeraDenovo(stcomb, minFoldParentOverAbundance=3.5, multithread=FALSE)
table(bim1)
sum(stcomb[,bim1])/sum(stcomb)

#######################################################################

## Extract Sample Names
sample_names <- read.csv("sample_data.csv", header=TRUE)
#sample.names1 <- sapply(strsplit(fns1, "-"), function(x) paste(x[1]))
#sample.names1 <- sapply(strsplit(sample.names1, "/"), function(x) paste(x[4]))
rownames(stcomb) <- sample_names$Sample_Name
sample.names1

## remove chimeras
seqtab <- removeBimeraDenovo(stcomb, method="consensus", multithread=FALSE)

# Write to disk
saveRDS(seqtab, "seqtab_final.rds") 
saveRDS(tax1, "tax_final.rds")

################### construct tree ####################################
##https://urldefense.proofpoint.com/v2/url?u=https-3A__f1000research.com_articles_5-2D1492_v2&d=DwIGAg&c=sJ6xIWYx-zLMB3EPkvcnVg&r=1ag44WqDjpGP-PgOfmXDVzD8kytJpCm5jahIwIWsSXE&m=eavk7tyFIigW4s74NueOiRORmtxLaCr0qqI15pm0FuI&s=BXKohgSfL78tGajE_hHqhSpElqTvr8Bensr-Lxuyx1Q&e= 
# 
# library(DECIPHER)
# 
# seqs <- getSequences(seqtab)
# names(seqs) <- seqs # This propagates to the tip labels of the tree
# alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
# 
# ## use phangorn to construct the tree
# 
# library(phangorn)
# 
# phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
# dm <- dist.ml(phang.align)
# treeNJ <- NJ(dm) # Note, tip order != sequence order
# fit = pml(treeNJ, data=phang.align)
# 
# ## negative edges length changed to 0!
# 
# fitGTR <- update(fit, k=4, inv=0.2)
# fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                     rearrangement = "stochastic", control = pml.control(trace = 0))
# 
# saveRDS(fitGTR, "fitGTR.RDS")
# detach("package:phangorn", unload=TRUE)