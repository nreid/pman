library(tidyverse)
library(Biostrings)     # Provides DNAString, DNAStringSet, etc
library(BSgenome)       # Provides getSeq()
library(GenomicRanges)  # Provides GRanges, etc
library(rtracklayer)    # Provides import() and export()


# This R script finds sbf1 cut sites, flanking sites and writes them to bed files

# get genome first:

# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/704/035/GCA_003704035.3_HU_Pman_2.1.3/GCA_003704035.3_HU_Pman_2.1.3_genomic.fna.gz

sbf1 <- "CCTGCAGG"
msei <- "TTAA"

ecoR1 <- "GAATTC"


# read in reference genome 
seqs <- readDNAStringSet("../genome/GCA_003704035.3_HU_Pman_2.1.3_genomic.fna.gz")

# get contig lengths
clen <- sapply(seqs,length)

# get total genome length
glen <- sum(clen)

# ignore small contigs
	# genome is 2.43gb (not what NCBI says though...)
seqs <- seqs[clen > 1e6]
clen <- clen[clen > 1e6]

# find sbf1 sites and make a GRanges object
rc <- vmatchPattern(ecoR1,seqs)
rc2 <- as(rc, "GRanges")

# find msei sites and get granges object
fc <- vmatchPattern(msei,seqs)
fc2 <- as(fc, "GRanges")

# get closest upstream and downstream msei cut sites to sbf1 cut sites
rightind <- precede(rc2,fc2,"first")
	NOrightNA <- !is.na(rightind)
leftind <- follow(rc2,fc2,"last")
	NOleftNA <- !is.na(leftind)

out <- matrix(nrow=length(rc2),ncol=2)

out[NOrightNA,1] <- distance(rc2[NOrightNA,],fc2[rightind[NOrightNA],])
out[NOleftNA,2] <- distance(rc2[NOleftNA,],fc2[leftind[NOleftNA],])

# `out` now contains all SbfI-MseI RAD fragment lengths in the reference genome. 
	# column 1 is the left side of SbfI, column 2 is the right. 

# get number of RAD sites with a given 
nsites <- (as.numeric(rowSums(out > 250 & out < 350) > 0) %>% table())[2]

# RAD site frequency:
glen/nsites






