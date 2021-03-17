library(tidyverse)
library(Biostrings)     # Provides DNAString,DNAStringSet,etc
library(BSgenome)       # Provides getSeq()
library(GenomicRanges)  # Provides GRanges,etc
library(rtracklayer)    # Provides import() and export()


# This R script finds sbf1 cut sites

# get genome first:

# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/704/035/GCA_003704035.3_HU_Pman_2.1.3/GCA_003704035.3_HU_Pman_2.1.3_genomic.fna.gz

# read in reference genome 
seqs <- readDNAStringSet("../genome/GCA_003704035.3_HU_Pman_2.1.3_genomic.fna.gz")

# get contig lengths
clen <- nchar(seqs)

# get total genome length
glen <- sum(clen)

# ignore small contigs
	# genome is 2.43gb (not what NCBI says though...)
seqs <- seqs[clen > 1e6]
clen <- clen[clen > 1e6]

# restriction enzymes:

P7 <- list(
	sbfI = "CCTGCAGG",
	nsiI = "ATGCAT",
	pstI = "CTGCAG"
	)

RE <- list(
	aciI = "CCGC",
	# ageI
	aluI = "AGCT",
	# apaLI
	apeKI = c("GCAGC","GCTGC"),
	apoI = c("AAATTC","AAATTT","GAATTC","GAATTT"),
	aseI = "ATTAAT",
	bamHI = "GGATCC",
	# bbvCI
	bfaI = "CTAG",
	bfuCI = "GATC",
	# bgIII
	# bsaHI
	# bspDI
	# bstYI
	# cac8I
	# claI
	# csp6I
	# ddeI
	# dpnII
	# eaeI
	ecoRI = "GAATTC",
	# ecoRV
	# ecoT22I
	# haeIII
	# hinP1I
	# hindIII
	# hpaII
	# kpnI
	# mluCI
	mseI = "TTAA",
	# mslI
	mspI = "CCGG"
	# ncoI
	# ndeI
	# nheI
	# nlaIII
	# notI
	# nspI
	# rsaI
	# sacI
	# sau3AI
	# sexAI
	# sgrAI
	# speI
	# sphI
	# taqI
	# xbaI
	# xhoI
)

# a function to do the digest for a pair of enzymes
	# frequentcutter can have ambiguous sites, just supply a character vector with each resolution
	# returns 
		# a GRanges object with all the rarecutter sites
		# a corresponding matrix with left and right fragment lengths (NA if no adjacent frequentcutter site)
		# number of sites
		# frequency/bp of sites


digest <- function(rarecutter, frequentcutter, genomeseq, lowerbound=250, upperbound=450){

	# get contig lengths
	clen <- nchar(genomeseq)

	# get total genome length
	glen <- sum(clen)

	# find common sites and get granges object
	rc <- vmatchPattern(rarecutter, genomeseq)
	rc2 <- as(rc, "GRanges")
		
	# find common sites and get granges object
		# if multiple cutters are supplied (e.g. for a cutter with ambiguities), use second approach
	if(length(frequentcutter) == 1){
		fc <- vmatchPattern(frequentcutter, genomeseq)
		fc2 <- as(fc, "GRanges")
	}
	if(length(frequentcutter) > 1){
		fc <- sapply(frequentcutter,vmatchPattern,subject=genomeseq) %>% sapply(., as, Class="GRanges")
		fc2 <- fc[[1]]
		for(i in 2:length(fc)){
			fc2 <- c(fc2,fc[[i]])
		}
	}

	# get distances to closest upstream and downstream cut sites to rare cutter for both frequent and rare cutter

	# first get rare cutter distances
	rrightind <- precede(rc2, rc2, "first")
		rNOright <- !is.na(rrightind)
	rleftind <- follow(rc2, rc2, "last")
		rNOleft <- !is.na(rleftind)

	r2r <- matrix(nrow=length(rc2), ncol=2)
	r2r[rNOright, 1] <- distance(rc2[rNOright,],rc2[rrightind[rNOright],])
	r2r[rNOleft, 2] <- distance(rc2[rNOleft,],rc2[rleftind[rNOleft],])

	# then get frequent cutter distances
	frightind <- precede(rc2, fc2, "first")
		fNOright <- !is.na(frightind)
	fleftind <- follow(rc2, fc2, "last")
		fNOleft <- !is.na(fleftind)
	
	r2f <- matrix(nrow=length(rc2), ncol=2)
	r2f[fNOright, 1] <- distance(rc2[fNOright,],fc2[frightind[fNOright],])
	r2f[fNOleft, 2] <- distance(rc2[fNOleft,],fc2[fleftind[fNOleft],])
	
	# create "out" object, giving RAD fragment sizes to left and right for each fragment
	# if another rare cut site occurs closer than a frequent cut site, change the distance to NA
	out <- r2f
	out[r2f > r2r] <- NA


	# get number of RAD sites 
	nsites <- (as.numeric(rowSums(out > lowerbound & out < upperbound,na.rm=TRUE) > 0) %>% table())[2]
	# RAD site frequency:
	bpfreq <- glen/nsites

	return(list(rarecuttersites=rc2,fragementlengths=out,nsites=nsites,frequencyperbp=bpfreq))
}


# loop through REs
k <- 1
out <- vector("list",length(P7)*length(RE))
dat <- data.frame(rarecut=character(length=27),freqcut=character(length=27),nsites=numeric(length=27),freqperbp=numeric(length=27))

for(i in 1:length(P7)){
	for(j in 1:length(RE)){
		
		# out[[k]] <- digest(rarecutter=P7[[i]],frequentcutter=RE[[j]],genomeseq=seqs)
		dat[k,1] <- names(P7)[[i]]
		dat[k,2] <- names(RE)[[j]]
		dat[k,3] <- out[[k]][3]
		dat[k,4] <- out[[k]][4]
		k <- k+1
		print(c(i,j))
	}
}



