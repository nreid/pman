library(tidyverse)
library(vcfR)
library(ape)
library(beeswarm)
library(vroom)
library(RColorBrewer)

# to do
	# visualize mapping rate data

	# read in populations.sumstats.tsv.gz

	# dig into VCF
		# depth
		# heterozygosity
		# hdplot

##############################
# visualize mapping rate data
##############################

# these steps aggregate
	# 1) metadata
	# 2) mapping data
	# 3) barcode/pool data
	
# read in metadata table
meta <- read.table("../meta/metadata.txt",header=FALSE,sep="\t",quote="",comment.char="")
	colnames(meta) <- c("sample", "ID", "phenotype", "sample_in_pool","barcode","pool","barcode2", "spacer")
	meta[,1] <- paste(meta[,2],meta[,3],sep="_")

# read in sam stats table for raw aligned data
snraw <- read.table("../results/align_stats/SN.raw.txt",sep="\t") %>% t()
	snraw <- cbind(sample=str_extract(rownames(snraw),"(?<=^X)[^\\.]+"),data.frame(snraw))
	colnames(snraw) <- colnames(snraw) %>% str_replace(regex("\\.*$"),"")

	snraw <- snraw[,c(1,2,8,9,14,20,24,39)]

# read in sam stats table for deduped aligned data
snded <- read.table("../results/align_stats/SN.dedup.txt",sep="\t") %>% t()
	snded <- cbind(sample=str_extract(rownames(snded),"(?<=^X)[^\\.]+"),data.frame(snded))
	colnames(snded) <- colnames(snded) %>% str_replace(regex("\\.*$"),"")

	snded <- snded[,c(1,2,8,9,14,20,24,39)]


# add dedup mapping data to metadata tables
metaraw <- left_join(x=meta,y=snraw) 
metaded <- left_join(x=meta,y=snded) 

# edit column names to replace periods with underscores
colnames(metaraw) <-  gsub("\\.", "_",colnames(metaraw))
colnames(metaded) <-  gsub("\\.", "_",colnames(metaded))
 

# look at mapping rates for raw sequences
plot(metaraw$reads_mapped/metaraw$raw_total_sequences,col=factor(metaraw$pool))

# look at dedup/raw sequences fraction
	# duplication rates vary by pool, with pools 1 and 4 having ~5% higher duplication rates
plot(metaded$reads_mapped/metaraw$reads_mapped,col=factor(metaraw$pool))

# look at proper pairing rate for dedup sequences
plot(metaded$percentage_of_properly_paired_reads,col=factor(metaded$pool),pch=as.numeric(as.factor(metaded$phenotype)))

# look at overall mismatch rate between reads and ref genome
plot(metaded$mismatches/metaded$bases_mapped,col=factor(metaded$pool),pch=as.numeric(as.factor(metaded$phenotype)))

# there is a disturbing bimodality in the mismatch rate and the properly paired read mapping rate
	# a bunch of small spot individuals have a 60% higher mismatch rate
	# the same individuals have a 2-3% lower proper read pairing rate

	# this probably means those individuals have much higher non-pman ancestry, but why?


######################################
# read in populations.sumstats.tsv.gz
######################################

# gzipped in advance, not by stacks
sumstats <- read.table("../results/stacks/populations.sumstats.tsv",skip=1,header=TRUE,comment.char="",sep="\t")
	colnames(sumstats) <- colnames(sumstats) %>%
		str_replace(.,"^X..","") %>% 
		str_replace_all(.,"\\.","_")

hist(sumstats$Obs_Het)

######################################
# dig into VCF, run HDplot
######################################

# source McKinney et al 2017 hdplot function
	# https://github.com/gjmckinney/HDplot
source("hdplot.R")

vcf<-read.vcfR("../results/stacks/filtered.vcf.gz")

HDres <- HDplot(vcf)

# plot of heterozygosity against read ratio deviation
HDres %>% ggplot() + geom_point(aes(x=H,y=D),size=0.5)

# again but zoom in on Y
	HDres %>% ggplot() + geom_point(aes(x=H,y=D),size=0.5) + ylim(-10,10)


# plot of heterozygosity against read ratio
HDres %>% ggplot()+geom_point(aes(x=H,y=ratio),size=0.5)

# extract the depth of coverage for each individual genotype
dp<-extract.gt(vcf, element = "DP", mask = FALSE, as.numeric = TRUE, return.alleles = FALSE, 
                    IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)

# plot depth of each snp along the genome, colored by chromosome
plot(rowSums(dp),pch=20,cex=.2,col=factor(vcf@fix[,1]))
	# again but zoom in on the Y
	plot(rowSums(dp),pch=20,cex=.2,col=factor(vcf@fix[,1]),ylim=c(0,5000))

# plot a histogram of depth of coverage for each genotype
hist(rowSums(dp),breaks=1000,xlim=c(0,5000))

# plot total site depth as a function of heterozygosity
plot(HDres$H,rowSums(dp),ylim=c(0,20000),pch=20,cex=.4)


# look more closely at individuals---------------------------

# extract genotypes, recode as 0,1,2
genos<-extract.gt(vcf, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, 
                    IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
	genos[is.na(genos)]<-"./."
	genos <- gsub("\\|","/",genos)
	
	genos<-apply(genos,2,function(x) dplyr::recode(x,'0/0'=0,'1/1'=2,'0/1'=1,'1/0'=1,.default=NA_real_))

# create a matrix of distances between genotypes
gdist <- genos %>% t() %>% dist()

# plot an ugly neighbor-joining tree
	# it's not too hard to add color to the labels

nj(gdist) %>% plot(.,"unrooted")

# make an mds plot (similar to a PCA)
gdist <- genos %>% t() %>% dist()
mds <- cmdscale(gdist) 
	colnames(mds) <- c("mds1","mds2")
	mds <- data.frame(sample=rownames(mds),mds)
	mds <- left_join(mds,meta)
	
ggplot(mds,aes(x=mds1,y=mds2,color=phenotype)) +
    geom_point()



ihz <- colMeans(genos==1,na.rm=TRUE)
plot(ihz)

iaa <- colMeans(genos,na.rm=TRUE)
plot(iaa)

hoa <- colMeans(genos==2,na.rm=TRUE)
plot(hoa)


beeswarm(ihz ~ meta[colnames(genos),"region"],pch=20,ylab="individual heterozygosity",xlab="region")

plot(rowMeans(genos[,26:50],na.rm=TRUE)/2 - rowMeans(genos[,1:25],na.rm=TRUE)/2,pch=20,cex=.2,col=factor(vcf@fix[,1]))

plot(rowMeans(genos[,1:50 > 25 & iaa < 0.8],na.rm=TRUE)/2 - rowMeans(genos[,1:25],na.rm=TRUE)/2,pch=20,cex=.2,col=factor(vcf@fix[,1]))

plot(rowMeans(genos[,1:50 > 25 & iaa > 0.8],na.rm=TRUE)/2 - rowMeans(genos[,1:25],na.rm=TRUE)/2,pch=20,cex=.2,col=factor(vcf@fix[,1]))


#########################################
# how to create a site list for filtering
#########################################


# for use in vcftools: tab separated, two columns, SEQUENCE POSITION

	ex <- rowSums(dp) > 40000 

	vcf@fix[ex,1:2]

	# then use write.table() to write the file

# stacks wants a list of its locus IDs to use as a blacklist

	# this is pretty easy if you're using the sumstats.tsv file
		# grab a set of markers with, say > 5 failed HWE tests and write out the locus IDs

	ooHWE <- group_by(sumstats,Locus_ID,Chr,BP) %>% 
		summarize(pops_oo_hwe=sum(HWE_P_value < 0.005))

	# if you want to filter on information you extracted from the VCF, then you'll need to do a "join" operation


