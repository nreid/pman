# analysis of peromyscus maniculatus ddRAD data

## in silico digest of p. man genome

download genome for _P. maniculatus bairdii_ from NCBI: accession GCA_003704035.3. see script [01_get_genome.sh](/scripts/01_get_genome.sh). 

run _in silico_ digest in R. see script [find_cut_sites.R](/scripts/find_cut_sites.R)

### results

length of chromosome contigs is 2.43gb. 

restriction enzymes ecoRI + mseI size selected to 300-450bp produces 114,404 sites, at a frequency of 1 per 21,297bp. 