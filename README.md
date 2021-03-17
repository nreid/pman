# analysis of peromyscus maniculatus ddRAD data

## in silico digest of p. man genome

download genome for _P. maniculatus bairdii_ from NCBI: accession GCA_003704035.3. see script [01_get_genome.sh](/scripts/01_get_genome.sh). 

run _in silico_ digest in R. see script [find_cut_sites.R](/scripts/find_cut_sites.R)

### results

length of chromosome contigs is 2.43gb. 

table of RE pairs and frequency of ddRAD sites:

| rarecut | freqcut | nsites | freqperbp |
| --------|---------|--------|-----------
| sbfI  |  aciI |  12691 | 191982.038 |
| sbfI  |  aluI |  21484 | 113407.375 |
| sbfI  |  aseI |   7669 | 317700.358 |
| sbfI  | bamHI |   5132 | 474755.269 |
| sbfI  |  bfaI |  26797 |  90922.269 |
| sbfI  | bfuCI |  25702 |  94795.893 |
| sbfI  | ecoRI |   6292 | 387228.869 |
| sbfI  |  mseI |  23193 | 105050.836 |
| sbfI  |  mspI |  14504 | 167984.283 |
| nsiI  |  aciI |  95033 |  25637.874 |
| nsiI  |  aluI | 260301 |   9360.103 |
| nsiI  |  aseI | 106420 |  22894.607 |
| nsiI  | bamHI |  41755 |  58350.953 |
| nsiI  |  bfaI | 277080 |   8793.287 |
| nsiI  | bfuCI | 274608 |   8872.444 |
| nsiI  | ecoRI |  75582 |  32235.771 |
| nsiI  |  mseI | 193194 |  12611.386 |
| nsiI  |  mspI | 112054 |  21743.481 |
| pstI  |  aciI | 126578 |  19248.559 |
| pstI  |  aluI | 254282 |   9581.661 |
| pstI  |  aseI |  98219 |  24806.240 |
| pstI  | bamHI |  53159 |  45833.143 |
| pstI  |  bfaI | 301668 |   8076.574 |
| pstI  | bfuCI | 300711 |   8102.278 |
| pstI  | ecoRI |  73935 |  32953.865 |
| pstI  |  mseI | 264446 |   9213.390 |
| pstI  |  mspI | 146287 |  16655.233 |
 