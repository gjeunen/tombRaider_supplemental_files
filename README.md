# *tombRaider* supplemental files

## 1. Introduction

This document serves as a comprehensive guide and supplementary information for the bioinformatic and statistical analysis of the manuscript entitled "*tombRaider* - improved species and haplotype recovery from metabarcoding data through artefact and pseudogene exclusion" by Jeunen *et al*., 2024. This work is associated with the Marsden Fast-Start fund (MFP-UOO002116). This document is split up into three parts, one for each of the data sets reanalysed for the publication.

## 2. Supplement 2: mock community data analysis

## 3. Supplement 3: air eDNA data analysis

## 4. Supplement 4: salmon haplotype data analysis

### 4.1 Starting files

Separate ASV tables can be downloaded from the online Supplementary Files associated with the original publication [Weitemier et al., 2021](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.15811) using [this link](https://ir.library.oregonstate.edu/concern/datasets/mp48sm28w). Download the "Oncorhynchus_eDNA_Supplemental_Data.zip" file. and place it in the starting directory and unzip the document.

```{code-block} bash
unzip Salmon\ Trout\ eDNA\ Data.zip
```

Unzipping the document will have created the directory Oncorhynchus_eDNA_Supplemental_Data. The individual ASV tables for the three species incorporated in our analysis (*Oncorhynchus clarkii*, *O. kisutch*, and *O. tshawytscha*) and target gene (ND2) can be found in subdirectory `Oncorhynchus_eDNA_Supplemental_Data/Chimera_Removed_summaries/Denoise_summary_alpha_4/`.

The starting files needed for the analysis are:

1. clarkii_ND2_combine_chimera_removed2.denoise_summary.csv
2. kisutch_ND2_combine.denoise_summary.csv
3. tshawytscha_ND2_combine.denoise_summary.csv

Within each ASV table are the locations in which the specific species was encountered as columns. Locations (columns) for which the specific species was not detected are removed from the ASV table. ASVs are presented in the rows with identifiers "UniqXX", where "XX" is a number. While numbers are unique within a species, they are re-occurring between species. Hence, these identifiers will need to be reformatted. Finally, the ASV sequence is present in the last column of the ASV table.

For this analysis, we will reformat the separate ASV tables and generate a single ASV table (output file name: "asv_table.txt") and a single two-line ASV fasta file (output file name: "asvs.fasta") using the python script below.

```{code-block} python
#!/usr/bin/env Python3

import collections

fasta_dict = {}
table_dict = collections.defaultdict(dict)
location_list = []
seq_id_list = []

with open('clarkii_ND2_combine_chimera_removed2.denoise_summary.csv', 'r') as clarkii_in:
    for line in clarkii_in:
        if line.startswith('id'):
            locations = line.split(',')[1:-1]
            for loc in locations:
                location_list.append(loc)
        else:
            seq_id = '>clarkii_' + line.split(',')[0]
            sequence = line.split(',')[-1].rstrip('\n')
            read_count = line.split(',')[1:-1]
            fasta_dict[seq_id] = sequence
            seq_id_list.append(seq_id)
            for item in locations:
                table_dict[item][seq_id.lstrip('>')] = read_count[locations.index(item)]

with open('kisutch_ND2_combine.denoise_summary.csv', 'r') as kisutch_in:
    for line in kisutch_in:
        line = line.replace('"', '')
        if 'id' in line.split(',')[0]:
            locations = line.split(',')[1:-1]
            for loc in locations:
                if loc not in location_list:
                    location_list.append(loc)
        else:
            seq_id = '>kisutch_' + line.split(',')[0]
            sequence = line.split(',')[-1].rstrip('\n')
            read_count = line.split(',')[1:-1]
            fasta_dict[seq_id] = sequence
            seq_id_list.append(seq_id)
            for item in locations:
                table_dict[item][seq_id.lstrip('>')] = read_count[locations.index(item)]

with open('tshawytscha_ND2_combine.denoise_summary.csv', 'r') as tshawytscha_in:
    for line in tshawytscha_in:
        line = line.replace('"', '')
        if 'id' in line.split(',')[0]:
            locations = line.split(',')[1:-1]
            for loc in locations:
                if loc not in location_list:
                    location_list.append(loc)
        else:
            seq_id = '>tshawytscha_' + line.split(',')[0]
            sequence = line.split(',')[-1].rstrip('\n')
            read_count = line.split(',')[1:-1]
            fasta_dict[seq_id] = sequence
            seq_id_list.append(seq_id)
            for item in locations:
                table_dict[item][seq_id.lstrip('>')] = read_count[locations.index(item)]

with open('asvs.fasta', 'w') as asv_out:
    for k,v in fasta_dict.items():
        asv_out.write(k + '\n' + v + '\n')

with open('asv_table.txt', 'w') as table_out:
    headerstring = 'id\t' + "\t".join(table_dict.keys()) + '\n'
    table_out.write(headerstring)
    for item in seq_id_list:
        table_out.write(f'{item}' + '\t')
        for key in table_dict:
            try:
                table_out.write(table_dict[key][item.lstrip('>')] + '\t')
            except KeyError:
                table_out.write('0\t')
        table_out.write('\n')
```

### 4.2 Reference database

For assigning a taxonomic ID to each ASV, we will generate a highly curated reference database for the ND2 gene using CRABS *v* 0.1.8 [(Jeunen et al., 2022)](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13741). First, download the ND2 reference sequences from NCBI.

```{code-block} bash
crabs db_download -s ncbi -db nucleotide -q 'ND2[All Fields] AND (NADH2[All Fields] AND animals[filter])' -o ncbi_ND2.fasta -e gjeunen@gmail.com
```

Additionally, the NCBI taxonomy information needs to be downloaded to complete the reference database creation with CRABS.

```{code-block} bash
crabs db_download -s taxonomy
```

Second, extract amplicons from the downloaded reference sequences using an *in silico* PCR analysis.

```{code-block} bash
crabs insilico_pcr -i ncbi_ND2.fasta -o ncbi_ND2_insilico.fasta -f CAACMGCCGCAGCAATAATCCT -r ATTTTACGCAGTTGGGTTTGATTAAGCCC
```

Third, assign a taxonomic lineage to each sequence.

```{code-block} bash
crabs assign_tax -i ncbi_ND2_insilico.fasta -o ncbi_ND2_insilico_tax.tsv -a nucl_gb.accession2taxid -t nodes.dmp -n names.dmp -w yes
```

Fourth, dereplicate the reference database to reduce the file size and remove reduntant sequences.

```{code-block} bash
crabs dereplicate -i ncbi_ND2_insilico_tax.tsv -o ncbi_ND2_insilico_tax_derep.tsv -m uniq_species
```

Fifth, use various filtering parameters to retain only high quality references in the local database.

```{code-block} bash
crabs seq_cleanup -i ncbi_ND2_insilico_tax_derep.tsv -o ncbi_ND2_insilico_tax_derep_clean.tsv -e yes -s yes -n 0
```

### 4.3 Taxonomic assignment

### 4.4 *tombRaider*

Once the three input files ("asv_table.txt", "asvs.fasta", and "blast_taxonomy.txt") are generated, we will execute *tombRaider* to identify and remove artefacts using the `--use-accession-id` parameter and retain "biologically correct" salmon haplotypes for all three salmonid species simultaneously.

```{code-block} bash
tombRaider --method 'taxon-dependent co-occurrence' --frequency-input asv_table.txt --taxonomy-input bast_taxonomy.txt --sequence-input asvs.fasta --frequency-output asv_table_new.txt --taxonomy-output blast_taxonomy_new.txt --sequence-output asvs_new.fasta --occurrence-type abundance --detailed-log Supplement_4_tombRaider_log.txt --use-accession-id --similarity 95
```

The *tombRaider* log file can be found in the supplemental file named "Supplement_4_tombRaider_log.txt".

### 4.5 Phylogenetic tree

To assess the efficiency of *tombRaider* to identify true haplotypes, we will generate a phylogenetic tree that includes all ASVs, as well as all reference sequences for our three salmonid species. By visualising the tree for which ASVs were retained and removed, we can investigate the accuracy of *tombRaider*. To start, we need to extract all reference barcodes from the reference database that are assigned to either of the three salmonid species and combine those reference barcodes with our ASVs in a single file. We can do this using the following python script.

```{code-block} python
seq_list = []
species = ['Oncorhynchus clarkii', 'Oncorhynchus kisutch', 'Oncorhynchus tshawytscha']
with open('asvs.fasta', 'r') as infile:
    for line in infile:
        seq_list.append(line)
with open('ncbi_ND2_insilico_tax_derep_clean.tsv', 'r') as reffile:
    for line in reffile:
        for item in species:
            if item in line:
                seqID = line.split('\t')[0] + '_' + line.split('\t')[-2] + '\n'
                sequence = line.split('\t')[-1]
                seq_list.append(seqID)
                seq_list.append(sequence)
with open('O_all_ref_seqs_plus_ASVs_renamed_headers_v2.fasta', 'w') as outfile:
    for item in seq_list:
        outfile.write(item)
```

The Bayesian phylogenetic tree was generated using BEAST *v* 2.7.6 [(Bouckaert et al., 2019)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006650) after aligning the sequences using the `AlignSeqs` function of the DECIPHER *v* 2.28.0 R package [(Wright, 2020)](https://rnajournal.cshlp.org/content/26/5/531.short).

```{code-block} R
## load libraries
library(DECIPHER)
library(dada2)
library(Biostrings)
library(ape)

## read fasta, align, and write output as nexus
sequenceTable <- readDNAStringSet('O_all_ref_seqs_plus_ASVs_renamed_headers_v2.fasta')
seqs <- getSequences(sequenceTable)
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)
write.nexus.data(alignment, 'O_all_ref_seqs_plus_ASVs_aligned_v2.nex')
```

Phylogenetic tree construction was performed with a Markov chain Monte Carlo (MCMC) chain length of 108 iterations, sampling trees every 1000. Convergence of the MCMC chains and effective sample size was checked using TRACER v 1.7.2 [(Rambaut et al., 2018)](https://academic.oup.com/sysbio/article/67/5/901/4989127). The maximum credibility tree from the posterior sample of phylogenetic time-trees with a burn-in percentage of 85 % was identified through TreeAnnotator v 2.7.6 [(Bouckaert et al., 2019)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006650) and used for subsequent analyses. The final tree was read into R and visualised using the ggtree *v* 3.10.0 R package.

```{code-block} R
# load libraries
required.libraries <- c("dada2", "DECIPHER", "purrr", "ape", "picante", 
                        "pez", "phytools",
                        "vegan", "adephylo", 
                        "phylobase", "geiger", 
                        "mvMORPH", "OUwie", 
                        "hisse", "BAMMtools",
                        "phylosignal", "Biostrings",
                        "devtools","ggplot2", 
                        "kableExtra", "betapart", "gridExtra",
                        "reshape2", "ggtree", "car", "egg", "tidyverse", "dplyr",
                        "hrbrthemes", "readxl", "ggrepel", "pracma", "scales", "ggpubr", "lsmeans", "multcomp",
                        "phyloseq", "gplots", "tidytree", "ggridges")
lapply(required.libraries, require, character.only = TRUE)

# First, read the data into R.
sequenceTable <- readDNAStringSet('O_all_ref_seqs_plus_ASVs_renamed_headers_v2.fasta')
phyloTree <- read.nexus('O_all_ref_seqs_plus_ASVs_aligned_v2-tree.tree')

# Second, plot tree and figure out the node numbers for manipulation
ggtree(phyloTree) + geom_text(aes(label = node)) + geom_tiplab()

# Third, collapse nodes I don't further need (node 148 for all clarkii ref seqs without ASVs; node 214 for outgroup)
rawTree <- ggtree(phyloTree) + geom_text(aes(label = node)) + geom_tiplab()
collapsedTree <- rawTree %>% collapse(node = 148) +
  geom_point2(aes(subset=(node==148)), shape = 21, size = 5, fill = 'green')
collapsedTree <- collapsedTree %>% collapse(node = 214) +
  geom_point2(aes(subset=(node==214)), shape = 21, size = 5, fill = 'red')
collapsedTree
```

The final tree was annotated manually and combined with the haplotype networks (see **4.6 Haplotype networks**) to generate Figure 5 in the manuscript.

### 4.6 Haplotype networks

To generate the haplotype networks for each salmonid species independently, we will generate fasta files from the initial csv supplementary files using the python script below.

```{code-block} python
#!/usr/bin/env Python3

fasta_dict = {}
with open('clarkii_ND2_combine_chimera_removed2.denoise_summary.csv', 'r') as clarkii_in:
    for line in clarkii_in:
        if line.startswith('id'):
            continue
        else:
            seq_id = '>clarkii_' + line.split(',')[0]
            sequence = line.split(',')[-1].rstrip('\n')
            fasta_dict[seq_id] = sequence
with open('clarkii_asvs.fasta', 'w') as outfile:
    for k, v in fasta_dict.items():
        asv_out.write(k + '\n' + v + '\n')

fasta_dict = {}
with open('kisutch_ND2_combine.denoise_summary.csv', 'r') as kisutch_in:
    for line in kisutch_in:
        line = line.replace('"', '')
        if 'id' in line.split(',')[0]:
            continue
        else:
            seq_id = '>kisutch_' + line.split(',')[0]
            sequence = line.split(',')[-1].rstrip('\n')
            fasta_dict[seq_id] = sequence
with open('kisutch_asvs.fasta', 'w') as outfile:
    for k, v in fasta_dict.items():
        asv_out.write(k + '\n' + v + '\n')

fasta_dict = {}
with open('tshawytscha_ND2_combine.denoise_summary.csv', 'r') as tshawytscha_in:
    for line in tshawytscha_in:
        line = line.replace('"', '')
        if 'id' in line.split(',')[0]:
            continue
        else:
            seq_id = '>tshawytscha_' + line.split(',')[0]
            sequence = line.split(',')[-1].rstrip('\n')
            fasta_dict[seq_id] = sequence
with open('tshawytscha_asvs.fasta', 'w') as outfile:
    for k, v in fasta_dict.items():
        asv_out.write(k + '\n' + v + '\n')
```

The separate ASV files can be read into R to create individual haplotype networks using the pegas *v* 1.3 and adegenet *v* 2.1.10 R packages.

```{code-block} R
library(pegas)
library(adegenet)
dna <- fasta2DNAbin('kisutch_asvs.fasta')
h <- haplotype(dna)
d <- dist.dna(h, 'N')
nt <- rmst(d, quiet = TRUE)
x <- haploNet(h, d)
plot(x, threshold = 0, size = 2, scale.ratio = 2)
pdf("haploNet_plot_kisutch.pdf")
plot(x, threshold = 0, size = 2, scale.ratio = 2)
dev.off()

dna <- fasta2DNAbin('clarkii_asvs.fasta')
h <- haplotype(dna)
d <- dist.dna(h, 'N')
nt <- rmst(d, quiet = TRUE)
x <- haploNet(h, d)
plot(x, threshold = 0, size = 2)
pdf("haploNet_plot_clarkii.pdf")
plot(x, threshold = 0, size = 2, scale.ratio = 2)
dev.off()

dna <- fasta2DNAbin('tshawytscha_asvs.fasta')
h <- haplotype(dna)
d <- dist.dna(h, 'N')
nt <- rmst(d, quiet = TRUE)
x <- haploNet(h, d)
plot(x, threshold = 0)
pdf("haploNet_plot_tshawytscha.pdf")
plot(x, threshold = 0, size = 2, scale.ratio = 2)
dev.off()
```

Individual haplotype networks were manually annotated and combined with the phylogenetic tree to generate Figure 5 in the manuscript.
