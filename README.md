# *tombRaider* supplemental files

## 1. Introduction

This document serves as a comprehensive guide and supplementary information for the bioinformatic and statistical analysis of the manuscript entitled "*tombRaider* - improved species and haplotype recovery from metabarcoding data through artefact and pseudogene exclusion" by Jeunen *et al*., 2024. This work is associated with the Marsden Fast-Start fund (MFP-UOO002116). This document is split up into three parts, one for each of the data sets reanalysed for the publication.

## 2. Supplement 2: mock community data analysis

### 2.1 Starting files

The raw Illumina sequencing data for this project [(Braukmann et al., 2019)](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13008) has been deposited on NCBI's Short Read Archive (SRP158933). To download the raw fastq files, we can make use of the `fastq-dump` function from the SRA toolkit. The deposited sequences are split up into the different mock communities, so we will have to execute the `fastq-dump` function multiple times, once for each fastq file to download.

```{code-block} bash
fastq-dump SRR7759458
fastq-dump SRR7759449
fastq-dump SRR7759450
fastq-dump SRR7759473
fastq-dump SRR7759476
fastq-dump SRR7759475
fastq-dump SRR7759446
fastq-dump SRR7759443
fastq-dump SRR7759444
```

The `fastq-dump` downloads a single fastq file, even though the manuscript mentioned paired-end sequencing. Upon inspection using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), the forward (R1) and reverse (R2) reads are concatenated. We can use [VSEARCH *v* 2.16.0](https://github.com/torognes/vsearch) to split the file into an R1 and R2 fastq file.

```{code-block} bash
vsearch --fastq_filter SRR7759458_bulk_abdomen_replicate_1.fastq --fastq_stripright 301 --fastqout SRR7759458_bulk_abdomen_replicate_1_R1.fastq
vsearch --fastq_filter SRR7759458_bulk_abdomen_replicate_1.fastq --fastq_stripleft 301 --fastqout SRR7759458_bulk_abdomen_replicate_1_R2.fastq

vsearch --fastq_filter SRR7759449_bulk_abdomen_replicate_2.fastq --fastq_stripright 301 --fastqout SRR7759449_bulk_abdomen_replicate_2_R1.fastq
vsearch --fastq_filter SRR7759449_bulk_abdomen_replicate_2.fastq --fastq_stripleft 301 --fastqout SRR7759449_bulk_abdomen_replicate_2_R2.fastq

vsearch --fastq_filter SRR7759450_bulk_abdomen_replicate_3.fastq --fastq_stripright 301 --fastqout SRR7759450_bulk_abdomen_replicate_3_R1.fastq
vsearch --fastq_filter SRR7759450_bulk_abdomen_replicate_3.fastq --fastq_stripleft 301 --fastqout SRR7759450_bulk_abdomen_replicate_3_R2.fastq

vsearch --fastq_filter SRR7759473_bulk_leg_replicate_1.fastq --fastq_stripright 301 --fastqout SRR7759473_bulk_leg_replicate_1_R1.fastq
vsearch --fastq_filter SRR7759473_bulk_leg_replicate_1.fastq --fastq_stripleft 301 --fastqout SRR7759473_bulk_leg_replicate_1_R2.fastq

vsearch --fastq_filter SRR7759476_bulk_leg_replicate_2.fastq --fastq_stripright 301 --fastqout SRR7759476_bulk_leg_replicate_2_R1.fastq
vsearch --fastq_filter SRR7759476_bulk_leg_replicate_2.fastq --fastq_stripleft 301 --fastqout SRR7759476_bulk_leg_replicate_2_R2.fastq

vsearch --fastq_filter SRR7759475_bulk_leg_replicate_3.fastq --fastq_stripright 301 --fastqout SRR7759475_bulk_leg_replicate_3_R1.fastq
vsearch --fastq_filter SRR7759475_bulk_leg_replicate_3.fastq --fastq_stripleft 301 --fastqout SRR7759475_bulk_leg_replicate_3_R2.fastq

vsearch --fastq_filter SRR7759446_composite_leg_replicate_1.fastq --fastq_stripright 301 --fastqout SRR7759446_composite_leg_replicate_1_R1.fastq
vsearch --fastq_filter SRR7759446_composite_leg_replicate_1.fastq --fastq_stripleft 301 --fastqout SRR7759446_composite_leg_replicate_1_R2.fastq

vsearch --fastq_filter SRR7759443_composite_leg_replicate_2.fastq --fastq_stripright 301 --fastqout SRR7759443_composite_leg_replicate_2_R1.fastq
vsearch --fastq_filter SRR7759443_composite_leg_replicate_2.fastq --fastq_stripleft 301 --fastqout SRR7759443_composite_leg_replicate_2_R2.fastq

vsearch --fastq_filter SRR7759444_composite_leg_replicate_3.fastq --fastq_stripright 301 --fastqout SRR7759444_composite_leg_replicate_3_R1.fastq
vsearch --fastq_filter SRR7759444_composite_leg_replicate_3.fastq --fastq_stripleft 301 --fastqout SRR7759444_composite_leg_replicate_3_R2.fastq
```

### 2.2 Bioinformatic analysis

To start the bioinformatic analysis, we merge the forward and reverse reads to generate the amplicons.

```{code-block} bash
vsearch --fastq_mergepairs SRR7759458_bulk_abdomen_replicate_1_R1.fastq --reverse SRR7759458_bulk_abdomen_replicate_1_R2.fastq --fastqout SRR7759458_bulk_abdomen_replicate_1_merged.fastq

vsearch --fastq_mergepairs SRR7759449_bulk_abdomen_replicate_2_R1.fastq --reverse SRR7759449_bulk_abdomen_replicate_2_R2.fastq --fastqout SRR7759449_bulk_abdomen_replicate_2_merged.fastq

vsearch --fastq_mergepairs SRR7759450_bulk_abdomen_replicate_3_R1.fastq --reverse SRR7759450_bulk_abdomen_replicate_3_R2.fastq --fastqout SRR7759450_bulk_abdomen_replicate_3_merged.fastq

vsearch --fastq_mergepairs SRR7759473_bulk_leg_replicate_1_R1.fastq --reverse SRR7759473_bulk_leg_replicate_1_R2.fastq --fastqout SRR7759473_bulk_leg_replicate_1_merged.fastq

vsearch --fastq_mergepairs SRR7759476_bulk_leg_replicate_2_R1.fastq --reverse SRR7759476_bulk_leg_replicate_2_R2.fastq --fastqout SRR7759476_bulk_leg_replicate_2_merged.fastq

vsearch --fastq_mergepairs SRR7759475_bulk_leg_replicate_3_R1.fastq --reverse SRR7759475_bulk_leg_replicate_3_R2.fastq --fastqout SRR7759475_bulk_leg_replicate_3_merged.fastq

vsearch --fastq_mergepairs SRR7759446_composite_leg_replicate_1_R1.fastq --reverse SRR7759446_composite_leg_replicate_1_R2.fastq --fastqout SRR7759446_composite_leg_replicate_1_merged.fastq

vsearch --fastq_mergepairs SRR7759443_composite_leg_replicate_2_R1.fastq --reverse SRR7759443_composite_leg_replicate_2_R2.fastq --fastqout SRR7759443_composite_leg_replicate_2_merged.fastq

vsearch --fastq_mergepairs SRR7759444_composite_leg_replicate_3_R1.fastq --reverse SRR7759444_composite_leg_replicate_3_R2.fastq --fastqout SRR7759444_composite_leg_replicate_3_merged.fastq
```

For ease, move files to a new directory.

```{code-block} bash
mkdir merged
cp *merged.fastq merged/
cd merged/
```

Thus far, reads still contain the primer sequences. They can now be removed using [cutadapt *v* 4.4](https://cutadapt.readthedocs.io/en/stable/).

```{code-block} bash
mkdir demux
cutadapt SRR7759458_bulk_abdomen_replicate_1_merged.fastq -g 'GCTTTCCCACGAATAAATAATA...TGATTTTTTGGWCAYCCWGAAGTTTA' -o demux/SRR7759458_bulk_abdomen_replicate_1_merged_demux.fastq --discard-untrimmed --no-indels -e 2 --cores=0

cutadapt SRR7759449_bulk_abdomen_replicate_2_merged.fastq -g 'GCTTTCCCACGAATAAATAATA...TGATTTTTTGGWCAYCCWGAAGTTTA' -o demux/SRR7759449_bulk_abdomen_replicate_2_merged_demux.fastq --discard-untrimmed --no-indels -e 2 --cores=0

cutadapt SRR7759450_bulk_abdomen_replicate_3_merged.fastq -g 'GCTTTCCCACGAATAAATAATA...TGATTTTTTGGWCAYCCWGAAGTTTA' -o demux/SRR7759450_bulk_abdomen_replicate_3_merged_demux.fastq --discard-untrimmed --no-indels -e 2 --cores=0

cutadapt SRR7759473_bulk_leg_replicate_1_merged.fastq -g 'GCTTTCCCACGAATAAATAATA...TGATTTTTTGGWCAYCCWGAAGTTTA' -o demux/SRR7759473_bulk_leg_replicate_1_merged_demux.fastq --discard-untrimmed --no-indels -e 2 --cores=0

cutadapt SRR7759476_bulk_leg_replicate_2_merged.fastq -g 'GCTTTCCCACGAATAAATAATA...TGATTTTTTGGWCAYCCWGAAGTTTA' -o demux/SRR7759476_bulk_leg_replicate_2_merged_demux.fastq --discard-untrimmed --no-indels -e 2 --cores=0

cutadapt SRR7759475_bulk_leg_replicate_3_merged.fastq -g 'GCTTTCCCACGAATAAATAATA...TGATTTTTTGGWCAYCCWGAAGTTTA' -o demux/SRR7759475_bulk_leg_replicate_3_merged_demux.fastq --discard-untrimmed --no-indels -e 2 --cores=0

cutadapt SRR7759446_composite_leg_replicate_1_merged.fastq -g 'GCTTTCCCACGAATAAATAATA...TGATTTTTTGGWCAYCCWGAAGTTTA' -o demux/SRR7759446_composite_leg_replicate_1_merged_demux.fastq --discard-untrimmed --no-indels -e 2 --cores=0

cutadapt SRR7759443_composite_leg_replicate_2_merged.fastq -g 'GCTTTCCCACGAATAAATAATA...TGATTTTTTGGWCAYCCWGAAGTTTA' -o demux/SRR7759443_composite_leg_replicate_2_merged_demux.fastq --discard-untrimmed --no-indels -e 2 --cores=0

cutadapt SRR7759444_composite_leg_replicate_3_merged.fastq -g 'GCTTTCCCACGAATAAATAATA...TGATTTTTTGGWCAYCCWGAAGTTTA' -o demux/SRR7759444_composite_leg_replicate_3_merged_demux.fastq --discard-untrimmed --no-indels -e 2 --cores=0
```

Once the primer sequences are removed, we will further clean the data based on minimum length, maximum length, and number of ambiguous basecalls.

```{code-block} bash
cd demux/
mkdir qual
for fq in *.fastq
do
echo "\n\n\nAnalysing: ${fq}"
fasta=${fq/.fastq/.fasta} 
vsearch --fastq_filter ${fq} --fastq_maxee 1.0 --fastq_minlen 400 --fastq_maxlen 415 --fastq_maxns 0 --fastqout qual/${fq} --fastaout qual/${fasta} --relabel ${fq/.fastq/}.
done
```

To save computational resources, dereplicate reads.

```{code-block} bash
cd qual/
mkdir derep
cat *.fasta > derep/combined.fasta
cat *.fastq > derep/combined.fastq
cd derep/
vsearch --derep_fulllength combined.fasta --sizeout --relabel uniq. --output uniques.fasta
```

From the dereplicated fasta file, find all ASVs.

```{code-block} bash
vsearch --cluster_unoise uniques.fasta --sizein --sizeout --relabel denoised. --centroids denoised.fasta
```

VSEARCH does not automatically remove chimeric sequences during denoising, so let's remove those sequences now.

```{code-block} bash
vsearch --uchime3_denovo denoised.fasta --sizein --nonchimeras asv.fasta --relabel asv.
```

Now that we have the final list of ASVs, we can generate the ASV table.

```{code-block} bash
vsearch --usearch_global combined.fasta --db asv.fasta --strand plus --id 0.97 --otutabout asv_table.txt
```

### 2.3 Taxonomic assignment

Since the amplicons are within the COI gene region, we can use BOLDigger [Buchner and Leese, 2020](https://mbmg.pensoft.net/article/53535/) to assign a taxonomic ID to each ASV.

```{code-block} bash
mkdir BOLD
boldigger-cline ie_coi username password asv.fasta BOLD/
```

This will generate a file named "BOLDResults_asv_part_1.txt". The following python script can be used to reformat this file into a standard taxonomy output file, which will be used during subsequent analyses.

```{code-block} python
#! /usr/bin/env python3

## import modules
import collections

boldDict = collections.defaultdict(dict)
speciesListDict = collections.defaultdict(list)

seqID = 0
with open('BOLDResults_zotus2line_part_1.txt', 'r') as boldinfile:
    next(boldinfile)
    for line in boldinfile:
        line = line.rstrip('\n')
        lastSeqID = seqID
        seqID = line.split('\t')[0].lstrip('>')
        if seqID == '':
            seqID = lastSeqID
        taxID = ", ".join([item for item in line.split('\t')[1:8] if item != ''])
        similarity = line.split('\t')[8]
        if similarity != '':
            similarity = float(similarity)
        status = line.split('\t')[9]
        processID = line.split('\t')[10]
        if seqID not in boldDict:
            boldDict[seqID]['taxInfo'] = taxID
            boldDict[seqID]['score'] = similarity
            boldDict[seqID]['status'] = status
            boldDict[seqID]['processID'] = processID
        else:
            if boldDict[seqID]['score'] <= similarity and boldDict[seqID]['taxInfo'] != taxID and taxID not in speciesListDict[seqID]:
                speciesListDict[seqID].append(taxID)

zotuList = []
with open('asv_table.txt', 'r') as zotufile:
    for line in zotufile:
        seqID = line.split('\t')[0]
        if seqID != 'ID':
            zotuList.append(seqID)

boldFinalDict = collections.defaultdict(dict)
for item in zotuList:
    if item not in boldDict:
        boldFinalDict[item]['taxInfo'] = 'NA'
        boldFinalDict[item]['score'] = 'NA'
        boldFinalDict[item]['status'] = 'NA'
        boldFinalDict[item]['processID'] = 'NA'
    else:
        boldFinalDict[item] = boldDict[item]

with open('boldFormattedNew.txt', 'w') as outfile:
    outfile.write(f'#OTU ID\ttaxInfo\tsimilarity\tstatus\tprocessID\totherID\n')
    for item in boldFinalDict:
        outfile.write(f'{item}\t{boldFinalDict[item]["taxInfo"]}\t{boldFinalDict[item]["score"]}\t{boldFinalDict[item]["status"]}\t{boldFinalDict[item]["processID"]}\t{"; ".join(speciesListDict[item])}\n')
```

### 2.4 *tombRaider*

Now that we have generated the three input files needed for *tombRaider*, we can create all the separate ASV tables and ASV sequence files for each data curation method, including:

1. unfiltered (no data curation post denoising): these will be the files we have generated during the bioinformatic analysis
2. taxon-dependent co-occurrence merging (criteria: taxonomic ID, sequence similarity, co-occurrence pattern)
3. taxonomic ID merging (criteria: taxonomic ID)
4. taxon-independent co-occurrence merging (criteria: sequence similarity, co-occurrence pattern; i.e., LULU)
5. abundance filtering (will be done outside of *tombRaider*)
6. mock community (files can be downloaded from the Supplementary Information associated with the original publication)

#### 2.4.1 unfiltered

Data files created during bioinformatic analysis.

#### 2.4.2 taxon-dependent co-occurrence merging

```{code-block} bash
tombRaider --method 'taxon-dependent co-occurrence' --frequency-input zotutable.txt --frequency-output tombRaider_taxon_dependent_BOLD/zotutableNew.txt --sequence-input zotus_mock_community.fasta --sequence-output tombRaider_taxon_dependent_BOLD/zotus_mock_communityNew.fasta --bold-input boldigger/BOLDResults_zotus2line_part_1.txt --bold-output tombRaider_taxon_dependent_BOLD/bold-outputNew.txt --occurrence-type abundance --count 1 --sort 'total read count' --condensed-log tombRaider_taxon_dependent_BOLD/condensedLog.txt --detailed-log tombRaider_taxon_dependent_BOLD/detailedLog.txt --bold-format complete --similarity 88
```

#### 2.4.3 taxonomic ID merging

```{code-block} bash
tombRaider --method 'taxon-dependent merging' --frequency-input zotutable.txt --frequency-output tombRaider_taxon_merging/zotutableNew.txt --sequence-input zotus_mock_community.fasta --sequence-output tombRaider_taxon_merging/zotus_mock_communityNew.fasta --bold-input boldigger/BOLDResults_zotus2line_part_1.txt --bold-output tombRaider_taxon_merging/bold-outputNew.txt --occurrence-type abundance --count 0 --sort 'total read count' --condensed-log tombRaider_taxon_merging/condensedLog.txt --detailed-log tombRaider_taxon_merging/detailedLog.txt --bold-format complete
```

#### 2.4.4 taxon-independent co-occurrence merging

```{code-block} bash
tombRaider --method 'taxon-independent co-occurrence' --frequency-input zotutable.txt --frequency-output tombRaider_taxon_independent_BOLD/zotutableNew.txt --sequence-input zotus_mock_community.fasta --sequence-output tombRaider_taxon_independent_BOLD/zotus_mock_communityNew.fasta --bold-input boldigger/BOLDResults_zotus2line_part_1.txt --bold-output tombRaider_taxon_independent_BOLD/bold-outputNew.txt --occurrence-type abundance --count 1 --sort 'total read count' --condensed-log tombRaider_taxon_independent_BOLD/condensedLog.txt --detailed-log tombRaider_taxon_independent_BOLD/detailedLog.txt --bold-format complete --similarity 84
```

#### 2.4.5 abundance filtering

Abundance filtering is currently not incorporated in *tombRaider*, as many different approaches exist. Hence, we will create these files in excel. Open the ASV table and set reads to 0 if it does not contain at least 0.01% of the total reads within the sample. Finally, remove all ASVs that are dropped and update the ASV sequence file accordingly.

#### 2.4.6 mock community

The mock community information can be found in the Supplementary Information files associated with the original publication (file name: "men13008-sup-0009-TableS1.xlsx").

### 2.5 Statistical analysis

Below is the R code used for all statistical analyses and the creation of Figure 3 in the manuscript.

```{code-block} R
###################################
## 0. set up working environment ##
###################################
## set working directory
## load libraries
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

############################
## 1. GENERATE PHYLO TREE ##
############################
# From the bioinformatic analysis, we ended up with 6 ASV tables from the different filtering techniques,
# including the (i) raw table, (ii) taxon merge table, (iii) taxon independent co-occurrence table, 
# (iv) taxon dependent co-occurrence table, (v) abundance table, and (vi) mock community table.
# We also ended up with 2 ASV fasta files, including one with the sequences from the metabarcoding data set and
# another one where we supplemented the missed species from the mock community 
# (sequences obtained from BOLD, in silico PCR analysis to retain the amplicon).
# To generate the phylogenetic tree, we will use the ASV fasta file with the supplemented sequences from the missed mock community species.

# First, we need to align the sequences. We will do this in R using the DECIPHER package.
# a) load sequences into R using Biostrings and format appropriately using dada2
sequenceTable <- readDNAStringSet('../phyloTree2/zotusWithMissing.fasta')
seqs <- getSequences(sequenceTable)
# b) align sequences using DECIPHER
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)
# c) visualise alignment to check it worked properly
alignment[[1]]
# d) export alignment to nexus format using ape
write.nexus.data(alignment, '../phyloTree2/zotusWithMissingAlign.nex')

# Second, we need to generate the phylogenetic tree, which we will do in BEAST2.
# The next steps are not conducted in R.
# a) generate the .xml file using BEAUTi 2. 
#       a.1) Import nexus alignment
#       a.2) Set substitution model to HKY
#       a.3) Set starting tree to Cluster Tree and cluster type to UPGMA
#       a.4) Set Chain length to 10^8 and log every 10,000 trees
#       a.5) Keep remaining settings on default and export .xml file
# b) run BEAST 2 to generate phylogenetic trees.
# c) after BEAST 2 completes, check log file in Tracer.
# d) find best supported tree using TreeAnnotator and export to .tree file.
# e) visualise tree quality using FigTree.

# We have now created the phylogenetic tree for this data set.

###########################
## 2. COMBINE OTU TABLES ##
###########################
# After generating the phylogenetic tree, we need to combine all 6 ASV tables, so that we can generate the figures
# and compare the different filtering approaches.
# First, we need to read the data into R
rawTable <- read.table('../zotuTable.txt', header = TRUE, row.names = 1, sep = '\t', check.names = FALSE, comment.char = '')
mergeTable <- read.table('../tombRaider_taxon_merging/zotuTableNew.txt', header = TRUE, row.names = 1, sep = '\t', check.names = FALSE, comment.char = '')
independentTable <- read.table('../tombRaider_taxon_independent/zotutableTaxonIndependent.txt', header = TRUE, row.names = 1, sep = '\t', check.names = FALSE, comment.char = '')
dependentTable <- read.table('../tombRaider_taxon_dependent_BOLD/zotuTableNew.txt', header = TRUE, row.names = 1, sep = '\t', check.names = FALSE, comment.char = '')
abundanceTable <- read.table('../zotuTableAbundance.txt', header = TRUE, row.names = 1, sep = '\t', check.names = FALSE, comment.char = '')
mockTable<- read.table('../zotuTableMock.txt', header = TRUE, row.names = 1, sep = '\t', check.names = FALSE, comment.char = '')

# Second, we will change header names for ease
rawNames <- c("raw_CL2", "raw_CL3", "raw_CL1", "raw_BA2", "raw_BA3", "raw_BA1", "raw_BL1", "raw_BL3", "raw_BL2")
colnames(rawTable) <- rawNames
mergeNames <- c("merge_CL2", "merge_CL3", "merge_CL1", "merge_BA2", "merge_BA3", "merge_BA1", "merge_BL1", "merge_BL3", "merge_BL2")
colnames(mergeTable) <- mergeNames
independentNames <- c("independent_CL2", "independent_CL3", "independent_CL1", "independent_BA2", "independent_BA3", "independent_BA1", "independent_BL1", "independent_BL3", "independent_BL2")
colnames(independentTable) <- independentNames
dependentNames <- c("dependent_CL2", "dependent_CL3", "dependent_CL1", "dependent_BA2", "dependent_BA3", "dependent_BA1", "dependent_BL1", "dependent_BL3", "dependent_BL2")
colnames(dependentTable) <- dependentNames
abundanceNames <- c("abundance_CL2", "abundance_CL3", "abundance_CL1", "abundance_BA2", "abundance_BA3", "abundance_BA1", "abundance_BL1", "abundance_BL3", "abundance_BL2")
colnames(abundanceTable) <- abundanceNames
mockNames <- c("mock_CL2", "mock_CL3", "mock_CL1", "mock_BA2", "mock_BA3", "mock_BA1", "mock_BL1", "mock_BL3", "mock_BL2")
colnames(mockTable) <- mockNames

# Third, convert row names to a column in each dataframe
df_list <- list(rawTable, abundanceTable, independentTable, mergeTable, dependentTable, mockTable)
df_list <- lapply(df_list, function(df) {
  df$category <- rownames(df)
  return(df)
})

# Fourth, perform a full join on the new column to combine all dataframes using tidyverse
merged_df <- df_list %>%
  reduce(full_join, by = 'category')

# Fifth, set row names back, remove the extra column, set NA to 0, and transform data to pa
rownames(merged_df) <- merged_df$category
merged_df <- select(merged_df, -category)
merged_df[is.na(merged_df)] <- 0
merged_df[merged_df > 1] <- 1

# Sixth, write new dataframe to txt file
write.table(merged_df, file = 'zotuTableCombinedAll.txt', sep = '\t', row.names = TRUE)

#############################
## 3. FIGURE 1A PHYLO TREE ##
#############################
# Now that we have all files generated, we can start creating the figures. 
# The first subplot is the phylogenetic tree.
# First, read the data into R.
phyloTreeLongTable <- read.table('zotuTableCombinedAll.txt', header = TRUE, row.names = 1, sep = '\t', check.names = FALSE, comment.char = '')
phyloTree <- read.nexus('../phyloTree2/zotusWithMissing-tree.tree')

# Before creating the tree figure, we need to reformat the table.
# We need to group the samples by filter method, set 0 to NA, and change detection to header name for the visualization.
# Second, reformat phyloTreeLongTable according to figure specifications
column_names <- unique(gsub("^(.*?)_.*", "\\1", names(phyloTreeLongTable)))
phyloTreeShortTable <- data.frame(matrix(NA, nrow = nrow(phyloTreeLongTable), ncol = length(column_names)))
names(phyloTreeShortTable) <- column_names
for (col_name in column_names) {
  cols_to_combine <- grep(paste0("^", col_name, "_"), names(phyloTreeLongTable), value = TRUE)
  phyloTreeShortTable[[col_name]] <- rowSums(phyloTreeLongTable[cols_to_combine])
}
rownames(phyloTreeShortTable) <- rownames(phyloTreeLongTable)
phyloTreeShortTable[phyloTreeShortTable == 0] <- NA
for (col in column_names) {
  phyloTreeShortTable[[col]] <- ifelse(phyloTreeShortTable[[col]] > 0, col, phyloTreeShortTable[[col]])
}

# Third, reorder columns
phyloTreeShortTable <- phyloTreeShortTable %>%
  select(raw, abundance, independent, merge, dependent, mock, everything())

# Fourth, set colours for the figure
sample_colors <- c("raw" = "#c4d8e1", "merge" = "#BFB8DA", "independent" = "#4e6c82", "dependent" = "#f2d379", "abundance" = "#90adbf", "mock" = "grey30")

# Fifth, generate phylogenetic tree figure
circularTree <-ggtree(phyloTree, layout = "fan", open.angle = 20, size = 0.05, color = 'grey30')
circularTree <- rotate_tree(circularTree, 280)
phyloTreePlot <- gheatmap(circularTree, phyloTreeShortTable, width = 0.5, offset = -0.018, color = NA, hjust = 0.5, font.size = 3, custom_column_labels = c('raw', 'abundance', 'independent', 'merge', 'dependent', 'mock')) + 
  scale_fill_manual(values = sample_colors, na.value = NA) + 
  theme(legend.position = 'none') 
phyloTreePlot

#################################
## 4. FIGURE 1B COUNT BAR PLOT ##
#################################
# For this figure, we will use a file generated in excel, whereby we counted the number of instances for each type, including
# expected = sequences in the mock community,
# missing = sequences in the mock community but not detected,
# contaminant = sequences that are not in the mock community but match 100% to another species,
# artefact = sequences not in the mock community and not matching 100% to another species
# First, read the data into R
countTable <- read.table('count_bar_plot.txt', header = TRUE, row.names = 1, sep = '\t', check.names = FALSE, comment.char = '')

# Second, reshape data for plotting
countTable$category <- rownames(countTable)
countTable <- tidyr::pivot_longer(countTable, cols = -category, names_to = "Variable", values_to = "Value")

# Third, reorder variables for plotting
countTable$Variable <- factor(countTable$Variable, levels = c("raw", "abundance", "independent", "merge", "dependent", "mock"))
countTable$category <- factor(countTable$category, levels = c("Artefact", "Contaminant", "Missing", "Expected"))

# Fourth, set colours for plotting
type_colours <- c("Artefact" = "#ce7c7d", "Contaminant" = "#e7dfbb", "Missing" = "white", "Expected" = "grey40")

# Fifth, plot using ggplot2
barPlot <- ggplot(countTable, aes(x = Value, y = Variable, fill = category)) +
  geom_bar(stat = "identity") +
  labs(title = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_blank(),
        legend.position = 'none') +
  scale_fill_manual(values = type_colours, na.value = NA) + 
  scale_x_continuous(breaks = seq(0, 900, by = 450)) +
  expand_limits(x = c(0, 900))
barPlot

#################################
## 5. FIGURE 1C RIDGELINE PLOT ##
#################################
# For the ridgeline plot, we will again make use of a file generated in excel, whereby we listed the similarity score
# for each filtering method and split it up in kept sequences and filtered out sequences.
# First, read the data into R.
violinTable <- read.table('violin_plot_data.txt', header = TRUE, row.names = 1, sep = '\t', check.names = FALSE, comment.char = '')
violinMetadata <- read.table('violin_plot_metadata.txt', header = TRUE, row.names = 1, sep = '\t', check.names = FALSE, comment.char = '')

# Second, change 0 to NA
violinTable[violinTable == 0] <- NA

# Third, reshape to long format
violinTableLong <- gather(violinTable, key = "variable", value = "value")

# Fourth, join table with metadata
violinTableLong <- merge(violinTableLong, violinMetadata, by.x = "variable", by.y = "row.names")

# Fifth, reorder variables for plotting
violinTableLong$group <- factor(violinTableLong$group, levels = c("raw", "abundance", "independent", "merge", "dependent", "mock"))

# Sixth, plot using ggplot2
ridgePlot <- ggplot(violinTableLong, aes(x = value, y = group, fill = type, point_color = type, color = type)) +
  geom_density_ridges(
    jittered_points = TRUE, scale = 0.95, rel_min_height = 0.01, point_shape = '|', point_size = 3, size = 0.25, position = position_points_jitter(height = 0)
  ) +
  scale_x_continuous(expand = c(0, 0), name = "percent identity (%)") +
  scale_fill_manual(values = c("#e7dfbb50", "#69b3a250"), labels = c("discarded", "included")) +
  scale_color_manual(values = c("#e7dfbb", "#69b3a2"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#e7dfbb", "#69b3a2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(
    override.aes = list(
      fill = c("#e7dfbbA0", "#69b3a2A0"),
      color = NA, point_color = NA))
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_blank(), legend.position = 'none',
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  scale_x_continuous(breaks = seq(70, 100, by = 15)) +
  expand_limits(x = c(70, 100))
ridgePlot

###############################################
## 6. FIGURE 1D - 1F ALPHA DIVERSITY BOXPLOT ##
###############################################
# For the alpha diversity analysis, we will plot species richness and Faith's PD.
alphaTable <- read.table('zotuTableCombinedAll.txt', header = TRUE, row.names = 1, sep = '\t', check.names = FALSE, comment.char = '')
phyloTree <- read.nexus("../phyloTree2/zotusWithMissing-tree.tree")
alphaMetaTable <- read.table('sampleInfo.txt', header = TRUE, sep = '\t', check.names = FALSE, comment.char = '')

# First, prepare data for analysis
is.binary(phyloTree)
is.ultrametric(phyloTree)
alphaTable.T <- t(alphaTable)
all(sort(phyloTree$tip.label) == sort(colnames(alphaTable.T)))
all(phyloTree$tip.label == colnames(alphaTable.T))
phyloTree.clean <- match.phylo.comm(phy = phyloTree, comm = alphaTable.T)$phy
alphaTable.clean <- match.phylo.comm(phy = phyloTree, comm = alphaTable.T)$comm
all(phyloTree.clean$tip.label == colnames(alphaTable.clean))
plot(phyloTree.clean, cex = 0.4)

# Second, calculate Faith's PD and species richness (SR)
alpha.PD <- pd(samp = alphaTable.clean, tree = phyloTree.clean, include.root = FALSE)
cor.test(alpha.PD$PD, alpha.PD$SR)
plot(alpha.PD$PD, alpha.PD$SR, xlab = 'Phylogenetic Diversity', ylab = 'Species Richness', pch = 16)

# Third, combine metadata with Faith's PD and SR
alpha.PD$ID <- rownames(alpha.PD)
alpha.PD.meta <- merge(alpha.PD, alphaMetaTable, by = 'ID')

# Fourth, plot boxplot for PD and SR
alpha.PD.meta$type <- factor(alpha.PD.meta$type, levels = c("raw", "abundance", "independent", "merge", "dependent", "mock"))
sample_colors <- c("raw" = "#c4d8e1", "merge" = "#BFB8DA", "independent" = "#4e6c82", "dependent" = "#f2d379", "abundance" = "#90adbf", "mock" = "grey30")
boxplot.PD <- ggplot(alpha.PD.meta, aes(x = type, y = PD, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_point(shape = 1, position = position_jitterdodge(jitter.width = 0.1), size = 2, colour = 'black') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = 'none', axis.title.x = element_blank()) +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(breaks = seq(15, 45, by = 15)) +
  expand_limits(y = c(15, 45))
boxplot.PD

boxplot.SR <- ggplot(alpha.PD.meta, aes(x = SR, y = type, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_point(shape = 1, position = position_jitterdodge(jitter.width = 0.1), size = 2, colour = 'black') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = 'none', axis.title.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  scale_fill_manual(values = sample_colors) +
  scale_x_continuous(breaks = seq(100, 700, by = 300)) +
  expand_limits(x = c(100, 700))
boxplot.SR

# Fifth, run ANOVA for PD and SR
anova.PD <- aov(PD ~ as.factor(type), data = alpha.PD.meta)
summary(anova.PD)
TukeyHSD(anova.PD)

anova.SR <- aov(SR ~ as.factor(type), data = alpha.PD.meta)
summary(anova.SR)
TukeyHSD(anova.SR)

# Sixth, check assumptions of ANOVA
# assumption 1: normal distribution of residuals
hist(anova.PD$residuals, main = 'Histogram of PD residuals', xlab = 'Residuals')
hist(anova.SR$residuals, main = 'Histogram of SR residuals', xlab = 'Residuals')
# assumption 2: homogeneity of variance
leveneTest(PD ~ as.factor(type), data = alpha.PD.meta)
leveneTest(SR ~ as.factor(type), data = alpha.PD.meta)

# Seventh, significant levene's test, so need to run Welch's ANOVA instead
oneway.test(PD ~ as.factor(type), data = alpha.PD.meta, var.equal = FALSE)
oneway.test(SR ~ as.factor(type), data = alpha.PD.meta, var.equal = FALSE)
modelPD <- lm(PD ~ as.factor(type), data= alpha.PD.meta)
modelSR <- lm(SR ~ as.factor(type), data= alpha.PD.meta)
leastsquarePD <- lsmeans(modelPD, pairwise ~ type, adjust = 'tukey')
cld(leastsquarePD, alpha = 0.05, Letters = letters, adjust = 'tukey')
leastsquareSR <- lsmeans(modelSR, pairwise ~ type, adjust = 'tukey')
cld(leastsquareSR, alpha = 0.05, Letters = letters, adjust = 'tukey')

#################################################
## 7. FIGURE 1E - 1G BETA DIVERSITY ORDINATION ##
#################################################
# For the beta diversity analysis, we will both analyse Taxonomic Diversity (TD) and Phylogenetic Diversity (PD)
# First, read the data into R
betaTable <- as.data.frame(t(read.table('zotuTableCombinedAll.txt', header = TRUE, row.names = 1, sep = '\t', check.names = FALSE, comment.char = '')))
phyloTree <- read.nexus("../phyloTree2/zotusWithMissing-tree.tree")
betaMetaTable <- read.table('sampleInfo.txt', header = TRUE, sep = '\t', check.names = FALSE, comment.char = '')

# Second, import data into phyloseq
OTU <- otu_table(betaTable, taxa_are_rows = FALSE)
betaMetaTable$type <- as.factor(betaMetaTable$type)
META <- sample_data(betaMetaTable)
rownames(META) <- META$ID
physeq <- merge_phyloseq(OTU, META, phyloTree)

# Third, generate ordination plot TD
sample_colors <- c("raw" = "#c4d8e1", "merge" = "#BFB8DA", "independent" = "#4e6c82", "dependent" = "#f2d379", "abundance" = "#90adbf", "mock" = "grey30")
TD.ord.PCoA <- ordinate(physeq, method = "PCoA", distance = 'jaccard', try = 100, trymax = 1000, k = 2)
TD.ord.PCoA.plot <- plot_ordination(physeq, TD.ord.PCoA, type = 'samples', color = 'type', shape = 'type') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = 'black', fill = NA), legend.position = 'none')  +
  geom_point(size = 5) +
  scale_color_manual(values = sample_colors)
TD.ord.PCoA.plot

# Fourth, generate ordination plot PD
PD.ord.PCoA <- ordinate(physeq, method = 'PCoA', distance = 'unifrac', weighted = FALSE)
PD.ord.PCoA.plot <- plot_ordination(physeq, PD.ord.PCoA, type = 'samples', color = 'type', shape = 'type') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = 'black', fill = NA), legend.position = 'none')  +
  geom_point(size = 5) +
  scale_color_manual(values = sample_colors)
PD.ord.PCoA.plot

ggarrange(phyloTreePlot, ggarrange(ggarrange(barPlot, ridgePlot, ncol = 2, labels = c('b', 'c')), ggarrange(boxplot.PD, PD.ord.PCoA.plot, ncol = 2, labels = c('d', 'e')), nrow = 2, heights = c(1, 1.5)), ncol = 2, labels = c('a', ''), widths = c(1, 2))
```

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
