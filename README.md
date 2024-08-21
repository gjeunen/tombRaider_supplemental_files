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
