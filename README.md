# Huanan Market Metagenomic Data and Analyses

This repository provides data, sequences, and analyses of the metagenomic detection of mammalian wildlife at the Huanan Market in early 2020.

**Article sections:**

- (R1):  SARS-CoV-2 genetic diversity linked to the Huanan market is consistent with market emergence  
- (R2):  Increased High SARS-CoV-2 positivity rate in and near a wildlife stall in the Huanan market   
- (R3):  Mammalian wildlife species detected in five SARS-CoV-2 positive samples from a wildlife stall and in other wildlife stalls   
- (R4):  Wildlife stalls and SARS-CoV-2 positive samples contain other mammalian viruses associated with the animal trade  
- (R5):  Reconstruction of mitochondrial genotypes of potential intermediate host species of SARS-CoV-2 within the Huanan market  

## Contents by folder
-  `metadata/`: sample and species metadata (R)   
-  `mitochondrial_mappings/`: Read mapping counts and covered bases (breadth) for animal mtDNA. (R3, R5)
-  `mitochondrial_genotypes/`: Animal mtDNA consensus genotypes and SNV VCFs. (R3, R5)  
-  `sars2_mappings/`: identification and quantification of SARS-CoV-2 reads. (R1, R2)  
-  `sars2_phylogenetics/`: SARS-CoV-2 phylogeny.  (R1)  
-  `viruses`:  identification and quantification of other viruses in the market. (R4)
-  `correlations/`: computations of correlations between viruses and animals (R3)  
-  `fastq_processing/`: code for sample (pre)processing (R)     
-  `figs/`: directory where figure outputs are stored (R)     
-  `figures_R/`: R code to generate figures (R)    
-  `summary_and_tables/`: Supplementary tables (R)  

## Contents by topic

### Sequences and genotypes of animals and viruses found at the market

#### Animals

-  *Mammalian mtDNA consensus genomes*: `./mitochondrial_genotypes/lowconfidence_50_percent_mt_genomes.fasta` and `mediumconfidence_50_percent_mt_genomes.fasta`.
-  *Script to generate mtDNA consensus genomes*: `./mitochondrial_genotypes/mt_consensus_genomes.py`.
-  *BAMs and SNVs of animal mtDNA mappings*: `./mitochondrial_genotypes/bams_and_snvs/`.
-  *mtDNA dereplicated databases*: `./mitochondrial_mappings/mitochondria_dereplicated93_final.fasta.gz` and `./mitochondrial_mappings/mitochondria_dereplicated98.fasta.gz`.
-  *Scripts to count mapped reads and covered bases*: `./mitochondrial_mappings/count_reads98.py` and `./mitochondrial_mappings/count_reads93.py`.
-  *Mapping mtDNA counts and covered bases*: `./mitochondrial_mappings/mitochondrial_metazoa_counts_93.tsv` and `./mitochondrial_mappings/mitochondrial_metazoa_coveredbases_93.tsv`.

#### Viruses
-  *Viral sequences (bamboo rat CoV, raccoon dog amdovirus, civet kobuvirus, canine coronaviruses)*: `./viruses/market_mammalian_viruses.fna`  
-  *Influenza H9N2 fragmented genomes reconstructed from the market and trees:* `./viruses/Influenza/*.fasta` and `./viruses/Influenza/HA_tree.newick`.

### Sample Metadata

`./metadata/Liu_etal_2023_with_sequencing.csv`: Full metadata of all samples, linking samples to their primary metagenomic (untargeted) sequencing runs.

`./metadata/Sequencing_run_info.tsv`: Information on the 184 sequencing runs in Liu et al.

`./metadata/Liu_PCR_results.csv`: PCR CT values reported by Liu et al.

## Contents by Figure

### Figure 1

-  `sars2_phylogenetics/20210903.WHO.masked.fasta`: Sequences used in Fig 1A tree.  
-  `sars2_phylogenetics/HSM_sequences/*.fa`: consensus FASTA sequences for HSM SARS-CoV-2 genomes.   
-  `sars2_phylogenetics/combined_tMRCA.csv`: tMRCA results (Fig 1C)  
-  `R/Fig1C.R`: Script to plot tMRCA distributions (Fig 1C)  

### Figure 2: SARS2 abundances

-  `sars2_mappings/count_reads_sars2.py`: Script used to count the number of SARS-CoV-2 reads per sample.    
-  `sars2_mappings/sars2_reads_post_trimming.tsv`: counts of SARS-CoV-2 reads in each sample.  
- `R/Fig2ABC.R` and `R/Fig2DEFGH.R`: Scripts to plot the different panels of Figure 2.  

### Figure 2/3: mtDNA abundances

-  `mitochondrial_mappings/clustering98.py`: Script for first clustering the MT reference database at 98% identity.  
-  `mitochondrial_mappings/clustering93.py`: Script for clustering the MT reference database at 93% identity.  
-  `mitochondrial_mappings/count_reads98.py`: Script for counting mtDNA reads after mapping to the 98% clustered database.  
-  `mitochondrial_mappings/count_reads93.py`: Script for counting mtDNA reads after mapping to the 93% clustered database (final species-level assignments).  
-  `mitochondrial_mappings/mitochondria_dereplicated93_final.fasta.gz`: the final 93% ANI dereplicated database used in this study.  
-  `mitochondrial_mappings/mitochondrial_metazoa_counts_93.tsv`: Read counts for animal mtDNA.  
-  `mitochondrial_mappings/mitochondrial_metazoa_coveredbases_93.tsv`: Covered bases count for animal mtDNA  
-  `mitochondrial_mappings/species_descriptions_with_common_name.csv`: Descriptions of species detected in this study.  
-  `R/Fig3A.R` and `R/FigBCDE.R`: Scripts to plot the different panels of Figure 3.  

### Figure 4: viruses

-  `viruses/Influenza/*.fasta`: consensus Influenza H9N2 sequences identified in this study.  
-  `viruses/Influenza/HA_tree.newick`: Influenza H9N2 HA tree constructed.  
-  `viruses/count_reads_viral.py`: Script for filtering and counting reads mapped to the  viral database.  
-  `viruses/filtered_viral_counts_97_95_20_200.tsv`: counts of mapped viral reads.  
-  `viruses/market_mammalian_viruses.fna`: FASTA files of the consensus genomes assembled for mammalian viruses from the market.  
-  `viruses/plant_insect_phage_exclude.txt`: List of viruses to ignore due to being phage / plant / insect viruses.  
-  `viruses/viral_metadata.tsv`: Metadata on various viruses.  
-  `viruses/cluster_viral_references.py`: Script to cluster the viral reference database.  
-  `viruses/get_viral_consensus.py`: Script to reconstruct viral consensus genomes from multiple samples.  
-  `viruses/viral_db.list.txt`: a list of all viral sequences included in the dereplicated mapping database.   
-  `R/Fig4ABC.R`: Script to plot panels A-C of Figure 4.  

### Figure 5: mtDNA phylogenetics

-  `mitochondrial_genotypes/bams_and_snvs`: All mapped mammalian mtDNA BAMs from market samples.   
-  `mitochondrial_genotypes/MT_Maps.ipynb`: Script to create Figure 5a mtDNA genome maps.  
-  `mitochondrial_genotypes/mt_consensus_genomes.py`: Script to reconstruct mtDNA consensus genomes.  
-  `mitochondrial_genotypes/lowconfidence_50_percent_mt_genomes.fasta`: MT genomes obtained in this study, requiring just 1x read coverage.  
-  `mitochondrial_genotypes/mediumconfidence_50_percent_mt_genomes.fasta`: MT genomes obtained in this study, requiring 3x read coverage to call a base.  
-  `mitochondrial_genotypes/mtDNA_RD_phylogenetics.clean.ipynb`: Raccoon dog phylogenetics analysis notebook.  
-  `mitochondrial_genotypes/raccoon_dog_cytB_all.duplicated.fasta`: Raccon dog reference cytB alignment.  

## Supplementary: rRNA analysis

- `mitochondrial_mappings/rRNA/count_reads93_rRNAs.py`: Script to count reads excluding rRNA regions of the mitochondria.
-  `mitochondrial_mappings/rRNA/mt93_counts_rrna.tsv`: Read counts ignoring rRNA regions.  
-  `mitochondrial_mappings/rRNA/rRNA_analysis.ipynb`: Script to identify rRNA regions.  
-  `mitochondrial_mappings/rRNA/rrna_positions.tsv`: positions of rRNA genes in the MT genomes.  

## Supplementary: Correlations

-  `correlations/Correlation.ipynb`: Script to calculate correlations between animal sequence read abundances and SARS-CoV-2.

### He et al. 2022 Cell data: Animal data from other markets

`HeCell2022/`: location of these files  
  -  `He2022_metadata.csv`: Metadata about all of the files from that study  
  -  `cell2022_lowconfidence_mitochondria.fasta`: FASTA file of consensus mitochondrial genomes with at least 1 read covering each position.  
  -  `cell2022_mediumconfidence_mitochondria.fasta`: FASTA file of consensus mitochondrial genomes with at least 3 reads covering each position.  
  -  `scripts/mt_consensus_genomes_He2022.py`: script used to generate Mitochondrial genomes above   
  -  `game_animals.xlsx`: Additional metadata from He et al. 2022.  

### FASTQ processing

- `fastq_processing/run_bbduk.sh`: BBDuk commands used for quality filtering.
- `fastq_processing/frequent_adapter_summary.tsv` - Adapter and quality filtering summary per file.
`data/fastq_processing/`: data related to FASTQ files  
  -  ` all_adapters.fa `: list of adapters used for trimming  
  -  ` frequent_adapter_summary.tsv`: results of adapter trimming, 5 most common adapters per library.  
  -  ` run_bbduk.sh`: List of BBDuk commands run for each file for trimming.  
  -  ` bbduk_results/*.txt`: Full adapter trimming output for each file.  
