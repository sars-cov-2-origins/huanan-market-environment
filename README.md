# huanan-env-paper-private
For data and figures working towards the final paper

### Sequences: quick links

- *Viral sequences (bamboo rat CoV, raccoon dog amdovirus, civet kobuvirus)*: `./viruses/market_mammalian_viruses.fna`
- *Influenza H9N2 fragmented genome reconstructed from the market:* `./viruses/Influenza.fasta`
- *mtDNA consensus genomes*: `./MT_phylogenetics/lowconfidence_50_percent_mt_genomes.fasta`
- *BAMs of mtDNA mappings*: `./MT_phylogenetics/bams_and_snvs/`


### Sample Metadata


`./metadata/Liu_etal_2023_market_samples.tsv` - Full metadata of all samples, including estimated geolocations.

`./metadata/Liu_etal_2023_with_sequencing.csv` - Full metadata of all samples, linking samples to their primary metagenomic (untargeted) sequencing runs.

`./metadata/Sequencing_run_info.tsv` - Information on the 184 sequencing runes in Liu et al.

`./metadata/Liu_PCR_results.csv` - PCR CT values reported by Liu et al.

## By Figure

### Figure 1

`SARS2/20210903.WHO.masked.fasta` - Sequences used in Fig 1 tree.

`SARS2/Fig1_tMRCA.ipynb` - plotting the tMRCA estimates.

`SARS2/HSM_sequences/*.fa` - consensus FASTA sequences for HSM SARS-CoV-2 genomes.

### Figure 2: SARS2 abundances

`./SARS2/count_reads_sars2.py` - Script used to count the number of SARS-CoV-2 reads per sample.

`./SARS2/sars2_reads_post_trimming.tsv` - counts of SARS-CoV-2 reads in each sample.

`./R/Fig2ABC.R` and `./R/Fig2DEFGH.R` - Scripts to plot the different panels of Figure 2.

### Figure 2/3: mtDNA abundances

`./mtDNA/clustering98.py` - Script for first clustering the MT reference database at 98% identity.

`./mtDNA/clustering93.py` - Script for clustering the MT reference database at 93% identity.

`./mtDNA/count_reads98.py` - Script for counting mtDNA reads after mapping to the 98% clustered database.

`./mtDNA/count_reads93.py` - Script for counting mtDNA reads after mapping to the 93% clustered database (final species-level assignments).

`./mtDNA/mitochondria_dereplicated93_final.fasta.gz` - the final 93% ANI dereplicated database used in this study.

`./mtDNA/mitochondrial_metazoa_counts_93.tsv` - Read counts for animal mtDNA.

`./mtDNA/mitochondrial_metazoa_coveredbases_93.tsv` - Covered bases count for animal mtDNA

`./mtDNA/species_descriptions_with_common_name.csv` - Descriptions of species detected in this study.

`./R/Fig3A` and `./R/FigBCDE.R` - Scripts to plot the different panels of Figure 3.

### Figure 4: viruses

`./viruses/Influenza.fasta` - consensus Influenza H9N2 sequences identified in this study.

`./viruses/count_reads_viral.py` - Script for filtering and counting reads mapped to the viral database.

`./viruses/filtered_viral_counts_97_95_20_200.tsv` - counts of mapped viral reads.

`./viruses/market_mammalian_viruses.fna` - FASTA files of the consensus genomes assembled for mammalian viruses from the market.

`./viruses/plant_insect_phage_exclude.txt` - List of viruses to ignore due to being phage / plant / insect viruses.

`./viruses/viral_metadata.tsv` - Metadata on various viruses.

`./viruses/cluster_viral_references.py` - Script to cluster the viral reference database.

`./viruses/get_viral_consensus.py` - Script to reconstruct viral consensus genomes from multiple samples.

`./viruses/viral_db.list.txt` - a list of all viral sequences included in the dereplicated mapping database.

`./R/Fig4ABC.R` - Script to plot panels A-C of Figure 4.

### Figure 5: mtDNA phylogenetics

`./MT_phylogenetics/bams_and_snvs` - All mapped mammalian mtDNA BAMs from market samples.

`./MT_phylogenetics/MT_Maps.ipynb` - Script to create Figure 5a mtDNA genome maps.

`./MT_phylogenetics/mt_consensus_genomes.py` - Script to reconstruct mtDNA consensus genomes.

`./MT_phylogenetics/lowconfidence_50_percent_mt_genomes.fasta` - MT genomes obtained in this study, requiring just 1x read coverage.

`./MT_phylogenetics/mediumconfidence_50_percent_mt_genomes.fasta` - MT genomes obtained in this study, requiring 3x read coverage to call a base.

`./MT_phylogenetics/mtDNA_RD_phylogenetics.clean.ipynb` - Raccoon dog phylogenetics analysis notebook.

`./MT_phylogenetics/raccoon_dog_cytB_all.duplicated.fasta` - Raccon dog reference cytB alignment.


## Supplementary: rRNA analysis

`./mtDNA/rRNA/count_reads93_rRNAs.py` - Script to count mtDNA reads ignoring rRNA regions.

`./mtDNA/rRNA/mt93_counts_rrna.tsv` - Read counts ignoring rRNA regions.

`./mtDNA/rRNA/rRNA_analysis.ipynb` - Script to identify rRNA regions.

`./mtDNA/rRNA/rrna_positions.tsv` - positions of rRNA genes in the MT genomes.

## Supplementary: Correlations

`./correlations/Correlation.ipynb` - Script to calculate correlations between animal sequence read abundances and SARS-CoV-2.


### He et al. 2022 Cell data: Animal data from other markets

`./HeCell2022/` - location of these files

`He2022_metadata.csv` - Metadata about all of the files from that study

`cell2022_lowconfidence_mitochondria.fasta` - FASTA file of consensus mitochondrial genomes with at least 1 read covering each position.

`cell2022_mediumconfidence_mitochondria.fasta` - FASTA file of consensus mitochondrial genomes with at least 3 reads covering each position.

`./scripts/mt_consensus_genomes_He2022.py` - script used to generate Mitochondrial genomes above

`game_animals.xlsx` - Additional metadata from He et al. 2022.


### FASTQ processing

`./data/fastq_processing/` - data related to FASTQ files

` all_adapters.fa ` - list of adapters used for trimming

` frequent_adapter_summary.tsv` - results of adapter trimming, 5 most common adapters per library.

` run_bbduk.sh` - List of BBDuk commands run for each file for trimming.

` ./bbduk_results/*.txt` - Full adapter trimming output for each file.
