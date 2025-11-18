
conda activate qiime2-amplicon-2024.10


# Import data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path SG_manifest_ITS2.txt \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
  
# Trim with CUTADAPT. ITS2 data needs more detailed trimming than what DADA2 allows due to variable fragment length. When fragment is shorter than read length, read-through primers on 3' ends must also be found and removed.
qiime cutadapt trim-paired \
	--i-demultiplexed-sequences paired-end-demux.qza \
	--o-trimmed-sequences trimmed-seqs.qza \
	--p-front-f GCATCGATGAAGAACGCAG \
	--p-adapter-f GCATATCAATAAGCGSAGGA \
	--p-front-r TCCTSCGCTTATTGATATGC \
	--p-adapter-r CTGCGTTCTTCATCGATGC \
	--p-error-rate 0.2 \
	--p-discard-untrimmed TRUE \
	--p-quality-cutoff-3end 20 \
	--p-cores 4
	
# Get summaries
qiime demux summarize \
  --i-data trimmed-seqs.qza \
  --o-visualization trimmed-seqs.qzv
  
  
# DADA2 denoising. Leave trimming options at 0, cutadapt already does that
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed-seqs.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-max-ee-r 4 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats dada2-stats.qza \
  --p-n-threads 4

# Reads passed through DADA2  
qiime metadata tabulate \
  --m-input-file dada2-stats.qza \
  --o-visualization dada2-stats.qzv

qiime metadata tabulate \
  --m-input-file rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file table.qza \
  --o-visualization table.qzv

# Decontam
qiime quality-control decontam-identify \
  --i-table dada2-asv-table-maxEEr4.qza \
  --m-metadata-file SOILGUT_cleaned_metadata_ITS2.txt \
  --p-method 'prevalence' \
  --p-prev-control-column is.neg \
  --p-prev-control-indicator "TRUE" \
  --o-decontam-scores decontam_table_ITS2.qza
 
qiime quality-control decontam-score-viz \
	--i-decontam-scores decontam_table_ITS2.qza \
	--i-table dada2-asv-table-maxEEr4.qza \
	--o-visualization decontam_viz.qzv

qiime quality-control decontam-remove \
	--i-decontam-scores decontam_table_ITS2.qza \
	--i-table dada2-asv-table-maxEEr4.qza \
	--o-filtered-table feature-table.qza \

# VSEARCH
# cluster the sequences 97%
qiime vsearch cluster-features-de-novo \
--i-sequences rep-seqs.qza \
--i-table table.qza  \
--p-perc-identity 0.97 \
--o-clustered-table table-dada2-0_97_clust.qza \
--o-clustered-sequences rep-seqs-dada2-0_97_clust.qza

# Sequences  
qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2-0_97_clust.qza \
  --o-visualization rep-seqs-dada2-0_97_clust.qzv

# Skip phylogenetic tree for ITS2 data, the region is not suitable.

# Classifications. Download a pre-trained QIIME2-format classifier at https://github.com/colinbrislawn/unite-train/releases
qiime feature-classifier classify-sklearn \
  --i-classifier unite_ver10_dynamic_s_all_04.04.2024-Q2-2024.2.qza \
  --i-reads rep-seqs-dada2-0_97_clust.qza \
  --o-classification dada2-taxonomy-pretrained.qza \
  --p-n-jobs 4
  
# Examine
qiime metadata tabulate \
  --m-input-file dada2-taxonomy-pretrained.qza \
  --o-visualization dada2-taxonomy-pretrained.qzv

  


# Export data for analysis in R
# You can either export data from artifacts to be loaded in R or import artifacts directly in R using r library 'qiime2R'
qiime tools export --input-path taxonomy.qza --output-path exported/
#qiime tools export --input-path rooted-tree.qza --output-path exported/
qiime tools export --input-path rep-seqs.qza --output-path exported/
qiime tools export --input-path table.qza --output-path exported/
biom convert -i exported/feature-table.biom -o exported/feature-table.tsv --to-tsv


# EXTRA:
# Below are some analysis steps if you wish to do them in QIIME2 environment.


# Filter to k__Fungi
qiime taxa filter-table \
  --i-table table_97_clust.qza \
  --i-taxonomy taxonomy_97_clust.qza \
  --p-include "k__Fungi" \
  --o-filtered-table table-fungi.qza


# Core metrics
qiime diversity core-metrics \
  --i-table table.qza \
  --p-sampling-depth 8251 \
  --m-metadata-file SOIL2GUT_metadata_ITS2.tsv \
  --output-dir core-metrics-results


# Longitudinal pairwise distances (First distances)
cd core-metrics-results/
qiime longitudinal pairwise-distances \
  --i-distance-matrix jaccard_distance_matrix.qza \
  --m-metadata-file ../SOIL2GUT_metadata_ITS2.tsv \
  --p-group-column Soil_type \
  --p-state-column Sampling_timepoint \
  --p-state-1 "T1" \
  --p-state-2 "T2" \
  --p-individual-id-column Pup_ID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-jaccard.qzv 
 
qiime longitudinal pairwise-distances \
  --i-distance-matrix bray_curtis_distance_matrix.qza \
  --m-metadata-file ../SOIL2GUT_metadata_ITS2.tsv \
  --p-group-column Soil_type \
  --p-state-column Sampling_timepoint \
  --p-state-1 "T1" \
  --p-state-2 "T2" \
  --p-individual-id-column Pup_ID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-bray_curtis.qzv 

# FUNguild
python FUNGuild/Guilds_v1.1.py -otu fungi_zeroes.txt -db fungi -m -u
