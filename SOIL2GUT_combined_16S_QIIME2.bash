
# This pipeline preprocesses SOIL2GUT JYU (faeces+ileum+soil) and FEM (skin) Miseq runs separately at first, then combines them before 97% clustering.
# Separate runs (especially on separate machines) should be denoised with their own run-specific error profiles. Use exact same settings for CUTADAPT and DADA2 otherwise, decontam and combine with merge plugin.1

# Process JYU dataset.
module load qiime2/2024.10-amplicon

PROJECT_DIR=/scratch/project_2003826/SOIL2GUT_16S
MANIFEST=$PROJECT_DIR/manifest.txt
METADATA=$PROJECT_DIR/SOIL2GUT_metadata_mastersheet.txt
CLASSIFIER=$PROJECT_DIR/2024.09.backbone.full-length.nb.sklearn-1.4.2.qza
cd $PROJECT_DIR

# Import
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $MANIFEST \
  --output-path $PROJECT_DIR/demux-paired-end.qza \
  --input-format PairedEndFastqManifestPhred33V2

# Quality summary
qiime demux summarize \
  --i-data $PROJECT_DIR/demux-paired-end.qza \
  --o-visualization $PROJECT_DIR/demux-summary.qzv

# Primer trimming with cutadapt
# 451F: 		CCTAYGGGRBGCASCAG
# 805R (JYU):	GACTACNVGGGTWTCTAATCC
# 806R (FEM):   GGACTACNVGGGTWTCTAAT
# To combine JYU (ileum, soil, feces) and FEM (skin) datasets, use overlap for R for both datasets:
# 805R (JYU).....GACTACNVGGGTWTCTAATCC..    
# 806R (FEM)....GGACTACNVGGGTWTCTAAT....              
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences $PROJECT_DIR/demux-paired-end.qza \
  --p-front-f CCTAYGGGRBGCASCAG \
  --p-front-r GACTACNVGGGTWTCTAAT \
  --p-discard-untrimmed \
  --o-trimmed-sequences $PROJECT_DIR/demux-trimmed.qza \
  --p-cores 4

# Quality summary
qiime demux summarize \
  --i-data $PROJECT_DIR/demux-trimmed.qza \
  --o-visualization $PROJECT_DIR/demux-trimmed-summary.qzv


# DADA2 (examine demux-trimmed-summary.qzv and adjust trunc-len)
# Do this for JYU and FEM datasets separately. Each individual Miseq run should use their own error profile.
# 19.9.2025 update: you must also use the same dada2 settings to make ASVs comparable between datasets.
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $PROJECT_DIR/demux-trimmed.qza \
  --p-trunc-len-f 280 \
  --p-trunc-len-r 240 \
  --o-table $PROJECT_DIR/table.qza \
  --o-representative-sequences $PROJECT_DIR/rep-seqs.qza \
  --o-denoising-stats $PROJECT_DIR/denoising-stats.qza \
  --p-n-threads 4

# Quality summary
qiime metadata tabulate \
  --m-input-file $PROJECT_DIR/denoising-stats.qza \
  --o-visualization $PROJECT_DIR/denoising-stats.qzv


# Decontam (Prevalence-based)
qiime quality-control decontam-identify \
  --i-table $PROJECT_DIR/table.qza \
  --m-metadata-file $METADATA \
  --p-method 'prevalence' \
  --p-prev-control-column is.neg \
  --p-prev-control-indicator "TRUE" \
  --o-decontam-scores $PROJECT_DIR/decontam-scores.qza
 
qiime quality-control decontam-score-viz \
	--i-decontam-scores $PROJECT_DIR/decontam-scores.qza \
	--i-table $PROJECT_DIR/table.qza \
	--o-visualization $PROJECT_DIR/decontam-visual.qzv

qiime quality-control decontam-remove \
	--i-decontam-scores $PROJECT_DIR/decontam-scores.qza \
	--i-table $PROJECT_DIR/table.qza \
	--i-rep-seqs $PROJECT_DIR/rep-seqs.qza \
	--p-threshold 0.1 \
	--o-filtered-table $PROJECT_DIR/table-decontam.qza \
	--o-filtered-rep-seqs $PROJECT_DIR/rep-seqs-decontam.qza

####################

#Repeat for skin dataset (FEM). Each individual Miseq run should use their own DADA2 error profile.

PROJECT_DIR=/scratch/project_2003826/SOIL2GUT_16S
MANIFEST=$PROJECT_DIR/manifest.txt
METADATA=$PROJECT_DIR/SOIL2GUT_metadata_ileum-fecal-soil_FULL_nodupe.txt
CLASSIFIER=$PROJECT_DIR/2024.09.backbone.full-length.nb.sklearn-1.4.2.qza
cd $PROJECT_DIR

# Import
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $MANIFEST \
  --output-path $PROJECT_DIR/demux-paired-end.qza \
  --input-format PairedEndFastqManifestPhred33V2

# Quality summary
qiime demux summarize \
  --i-data $PROJECT_DIR/demux-paired-end.qza \
  --o-visualization $PROJECT_DIR/demux-summary.qzv

# Primer trimming with cutadapt
# Must use same sequences as the other dataset              
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences $PROJECT_DIR/demux-paired-end.qza \
  --p-front-f CCTAYGGGRBGCASCAG \
  --p-front-r GACTACNVGGGTWTCTAAT \
  --p-discard-untrimmed \
  --o-trimmed-sequences $PROJECT_DIR/demux-trimmed.qza \
  --p-cores 4

# Quality summary
qiime demux summarize \
  --i-data $PROJECT_DIR/demux-trimmed.qza \
  --o-visualization $PROJECT_DIR/demux-trimmed-summary.qzv


# DADA2 (examine demux-trimmed-summary.qzv and adjust trunc-len)
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $PROJECT_DIR/demux-trimmed.qza \
  --p-trunc-len-f 280 \
  --p-trunc-len-r 240 \
  --o-table $PROJECT_DIR/table.qza \
  --o-representative-sequences $PROJECT_DIR/rep-seqs.qza \
  --o-denoising-stats $PROJECT_DIR/denoising-stats.qza \
  --p-n-threads 4

# Quality summary
qiime metadata tabulate \
  --m-input-file $PROJECT_DIR/denoising-stats.qza \
  --o-visualization $PROJECT_DIR/denoising-stats.qzv


# Decontam (Prevalence-based)
qiime quality-control decontam-identify \
  --i-table $PROJECT_DIR/table.qza \
  --m-metadata-file $METADATA \
  --p-method 'prevalence' \
  --p-prev-control-column is.neg \
  --p-prev-control-indicator "TRUE" \
  --o-decontam-scores $PROJECT_DIR/decontam-scores.qza
 
qiime quality-control decontam-score-viz \
	--i-decontam-scores $PROJECT_DIR/decontam-scores.qza \
	--i-table $PROJECT_DIR/table.qza \
	--o-visualization $PROJECT_DIR/decontam-visual.qzv

qiime quality-control decontam-remove \
	--i-decontam-scores $PROJECT_DIR/decontam-scores.qza \
	--i-table $PROJECT_DIR/table.qza \
	--i-rep-seqs $PROJECT_DIR/rep-seqs.qza \
	--p-threshold 0.1 \
	--o-filtered-table $PROJECT_DIR/table-decontam.qza \
	--o-filtered-rep-seqs $PROJECT_DIR/rep-seqs-decontam.qza

mv table-decontam.qza table-decontam-skin.qza
mv rep-seqs-decontam.qza rep-seqs-decontam-skin.qza

# Combine JYU and FEM datasets
PROJECT_DIR_SKIN=/scratch/project_2003826/SOIL2GUT_SKIN
# Feature table
qiime feature-table merge \
  --i-tables $PROJECT_DIR/table-decontam.qza $PROJECT_DIR_SKIN/table-decontam-skin.qza \
  --p-overlap-method error_on_overlapping_sample \
  --o-merged-table merged-table.qza
# sequences
qiime feature-table merge-seqs \
    --i-data $PROJECT_DIR/rep-seqs-decontam.qza $PROJECT_DIR_SKIN/rep-seqs-decontam-skin.qza \
    --o-merged-data merged-seqs.qza

# VSEARCH clustering at 97% sequence similarity (species level) !!! OPTIONAL! For this dataset do ASVs and then tip_glom in phyloseq
#qiime vsearch cluster-features-de-novo \
#--i-sequences rep-seqs.qza \
#--i-table table.qza  \
#--p-perc-identity 0.97 \
#--o-clustered-table table-clust97.qza \
#--o-clustered-sequences rep-seqs-clust97.qza

# Phylogenetic tree construction
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $PROJECT_DIR/merged-seqs.qza \
  --p-n-threads 4 \
  --o-alignment $PROJECT_DIR/aligned-seqs.qza \
  --o-masked-alignment $PROJECT_DIR/masked-aligned-seqs.qza \
  --o-tree $PROJECT_DIR/unrooted-tree.qza \
  --o-rooted-tree $PROJECT_DIR/rooted-tree.qza

# Taxonomic classification (GREENGENES2). High memory required, submit as batch job.
qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER \
  --i-reads $PROJECT_DIR/merged-seqs.qza \
  --o-classification $PROJECT_DIR/taxonomy.qza

# Summarize outputs
qiime metadata tabulate \
  --m-input-file $PROJECT_DIR/taxonomy.qza \
  --o-visualization $PROJECT_DIR/taxonomy.qzv

qiime feature-table summarize \
  --i-table $PROJECT_DIR/merged-table.qza \
  --o-visualization $PROJECT_DIR/merged-table.qzv \
  --m-sample-metadata-file $METADATA

qiime feature-table tabulate-seqs \
  --i-data $PROJECT_DIR/merged-seqs.qza \
  --o-visualization $PROJECT_DIR/merged-seqs.qzv 


# Clean table from chloroplasts and mitochondria
qiime taxa filter-table \
  --i-table merged-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table merged-table-noMtCp.qza


# Export files
EXPORT_DIR=$PROJECT_DIR/exported
mkdir -p $EXPORT_DIR

qiime tools export \
  --input-path $PROJECT_DIR/table-decontam.qza \
  --output-path $EXPORT_DIR

biom convert -i $EXPORT_DIR/feature-table.biom -o $EXPORT_DIR/feature-table.tsv --to-tsv 

qiime tools export \
  --input-path $PROJECT_DIR/taxonomy.qza \
  --output-path $EXPORT_DIR

qiime tools export \
  --input-path $PROJECT_DIR/rep-seqs-decontam.qza \
  --output-path $EXPORT_DIR

qiime tools export \
  --input-path $PROJECT_DIR/rooted-tree.qza \
  --output-path $EXPORT_DIR

# Rename tree file
mv $EXPORT_DIR/tree.nwk $EXPORT_DIR/rooted-tree.nwk


# Extra:

# Clean data from chloroplasts and mitochondria
qiime taxa filter-table \
  --i-table table-decontam.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mt-chloro.qza

# Core metrics
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table-no-mt-chloro.qza \
  --p-sampling-depth 11232 \
  --m-metadata-file $METADATA \
  --output-dir core-metrics-results

# Longitudinal pairwise distances
cd core-metrics-results/
qiime longitudinal pairwise-distances \
  --i-distance-matrix bray_curtis_distance_matrix.qza \
  --m-metadata-file $METADATA \
  --p-group-column Soil_type \
  --p-state-column sampling_timepoint \
  --p-state-1 "T1" \
  --p-state-2 "T2" \
  --p-individual-id-column Pup_ID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-bray.qzv
  
qiime longitudinal pairwise-distances \
  --i-distance-matrix jaccard_distance_matrix.qza \
  --m-metadata-file $METADATA \
  --p-group-column Soil_type \
  --p-state-column sampling_timepoint \
  --p-state-1 "T1" \
  --p-state-2 "T2" \
  --p-individual-id-column Pup_ID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-jaccard.qzv 
  