
module load qiime2/2024.10-amplicon

PROJECT_DIR=/scratch/project_2003826/SOIL2GUT_16S
MANIFEST=$PROJECT_DIR/manifest_feces_soil.txt
METADATA=$PROJECT_DIR/SOIL2GUT_metadata_feces-soil.txt
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
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences $PROJECT_DIR/demux-paired-end.qza \
  --p-front-f CCTAYGGGRBGCASCAG \
  --p-front-r GGACTACNVGGGTWTCTAAT \
  --p-discard-untrimmed \
  --o-trimmed-sequences $PROJECT_DIR/demux-trimmed.qza \
  --p-cores 4

# Quality summary
qiime demux summarize \
  --i-data $PROJECT_DIR/demux-trimmed.qza \
  --o-visualization $PROJECT_DIR/demux-trimmed-summary.qzv


# DADA2 (adjust trunc-len as per demux-trimmed-summary.qzv)
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



# Phylogenetic tree construction
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $PROJECT_DIR/rep-seqs-decontam.qza \
  --p-n-threads 4 \
  --o-alignment $PROJECT_DIR/aligned-rep-seqs.qza \
  --o-masked-alignment $PROJECT_DIR/masked-aligned-rep-seqs.qza \
  --o-tree $PROJECT_DIR/unrooted-tree.qza \
  --o-rooted-tree $PROJECT_DIR/rooted-tree.qza

# Taxonomic classification (GREENGENES2)
qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER \
  --i-reads $PROJECT_DIR/rep-seqs-decontam.qza \
  --o-classification $PROJECT_DIR/taxonomy.qza

# Summarize outputs
qiime metadata tabulate \
  --m-input-file $PROJECT_DIR/taxonomy.qza \
  --o-visualization $PROJECT_DIR/taxonomy.qzv

qiime feature-table summarize \
  --i-table $PROJECT_DIR/table-decontam.qza \
  --o-visualization $PROJECT_DIR/table-decontam.qzv \
  --m-sample-metadata-file $METADATA

qiime feature-table tabulate-seqs \
  --i-data $PROJECT_DIR/rep-seqs-decontam.qza \
  --o-visualization $PROJECT_DIR/rep-seqs-decontam.qzv



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

# Longitudinal pairwise distances (First distances)
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
  