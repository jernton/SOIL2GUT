
module load qiime2/2024.10-amplicon

PROJECT_DIR=/scratch/project_2003826/SOIL2GUT_SKIN
MANIFEST=$PROJECT_DIR/manifest_skin.txt
METADATA=$PROJECT_DIR/metadata_skin.txt
CLASSIFIER=$PROJECT_DIR/2024.09.backbone.full-length.nb.sklearn-1.4.2.qza
cd $PROJECT_DIR

# 1. Import
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $MANIFEST \
  --output-path $PROJECT_DIR/demux-paired-end.qza \
  --input-format PairedEndFastqManifestPhred33V2

# 2. Quality summary
qiime demux summarize \
  --i-data $PROJECT_DIR/demux-paired-end.qza \
  --o-visualization $PROJECT_DIR/demux-summary.qzv

# 3. Primer trimming with cutadapt
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences $PROJECT_DIR/demux-paired-end.qza \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACNVGGGTWTCTAAT \
  --p-discard-untrimmed \
  --o-trimmed-sequences $PROJECT_DIR/demux-trimmed.qza \
  --p-cores 4

# 4. DADA2 (adjust trunc-len as per demux-summary.qzv)
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $PROJECT_DIR/demux-trimmed.qza \
  --p-trunc-len-f 280 \
  --p-trunc-len-r 240 \
  --o-table $PROJECT_DIR/table.qza \
  --o-representative-sequences $PROJECT_DIR/rep-seqs.qza \
  --o-denoising-stats $PROJECT_DIR/denoising-stats.qza \
  --p-n-threads 4


qiime metadata tabulate \
  --m-input-file $PROJECT_DIR/denoising-stats.qza \
  --o-visualization $PROJECT_DIR/denoising-stats.qzv


# 5. Decontam (Prevalence-based)
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


# 6. Phylogenetic tree construction
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $PROJECT_DIR/rep-seqs-decontam.qza \
  --p-n-threads 4 \
  --o-alignment $PROJECT_DIR/aligned-rep-seqs.qza \
  --o-masked-alignment $PROJECT_DIR/masked-aligned-rep-seqs.qza \
  --o-tree $PROJECT_DIR/unrooted-tree.qza \
  --o-rooted-tree $PROJECT_DIR/rooted-tree.qza

# 7. Taxonomic classification (GREENGENES2)
qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER \
  --i-reads $PROJECT_DIR/rep-seqs-decontam.qza \
  --o-classification $PROJECT_DIR/taxonomy.qza

# 8. Summarize outputs
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



# 9. Export files
EXPORT_DIR=$PROJECT_DIR/exported
mkdir -p $EXPORT_DIR

# Export feature table
qiime tools export \
  --input-path $PROJECT_DIR/table-decontam.qza \
  --output-path $EXPORT_DIR
  
biom convert -i $EXPORT_DIR/feature-table.biom -o $EXPORT_DIR/feature-table.tsv --to-tsv

# Export taxonomy
qiime tools export \
  --input-path $PROJECT_DIR/taxonomy.qza \
  --output-path $EXPORT_DIR

# Export sequences
qiime tools export \
  --input-path $PROJECT_DIR/rep-seqs-decontam.qza \
  --output-path $EXPORT_DIR

# Export tree
qiime tools export \
  --input-path $PROJECT_DIR/rooted-tree.qza \
  --output-path $EXPORT_DIR

# Rename tree file
mv $EXPORT_DIR/tree.nwk $EXPORT_DIR/rooted-tree.nwk


biom convert -i feature-table.biom -o feature-table.tsv --to-tsv
