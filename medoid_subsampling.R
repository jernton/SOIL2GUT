library(cluster)
library(tidyverse)
library(microeco)

# Unsupervised method for subsampling representative samples for RNA extraction from prior amplicon sequencing dataset.

# Calculate clusters and medoids (here 30), providing a distance matrix (here a jaccard matrix in microeco object dsFFecT2 from running section Data Preparation in faecal_metataxonomics_markdown.html) as input.
Med <- pam(dsFFecT2$beta_diversity$jaccard, 30, diss = TRUE)


# Mark selected subsamples (medoids) as "TRUE" in metadata
dat <- dsFFecT2$sample_table
dat$SampleID <- rownames(dat) # Need to add rownames as column for mutate script below to work
dat <- dat %>%
  mutate(Subsample = case_when(
    SampleID %in% Med$medoids ~ "TRUE",
    TRUE ~ "FALSE"
  ))

# Subset metadata
dat30 <- subset(dat, medoid30 == TRUE)

# Check subsample representation in treatment groups
sample_data(dat30) %>% 
  group_by(Soil_type) %>%
  tally()

# "Second best" samples
# Some subsamples may have sub-par RNA quality.
# To replace the these samples, we choose the "second best" sample from the same cluster, i.e. the sample closest to the medoid if n > 2, 
# or just the remaining sample if n = 2.
# # Convert microeco object to phyloseq and define data
library(phyloseq)
library(file2meco)
ps <- meco2phyloseq(mec)

# Get distances
# creates Jaccard matrix of pairs Var1 and Var2 are sample names
jc.m = phyloseq::distance(ps, method = "bray", binary = T) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols = c(everything(), -Var1), names_to = "Var2", values_to = "value") %>%
  filter(as.character(Var1) != as.character(Var2))

# get sample data

sd = data.frame(sample_data(ps)) %>% 
  select(SampleID, clustering) 

# check data and convert integer to character if needed
sapply(jc.m, class) 
sapply(sd, class)
sd$clustering <- sd$clustering %>% as.character()

# combine distances with sample data
# rename sample column name as Var1 to match with other matrix
colnames(sd) = c("Var1", "Type1")
sd$Var1 <- as.character(sd$Var1)
jc.sd = left_join(jc.m, sd, by = "Var1")
# rename again
colnames(sd) = c("Var2", "Type2")
jc.sd = left_join(jc.sd, sd, by = "Var2")

# Now extract rows where Type1 == Type 2 i.e. both samples belong to same cluster
sel_rows <- jc.sd[jc.sd$Type1 == jc.sd$Type2, ]

# Export the resulting data frame. Examine the table and select "next best" samples.
write.table(sel_rows, file = 'medoid_distances.txt', sep = "\t", row.names = FALSE)
