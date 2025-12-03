# clear environment
#rm(list=ls())

# own packages - in project libpath
.libPaths(c("/projappl/pwatts/rpackages_432", .libPaths()))
libpath <- .libPaths()[1]
# this command can be used to check that the folder is now visible:
.libPaths() 
# it should be first on the list

# set working directory
#setwd("/scratch/project_2003826/soil2gut/")
#for fungal analysis
setwd("/scratch/project_2003826/soil2gut/fungi")

# setup environment
library(phyloseq) #works
library(ggplot2)
library(ggpubr)
library(vegan)
library(dplyr)
library(hrbrthemes)
library(gcookbook)
library(tidyverse) #if needed?
library(devtools)
library(knitr)
library(rstatix)
library(ggsci)
library(ComplexHeatmap)
library(RColorBrewer)
library(igraph)
library(network)
library(sna)
library(tidyfst)
library(ggrepel)
library(psych)
#library(Biostrings)
library(remotes)
library(ggalluvial)
library(ggraph)
library(WGCNA)
library(ggnewscale)
library(pulsar)
library(patchwork)
library(EasyStat)
library(ggClusterNet)
library(SpiecEasi)
library(gameofthrones)
library(qiime2R)
library(gghalves)

#ps <- qza_to_phyloseq(
#  features="TABLE.qza",
#  taxonomy="TAXONOMY.qza",
#  tree="REP-SEQS-root-tree.qza",
#  metadata = "METADATA.txt"
#)

otu_table <- read.csv("SOILGUT_cleaned_ASVs_tipglom0.05.txt", sep = "\t", row.names = 1)
otu_table <- as.matrix(otu_table)
str(otu_table)

taxonomy <- read.csv("SOILGUT_cleaned_taxonomy_tipglom0.05.txt", sep = "\t", row.names = 1)
taxonomy <- as.matrix(taxonomy)
str(taxonomy)

metadata <- read.table("SOILGUT_cleaned_metadata_tipglom0.05.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
phy_tree <- read_tree("SOILGUT_cleaned_tree_rooted_tipglom0.05.tree")

OTU <- otu_table(otu_table, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy)
MET <- sample_data(metadata)

taxa_names(TAX)
taxa_names(OTU)
taxa_names(phy_tree)
MET[is.na(MET)] <- "NA"
MET

ps <- phyloseq(OTU, TAX, MET, phy_tree)

#subset data
#note that its - sample timepoint is t1, t2 - data cannot be merged....
ps.fecal = subset_samples(ps, Sample_type == "fecal_community")
ps.fecal
ps.start = subset_samples(ps.fecal, Sampling_timepoint == "Start")
ps.start
ps.end = subset_samples(ps.fecal, Sampling_timepoint == "End")
ps.end

ps.end.n = subset_samples(ps.end, Soil_type == "National park")
ps.end.n
ps.end.u = subset_samples(ps.end, Soil_type == "Urban forest")
ps.end.u
ps.end.c = subset_samples(ps.end, Soil_type == "Control")
ps.end.c

####
## NEW DATA FROM TONI

## TRY QIIME IMPORT
ps.fungi <- qza_to_phyloseq(
  features="table_97_clust.qza",
  taxonomy="taxonomy_97_clust.qza",
#  tree="REP-SEQS-root-tree.qza", # no tree for fungi
  metadata = "its-SOILGUT_cleaned_metadata_ITS2.txt"
)

ps.fungi

#otu_table.1 <- read.csv("its-feature-table.tsv", sep = "\t", row.names = 1)
#otu_table.1 <- as.matrix(otu_table.1)
#str(otu_table.1)

#taxonomy.1 <- read.csv("its-taxonomy.tsv", sep = "\t", row.names = 1)
#taxonomy.1 <- as.matrix(taxonomy.1)
#str(taxonomy.1)

#metadata.1 <- read.table("its-SOILGUT_cleaned_metadata_ITS2.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
#phy_tree.1 <- read_tree("its-tree.nwk")

#OTU1 <- otu_table(otu_table.1, taxa_are_rows = TRUE)
#TAX1 <- tax_table(taxonomy.1)
#MET1 <- sample_data(metadata.1)

#taxa_names(TAX1)
#taxa_names(OTU1)
#taxa_names(phy_tree.1)
#MET1[is.na(MET1)] <- "NA"
#MET1

#ps.fungi <- phyloseq(OTU1, TAX1, MET1)

#subset data
# remove plants - Viridiplantae_phy_Incertae_sedis
table(tax_table(ps.fungi)[, "Phylum"])

# Remove prefixes like k__, p__, etc., from taxonomy
tax <- tax_table(ps.fungi) %>% as.data.frame()
tax_clean <- apply(tax, 2, function(x) gsub("^[a-z]__","", x))

# Assign back the cleaned taxonomy
tax_table(ps.fungi) <- tax_table(as.matrix(tax_clean))
table(tax_table(ps.fungi)[, "Kingdom"])

ps.fungi.clean <- subset_taxa(ps.fungi, Kingdom == "Fungi")
ps.fungi.clean

ps.fungi.fecal = subset_samples(ps.fungi.clean, Sample_type == "fecal_community")
ps.fungi.fecal
ps.fungi.start = subset_samples(ps.fungi.fecal, Sampling_timepoint == "T1")
ps.fungi.start
ps.fungi.end = subset_samples(ps.fungi.fecal, Sampling_timepoint == "T2")
ps.fungi.end

ps.fungi.end.n = subset_samples(ps.fungi.end, Soil_type == "National park")
ps.fungi.end.n
ps.fungi.end.u = subset_samples(ps.fungi.end, Soil_type == "Urban forest")
ps.fungi.end.u
ps.fungi.end.c = subset_samples(ps.fungi.end, Soil_type == "Control")
ps.fungi.end.c

###
save.image(file = "soil2gut.fecal.bact.fungi.RData")
#load(file = "soil2gut.fecal.bact.fungi.RData")

### CODE CLEANED AND ANNOTATED
### SETUP & NETWORK CONSTRUCTION
### Build microbial co-occurrence networks for each group using ggClusterNet
### NOTE that result for bacteria is tab.r - here it has been named tab.fungi

tab.fungi = network.pip(
  ps = ps.fungi.end,
  N = 0, # number OTUs - 0 for all (=586 max for this ps object before subsetting)
  ra = 0.05,
  big = TRUE, #change if > 200 OTUs
  select_layout = TRUE,
  #layout_net = "model_maptree2",
  layout_net = "model_igraph",
  r.threshold = 0.6, #default is 0.6
  p.threshold = 0.05,
  #maxnode = 2,
  method = "pearson", #sparcc requires SpiecEasi
  label = FALSE,  #add node label
  #lab = "elements",
  group = "Soil_type", #needs group column
  fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE, # takes memory
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 500,
  R = 50,
  ncpus = 4
)

save.image(file = "soil2gut.fecal.bact.fungi.r06p05.RData")
#load(file = "soil2gut.fecal.bact.fungi.r06p05.RData")

### EXTRACT NETWORK VISUALISATION AND CORRELATION MATRICES
plot <- tab.fungi[[1]]
p0 <- plot[[1]]      # Network visualisation
p0.1 <- plot[[2]]    # ZiPi plot
p0.2 <- plot[[3]]    # Power law / random networks

dat <- tab.fungi[[2]]
cortab <- dat$net.cor.matrix$cortab

### NETWORK COMPARISON STATISTICS
dat.compare <- module.compare.net.pip(
  ps = NULL,
  corg = cortab,
  degree = TRUE,
  zipi = TRUE,
  r.threshold = 0.6,
  p.threshold = 0.05,
  method = "pearson",
  padj = FALSE,
  n = 4
)
res <- dat.compare[[1]]
head(res)

### CUSTOMISED NETWORK VISUALISATION
node <- dat$net.cor.matrix$node
edge <- dat$net.cor.matrix$edge

p <- ggplot() +
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = cor),
               data = edge, size = 0.05, alpha = 0.8) +
  geom_point(aes(x = X1, y = X2, fill = Phylum, size = igraph.degree),
             pch = 21, data = node, color = "grey80") +
  facet_wrap(. ~ label, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = c("black", "blue")) +
  scale_size(range = c(0.8, 5)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        panel.grid = element_blank())

print(p)

### CUSTOMISE networks
### add abundance information to node
library(tibble)

# transform to relative abundance if needed
ps.rel <- transform_sample_counts(ps.fungi.end, function(x) x / sum(x))

# get mean abundance per ASV/OTU
otu_abund <- otu_table(ps.rel) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%  # Use "ID" to match your node object
  mutate(mean_abundance = rowMeans(dplyr::select(., -ID)))

# combine with taxonomy
tax_data <- tax_table(ps.fungi.end) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID")

abund_tax <- dplyr::left_join(otu_abund, tax_data, by = "ID")

node <- dplyr::left_join(node, abund_tax, by = "ID") #if use of node fails later it may be because we have now changed node...

# identify top 4 most abundant phyla across all nodes ---
top_phyla <- node %>%
  dplyr::filter(!is.na(Phylum.y)) %>% #note that Phylum was renamed to Phylum.y
  dplyr::group_by(Phylum.y) %>%
  dplyr::summarise(total = sum(mean_abundance, na.rm = TRUE)) %>% #and colum called mean_abundance
  dplyr::arrange(desc(total)) %>%
  dplyr::slice_head(n = 4)

top_phyla_vec <- top_phyla$Phylum.y

#TONI's palatte - Control-National-Urban the palette is "Dark2" from package RColorBrewer, 
#or manually: Control = "#1b9e77", National = "#d95f02", Urban = "#7570b3"
#brewer.pal(6, "Dark2")
# "#1B9E77" "#D95F02" "#7570B3" "#E7298A"
# "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02"
#Control = "#1b9e77", National = "#d95f02", Urban = "#7570b3"

# define custom colors for top 4 phyla
custom_colors <- c(
  "Ascomycota"    = "#1B9E77",
  "Mucoromycota"  = "#D95F02",
  "Basidiomycota"        = "#7570B3",
  "Mortierellomycota"      = "#E7298A"
)
custom_colors <- custom_colors[names(custom_colors) %in% top_phyla_vec]
custom_colors <- c(custom_colors, Other = "grey80")  # Add "Other"

# assign phylum colors to node data
node$Phylum_colored <- ifelse(
  node$Phylum.y %in% names(custom_colors),
  node$Phylum.y,
  "Other"
)

# ensure facet order is Control, National park, Urban forest
# convert label to factor with desired order
desired_order <- c("Control", "National park", "Urban forest")
node$label <- factor(node$Group, levels = desired_order)
edge$label <- factor(edge$Group, levels = desired_order)

#node$label <- factor(node$Group, levels = c("Control", "National park", "Urban forest"))

# create customized network plot
# check labels for facet label
unique(node$label)
unique(edge$label)

p_custom <- ggplot() +
  geom_segment(aes(x = X1, y = Y1,
                   xend = X2, yend = Y2,
                   color = cor),
               data = edge, size = 0.4, alpha = 0.5) +
  geom_point(aes(x = X1, y = X2,
                 fill = Phylum_colored,
                 size = igraph.degree),
             pch = 21, data = node, color = "grey30") + #colour for the outline
  #facet_wrap(. ~ label, scales = "free_y", nrow = 1) +
  facet_wrap(. ~ label, scales = "fixed", nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_manual(values = custom_colors) +
# scale_colour_gradient2(low = "blue", mid = "grey80", high = "red", midpoint = 0) +
  scale_colour_manual(values = c("grey5", "black")) +
  scale_size(range = c(2, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(size = 14, face = "plain", color = "black", margin = margin(4, 4, 4, 4)),  # Larger font for facet labels
    strip.background = element_rect(fill = "gray90", color = "gray90", size = 2)  # Great background for facet labels
  )

# Display the plot
print(p_custom)

save.image(file = "soil2gut.fecal.bact.fungi.r06p05.1.RData")
#load(file = "soil2gut.fecal.bact.fungi.r06p05.1.RData")

### --- ZiPi PLOT ---
dat.z <- dat$zipi.data

# Define ZiPi role regions
x1 <- c(0, 0.62, 0, 0.62)
x2 <- c(0.62, 1, 0.62, 1)
y1 <- c(-Inf, 2.5, 2.5, -Inf)
y2 <- c(2.5, Inf, Inf, 2.5)
lab <- c("peripheral", "network hubs", "module hubs", "connectors")
roles.colors <- c("#E6E6FA", "#DCDCDC", "#F5FFFA", "#FAEBD7")

# Repeat for each group
tab <- data.frame(x1, y1, x2, y2, lab)
tab2 <- do.call(rbind, replicate(length(unique(dat.z$group)), tab, simplify = FALSE))

q <- ggplot() +
  geom_rect(data = tab2, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = lab)) +
  scale_fill_manual(values = roles.colors) +
  geom_point(data = dat.z, aes(x = p, y = z, color = module)) +
  ggrepel::geom_text_repel(data = dat.z, aes(x = p, y = z, label = label, color = module), size = 4) +
  facet_grid(. ~ group, scales = 'free') +
  theme_bw() +
  xlab("Participation coefficient") + ylab("Within-module connectivity z-score") +
  theme(strip.background = element_rect(fill = "white"))

print(q)

### --- COMPUTE GROUP-SPECIFIC NETWORK PROPERTIES (IGRAPH) ---
node_df <- dat$net.cor.matrix$node
edge_df <- dat$net.cor.matrix$edge
groups <- unique(node_df$Group)

get_net_props <- function(group_id) {
  nodes_group <- dplyr::filter(node_df, Group == group_id)
  edges_group <- dplyr::filter(edge_df, Group == group_id)
  g <- igraph::graph_from_data_frame(d = edges_group[, c("OTU_1", "OTU_2", "weight")],
                                     vertices = nodes_group[, c("elements", "Phylum")],
                                     directed = FALSE)
  props <- net_properties.4(g, n.hub = TRUE)
  props_df <- as.data.frame(t(props))
  props_df$Group <- group_id
  return(props_df)
}

net_props_df <- dplyr::bind_rows(lapply(groups, get_net_props))
print(net_props_df)

save.image(file = "soil2gut.fecal.bact.fungi.r06p05.2.RData")
#load(file = "soil2gut.fecal.bact.fungi.r06p05.2.RData")

### ALTERNATE CODE for ABOVE - saves g from igraph
get_net_props <- function(group_id) {
  nodes_group <- dplyr::filter(node_df, Group == group_id)
  edges_group <- dplyr::filter(edge_df, Group == group_id)
  
  # Create the graph from data
  g <- igraph::graph_from_data_frame(d = edges_group[, c("OTU_1", "OTU_2", "weight")],
                                     vertices = nodes_group[, c("elements", "Phylum")],
                                     directed = FALSE)
  
  # Compute network properties
  props <- net_properties.4(g, n.hub = TRUE)
  
  # Convert properties to a data frame
  props_df <- as.data.frame(t(props))
  
  # Add group information to the properties
  props_df$Group <- group_id
  
  # Return both the graph and the properties
  return(list(graph = g, properties = props_df))
}

# Apply the function to each group and store the results
results <- lapply(groups, get_net_props)

# Extract the network properties and graphs
net_props_df <- dplyr::bind_rows(lapply(results, function(x) x$properties))
graphs <- lapply(results, function(x) x$graph)

# View the network properties
print(net_props_df)

# Access the graphs for further analysis or visualization
print(graphs)

# Plot the graph for the first group
plot(graphs[[1]])


#############

### --- COMPUTE SAMPLE-LEVEL NETWORK PROPERTIES ---
dat.f2 <- NULL
for (group_name in names(cortab)) {
  pst <- ps.fungi.end %>%
    subset_samples.wt("Soil_type", group_name) %>%
    subset_taxa.wt("OTU", colnames(cortab[[group_name]])) %>%
    filter_taxa(function(x) sum(x) > 0, TRUE) %>%
    scale_micro("rela")
  dat.f <- netproperties.sample(pst = pst, cor = cortab[[group_name]])
  dat.f$Group <- group_name
  dat.f2 <- dplyr::bind_rows(dat.f2, dat.f)
}

### --- VISUALISE SAMPLE-LEVEL METRICS ---
dat.f2$Group <- as.factor(dat.f2$Group)
r <- ggplot(dat.f2, aes(x = Group, y = no.clusters, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5, color = "black", na.rm = TRUE) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "black", na.rm = TRUE) +
  labs(title = "", fill = "Group") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_blank())

print(r)

### SAVE WORKSPACE
save.image(file = "soil2gut.fecal.bact.fungi.r06p05.4.RData")
#load(file = "soil2gut.fecal.bact.fungi.r06p05.4.RData")











## VISUALISE KEYSTONE TAXA
# Assuming g is your igraph object
degree_centrality <- igraph::degree(g)
betweenness_centrality <- igraph::betweenness(g)
closeness_centrality <- igraph::closeness(g)

# Choose centrality measure to define keystone nodes, for example degree centrality
keystone_nodes <- names(degree_centrality)[degree_centrality >= quantile(degree_centrality, 0.95)]







































tab.fungi = network.pip(
  ps = ps.end,
  N = 0, # number OTUs - 0 for all (=586 max for this ps object before subsetting)
  ra = 0.05,
  big = TRUE, #change if > 200 OTUs
  select_layout = TRUE,
  #layout_net = "model_maptree2",
  layout_net = "model_igraph",
  r.threshold = 0.6, #default is 0.6
  p.threshold = 0.05,
  #maxnode = 2,
  method = "pearson", #sparcc requires SpiecEasi
  label = FALSE,  #add node label
  #lab = "elements",
  group = "Soil_type", #needs group column
  fill = "Phylum",
  size = "igraph.degree",
  zipi = FALSE, # takes memory
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 500,
  R = 50,
  ncpus = 4
)

# extract images
plot = tab.fungi[[1]]
# extract graph visualisation results
p0 = plot[[1]] #network
# zipi
#p0.1 = plot[[2]]
# power law distribution of random networks
p0.2 = plot[[3]] #same as plot

# extract correlation matrix
dat = tab.fungi[[2]]
cortab = dat$net.cor.matrix$cortab

# network saliancy comparison
dat.compare = module.compare.net.pip(
  ps = NULL,
  corg = cortab,
  degree = TRUE,
  zipi = FALSE,
  r.threshold= 0.6,
  p.threshold=0.05,
  method = "pearson",
  padj = F,
  n = 4)

res = dat.compare[[1]]
head(res)

## customise the plots
#dat = tab.fungi[[2]] #constructed above
node = dat$net.cor.matrix$node
edge = dat$net.cor.matrix$edge
head(edge)

#this is the same as p0 created above, but customisable
p <- ggplot() + geom_segment(aes(x = X1, y = Y1, 
                                 xend = X2, 
                                 yend = Y2, 
                                 color = cor #colour for correlations
),
data = edge, 
size = 0.05, 
alpha = 0.8) +
  geom_point(aes(X1, X2,
                 fill = Phylum,
                 size = igraph.degree),
             pch = 21, 
             data = node, 
             color = "grey80") +
  facet_wrap(.~ label, scales="free_y",  nrow = 1) +
  # geom_text_repel(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
  # geom_text(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
  scale_colour_manual(values = c("black", "blue")) + #colours the correlations
  # scale_fill_hue()+
  scale_size(range = c(0.8, 5)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)
  ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(legend.position = "none",  # Removes the legend
        legend.background = element_rect(colour = NA)) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

p

## network properties
# two groups are embedded in the data so manually create group specific igraph objects (requires igraph and dplyr)
# extract node and edge data
node_df <- tab.fungi[[2]]$net.cor.matrix$node
edge_df <- tab.fungi[[2]]$net.cor.matrix$edge

# split by group 
groups <- unique(node_df$Group)  

# create a function to build igraph and compute properties
get_net_props <- function(group_id) {
  nodes_group <- node_df %>% filter(Group == group_id)
  edges_group <- edge_df %>% filter(Group == group_id)
  
  # Build igraph
  g <- graph_from_data_frame(d = edges_group[, c("OTU_1", "OTU_2", "weight")],
                             vertices = nodes_group[, c("elements", "Phylum")],
                             directed = FALSE)
  
  # Compute network properties
  props <- net_properties.4(g, n.hub = TRUE)
  props_df <- as.data.frame(t(props))  # Convert to a single-row dataframe
  props_df$Group <- group_id
  return(props_df)
}

# apply to each group
net_props_list <- lapply(groups, get_net_props)

# combine into one data frame
net_props_df <- bind_rows(net_props_list)

# view result
print(net_props_df)

### modified the code in GGCLUSTERNET
#cortab <- tab.fungi[[2]]$net.cor.matrix$cortab # already have this

# check names
names(cortab)

#requires phyloseq, dplyr, ggClusterNet
# initialize an empty data frame to hold the results
dat.f2 <- NULL

# loop over each group
for (i in seq_along(names(cortab))) {
  group_name <- names(cortab)[i]
  # subset the phyloseq object by sample group
  pst <- ps.end %>%
    subset_samples.wt("Soil_type", group_name) %>%  # your grouping variable
    subset_taxa.wt("OTU", colnames(cortab[[group_name]])) %>%
    filter_taxa(function(x) sum(x) > 0, TRUE) %>%
    scale_micro("rela")  # normalize to relative abundance
  # compute sample-level network properties
  dat.f <- netproperties.sample(pst = pst, cor = cortab[[group_name]])
  # add group name to the results
  dat.f$Group <- group_name
  # combine with previous results
  if (is.null(dat.f2)) {
    dat.f2 <- dat.f
  } else {
    dat.f2 <- bind_rows(dat.f2, dat.f)
  }
}

#examine output
head(dat.f2)
summary(dat.f2)

# plot variables in the data
#check variables
colnames(dat.f2)
# check the data type of the columns you're working with
#str(dat.f2)
# check if the column is truly numeric
#is.numeric(dat.f2$"no.clusters")
# if it returns FALSE, check for non-numeric values
#unique(dat.f2$"no.clusters")
# convert all columns except 'Group' to numeric
#dat.f2_fixed <- dat.f2 %>%
#  mutate(across(-Group, ~ as.numeric(as.character(.))))
# convert 'Group' to a factor
dat.f2$Group <- as.factor(dat.f2$Group)

#plot
ggplot(dat.f2, aes(x = Group, y = no.clusters, fill = Group)) +
  geom_boxplot(outlier.shape = NA, 
               alpha = 0.7, 
               width = 0.5, 
               color = "black", 
               na.rm = TRUE) +  # boxplot with custom width and black borders
  geom_jitter(width = 0.2, 
              alpha = 0.6, 
              color = "black", 
              na.rm = TRUE) +  # add jittered points with custom transparency and black points
  labs(title = " ", fill = "Group") +  # title and legend label
  theme_minimal() +  # minimal theme
  theme(axis.title.x = element_blank(),  # remove x-axis title for cleaner look
        axis.title.y = element_text(size = 12, face = "bold"),  # bold y-axis title
        legend.position = "top",  # position the legend at the top
        plot.title = element_text(hjust = 0.5, 
                                  size = 14, 
                                  face = "bold"),  # center and style title
        panel.grid.major = element_blank(),  # remove major grid lines for cleaner look
        panel.grid.minor = element_blank())  # remove minor grid lines

# do groups differ (2 groups)
#wilcox.test(no.clusters ~ Group, data = dat.f2, na.rm = TRUE)

save.image(file = "soil2gut.fecal.r06.RData")
#load(file = "soil2gut.fecal.r06.RData")

##
tab.F = network.pip( ######
  ps = ps.fungi.end, ######
  N = 0, # number OTUs - 0 for all (=586 max for this ps object before subsetting)
  ra = 0.05,
  big = TRUE, #change if > 200 OTUs
  select_layout = TRUE,
  #layout_net = "model_maptree2",
  layout_net = "model_igraph",
  r.threshold = 0.6, #default is 0.6
  p.threshold = 0.05,
  #maxnode = 2,
  method = "pearson", #sparcc requires SpiecEasi
  label = FALSE,  #add node label
  #lab = "elements",
  group = "Soil_type", #needs group column
  fill = "Phylum",
  size = "igraph.degree",
  zipi = FALSE, # takes memory
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 500,
  R = 50,
  ncpus = 4
)

# extract images
plot.F = tab.F[[1]] ######
p0.F = plot.F[[1]] #network ######
p0.2.F = plot.F[[3]] #same as plot  ######

# extract correlation matrix
dat.F = tab.F[[2]]     ######
cortab.F = dat.F$net.cor.matrix$cortab  #####

# network saliancy comparison
dat.compare.F = module.compare.net.pip(
  ps = NULL,
  corg = cortab.F,   ######
  degree = TRUE,
  zipi = FALSE,
  r.threshold= 0.6,
  p.threshold=0.05,
  method = "pearson",
  padj = F,
  n = 4)

res.F = dat.compare.F[[1]]  #######
head(res.F)    ########

## customise the plots
#dat.F = tab.F[[2]] #constructed above   ######
node.F = dat.F$net.cor.matrix$node    #######
edge.F = dat.F$net.cor.matrix$edge    #######
head(edge.F)      #########

#this is the same as p0 created above, but customisable
p.F <- ggplot() + geom_segment(aes(x = X1, y = Y1,  ########
                                 xend = X2, 
                                 yend = Y2, 
                                 color = cor #colour for correlations
),
data = edge.F, #########
size = 0.05, 
alpha = 0.8) +
  geom_point(aes(X1, X2,
                 fill = Phylum,
                 size = igraph.degree),
             pch = 21, 
             data = node.F,  ########## 
             color = "grey80") +
  facet_wrap(.~ label, scales="free_y",  nrow = 1) +
  scale_colour_manual(values = c("black", "blue")) + #colours the correlations
  scale_size(range = c(0.8, 5)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)
  ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(legend.position = "none",  # Removes the legend
        legend.background = element_rect(colour = NA)) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())

p.F     #######

## network properties
# groups are embedded in the data so manually create group specific igraph objects (requires igraph and dplyr)
# extract node and edge data
node_df.F <- tab.F[[2]]$net.cor.matrix$node
edge_df.F <- tab.F[[2]]$net.cor.matrix$edge

# Split by group 
groups.F <- unique(node_df.F$Group)  

# Function to build igraph and compute properties
get_net_props.F <- function(group_id) {
  nodes_group.F <- node_df.F %>% filter(Group == group_id)
  edges_group.F <- edge_df.F %>% filter(Group == group_id)
  
  g.F <- graph_from_data_frame(d = edges_group.F[, c("OTU_1", "OTU_2", "weight")],
                               vertices = nodes_group.F[, c("elements", "Phylum")],
                               directed = FALSE)
  
  props.F <- net_properties.4(g.F, n.hub = TRUE)
  props_df.F <- as.data.frame(t(props.F))  # corrected here
  props_df.F$Group <- group_id
  return(props_df.F)
}

# Apply to each group
net_props_list.F <- lapply(groups.F, get_net_props.F)

# Combine into one dataframe
net_props_df.F <- bind_rows(net_props_list.F)

# View result
print(net_props_df.F)

## GGCLUSTER code to get network properties
### Modified version to analyze a second dataset (e.g., fungi or diet), using `.F` suffix
# Assume: cortab.F is your second correlation list
# Assume: ps.end.F is your second phyloseq object

# Initialize an empty data frame to hold the results
dat.f2.F <- NULL

# Loop over each group in the second dataset
for (i in seq_along(names(cortab.F))) {
  group_name.F <- names(cortab.F)[i]
  
  # Subset the phyloseq object by group
  pst.F <- ps.fungi.end %>%
    subset_samples.wt("Soil_type", group_name.F) %>%  # adjust if using a different grouping var
    subset_taxa.wt("OTU", colnames(cortab.F[[group_name.F]])) %>%
    filter_taxa(function(x) sum(x) > 0, TRUE) %>%
    scale_micro("rela")  # Normalize to relative abundance
  
  # Compute sample-level network properties
  dat.f.F <- netproperties.sample(pst = pst.F, cor = cortab.F[[group_name.F]])
  
  # Add group label
  dat.f.F$Group <- group_name.F
  
  # Combine results across groups
  if (is.null(dat.f2.F)) {
    dat.f2.F <- dat.f.F
  } else {
    dat.f2.F <- bind_rows(dat.f2.F, dat.f.F)
  }
}

# Examine output
head(dat.f2.F)
summary(dat.f2.F)

# Convert Group to factor
dat.f2.F$Group <- as.factor(dat.f2.F$Group)

# Plot (adjust color scheme if needed)
ggplot(dat.f2.F, aes(x = Group, y = `the.number.of.keystone.nodes`, fill = Group)) +
  geom_boxplot(outlier.shape = NA, 
               alpha = 0.7, 
               width = 0.5, 
               color = "black", 
               na.rm = TRUE) +
  geom_jitter(width = 0.2, 
              alpha = 0.6, 
              color = "black", 
              na.rm = TRUE) +
  labs(title = " ", fill = "Group") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#


################################# co-analysis networks
ps16s = ps.end %>% ggClusterNet::scale_micro()
psCOI = ps.fungi.end %>% ggClusterNet::scale_micro()

ps16s
psCOI

# Get the shared sample names
shared_samples <- intersect(sample_names(ps16s), sample_names(psCOI))

# Subset both phyloseq objects to those shared samples
ps16s.shared <- prune_samples(shared_samples, ps16s)
psCOI.shared  <- prune_samples(shared_samples, psCOI)

common_cols <- intersect(colnames(sample_data(ps16s)), colnames(sample_data(psCOI)))

sample_data(ps16s.shared) <- sample_data(ps16s.shared)[, common_cols]
sample_data(psCOI.shared) <- sample_data(psCOI.shared)[, common_cols]


#merge
ps.merge <- ggClusterNet::merge16S_ITS(ps16s = ps16s.shared,
                                       psITS = psCOI.shared,
                                       NITS = 3900,
                                       N16s = 1700)

map =  phyloseq::sample_data(ps.merge)
head(map)
map$Group = "one"
phyloseq::sample_data(ps.merge) <- map

## add env data
#data1 = env
#envRDA  = data.frame(row.names = env$ID,env[,-1])

#envRDA.s = vegan::decostand(envRDA,"hellinger")
#data1[,-1] = envRDA.s

#Gru = data.frame(ID = colnames(env)[-1],group = "env" )
#head(Gru)
result.bip <- ggClusterNet::corBionetwork(ps = ps.merge,
                                          N = 0,
                                          r.threshold = 0.6,
                                          p.threshold = 0.05,
                                          big = T,
                                          group = "zone", #separate networks for metadata
                                          #env = data1, # environment data
                                          #envGroup = Gru,# environment grouping factor
                                          #layout = "fruchtermanreingold",
                                          #path = "/scratch/project_2007483/combined_analysis/FINAL", # where to store data
                                          fill = "Phylum", 
                                          size = "igraph.degree", # date used to scale point sizes
                                          scale = TRUE, # relative abundance normalisation
                                          bio = TRUE, # bipartite network
                                          zipi = F, 
                                          step = 500, # number of random network samples
                                          width = 18,
                                          label = TRUE,
                                          height = 10)





























### modified the code in GGCLUSTERNET
#cortab <- tab.fungi[[2]]$net.cor.matrix$cortab # already have this

# check names
names(cortab)

#requires phyloseq, dplyr, ggClusterNet
# initialize an empty data frame to hold the results
dat.f2 <- NULL

# loop over each group
for (i in seq_along(names(cortab))) {
  group_name <- names(cortab)[i]
  # subset the phyloseq object by sample group
  pst <- ps.end %>%
    subset_samples.wt("Soil_type", group_name) %>%  # your grouping variable
    subset_taxa.wt("OTU", colnames(cortab[[group_name]])) %>%
    filter_taxa(function(x) sum(x) > 0, TRUE) %>%
    scale_micro("rela")  # normalize to relative abundance
  # compute sample-level network properties
  dat.f <- netproperties.sample(pst = pst, cor = cortab[[group_name]])
  # add group name to the results
  dat.f$Group <- group_name
  # combine with previous results
  if (is.null(dat.f2)) {
    dat.f2 <- dat.f
  } else {
    dat.f2 <- bind_rows(dat.f2, dat.f)
  }
}

#examine output
head(dat.f2)
summary(dat.f2)

# plot variables in the data
#check variables
colnames(dat.f2)
# check the data type of the columns you're working with
#str(dat.f2)
# check if the column is truly numeric
#is.numeric(dat.f2$"no.clusters")
# if it returns FALSE, check for non-numeric values
#unique(dat.f2$"no.clusters")
# convert all columns except 'Group' to numeric
#dat.f2_fixed <- dat.f2 %>%
#  mutate(across(-Group, ~ as.numeric(as.character(.))))
# convert 'Group' to a factor
dat.f2$Group <- as.factor(dat.f2$Group)

#plot
ggplot(dat.f2, aes(x = Group, y = no.clusters, fill = Group)) +
  geom_boxplot(outlier.shape = NA, 
               alpha = 0.7, 
               width = 0.5, 
               color = "black", 
               na.rm = TRUE) +  # boxplot with custom width and black borders
  geom_jitter(width = 0.2, 
              alpha = 0.6, 
              color = "black", 
              na.rm = TRUE) +  # add jittered points with custom transparency and black points
  labs(title = " ", fill = "Group") +  # title and legend label
  theme_minimal() +  # minimal theme
  theme(axis.title.x = element_blank(),  # remove x-axis title for cleaner look
        axis.title.y = element_text(size = 12, face = "bold"),  # bold y-axis title
        legend.position = "top",  # position the legend at the top
        plot.title = element_text(hjust = 0.5, 
                                  size = 14, 
                                  face = "bold"),  # center and style title
        panel.grid.major = element_blank(),  # remove major grid lines for cleaner look
        panel.grid.minor = element_blank())  # remove minor grid lines

# do groups differ (2 groups)
#wilcox.test(no.clusters ~ Group, data = dat.f2, na.rm = TRUE)

save.image(file = "soil2gut.fecal.r06.RData")
#load(file = "soil2gut.fecal.r06.RData")



library(igraph)
library(dplyr)

### Original dataset (e.g., Bacteria)
node_df <- tab.B[[2]]$net.cor.matrix$node
edge_df <- tab.B[[2]]$net.cor.matrix$edge
groups <- unique(node_df$Group)

get_net_props <- function(group_id) {
  nodes_group <- node_df %>% filter(Group == group_id)
  edges_group <- edge_df %>% filter(Group == group_id)
  
  g <- graph_from_data_frame(d = edges_group[, c("OTU_1", "OTU_2", "weight")],
                             vertices = nodes_group[, c("elements", "Phylum")],
                             directed = FALSE)
  
  props <- net_properties.4(g, n.hub = TRUE)
  props_df <- as.data.frame(t(props))
  props_df$Group <- group_id
  props_df$Dataset <- "Bacteria"
  return(props_df)
}

net_props_list <- lapply(groups, get_net_props)
net_props_df <- bind_rows(net_props_list)

### Second dataset (e.g., Diet or Fungi) with `.F` suffix

node_df.F <- tab.F[[2]]$net.cor.matrix$node
edge_df.F <- tab.F[[2]]$net.cor.matrix$edge
groups.F <- unique(node_df.F$Group)

get_net_props.F <- function(group_id) {
  nodes_group.F <- node_df.F %>% filter(Group == group_id)
  edges_group.F <- edge_df.F %>% filter(Group == group_id)
  
  g.F <- graph_from_data_frame(d = edges_group.F[, c("OTU_1", "OTU_2", "weight")],
                               vertices = nodes_group.F[, c("elements", "Phylum")],
                               directed = FALSE)
  
  props.F <- net_properties.4(g.F, n.hub = TRUE)
  props_df.F <- as.data.frame(t(props.F))
  props_df.F$Group <- group_id
  props_df.F$Dataset <- "Diet"  # or "Fungi", depending on what you're using
  return(props_df.F)
}

net_props_list.F <- lapply(groups.F, get_net_props.F)
net_props_df.F <- bind_rows(net_props_list.F)

### ✅ Optional: Combine both into one if needed (but keep them separate objects too)

net_props_all <- bind_rows(net_props_df, net_props_df.F)

# View the outputs
print(net_props_df)      # Bacteria
print(net_props_df.F)    # Diet / Fungi
print(net_props_all)     # Combined






##### REVISE TEXT ABOVE TO LOOP
library(igraph)
library(dplyr)

# Put your tabulated results into a named list
# For example:
tab_list <- list(Bacteria = tab.B, Fungi = tab.F)

# Create a function to calculate network properties for each dataset
get_net_props_from_tab <- function(tab_result, dataset_name) {
  # Extract node and edge data
  node_df <- tab_result[[2]]$net.cor.matrix$node
  edge_df <- tab_result[[2]]$net.cor.matrix$edge
  
  # Get groups (e.g., treatment groups)
  groups <- unique(node_df$Group)
  
  # Function to compute network properties for each group
  get_props_by_group <- function(group_id) {
    nodes_group <- node_df %>% filter(Group == group_id)
    edges_group <- edge_df %>% filter(Group == group_id)
    
    g <- graph_from_data_frame(d = edges_group[, c("OTU_1", "OTU_2", "weight")],
                               vertices = nodes_group[, c("elements", "Phylum")],
                               directed = FALSE)
    
    props <- net_properties.4(g, n.hub = TRUE)
    props_df <- as.data.frame(t(props))
    props_df$Group <- group_id
    props_df$Dataset <- dataset_name  # Add identifier for the dataset (bacteria/fungi/etc)
    return(props_df)
  }
  
  # Apply to all groups in the dataset
  net_props_list <- lapply(groups, get_props_by_group)
  net_props_df <- bind_rows(net_props_list)
  return(net_props_df)
}

# Run the function for all datasets in the list
all_net_props <- lapply(names(tab_list), function(name) {
  get_net_props_from_tab(tab_list[[name]], dataset_name = name)
})

# Combine everything into one big dataframe
net_props_combined <- bind_rows(all_net_props)

# View final result
print(net_props_combined)







######DUMP
tab.fungi <- network.pip(
  ps = ps.fungi.end,
  N = 0,
  ra = 0.05,
  big = TRUE,
  select_layout = TRUE,
  layout_net = "model_igraph",
  r.threshold = 0.6,
  p.threshold = 0.05,
  method = "pearson",
  label = FALSE,
  group = "Soil_type",
  fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 500,
  R = 50,
  ncpus = 4
)



# extract images
plot = tab.fungi[[1]]
# extract graph visualisation results
p0 = plot[[1]] #network
# zipi
p0.1 = plot[[2]]
# power law distribution of random networks
p0.2 = plot[[3]] #same as plot

# extract correlation matrix
dat = tab.fungi[[2]]
cortab = dat$net.cor.matrix$cortab

# network significance
dat.compare = module.compare.net.pip(
  ps = NULL,
  corg = cortab,
  degree = TRUE,
  zipi = TRUE,
  r.threshold= 0.6,
  p.threshold=0.05,
  method = "pearson",
  padj = F,
  n = 4)

res = dat.compare[[1]]
head(res)

## customise the plots
#dat = tab.fungi[[2]] #constructed above
node = dat$net.cor.matrix$node
edge = dat$net.cor.matrix$edge
head(edge)

#this is the same as p0 created above, but customisable
p <- ggplot() + geom_segment(aes(x = X1, y = Y1, 
                                 xend = X2, 
                                 yend = Y2, 
                                 color = cor), #colour for correlations
                             data = edge, 
                             size = 0.05, 
                             alpha = 0.8) +
  geom_point(aes(X1, X2,
                 fill = Phylum,
                 size = igraph.degree),
             pch = 21, 
             data = node, 
             color = "grey80") +
  facet_wrap(.~ label, scales="free_y",  nrow = 1) +
  # geom_text_repel(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
  # geom_text(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
  scale_colour_manual(values = c("black", "blue")) + #colours the correlations
  # scale_fill_hue()+
  scale_size(range = c(0.8, 5)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(legend.position = "none",  # Removes the legend
        legend.background = element_rect(colour = NA)) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())

p

## customise ZIPI
dat.z = dat$zipi.data
head(dat.z)

x1 <- c(0, 0.62, 0, 0.62)
x2 <- c(0.62, 1, 0.62, 1)
y1 <- c(-Inf, 2.5, 2.5, -Inf)
y2 <- c(2.5, Inf, Inf, 2.5)

lab <- c("peripheral",'network hubs','module hubs','connectors')
roles.colors <- c("#E6E6FA","#DCDCDC","#F5FFFA", "#FAEBD7")
tab = data.frame(x1 = x1, y1 = y1, x2 = x2, y2 = y2, lab = lab)
tem = dat.z$group %>% unique() %>% length()
for (i in 1:tem) {
  if (i == 1) {
    tab2 = tab
  } else{
    tab2 = rbind(tab2, tab)
  }
}


q <- ggplot() +
  geom_rect(data = tab2,
            mapping = aes(xmin = x1,
                          xmax = x2,
                          ymin = y1,
                          ymax = y2,
                          fill = lab))+
  guides(fill=guide_legend(title = "topological roles")) +
  scale_fill_manual(values = roles.colors)+
  geom_point(data = dat.z, aes(x = p, y = z, color = module)) + theme_bw()+
  guides(color= F) +
  ggrepel::geom_text_repel(data = dat.z,
                           aes(x = p, y = z,
                               color = module, label = label), size = 4)+
  # facet_wrap(.~group) +
  facet_grid(.~ group, scale = 'free') +
  theme(strip.background = element_rect(fill = "white"))+
  xlab("participation coefficient") + ylab("within-module connectivity z-score")

q

## network properties
# groups are embedded in the data so manually create group-specific igraph objects (requires igraph and dplyr)
# extract node and edge data
node_df <- tab.fungi[[2]]$net.cor.matrix$node
edge_df <- tab.fungi[[2]]$net.cor.matrix$edge

# split by group 
groups <- unique(node_df$Group)  

# create a function to build igraph and compute properties
get_net_props <- function(group_id) {
  nodes_group <- node_df %>% filter(Group == group_id)
  edges_group <- edge_df %>% filter(Group == group_id)
  
  # Build igraph
  g <- graph_from_data_frame(d = edges_group[, c("OTU_1", "OTU_2", "weight")],
                             vertices = nodes_group[, c("elements", "Phylum")],
                             directed = FALSE)
  
  # Compute network properties
  props <- net_properties.4(g, n.hub = TRUE)
  props_df <- as.data.frame(t(props))  # Convert to a single-row dataframe
  props_df$Group <- group_id
  return(props_df)
}

# apply to each group
net_props_list <- lapply(groups, get_net_props)

# combine into one data frame
net_props_df <- bind_rows(net_props_list)

# view result
print(net_props_df)

## modified code in GGCLUSTERNET to get network properties
#cortab <- tab.fungi[[2]]$net.cor.matrix$cortab # already have this
# check names
names(cortab)

# requires phyloseq, dplyr, ggClusterNet
# initialize an empty data frame to hold the results
dat.f2 <- NULL

# loop over each group
for (i in seq_along(names(cortab))) {
  group_name <- names(cortab)[i]
  # subset the phyloseq object by sample group
  pst <- ps.fungi.end %>%
    subset_samples.wt("Soil_type", group_name) %>%  # grouping variable
    subset_taxa.wt("OTU", colnames(cortab[[group_name]])) %>%
    filter_taxa(function(x) sum(x) > 0, TRUE) %>%
    scale_micro("rela")  # normalize to relative abundance
  # compute sample-level network properties
  dat.f <- netproperties.sample(pst = pst, cor = cortab[[group_name]])
  # add group name to the results
  dat.f$Group <- group_name
  # combine with previous results
  if (is.null(dat.f2)) {
    dat.f2 <- dat.f
  } else {
    dat.f2 <- bind_rows(dat.f2, dat.f)
  }
}

#examine output
head(dat.f2)
summary(dat.f2)

# plot variables in the data
colnames(dat.f2)
dat.f2$Group <- as.factor(dat.f2$Group)

r <- ggplot(dat.f2, aes(x = Group, y = no.clusters, fill = Group)) +
  geom_boxplot(outlier.shape = NA, 
               alpha = 0.7, 
               width = 0.5, 
               color = "black", 
               na.rm = TRUE) +  # boxplot with custom width and black borders
  geom_jitter(width = 0.2, 
              alpha = 0.6, 
              color = "black", 
              na.rm = TRUE) +  # add jittered points with custom transparency and black points
  labs(title = " ", fill = "Group") +  # title and legend label
  theme_minimal() +  # minimal theme
  theme(axis.title.x = element_blank(),  # remove x-axis title for cleaner look
        axis.title.y = element_text(size = 12, face = "bold"),  # bold y-axis title
        legend.position = "top",  # position the legend at the top
        plot.title = element_text(hjust = 0.5, 
                                  size = 14, 
                                  face = "bold"),  # center and style title
        panel.grid.major = element_blank(),  # remove major grid lines for cleaner look
        panel.grid.minor = element_blank())  # remove minor grid lines

r

# do groups differ (2 groups)
#wilcox.test(no.clusters ~ Group, data = dat.f2, na.rm = TRUE)

save.image(file = "soil2gut.fecal.bact.fungi.con.r06p05.1.RData")
#load(file = "soil2gut.fecal.bact.fungi.con.r06p05.RData")

## compute sample-based network properties over multiple networks
for (i in 1:length(names(cortab))) {
  pst = ps.fungi.end %>%   ##pst has been calculated before
    subset_samples.wt("Soil_type", names(cortab)[i]) %>%
    subset_taxa.wt("OTU", colnames(cortab[[names(cortab)[i]]])) %>%
    filter_taxa(function(x) sum(x) > 0, TRUE) %>%
    scale_micro("rela")
  
  dat.f = netproperties.sample(pst = pst, cor = cortab[[names(cortab)[i]]])  ##pst used here
  head(dat.f)
  if (i == 1) {
    dat.f2 = dat.f  ##and datf2 created previously also
  }else{
    dat.f2 = rbind(dat.f2, dat.f)
  }
  
}


