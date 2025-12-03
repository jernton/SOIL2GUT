library(tidyverse) # For utilities
library(file2meco) # meco<->physeq data conversions
library(ggbreak) # Axis breaks for ggplot
library(glue) # String manipulation for filenames
library(microViz) # fix taxonomy table ranks
library(scales) # ggplot text wrap functionality
library(phyloseq)
library(ANCOMBC)


# ANCOMBC
## Define taxonomy function
Define a taxonomy editing function to fix overlap issues before plotting heatmaps.
# Function to fill missing fields and prevent repetitive "_unclassified" suffixes
# Function to fill missing taxonomy fields and clean unclassified labels
# Apply this function to taxonomy tables during results extraction below 
fix_taxonomy_row <- function(taxa_row) {
  new_row <- taxa_row
  last_valid <- NA

  for (i in seq_along(taxa_row)) {
    current <- taxa_row[i]

    # Check if the field is missing: pattern "__" with no classification
    if (grepl("^.__$", current)) {
      if (!is.na(last_valid)) {
        new_row[i] <- paste0(last_valid, "_unclassified")
      } else {
        new_row[i] <- "Unclassified"
      }
    } else {
      last_valid <- current
    }
  }

  # Remove repeated _unclassified suffixes (e.g. "_unclassified_unclassified")
  new_row <- gsub("(_unclassified)+", "_unclassified", new_row)

  return(new_row)
}



# ANCOMBC Fungi

## Model

#Convert to phyloseq so I can run that object externally in ANCOMBCII.
#ATTN: use non-normalized data for DA analyses! Use lib_cut in ancombc to
#filter out low libsize samples instead

ps <- meco2phyloseq(dsFFecT2)
tax_table(ps) <- tax_fix(ps)
head(tax_table(ps))


rank_level <- "Species" # change accordingly between Family, Genus, Species level to set ancombc2() function and output filenames

set.seed(42)
output = ancombc2(data = ps,
                  fix_formula = "Soil_type",
                  rank = rank_level,
                  prv_cut = 0.1, lib_cut = 8000,
                  group = "Soil_type", struc_zero = TRUE, neg_lb = TRUE, # neg_lb=TRUE recommended when sample size >30
                  alpha = 0.05, verbose = FALSE,
                  global = TRUE, pairwise = TRUE,
                  dunnet = FALSE)



# Create taxa object for attaching taxa info to results
taxa.df <- dsF$tax_table %>% as.data.frame() %>%  rownames_to_column("ASV")


# Apply the taxonomy fixing function row-wise to taxa.df to avoid overlap of unclassified taxa in the figure
taxa.df.fixed <- as.data.frame(t(apply(taxa.df, 1, fix_taxonomy_row)))

taxa.df.fixed$Species <- make.unique(taxa.df.fixed$Species, sep = "_")


# Replace original
taxa.df <- taxa.df.fixed

# Extract result pairwise result.
res_pair = output$res_pair %>%
    mutate_if(is.numeric, function(x) round(x, 3)) # round some decimals
colnames(res_pair)[1] <- "ASV" # change ASV column name to match with taxa.df
res_pair <- left_join(
    x = res_pair,
    y = taxa.df,
    by = "ASV" # Optional: Attach taxonomy data to the table, but this is not necessary.
  )
write.table(res_pair, file = glue("ANCOMres/ANCOMBC2_fungi_97clust_pair_{rank_level}.txt"), sep = "\t", row.names = FALSE)

```

## Pairwise heatmap

```{r}
rank_level <- "Species" # Will affect output filenames
res_pair = read.table(glue("ANCOMres/ANCOMBC2_fungi_97clust_pair_{rank_level}.txt"), sep = "\t", header=TRUE)
colnames(res_pair)
colnames(res_pair)[1] <- "taxon" # Set the name vector to be used in the heatmap by adjusting the col number 

df_fig_pair1 = res_pair %>%
  dplyr::filter(diff_Soil_typePark == 1 |
                  diff_Soil_typeUrban == 1 | 
                  diff_Soil_typeUrban_Soil_typePark == 1) %>%
  dplyr::mutate(lfc1 = ifelse(diff_Soil_typePark == 1, 
                              round(lfc_Soil_typePark, 2), 0),
                lfc2 = ifelse(diff_Soil_typeUrban == 1, 
                              round(lfc_Soil_typeUrban, 2), 0),
                lfc3 = ifelse(diff_Soil_typeUrban_Soil_typePark == 1, 
                              round(lfc_Soil_typeUrban_Soil_typePark, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)



df_fig_pair1$group = recode(df_fig_pair1$group, 
                           `lfc1` = "Park - Control",
                           `lfc2` = "Urban - Control",
                           `lfc3` = "Urban - Park")
df_fig_pair1$group = factor(df_fig_pair1$group, 
                           levels = c("Park - Control",
                                      "Urban - Control", 
                                      "Urban - Park"))


lo = floor(min(df_fig_pair1$value))
up = ceiling(max(df_fig_pair1$value))
fig_pair = df_fig_pair1 %>%
  #ggplot(aes(x = group, y = taxon, fill = value)) + 
  ggplot(aes(x = group, y = factor(taxon, rev(levels(factor(taxon)))), fill = value)) + #code to make names into factor and reverse them on y axis
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#1b9e77" , high = "#d95f02", mid = "white",                        
                       na.value = "white", midpoint = 0, limit = c(lo, up),
                       name = NULL) +
  # geom_text(aes(group, taxon, label = value, color = 'black'), size = 4) + #changed this parameter as colour not defined in sensitivity analysis
  geom_text(aes(group, taxon, label = round(value, 2), color = 'black'), size = 5) + 
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = " ") +
  coord_fixed() + #should keep cells as squares - but maybe overridden by aspect ratio
  scale_x_discrete(labels = label_wrap(10)) +
  #scale_y_reverse() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

fig_pair.1 <- fig_pair +
  theme(
    aspect.ratio = 2/1, 
    legend.title = element_text(color = "black", size = 12, vjust = 0.5),
    legend.text = element_text(size = 14),
    legend.position = "right",
    panel.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.title.x = element_text(color = "black", size = 16, vjust = -0.5),
    axis.title.y = element_text(color = "black", size = 16, vjust = 0.5),
    axis.text.x = element_text(angle = 0, size = 14, vjust = 0.5),
    axis.text.y = element_text(angle = 360, size = 14, vjust = 0.5),
    plot.margin = ggplot2::margin(t = 0.0, r = 0.0, b = 0.0, l = 0.0, unit = "cm")
  )
fig_pair.1
ggsave(glue("Figure_ANCOMBC2_heatmap_fungi_pairwise_{rank_level}.png"), width=6, height = 4, units = "in")


```