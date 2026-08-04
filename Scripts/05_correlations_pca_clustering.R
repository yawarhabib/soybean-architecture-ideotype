###############################################################################################
###               Multi-trait selection defines a plant-architecture ideotype               ###
###                           for soybean in the Russian Far East                           ###
###                                                                                         ###
###                                       ----------                                        ###
###           Script 05: Trait correlations, principal components and clustering            ###
###                                                                                         ###
###                            R script by Yawar Habib, Skoltech                            ###
###           Center for Bio and Medical Technologies, ORCID 0009-0008-1812-3910            ###
###############################################################################################

######### PURPOSE #############################################################################
## Examine the relationships among the six architecture traits and the multivariate structure
## of the panel. Phenotypic correlations are computed among genotype means across the three
## years and genotypic-value correlations among the genotypic BLUPs. The BLUPs are then
## subjected to principal component analysis on the correlation matrix, and the genotypes are
## grouped by hierarchical clustering of their Euclidean distances using Ward's
## minimum-variance method.

######### DATA ################################################################################
## data/Architecture_63genotypes_2018_2020_clean.csv  (see script 01)
##
## The genotypic BLUPs are recomputed here from the same mixed model used in script 04, so
## this script is self-contained and does not depend on script 04 having been run.

######### OUTPUT ##############################################################################
## Figure 5       output/figures/Figure_correlation_panel.tiff
## Figure 6       output/figures/Figure_PCA_panel.tiff
## Supplementary  output/tables/Supplementary_phenotypic_correlation.csv
## Supplementary  output/tables/Supplementary_genotypic_correlation.csv
## Supplementary  output/tables/Table_cluster_profiles.csv
## Supplementary  output/figures/Supplementary_cluster_validation.tiff

######### PREPARATION OF THE WORKING INTERFACE IN R ###########################################
### I. Set working directory
# Set the working directory to the root of this repository, i.e. the folder that contains
# the data/, Scripts/ and output/ folders.
# On RStudio: tab 'Session' -> Set Working Directory -> Choose Directory.

### II. R packages needed for the analysis
# dplyr, tidyr, lme4, GGally, factoextra, cluster, ggplot2, ggrepel, patchwork

### III. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(lme4)
library(GGally)
library(factoextra)
library(cluster)
library(ggplot2)
library(ggrepel)
library(patchwork)

dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

Architecture_data <- read.csv("data/Architecture_63genotypes_2018_2020_clean.csv", header = TRUE) %>%
  mutate(Year = as.factor(Year), Replication = as.factor(Replication), Genotype = as.factor(Genotype))

######### CANONICAL TRAIT DEFINITIONS #########################################################
architecture_traits <- c("Height", "First_pod_height", "No_branches",
                         "No_ProdNodes_MainStem", "No_ProdNodes_LateralBranches",
                         "Vegetation_period")
trait_abbreviations <- c(Height="PH", First_pod_height="FPH", No_branches="NB",
                         No_ProdNodes_MainStem="PNMS", No_ProdNodes_LateralBranches="PNLB",
                         Vegetation_period="DM")

######### GENOTYPIC VALUES (BLUPS) ############################################################
## 1. BLUPs (genotypic values) — shared basis for everything below
blup_values <- sapply(architecture_traits, function(trait_name) {
  m <- lmer(as.formula(paste0(trait_name, " ~ (1|Genotype) + (1|Year) + (1|Genotype:Year)")),
            data = Architecture_data)
  ranef(m)$Genotype[, "(Intercept)"]
})
blup_values <- as.data.frame(blup_values)
colnames(blup_values) <- architecture_traits
rownames(blup_values) <- levels(Architecture_data$Genotype)

######### TRAIT CORRELATIONS ##################################################################
## 2. Correlations: phenotypic (genotype means) + genotypic-value (BLUPs)
genotype_means <- Architecture_data %>%
  group_by(Genotype) %>%
  summarise(across(all_of(architecture_traits), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

phenotypic_correlation        <- cor(genotype_means[architecture_traits], use = "complete.obs")
genotypic_value_correlation   <- cor(blup_values, use = "complete.obs")

cat("\n--- Phenotypic correlation (genotype means) ---\n"); print(round(phenotypic_correlation, 2))
cat("\n--- Correlation of genotypic values (BLUPs) ---\n");  print(round(genotypic_value_correlation, 2))

## abbreviated correlation tables (for supplementary)
phenotypic_correlation_table <- phenotypic_correlation
dimnames(phenotypic_correlation_table) <- list(trait_abbreviations[rownames(phenotypic_correlation)],
                                               trait_abbreviations[colnames(phenotypic_correlation)])
genotypic_correlation_table <- genotypic_value_correlation
dimnames(genotypic_correlation_table) <- list(trait_abbreviations[rownames(genotypic_value_correlation)],
                                              trait_abbreviations[colnames(genotypic_value_correlation)])

######### FIGURE 5 - CORRELATIONS AMONG GENOTYPIC VALUES ######################################
## correlation panel figure (ggpairs on BLUPs)
blup_plot_data <- blup_values
colnames(blup_plot_data) <- trait_abbreviations[colnames(blup_plot_data)]

correlation_panel_figure <- ggpairs(
  blup_plot_data,
  upper = list(continuous = wrap("cor", size = 4, color = "#9B2226")),
  lower = list(continuous = wrap("points", alpha = 0.55, size = 1.3, colour = "#C75D5D")),
  diag  = list(continuous = wrap("densityDiag", fill = "#C75D5D", colour = "#9B2226", alpha = 0.7)),
  axisLabels = "show"
) +
  theme_classic(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        strip.background = element_rect(colour = "grey8", linewidth = 0.3),
        strip.text = element_text(size = 11, colour = "grey8"),
        panel.border = element_rect(colour = "grey8", fill = NA, linewidth = 0.4),
        axis.text = element_text(size = 7, colour = "grey50"))

######### PRINCIPAL COMPONENT ANALYSIS ########################################################
## 3. PCA on BLUPs (scaled)
blup_display <- blup_values
colnames(blup_display) <- trait_abbreviations[colnames(blup_display)]
architecture_pca_display <- prcomp(blup_display, center = TRUE, scale. = TRUE)

variance_explained <- summary(architecture_pca_display)$importance
cat("\n--- Variance explained ---\n"); print(round(variance_explained, 3))
cat("\n--- Loadings ---\n"); print(round(architecture_pca_display$rotation, 3))

trait_groups <- factor(
  c(PH="Main-axis", FPH="Main-axis", PNMS="Main-axis",
    NB="Branching", PNLB="Branching", DM="Phenology")[colnames(blup_display)],
  levels = c("Branching", "Main-axis", "Phenology")
)

## variable coords + scaled individual scores
var_data   <- facto_summarize(architecture_pca_display, element = "var", axes = c(1, 2))
ind_scores <- as.data.frame(architecture_pca_display$x[, 1:2]); colnames(ind_scores) <- c("Dim.1", "Dim.2")
scale_ind  <- 0.9 / max(sqrt(ind_scores$Dim.1^2 + ind_scores$Dim.2^2))
ind_scores$Dim.1 <- ind_scores$Dim.1 * scale_ind
ind_scores$Dim.2 <- ind_scores$Dim.2 * scale_ind

pc1 <- round(variance_explained[2, 1] * 100, 1)
pc2 <- round(variance_explained[2, 2] * 100, 1)

######### HIERARCHICAL CLUSTERING #############################################################
## 4. Clustering (Ward.D2 on standardized BLUPs), k = 2
## k = 2 is the number supported by the silhouette and gap-statistic criteria, both of which
## are computed in the CLUSTER VALIDATION section below.
blup_scaled <- scale(blup_values)
distance_matrix <- dist(blup_scaled, method = "euclidean")
hierarchical_cluster <- hclust(distance_matrix, method = "ward.D2")
cluster_assignment <- cutree(hierarchical_cluster, k = 2)

cat("\n--- Cluster sizes ---\n"); print(table(cluster_assignment))

## raw-trait cluster profiles (interpretable description)
cluster_raw_profile <- Architecture_data %>%
  group_by(Genotype) %>%
  summarise(across(all_of(architecture_traits), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(Cluster = factor(cluster_assignment)) %>%
  group_by(Cluster) %>%
  summarise(across(all_of(architecture_traits), ~ round(mean(.x), 1)), .groups = "drop")
cat("\n--- Cluster means (raw trait values) ---\n"); print(as.data.frame(cluster_raw_profile))

## attach cluster (1 = Compact, 2 = Large) to the PCA individual scores
ind_scores$Cluster <- factor(cluster_assignment, labels = c("Compact", "Large"))

######### FIGURE 6 - PRINCIPAL COMPONENTS AND CLUSTERING ######################################
## PCA panel: (A) scree + (B) clustered biplot
## (A) scree — custom geom_col (avoids fviz_eig warning)
scree_data <- data.frame(
  PC = factor(paste0("PC", seq_len(ncol(variance_explained))), levels = paste0("PC", seq_len(ncol(variance_explained)))),
  Percent = as.numeric(variance_explained[2, ] * 100)
)
scree_plot <- ggplot(scree_data, aes(PC, Percent, group = 1)) +
  geom_col(fill = "#C75D5D", colour = "#9B2226", width = 0.7) +
  geom_line(colour = "#9B2226", linewidth = 0.5) +
  geom_point(colour = "#9B2226", size = 1.8) +
  geom_text(aes(label = sprintf("%.1f%%", Percent)), vjust = -0.6, size = 3.5, colour = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Principal component", y = "Variance explained (%)") +
  theme_classic(base_size = 12) + theme(axis.text = element_text(colour = "black"))

## (B) clustered biplot — grouped trait arrows + genotype points colored by cluster
pca_biplot_clustered <- fviz_pca_var(
  architecture_pca_display, col.var = trait_groups,
  palette = c("#C1466B", "#1F6F8B", "#E09F3E"), geom.var = "arrow", arrowsize = 0.8
) +
  geom_point(data = ind_scores, aes(x = Dim.1, y = Dim.2, fill = Cluster),
             shape = 21, colour = "grey30", size = 2, alpha = 0.85, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Compact" = "#E0C766", "Large" = "#7E5BB0")) +
  geom_text_repel(data = var_data, aes(x = Dim.1, y = Dim.2, label = name),
                  colour = "black", fontface = "bold", size = 4.5, inherit.aes = FALSE,
                  box.padding = 0.4, min.segment.length = 0) +
  labs(title = NULL, colour = "Trait group", fill = "Architectural group",
       x = paste0("PC1 (", pc1, "%)  —  architectural size"),
       y = paste0("PC2 (", pc2, "%)  —  branching vs. main-axis")) +
  theme_classic(base_size = 12) +
  ## the two legends are stacked rather than placed side by side: laid out horizontally
  ## they are wider than the 12-inch canvas and the "Large" label is clipped off the edge
  theme(axis.text = element_text(colour = "black"), legend.position = "bottom",
        legend.box = "vertical", legend.box.just = "left",
        legend.margin = margin(t = 0, b = 0))

pca_panel <- scree_plot + pca_biplot_clustered +
  plot_layout(widths = c(0.85, 1.15)) + plot_annotation(tag_levels = "A")

######### SUPPLEMENTARY FIGURE - CLUSTER VALIDATION ###########################################
## 6. Cluster validation (supplementary): silhouette + gap
## hcut() uses Ward's minimum-variance method on Euclidean distances by default, the same
## criteria as the clustering above.
##
## Both panels are drawn from the criterion values rather than through fviz_nbclust() and
## fviz_gap_stat() directly, because those functions annotate the plot with a dashed line at
## the k chosen by their own decision rule. Here the silhouette is effectively flat between
## k = 2 and k = 4 and the gap statistic does not support subdivision, so neither rule picks
## the k used in the analysis and the annotation would assert an optimum that is not there.
## The curves themselves are unchanged, and both panels run to k = 10.
k_max <- 10

silhouette_data <- fviz_nbclust(blup_scaled, FUN = hcut, method = "silhouette",
                                k.max = k_max)$data %>%
  mutate(k = as.integer(as.character(clusters)))

cat("\n--- Average silhouette width by k ---\n")
print(data.frame(k = silhouette_data$k, silhouette = round(silhouette_data$y, 4)), row.names = FALSE)

p_sil <- ggplot(silhouette_data, aes(k, y)) +
  geom_line(colour = "steelblue", linewidth = 0.6) +
  geom_point(colour = "steelblue", size = 2) +
  scale_x_continuous(breaks = 1:k_max, limits = c(0.7, k_max + 0.3)) +
  labs(title = "Silhouette method", x = "Number of clusters k", y = "Average silhouette width") +
  theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 12))

## the gap statistic is resampled, so the seed is fixed to make the figure reproducible
set.seed(123)
gap_stat <- clusGap(blup_scaled, FUN = hcut, K.max = k_max, B = 100)

gap_data <- as.data.frame(gap_stat$Tab) %>% mutate(k = seq_len(n()))

cat("\n--- Gap statistic by k ---\n")
print(data.frame(k = gap_data$k, gap = round(gap_data$gap, 4), SE = round(gap_data$SE.sim, 4)),
      row.names = FALSE)

p_gap <- ggplot(gap_data, aes(k, gap)) +
  geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.2, colour = "steelblue") +
  geom_line(colour = "steelblue", linewidth = 0.6) +
  geom_point(colour = "steelblue", size = 2) +
  scale_x_continuous(breaks = 1:k_max, limits = c(0.7, k_max + 0.3)) +
  labs(title = "Gap statistic", x = "Number of clusters k", y = "Gap statistic (k)") +
  theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 12))

cluster_validation_panel <- p_sil + p_gap + plot_layout(ncol = 2) + plot_annotation(tag_levels = "A")

######### SAVE OUTPUT #########################################################################
## figures
ggsave("output/figures/Figure_correlation_panel.tiff", correlation_panel_figure,
       width = 9, height = 8, units = "in", dpi = 600, device = "tiff", compression = "lzw")
ggsave("output/figures/Figure_PCA_panel.tiff", pca_panel,
       width = 12, height = 6, units = "in", dpi = 600, device = "tiff", compression = "lzw")
ggsave("output/figures/Supplementary_cluster_validation.tiff", cluster_validation_panel,
       width = 9, height = 4.5, units = "in", dpi = 600, device = "tiff", compression = "lzw")

## tables
write.csv(round(phenotypic_correlation_table, 2),
          "output/tables/Supplementary_phenotypic_correlation.csv", row.names = TRUE)
write.csv(round(genotypic_correlation_table, 2),
          "output/tables/Supplementary_genotypic_correlation.csv",  row.names = TRUE)
write.csv(cluster_raw_profile, "output/tables/Table_cluster_profiles.csv", row.names = FALSE)

cat("\n=== Correlations, PCA and clustering: all figures and tables saved to output/ ===\n")
