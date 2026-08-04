###############################################################################################
###               Multi-trait selection defines a plant-architecture ideotype               ###
###                           for soybean in the Russian Far East                           ###
###                                                                                         ###
###                                       ----------                                        ###
###                Script 06: Multi-trait selection - MGIDI and MTSI                        ###
###                                                                                         ###
###                            R script by Yawar Habib, Skoltech                            ###
###           Center for Bio and Medical Technologies, ORCID 0009-0008-1812-3910            ###
###############################################################################################

######### PURPOSE #############################################################################
## Select genotypes towards the plant-architecture ideotype with the multi-trait
## genotype-ideotype distance index (MGIDI), and identify those combining that proximity with
## year-to-year stability using the multi-trait stability index (MTSI).
##
## Ideotype    : c("h", "h", "h", "h", "h", "l") in the order PH, FPH, NB, PNMS, PNLB, DM.
##               Increase is sought for the five architectural traits and decrease for days
##               to maturity.
## MGIDI weights: 0.6 on the five size traits, 1 on days to maturity. The sensitivity of the
##               selection differentials to this choice is swept from 0.5 to 1.0 below.
## MTSI        : equal weighting across traits; mtsi() takes no weights argument.
## Selection intensity: 15% for both indices.

######### DATA ################################################################################
## data/Architecture_63genotypes_2018_2020_clean.csv  (see script 01)
##
## The genotypic BLUPs used for the ideotype profile in Figure 9B are computed here from the
## same mixed model as script 04, so this script does not read any table produced elsewhere.

######### OUTPUT ##############################################################################
## Table 4        output/tables/Table_MGIDI_factor_analysis.csv
## Table 5        output/tables/Table_MGIDI_selection_differential.csv
## Figure 8       output/figures/Figure_MGIDI_ranking.tiff
## Figure 9       output/figures/Figure_MGIDI_MTSI_ideotype.tiff
## Supplementary  output/tables/Supplementary_MGIDI_weight_sweep.csv
## Supplementary  output/tables/Supplementary_MGIDI_full_ranking.csv
## Supplementary  output/tables/Supplementary_MTSI_factor_analysis.csv
## Supplementary  output/tables/Supplementary_MTSI_selection_differential.csv
## Supplementary  output/tables/Supplementary_MTSI_full_ranking.csv
## Supplementary  output/tables/Table_MGIDI_MTSI_comparison.csv
## Supplementary  output/figures/Supplementary_MGIDI_contribution_bars.tiff

######### PREPARATION OF THE WORKING INTERFACE IN R ###########################################
### I. Set working directory
# Set the working directory to the root of this repository, i.e. the folder that contains
# the data/, Scripts/ and output/ folders.
# On RStudio: tab 'Session' -> Set Working Directory -> Choose Directory.

### II. R packages needed for the analysis
# dplyr, tidyr, lme4, metan, ggplot2, ggrepel, patchwork

### III. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(lme4)
library(metan)
library(ggplot2)
library(ggrepel)
library(patchwork)

dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

######### GENOTYPE LABELS #####################################################################
## ---- genotype label formatting (data-file names -> manuscript names) ----
## Applied to the raw data BEFORE the model is fitted, so every downstream object
## (model, sel_gen, sel_dif, MGIDI ranking, the BLUPs and all plots) carries the clean names.
## The forms below are those of Supplementary Table S3, the accession list, so that every
## label in Figures 8 and 9 matches the name under which the accession is catalogued.
pretty_geno <- function(x) {
  x <- gsub("_", " ", x)                       # underscores become spaces
  x <- sub("^No 700$",      "No. 700",      x)  # the "No." abbreviation keeps its point
  x <- sub("^No 368 Bila$", "No. 368 Bila", x)
  x <- sub("^No 53 83$",    "No. 53-83",    x)
  x <- sub("^907 66$",      "907-66",       x)  # numeric line codes take hyphens
  x <- sub("^1013 1 18$",   "1013-1-18",    x)
  x <- sub("^1218 4 4$",    "1218-4-4",     x)
  x <- sub("^1037 2$",      "1037-2",       x)
  x <- sub("^A 5 15$",      "A-5-15",       x)  # breeding-line codes likewise
  x <- sub("^D 5 15$",      "D-5-15",       x)
  x <- sub("^F 5 15$",      "F-5-15",       x)
  x <- sub("^L 74 10$",     "L-74-10",      x)
  x <- sub("^N m 4888$",    "N-m-4888",     x)
  x <- sub("^Gaterslebener stamm 54 4$", "Gaterslebener stamm 54-4", x)   # distinct accession
  x
}

Architecture_data <- read.csv("data/Architecture_63genotypes_2018_2020_clean.csv", header = TRUE) %>%
  mutate(Year = as.factor(Year), Replication = as.factor(Replication),
         Genotype = as.factor(pretty_geno(Genotype)))

architecture_traits <- c("Height", "First_pod_height", "No_branches",
                         "No_ProdNodes_MainStem", "No_ProdNodes_LateralBranches",
                         "Vegetation_period")

trait_abbreviations <- c(
  Height = "PH", First_pod_height = "FPH", No_branches = "NB",
  No_ProdNodes_MainStem = "PNMS", No_ProdNodes_LateralBranches = "PNLB",
  Vegetation_period = "DM"
)

######### MGIDI - MULTI-TRAIT GENOTYPE-IDEOTYPE DISTANCE INDEX ################################
## ideotype/weight order: PH, FPH, NB, PNMS, PNLB, DM
mgidi_model <- gamem_met(
  .data = Architecture_data, env = Year, gen = Genotype, rep = Replication,
  resp = c(Height, First_pod_height, No_branches,
           No_ProdNodes_MainStem, No_ProdNodes_LateralBranches, Vegetation_period),
  random = "gen", verbose = TRUE
)

mgidi_final <- mgidi(
  mgidi_model,
  ideotype = c("h", "h", "h", "h", "h", "l"),   # all higher except DM
  weights  = c(0.6, 0.6, 0.6, 0.6, 0.6, 1),     # 0.6 balance point; DM full
  SI = 15
)

cat("\n--- MGIDI selected genotypes ---\n"); print(mgidi_final$sel_gen)
cat("\n--- MGIDI selection differential ---\n"); print(mgidi_final$sel_dif)

## ---- MGIDI tables ----
mgidi_fa_table <- mgidi_final$FA %>%
  transmute(Trait = trait_abbreviations[VAR],
            FA1 = round(FA1, 2), FA2 = round(FA2, 2),
            Communality = round(Communality, 2), Uniqueness = round(Uniquenesses, 2))

mgidi_seldiff_table <- mgidi_final$sel_dif %>%
  transmute(Trait = trait_abbreviations[VAR], Factor,
            Xo = round(Xo, 2), Xs = round(Xs, 2),
            SD = round(SD, 2), `SD%` = round(SDperc, 1),
            Desired = ifelse(sense == "increase", "increase", "decrease")) %>%
  arrange(Factor, Trait)

mgidi_full_ranking <- mgidi_final$MGIDI %>%
  mutate(Rank = row_number(), MGIDI = round(MGIDI, 3)) %>%
  select(Rank, Genotype, MGIDI)

######### MGIDI WEIGHT SENSITIVITY ############################################################
## The size traits are positively correlated with later maturity, so the weighting placed on
## them determines whether the selected group improves in size or in earliness. The full
## sweep is reported so the choice of 0.6 can be judged rather than taken on trust.
size_weight_levels <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
weight_sweep_table <- lapply(size_weight_levels, function(w) {
  sd_tab <- mgidi(mgidi_model, ideotype = c("h","h","h","h","h","l"),
                  weights = c(w, w, w, w, w, 1), SI = 15, verbose = FALSE)$sel_dif
  data.frame(
    Size_weight = w,
    DM   = round(sd_tab$SDperc[sd_tab$VAR == "Vegetation_period"], 1),
    PH   = round(sd_tab$SDperc[sd_tab$VAR == "Height"], 1),
    FPH  = round(sd_tab$SDperc[sd_tab$VAR == "First_pod_height"], 1),
    PNMS = round(sd_tab$SDperc[sd_tab$VAR == "No_ProdNodes_MainStem"], 1),
    NB   = round(sd_tab$SDperc[sd_tab$VAR == "No_branches"], 1),
    PNLB = round(sd_tab$SDperc[sd_tab$VAR == "No_ProdNodes_LateralBranches"], 1)
  )
}) %>% bind_rows()

weight_sweep_table

######### FIGURE 8 - MGIDI RANKING AND TRAIT CONTRIBUTION #####################################
mgidi_rank_plot <- plot(
  mgidi_final, SI = 15, radar = TRUE, arrange.label = TRUE,
  size.point = 3.5, size.text = 3.2, col.sel = "forestgreen", col.nonsel = "red3"
) +
  coord_polar(clip = "off") +                                  # let labels run past the panel
  labs(x = NULL, y = NULL) +                                   # axis titles removed
  theme(axis.title  = element_blank(),                         # "Genotypes" and the index title
        axis.text.y = element_blank(),                         # radial tick numbers 1, 2, 3
        axis.ticks  = element_blank(),
        plot.margin = margin(24, 52, 16, 58))   # extra on the left: longest label sits at 9 o'clock

mgidi_contribution_radar <- suppressWarnings(
  plot(mgidi_final, type = "contribution", genotypes = "selected", size.text = 3.6) +
    scale_colour_manual(values = c("#C75D5D", "#5A8FAE")) +
    scale_fill_manual(values  = c("#C75D5D", "#5A8FAE")) +
    coord_polar(clip = "off") +
    labs(title = NULL) +                                       # caption carries this
    theme(plot.margin = margin(18, 40, 18, 40))
)

mgidi_contribution_bars <- suppressWarnings(
  plot(mgidi_final, type = "contribution", genotypes = "selected", radar = FALSE, size.text = 4) +
    scale_colour_manual(values = c("#C75D5D", "#5A8FAE")) +
    scale_fill_manual(values  = c("#C75D5D", "#5A8FAE")) +
    labs(title = NULL, fill = NULL, colour = NULL) +   # drop the metan default title and
    theme_classic(base_size = 11) +                    # the "name" legend heading, as in the radar
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
          legend.position = "bottom")
)

mgidi_panel <- (mgidi_rank_plot | mgidi_contribution_radar) + plot_annotation(tag_levels = "A")

######### MTSI - MULTI-TRAIT STABILITY INDEX ##################################################
## waasb() defaults to random = "gen", that is, genotype random and year fixed
mtsi_model <- waasb(
  Architecture_data, env = Year, gen = Genotype, rep = Replication,
  resp = c(Height, First_pod_height, No_branches,
           No_ProdNodes_MainStem, No_ProdNodes_LateralBranches, Vegetation_period),
  mresp = c("h", "h", "h", "h", "h", "l"), verbose = TRUE
)

mtsi_index <- mtsi(mtsi_model, index = "waasby", SI = 15)

cat("\n--- MTSI selected genotypes ---\n"); print(mtsi_index$sel_gen)

## ---- MTSI tables ----
mtsi_fa_table <- mtsi_index$FA %>%
  transmute(Trait = trait_abbreviations[VAR],
            FA1 = round(FA1, 2), FA2 = round(FA2, 2), FA3 = round(FA3, 2),
            Communality = round(Communality, 2), Uniqueness = round(Uniquenesses, 2))

mtsi_seldiff_table <- mtsi_index$sel_dif_trait %>%
  transmute(Trait = trait_abbreviations[VAR], Factor,
            Xo = round(Xo, 2), Xs = round(Xs, 2),
            `SD%` = round(SDperc, 1), `SG%` = round(SGperc, 1),
            Desired = ifelse(sense == "increase", "increase", "decrease"))

mtsi_full_ranking <- mtsi_index$MTSI %>%
  transmute(Genotype, MTSI = round(MTSI, 3))

######### MGIDI AGAINST MTSI ##################################################################
mgidi_selected <- mgidi_final$sel_gen
mtsi_selected  <- mtsi_index$sel_gen

comparison_data <- mgidi_full_ranking %>%
  select(Genotype, MGIDI) %>%
  left_join(mtsi_full_ranking, by = "Genotype") %>%
  mutate(
    sel_MGIDI = Genotype %in% mgidi_selected,
    sel_MTSI  = Genotype %in% mtsi_selected,
    Selection = case_when(
      sel_MGIDI & sel_MTSI ~ "Both",
      sel_MGIDI            ~ "MGIDI only",
      sel_MTSI             ~ "MTSI only",
      TRUE                 ~ "Neither"
    )
  )

cat("\n--- Selected by BOTH ---\n")
print(comparison_data$Genotype[comparison_data$Selection == "Both"])
cat("\n--- Counts ---\n"); print(table(comparison_data$Selection))

######### GENOTYPIC VALUES FOR THE IDEOTYPE PROFILE ###########################################
## Genotypic BLUPs, 63 x 6, from the same mixed model as script 04. Used for panel B of
## Figure 9, where each trait is scaled 0-1 across all 63 genotypes and oriented so that the
## outer edge of the radar is always the favourable direction.
blup_values <- sapply(architecture_traits, function(trait_name) {
  m <- lmer(as.formula(paste0(trait_name, " ~ (1|Genotype) + (1|Year) + (1|Genotype:Year)")),
            data = Architecture_data)
  ranef(m)$Genotype[, "(Intercept)"]
})
blup <- data.frame(Genotype = levels(Architecture_data$Genotype), blup_values,
                   check.names = FALSE)
colnames(blup)[-1] <- trait_abbreviations[architecture_traits]

######### FIGURE 9 - MULTI-TRAIT SELECTION ####################################################
## Panel A  <- comparison_data (Genotype, MGIDI, MTSI, sel_MGIDI, sel_MTSI, Selection)
## Panel B  <- genotypic BLUPs, 63 x 6

## ---- inputs ----
## Panel A is drawn on ranks rather than on the raw index values. The two indices are both
## strongly compressed at the favourable end (the nine MGIDI-selected genotypes span only
## 0.14 to 0.93), so on the raw scale the selected genotypes overplot in one corner. Ranks
## are evenly spaced, and both axes are reversed so that the favourable direction is up and
## to the right for both indices. The index values themselves are in the full rankings
## deposited as supplementary tables.
cmp <- comparison_data %>%
  mutate(rank_mgidi = rank(MGIDI), rank_mtsi = rank(MTSI))

dual        <- cmp$Genotype[cmp$Selection == "Both"]      # the four dual-selected
trait_order <- c("PH","FPH","NB","PNMS","PNLB","DM")
desired_low <- "DM"                                       # lower BLUP = favourable

## ---- panel A palette ----
col_both <- "#009E73"; col_mgidi <- "#CC79A7"; col_mtsi <- "#0072B2"; col_none <- "grey80"
mgidi_cut <- max(cmp$rank_mgidi[cmp$sel_MGIDI])           # last rank retained by MGIDI
mtsi_cut  <- max(cmp$rank_mtsi[cmp$sel_MTSI])             # last rank retained by MTSI

pA <- ggplot(cmp, aes(rank_mgidi, rank_mtsi)) +
  annotate("rect", xmin = 0.5, xmax = mgidi_cut + 0.5, ymin = 0.5, ymax = mtsi_cut + 0.5,
           fill = col_both, alpha = 0.10) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey80", linewidth = 0.4) +       # equal rank on both indices
  geom_vline(xintercept = mgidi_cut + 0.5, linetype = "dotted", colour = "grey55") +
  geom_hline(yintercept = mtsi_cut + 0.5,  linetype = "dotted", colour = "grey55") +
  geom_point(data = filter(cmp, Selection == "Neither"),
             aes(colour = Selection), size = 2.2, alpha = 0.85) +
  geom_point(data = filter(cmp, Selection != "Neither"),
             aes(colour = Selection), size = 3.8, alpha = 0.95) +
  scale_colour_manual(values = c("Both" = col_both, "MGIDI only" = col_mgidi,
                                 "MTSI only" = col_mtsi, "Neither" = col_none),
                      breaks = c("Both", "MGIDI only", "MTSI only", "Neither"),
                      labels = c("Ideal + stable", "Ideal (MGIDI)", "Stable (MTSI)",
                                 "Not selected")) +
  geom_text_repel(data = filter(cmp, Selection != "Neither"),
                  aes(label = Genotype, colour = Selection),
                  size = 3, fontface = "bold", box.padding = 0.55, min.segment.length = 0,
                  max.overlaps = Inf, segment.colour = "grey70", show.legend = FALSE, seed = 42) +
  scale_x_reverse(breaks = c(1, seq(10, 60, 10))) +
  scale_y_reverse(breaks = c(1, seq(10, 60, 10))) +
  labs(x = "MGIDI rank (1 = closest to the ideotype)",
       y = "MTSI rank (1 = most stable)", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text = element_text(colour = "black"), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size = 3.5)))

## ---- panel B: min-max scale each trait to 0-1 across all 63, oriented to desired direction ----
sc <- blup
for (t in trait_order) {
  v <- blup[[t]]; r <- (v - min(v)) / (max(v) - min(v))
  if (t == desired_low) r <- 1 - r          # reverse DM so outer = earlier = favourable
  sc[[t]] <- r
}
vals <- sc %>% filter(Genotype %in% dual) %>%
  pivot_longer(all_of(trait_order), names_to = "Trait", values_to = "v") %>%
  mutate(Trait = factor(Trait, levels = trait_order),
         Genotype = factor(Genotype, levels = dual))

## radar geometry (6 axes from top, clockwise)
n <- length(trait_order); ang <- pi/2 - 2*pi*(seq_len(n) - 1)/n; names(ang) <- trait_order
poly <- vals %>% mutate(a = ang[as.character(Trait)], x = v*cos(a), y = v*sin(a))
poly_closed <- bind_rows(poly, poly %>% group_by(Genotype) %>% slice(1) %>% ungroup()) %>%
  arrange(Genotype)

gridhex <- do.call(rbind, lapply(seq(0.2, 1, 0.2), function(L) {
  d <- data.frame(level = L, x = L*cos(ang[trait_order]), y = L*sin(ang[trait_order])); rbind(d, d[1, ])
}))
spokes  <- data.frame(x = cos(ang[trait_order]), y = sin(ang[trait_order]))
dir_lab <- c(PH="PH ↑", FPH="FPH ↑", NB="NB ↑",
             PNMS="PNMS ↑", PNLB="PNLB ↑", DM="DM ↓")
axlab   <- data.frame(lx = 1.15*cos(ang[trait_order]), ly = 1.15*sin(ang[trait_order]),
                      lab = dir_lab[trait_order])
ticklab <- data.frame(x = 0.02, y = seq(0.2, 1, 0.2), L = seq(0.2, 1, 0.2))

## panel B palette — distinct from panel A's category colours (genotype identity, not category)
pal4 <- c("Gaterslebener stamm" = "#0072B2", "VNIIS 1" = "#E69F00",
          "ISZ 16" = "#CC79A7", "No. 700" = "#009E73")

pB <- ggplot() +
  geom_polygon(data = gridhex, aes(x, y, group = level), fill = NA, colour = "grey85", linewidth = 0.4) +
  geom_segment(data = spokes, aes(x = 0, y = 0, xend = x, yend = y), colour = "grey85", linewidth = 0.4) +
  geom_polygon(data = poly_closed, aes(x, y, fill = Genotype, colour = Genotype),
               alpha = 0.08, linewidth = 1.1) +
  geom_point(data = poly, aes(x, y, colour = Genotype), size = 2.3) +
  geom_text(data = axlab, aes(lx, ly, label = lab), fontface = "bold", size = 3.5) +
  geom_text(data = ticklab, aes(x, y, label = sprintf("%.1f", L)), size = 2.5, colour = "grey45", hjust = 0) +
  scale_fill_manual(values = pal4) +
  scale_colour_manual(values = pal4) +
  coord_equal(clip = "off") + labs(fill = NULL, colour = NULL) +
  theme_void(base_size = 12) +
  theme(legend.position = "bottom", plot.margin = margin(6, 16, 6, 16),
        legend.text = element_text(size = 9))

## ---- compose ----
combined <- pA + pB + plot_layout(widths = c(1.05, 1)) +
  plot_annotation(tag_levels = "A")

######### SAVE OUTPUT #########################################################################
## --- figures (600-dpi TIFF lzw) ---
ggsave("output/figures/Figure_MGIDI_ranking.tiff", mgidi_panel,
       width = 16, height = 8.5, units = "in", dpi = 600, device = "tiff", compression = "lzw")
ggsave("output/figures/Figure_MGIDI_MTSI_ideotype.tiff", combined,
       width = 13.5, height = 7.5, units = "in", dpi = 600, device = "tiff", compression = "lzw")
ggsave("output/figures/Supplementary_MGIDI_contribution_bars.tiff", mgidi_contribution_bars,
       width = 9, height = 6, units = "in", dpi = 600, device = "tiff", compression = "lzw")

## --- tables (csv) ---
write.csv(mgidi_fa_table,      "output/tables/Table_MGIDI_factor_analysis.csv",              row.names = FALSE)
write.csv(mgidi_seldiff_table, "output/tables/Table_MGIDI_selection_differential.csv",       row.names = FALSE)
write.csv(weight_sweep_table,  "output/tables/Supplementary_MGIDI_weight_sweep.csv",         row.names = FALSE)
write.csv(mgidi_full_ranking,  "output/tables/Supplementary_MGIDI_full_ranking.csv",         row.names = FALSE)
write.csv(mtsi_fa_table,       "output/tables/Supplementary_MTSI_factor_analysis.csv",       row.names = FALSE)
write.csv(mtsi_seldiff_table,  "output/tables/Supplementary_MTSI_selection_differential.csv", row.names = FALSE)
write.csv(mtsi_full_ranking,   "output/tables/Supplementary_MTSI_full_ranking.csv",          row.names = FALSE)
write.csv(comparison_data,     "output/tables/Table_MGIDI_MTSI_comparison.csv",              row.names = FALSE)

cat("\n=== MGIDI and MTSI: all figures and tables saved to output/ ===\n")
