###############################################################################################
###               Multi-trait selection defines a plant-architecture ideotype               ###
###                           for soybean in the Russian Far East                           ###
###                                                                                         ###
###                                       ----------                                        ###
###             Script 04: Genetic parameters, variance components and BLUPs                ###
###                                                                                         ###
###                            R script by Yawar Habib, Skoltech                            ###
###           Center for Bio and Medical Technologies, ORCID 0009-0008-1812-3910            ###
###############################################################################################

######### PURPOSE #############################################################################
## Partition the phenotypic variance of each architecture trait into its genetic and
## non-genetic components, estimate broad-sense heritability on a genotype-mean basis, the
## genotypic and phenotypic coefficients of variation, the expected genetic advance and the
## selective accuracy, and predict the genotypic values (BLUPs) of the 63 accessions.
## The genotypic predictions are then cross-validated by leaving out one year at a time.
##
## Model: trait ~ (1|Genotype) + (1|Year) + (1|Genotype:Year) + (1|Year:Replication)
## Equation (1) of the manuscript with every effect random. The replication-within-year
## variance is estimated at or near zero for most traits and the fit is singular, so the
## reported parameters come from the reduced model without that term. Both fits are shown
## below and the two agree to the precision reported.

######### DATA ################################################################################
## data/Architecture_63genotypes_2018_2020_clean.csv  (see script 01)

######### OUTPUT ##############################################################################
## Table 3        output/tables/Table_genetic_parameters_architecture.csv
## Figure 3       output/figures/Figure_variance_partition_and_CV.tiff
## Figure 4       output/figures/Figure_BLUP_ranking.tiff
## Supplementary  output/tables/Supplementary_variance_components_and_cross_validation.csv
## Supplementary  output/tables/Supplementary_BLUP_full_ranking.csv
## Supplementary  output/figures/Supplementary_Figure_H2_GAM.tiff

######### PREPARATION OF THE WORKING INTERFACE IN R ###########################################
### I. Set working directory
# Set the working directory to the root of this repository, i.e. the folder that contains
# the data/, Scripts/ and output/ folders.
# On RStudio: tab 'Session' -> Set Working Directory -> Choose Directory.

### II. R packages needed for the analysis
# dplyr, tidyr, lme4, lmerTest, ggplot2, ggrepel, patchwork

### III. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)   # ranova() for LRT on random effects
library(ggplot2)
library(ggrepel)
library(patchwork)

dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

## ---- genotype label formatting (data-file names -> manuscript names) ----
## The forms below are those of Supplementary Table S3, the accession list, so that the
## deposited BLUP table names each accession as it is catalogued there. Applied to the raw
## data before the models are fitted; relabelling does not affect any estimate.
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

######### CANONICAL TRAIT DEFINITIONS #########################################################
architecture_traits <- c("Height", "First_pod_height", "No_branches",
                         "No_ProdNodes_MainStem", "No_ProdNodes_LateralBranches",
                         "Vegetation_period")

trait_labels <- c(
  Height = "Plant height", First_pod_height = "First pod height",
  No_branches = "Number of branches", No_ProdNodes_MainStem = "Productive nodes on main stem",
  No_ProdNodes_LateralBranches = "Productive nodes on lateral branches",
  Vegetation_period = "Days to maturity"
)
trait_abbreviations <- c(
  Height = "PH", First_pod_height = "FPH", No_branches = "NB",
  No_ProdNodes_MainStem = "PNMS", No_ProdNodes_LateralBranches = "PNLB",
  Vegetation_period = "DM"
)

######### MIXED MODELS AND DESIGN CONSTANTS ###################################################
## helper: fit the (reduced) mixed model for one trait
fit_trait_model <- function(trait_name) {
  lmer(as.formula(paste0(trait_name, " ~ (1|Genotype) + (1|Year) + (1|Genotype:Year)")),
       data = Architecture_data)
}

## helper: fit the full model of Equation (1), replication within year included
fit_trait_model_full <- function(trait_name) {
  lmer(as.formula(paste0(trait_name,
       " ~ (1|Genotype) + (1|Year) + (1|Genotype:Year) + (1|Year:Replication)")),
       data = Architecture_data)
}

## the six models are fitted once here and reused throughout the script
fitted_models <- lapply(architecture_traits, fit_trait_model)
names(fitted_models) <- architecture_traits

## the full model of Equation (1), shown for comparison. The replication-within-year
## variance is zero or near zero for every trait and lme4 reports a singular fit, which is
## why the reduced model above is used for the reported parameters.
replication_variance_check <- lapply(architecture_traits, function(trait_name) {
  vc_full <- as.data.frame(VarCorr(fit_trait_model_full(trait_name)))
  vc_red  <- as.data.frame(VarCorr(fitted_models[[trait_name]]))
  data.frame(
    Trait      = trait_abbreviations[trait_name],
    Vrep_full  = round(vc_full$vcov[vc_full$grp == "Year:Replication"], 4),
    Vg_full    = round(vc_full$vcov[vc_full$grp == "Genotype"], 3),
    Vg_reduced = round(vc_red$vcov[vc_red$grp == "Genotype"], 3),
    row.names  = NULL
  )
}) %>% bind_rows()
print(replication_variance_check)

## 1. Design constants (years, harmonic-mean reps)
number_of_years <- nlevels(Architecture_data$Year)            # 3

reps_per_cell <- Architecture_data %>% count(Genotype, Year, name = "reps") %>% pull(reps)
number_of_replicates <- length(reps_per_cell) / sum(1 / reps_per_cell)   # harmonic mean ~ 2.503

selection_intensity_k <- 2.06                                  # k at 5% selection intensity

######### GENETIC PARAMETERS ##################################################################
## 2. Genetic parameters table (H2, GCV, PCV, GA, GAM, accuracy)
genetic_parameters <- lapply(architecture_traits, function(trait_name) {
  m  <- fitted_models[[trait_name]]
  vc <- as.data.frame(VarCorr(m))
  Vg  <- vc$vcov[vc$grp == "Genotype"]
  Vgy <- vc$vcov[vc$grp == "Genotype:Year"]
  Vy  <- vc$vcov[vc$grp == "Year"]
  Ve  <- vc$vcov[vc$grp == "Residual"]
  grand_mean <- mean(Architecture_data[[trait_name]], na.rm = TRUE)

  Vp  <- Vg + Vgy / number_of_years + Ve / (number_of_years * number_of_replicates)
  H2  <- Vg / Vp
  GCV <- sqrt(Vg) / grand_mean * 100
  PCV <- sqrt(Vp) / grand_mean * 100
  GA  <- selection_intensity_k * sqrt(Vp) * H2
  GAM <- GA / grand_mean * 100

  data.frame(
    Trait = trait_labels[trait_name], Mean = round(grand_mean, 2),
    Vg = round(Vg, 2), Vgy = round(Vgy, 2), Ve = round(Ve, 2),
    H2 = round(H2, 2), GCV = round(GCV, 2), PCV = round(PCV, 2),
    GCV_PCV = round(GCV / PCV, 2), GA = round(GA, 2), GAM = round(GAM, 2),
    Accuracy = round(sqrt(H2), 2)
  )
}) %>% bind_rows()
rownames(genetic_parameters) <- NULL

genetic_parameters_table <- genetic_parameters %>%
  mutate(Trait = trait_abbreviations[match(Trait, trait_labels)]) %>%
  transmute(Trait, Mean, H2, `GCV (%)` = GCV, `PCV (%)` = PCV,
            `GCV/PCV` = GCV_PCV, GA, `GAM (%)` = GAM, Accuracy)
print(genetic_parameters_table)

######### VARIANCE COMPONENTS AND LIKELIHOOD RATIO TESTS ######################################
## 3. Variance components + LRT significance (genotype term)
variance_components_table <- lapply(architecture_traits, function(trait_name) {
  m  <- fitted_models[[trait_name]]
  vc <- as.data.frame(VarCorr(m))
  p_geno <- as.data.frame(ranova(m))["(1 | Genotype)", "Pr(>Chisq)"]
  star <- ifelse(p_geno < 0.001, "***", ifelse(p_geno < 0.01, "**", ifelse(p_geno < 0.05, "*", "")))
  data.frame(
    Trait = trait_abbreviations[trait_name],
    Vg = paste0(round(vc$vcov[vc$grp == "Genotype"], 2), star),
    Vgy = round(vc$vcov[vc$grp == "Genotype:Year"], 2),
    Ve  = round(vc$vcov[vc$grp == "Residual"], 2)
  )
}) %>% bind_rows()
print(variance_components_table)

######### BLUPS OF THE GENOTYPIC EFFECTS ######################################################
## 5. BLUPs: extract, rank, flag extremes, figure
blup_data <- lapply(architecture_traits, function(trait_name) {
  re <- ranef(fitted_models[[trait_name]])$Genotype
  data.frame(Trait = trait_abbreviations[trait_name],
             Genotype = rownames(re), BLUP = re[, "(Intercept)"])
}) %>% bind_rows()

n_extreme <- ceiling(0.05 * nlevels(Architecture_data$Genotype))   # 4 of 63

# full trait names with units (heights in cm, maturity in days; branch/node traits are counts)
trait_full <- c("Plant height (cm)", "First pod height (cm)", "Number of branches",
                "Productive nodes on main stem", "Productive nodes on lateral branches",
                "Days to maturity (days)")
names(trait_full) <- trait_abbreviations[architecture_traits]

## the deposited BLUP table keeps the trait abbreviations as column names, so that
## script 06 can read it directly
blup_full_table <- blup_data %>%
  mutate(Trait = factor(Trait, levels = trait_abbreviations[architecture_traits])) %>%
  select(Genotype, Trait, BLUP) %>% mutate(BLUP = round(BLUP, 2)) %>%
  pivot_wider(names_from = Trait, values_from = BLUP) %>% arrange(Genotype)

blup_data <- blup_data %>%
  group_by(Trait) %>% arrange(BLUP, .by_group = TRUE) %>%
  mutate(rank_low = row_number(), rank_high = row_number(desc(BLUP)),
         Extreme = case_when(rank_high <= n_extreme ~ "Top 5%",
                             rank_low  <= n_extreme ~ "Bottom 5%", TRUE ~ "Middle"),
         order_index = row_number()) %>% ungroup() %>%
  mutate(Trait = factor(trait_full[as.character(Trait)], levels = trait_full),
         Extreme = factor(Extreme, levels = c("Top 5%", "Middle", "Bottom 5%")))

######### LEAVE-ONE-YEAR-OUT CROSS-VALIDATION #################################################
## Each year in turn is held out. The model is refitted on the remaining two years and the
## genotypic BLUPs from that reduced fit are correlated with the genotype means observed in
## the held-out year. Predictive ability is the mean of the three fold correlations.
##
## A correlation between a predicted genotypic value and an observed phenotypic mean is
## attenuated by the noise in that single-year mean, so it is divided by the square root of
## the single-year heritability to give predictive accuracy, which is on the same scale as
## the selective accuracy reported in Table 3:
##
##      h2 (single year) = Vg / (Vg + Vgy + Ve / r)
##      predictive accuracy = predictive ability / sqrt(h2 single year)
##
## The correlation between the fold BLUPs and the full-data BLUPs is also reported, as a
## measure of how much the ranking depends on any one year.
years_for_cv <- levels(Architecture_data$Year)

cross_validation <- lapply(architecture_traits, function(trait_name) {
  vc <- as.data.frame(VarCorr(fitted_models[[trait_name]]))
  Vg  <- vc$vcov[vc$grp == "Genotype"]
  Vgy <- vc$vcov[vc$grp == "Genotype:Year"]
  Ve  <- vc$vcov[vc$grp == "Residual"]
  h2_single_year <- Vg / (Vg + Vgy + Ve / number_of_replicates)

  blup_full <- ranef(fitted_models[[trait_name]])$Genotype[, "(Intercept)"]
  names(blup_full) <- rownames(ranef(fitted_models[[trait_name]])$Genotype)

  folds <- sapply(years_for_cv, function(held_out_year) {
    training_data <- droplevels(Architecture_data[Architecture_data$Year != held_out_year, ])
    m_train <- lmer(as.formula(paste0(trait_name, " ~ (1|Genotype) + (1|Year) + (1|Genotype:Year)")),
                    data = training_data)
    blup_train <- ranef(m_train)$Genotype[, "(Intercept)"]
    names(blup_train) <- rownames(ranef(m_train)$Genotype)

    observed_mean <- Architecture_data[Architecture_data$Year == held_out_year, ] %>%
      group_by(Genotype) %>%
      summarise(observed = mean(.data[[trait_name]], na.rm = TRUE), .groups = "drop")
    g <- as.character(observed_mean$Genotype)

    c(ability   = cor(blup_train[g], observed_mean$observed),
      vs_full   = cor(blup_train[g], blup_full[g]))
  })

  predictive_ability <- mean(folds["ability", ])

  data.frame(
    Trait                 = trait_abbreviations[trait_name],
    `r 2018`              = round(folds["ability", "2018"], 3),
    `r 2019`              = round(folds["ability", "2019"], 3),
    `r 2020`              = round(folds["ability", "2020"], 3),
    `Predictive ability`  = round(predictive_ability, 3),
    `h2 single year`      = round(h2_single_year, 3),
    `Predictive accuracy` = round(predictive_ability / sqrt(h2_single_year), 3),
    `BLUP vs full min`    = round(min(folds["vs_full", ]), 3),
    `BLUP vs full max`    = round(max(folds["vs_full", ]), 3),
    check.names = FALSE
  )
}) %>% bind_rows()
rownames(cross_validation) <- NULL
print(cross_validation)

## The variance components and the cross-validation are reported as a single supplementary
## table, keyed by trait. The single-year heritability and the fold-versus-full BLUP
## correlations are printed above as diagnostics but are not part of that table.
variance_and_cross_validation <- variance_components_table %>%
  left_join(cross_validation %>%
              select(Trait, `r 2018`, `r 2019`, `r 2020`,
                     `Predictive ability`, `Predictive accuracy`),
            by = "Trait")
print(variance_and_cross_validation)

######### FIGURE 4 - BLUP RANKING #############################################################
blup_figure <- ggplot(blup_data, aes(order_index, BLUP, colour = Extreme)) +
  geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey70") +
  geom_segment(aes(xend = order_index, yend = 0), linewidth = 0.4) +   # lollipop stems
  geom_point(size = 2.0) +
  facet_wrap(~ Trait, scales = "free", ncol = 3) +
  scale_colour_manual(values = c("Top 5%" = "#2C6E91", "Middle" = "grey75", "Bottom 5%" = "#C0392B")) +
  labs(x = "Genotypes (ranked by BLUP)", y = "BLUP (genotypic deviation from population mean)", colour = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = "bold"))

######### SUPPLEMENTARY FIGURE - HERITABILITY AGAINST GENETIC ADVANCE #########################
## 6. H2 vs GAM figure  (builds h2_gam_data from genetic_parameters)
trait_groups <- c(PH = "Main-axis / structural", FPH = "Main-axis / structural",
                  PNMS = "Main-axis / structural", NB = "Branching", PNLB = "Branching",
                  DM = "Phenology")

h2_gam_data <- genetic_parameters %>%
  mutate(Trait = trait_abbreviations[match(Trait, trait_labels)],
         Category = trait_groups[Trait])

h2_gam_figure <- ggplot(h2_gam_data, aes(H2, GAM, colour = Category)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(aes(label = Trait), size = 4.6, fontface = "bold",
                  box.padding = 0.8, point.padding = 0.5, segment.colour = "grey70",
                  min.segment.length = 0, show.legend = FALSE) +
  scale_colour_manual(values = c("Branching" = "#C1466B", "Main-axis / structural" = "#1F6F8B",
                                 "Phenology" = "#E09F3E"), name = "Trait group") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(x = expression("Broad-sense heritability ("*italic(H)^2*")"),
       y = "Genetic advance as % of mean (GAM)") +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(colour = "black"), legend.position = "right",
        legend.title = element_text(face = "bold"))

######### FIGURE 3 - VARIANCE PARTITION AND COEFFICIENTS OF VARIATION #########################
## Panel A: variance partition (% of total). Panel B: GCV against PCV.
## All values are derived live from the models fitted above, nothing is hardcoded and no
## model is refitted.
## ---- one pass over the fitted models: variance components + coefficients of variation ----
gen_params <- lapply(architecture_traits, function(trait_name) {
  m  <- fitted_models[[trait_name]]
  vc <- as.data.frame(VarCorr(m))

  vg  <- vc$vcov[vc$grp == "Genotype"]
  vge <- vc$vcov[vc$grp == "Genotype:Year"]
  vy  <- vc$vcov[vc$grp == "Year"]
  ve  <- vc$vcov[vc$grp == "Residual"]

  xbar <- mean(Architecture_data[[trait_name]], na.rm = TRUE)   # trait mean

  ## phenotypic variance, genotype-mean basis
  vp_mean <- vg + vge / number_of_years + ve / (number_of_years * number_of_replicates)

  data.frame(
    Trait = trait_abbreviations[trait_name],
    Genotype = vg, `Genotype x Year` = vge, Year = vy, Residual = ve,
    GCV = sqrt(vg)      / xbar * 100,
    PCV = sqrt(vp_mean) / xbar * 100,
    check.names = FALSE
  )
}) %>% bind_rows()

## ---- panel A: variance partition (% of total) ----------------------------
## largest-remainder rounding for the printed segment labels. Rounding each share on its own
## does not always give 100: days to maturity is 67.45 / 12.91 / 4.23 / 15.41, which rounds to
## 67 + 13 + 4 + 15 = 99. Flooring each share and giving the remaining whole units to the
## largest fractional parts keeps every label within one point of its exact value and makes
## the four labels sum to 100 for every trait.
largest_remainder <- function(p) {
  base  <- floor(p)
  short <- round(sum(p)) - sum(base)
  if (short > 0) {
    take <- order(p - base, decreasing = TRUE)[seq_len(short)]
    base[take] <- base[take] + 1
  }
  base
}

long <- gen_params %>%
  select(Trait, Genotype, `Genotype x Year`, Year, Residual) %>%
  pivot_longer(-Trait, names_to = "Source", values_to = "Variance") %>%
  group_by(Trait) %>%
  mutate(Percent      = Variance / sum(Variance) * 100,
         PercentLabel = largest_remainder(Percent)) %>%
  ungroup()

g_order <- long %>% filter(Source == "Genotype") %>% arrange(desc(Percent)) %>% pull(Trait)

long <- long %>%
  mutate(Trait  = factor(Trait, levels = g_order),
         Source = factor(Source, levels = c("Residual", "Year", "Genotype x Year", "Genotype")))

## palette B
pal <- c("Genotype" = "#59A14F", "Genotype x Year" = "#8CB3D9",
         "Year" = "#C98BAA", "Residual" = "#D4D4CE")   # Year -> dusty rose

seg_lab <- long %>%
  mutate(lab = ifelse(Percent >= 4, sprintf("%d%%", PercentLabel), ""),
         tc  = ifelse(Source %in% c("Year", "Residual"), "grey25", "white"))  # dark text on pale fills

theme_fig <- theme_classic(base_size = 12) +
  theme(axis.text = element_text(colour = "black"),
        axis.title.y = element_text(margin = margin(r = 6)),
        legend.position = "bottom", legend.key.size = unit(0.9, "lines"),
        plot.margin = margin(8, 10, 6, 8))

pA <- ggplot(long, aes(Trait, Percent, fill = Source)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.6) +
  geom_text(data = seg_lab, aes(label = lab, colour = tc),
            position = position_stack(vjust = 0.5),
            size = 3.4, fontface = "bold", show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_manual(values = pal, breaks = c("Genotype", "Genotype x Year", "Year", "Residual")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Proportion of total variance (%)", fill = NULL) +
  theme_fig + guides(fill = guide_legend(nrow = 2, byrow = TRUE))

## ---- panel B: GCV vs PCV (derived above, not hardcoded) ------------------
cv <- gen_params %>%
  select(Trait, GCV, PCV) %>%
  mutate(Trait = factor(Trait, levels = g_order)) %>%
  pivot_longer(-Trait, names_to = "CV", values_to = "value") %>%
  mutate(CV = factor(CV, levels = c("GCV", "PCV")))

pB <- ggplot(cv, aes(Trait, value, fill = CV)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64,
           colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", value)),
            position = position_dodge(width = 0.72), vjust = -0.45,
            size = 2.9, fontface = "bold", colour = "grey25") +
  scale_fill_manual(values = c("GCV" = "#59A14F", "PCV" = "#C98BAA"),
                    breaks = c("GCV", "PCV"),
                    labels = c("Genotypic (GCV)", "Phenotypic (PCV)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14)), breaks = seq(0, 50, 10)) +
  labs(x = NULL, y = "Coefficient of variation (%)", fill = NULL) +
  theme_fig

## ---- compose ------------------------------------------------------------
combined <- pA + pB + plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A")

######### SAVE OUTPUT #########################################################################
ggsave("output/figures/Figure_variance_partition_and_CV.tiff", combined,
       width = 12, height = 6, units = "in", dpi = 600, device = "tiff", compression = "lzw")
ggsave("output/figures/Figure_BLUP_ranking.tiff", blup_figure,
       width = 11, height = 7, units = "in", dpi = 600, device = "tiff", compression = "lzw")
ggsave("output/figures/Supplementary_Figure_H2_GAM.tiff", h2_gam_figure,
       width = 8, height = 6, units = "in", dpi = 600, device = "tiff", compression = "lzw")

write.csv(genetic_parameters_table,  "output/tables/Table_genetic_parameters_architecture.csv", row.names = FALSE)
write.csv(variance_and_cross_validation,
          "output/tables/Supplementary_variance_components_and_cross_validation.csv", row.names = FALSE)
write.csv(blup_full_table,           "output/tables/Supplementary_BLUP_full_ranking.csv",       row.names = FALSE)

cat("\n=== Genetic parameters: all figures and tables saved to output/ ===\n")
