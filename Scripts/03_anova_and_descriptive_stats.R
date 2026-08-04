###############################################################################################
###               Multi-trait selection defines a plant-architecture ideotype               ###
###                           for soybean in the Russian Far East                           ###
###                                                                                         ###
###                                       ----------                                        ###
###                Script 03: Descriptive statistics and analysis of variance               ###
###                                                                                         ###
###                            R script by Yawar Habib, Skoltech                            ###
###           Center for Bio and Medical Technologies, ORCID 0009-0008-1812-3910            ###
###############################################################################################

######### PURPOSE #############################################################################
## Descriptive statistics for the six architecture traits, the combined analysis of variance
## across the three years, the model assumption diagnostics, and the year-to-year comparison
## of trait distributions.
##
## Model: trait ~ Year + Year:Replication + Genotype + Year:Genotype
## Type III sums of squares under the sum-to-zero constraint, because the replication is
## unbalanced (two to three replicates per genotype-year cell) and the genotype x year
## interaction is significant for every trait.

######### DATA ################################################################################
## data/Architecture_63genotypes_2018_2020_clean.csv  (see script 01)

######### OUTPUT ##############################################################################
## Table 1        output/tables/Table_descriptive_statistics_architecture.csv
## Table 2        output/tables/Table_ANOVA_architecture.csv
## Figure 2       output/figures/Figure_trait_distribution_across_years.tiff
## Supplementary  output/figures/Supplementary_QQ_plots.tiff
## Supplementary  output/figures/Supplementary_Residual_vs_Fitted.tiff

######### PREPARATION OF THE WORKING INTERFACE IN R ###########################################
### I. Set working directory
# Set the working directory to the root of this repository, i.e. the folder that contains
# the data/, Scripts/ and output/ folders.
# On RStudio: tab 'Session' -> Set Working Directory -> Choose Directory.

### II. R packages needed for the analysis
# dplyr, tidyr, psych, car, emmeans, ggplot2, ggpubr

### III. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(psych)
library(car)
library(emmeans)
library(ggplot2)
library(ggpubr)

dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

Architecture_data <- read.csv("data/Architecture_63genotypes_2018_2020_clean.csv", header = TRUE) %>%
  mutate(
    Year        = as.factor(Year),
    Replication = as.factor(Replication),
    Genotype    = as.factor(Genotype)
  )

######### CANONICAL TRAIT DEFINITIONS #########################################################
## six traits, DF removed
architecture_traits <- c("Height", "First_pod_height", "No_branches",
                         "No_ProdNodes_MainStem", "No_ProdNodes_LateralBranches",
                         "Vegetation_period")

trait_labels <- c(
  Height                       = "Plant height",
  First_pod_height             = "First pod height",
  No_branches                  = "Number of branches",
  No_ProdNodes_MainStem        = "Productive nodes on main stem",
  No_ProdNodes_LateralBranches = "Productive nodes on lateral branches",
  Vegetation_period            = "Days to maturity"
)

trait_abbreviations <- c(
  Height                       = "PH",
  First_pod_height             = "FPH",
  No_branches                  = "NB",
  No_ProdNodes_MainStem        = "PNMS",
  No_ProdNodes_LateralBranches = "PNLB",
  Vegetation_period            = "DM"
)

str(Architecture_data)

######### DESCRIPTIVE STATISTICS ##############################################################
descriptive_table_for_word <- describe(Architecture_data[, architecture_traits]) %>%
  as.data.frame() %>%
  mutate(
    Trait    = trait_labels[architecture_traits],
    `CV (%)` = round(sd / mean * 100, 2)
  ) %>%
  transmute(
    Trait,
    Mean      = round(mean, 2),
    SD        = round(sd, 2),
    Min       = round(min, 2),
    Max       = round(max, 2),
    Skewness  = round(skew, 2),
    Kurtosis  = round(kurtosis, 2),
    `CV (%)`
  )

descriptive_table_for_word

write.csv(descriptive_table_for_word,
          "output/tables/Table_descriptive_statistics_architecture.csv", row.names = FALSE)

######### COMBINED ANALYSIS OF VARIANCE #######################################################
## Type III sums of squares are only meaningful under sum-to-zero contrasts, so
## the contrast option is set before the models are fitted and restored after.
previous_contrasts <- options(contrasts = c("contr.sum", "contr.poly"))

## store the fitted lm models (needed for both ANOVA tables AND residual diagnostics)
fitted_models <- lapply(architecture_traits, function(trait_name) {
  model_formula <- as.formula(
    paste0(trait_name, " ~ Year + Year:Replication + Genotype + Year:Genotype")
  )
  lm(model_formula, data = Architecture_data)
})
names(fitted_models) <- architecture_traits

combined_anova_results <- lapply(fitted_models, function(m) Anova(m, type = "III"))

options(previous_contrasts)

combined_anova_results[["Height"]]

######### ANOVA SUMMARY TABLE #################################################################
add_stars <- function(p) ifelse(p < 0.001, "***",
                                ifelse(p < 0.01,  "**",
                                       ifelse(p < 0.05,  "*", "")))

source_terms <- c("Year", "Year:Replication", "Genotype", "Year:Genotype")

anova_table_for_word <- lapply(architecture_traits, function(trait_name) {
  tab <- combined_anova_results[[trait_name]]
  ms  <- tab[source_terms, "Sum Sq"] / tab[source_terms, "Df"]
  p   <- tab[source_terms, "Pr(>F)"]
  setNames(paste0(round(ms, 2), add_stars(p)), source_terms)
}) %>%
  bind_rows() %>%
  mutate(Trait = trait_labels[architecture_traits], .before = 1) %>%
  rename(`Replication (Year)` = `Year:Replication`,
         `Genotype x Year`    = `Year:Genotype`)

anova_table_for_word

write.csv(anova_table_for_word,
          "output/tables/Table_ANOVA_architecture.csv", row.names = FALSE)

######### MODEL ASSUMPTION DIAGNOSTICS ########################################################
## (a) Normality of residuals — Shapiro-Wilk per trait + residual skew/kurtosis
normality_check <- lapply(architecture_traits, function(trait_name) {
  res <- residuals(fitted_models[[trait_name]])
  sw  <- shapiro.test(res)
  data.frame(
    Trait            = trait_labels[trait_name],
    Shapiro_W        = round(sw$statistic, 3),
    Shapiro_p        = signif(sw$p.value, 3),
    Residual_skew    = round(psych::skew(res), 2),
    Residual_kurtosis= round(psych::kurtosi(res), 2)
  )
}) %>% bind_rows()
rownames(normality_check) <- NULL
normality_check

## (b) Homogeneity of variance — Levene's test across years (per trait)
homogeneity_check <- lapply(architecture_traits, function(trait_name) {
  lev <- leveneTest(Architecture_data[[trait_name]] ~ Architecture_data$Year)
  data.frame(
    Trait      = trait_labels[trait_name],
    Levene_F   = round(lev$`F value`[1], 3),
    Levene_p   = signif(lev$`Pr(>F)`[1], 3)
  )
}) %>% bind_rows()
rownames(homogeneity_check) <- NULL
homogeneity_check

## extract residuals + fitted values from the six stored models into one long frame
diagnostic_data <- lapply(architecture_traits, function(trait_name) {
  m <- fitted_models[[trait_name]]
  data.frame(
    Trait    = trait_abbreviations[trait_name],
    Residual = residuals(m),
    Fitted   = fitted(m)
  )
}) %>%
  bind_rows() %>%
  mutate(Trait = factor(Trait, levels = trait_abbreviations[architecture_traits]))

## ---- Supplementary Figure: Q-Q plots (normality of residuals) ----
qq_plot_figure <- ggplot(diagnostic_data, aes(sample = Residual)) +
  stat_qq(size = 0.9, colour = "#2C6E91") +
  stat_qq_line(colour = "#C0392B", linewidth = 0.7) +
  facet_wrap(~ Trait, scales = "free", ncol = 3) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles (residuals)") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

qq_plot_figure

## ---- Supplementary Figure: Residual vs Fitted (homogeneity / pattern) ----
resid_fitted_figure <- ggplot(diagnostic_data, aes(x = Fitted, y = Residual)) +
  geom_point(size = 0.9, colour = "#2C6E91", alpha = 0.65) +
  geom_hline(yintercept = 0, colour = "#C0392B", linewidth = 0.7) +
  geom_smooth(method = "loess", se = FALSE, colour = "#E67E22", linewidth = 0.7) +
  facet_wrap(~ Trait, scales = "free", ncol = 3) +
  labs(x = "Fitted values", y = "Residuals") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

resid_fitted_figure

## ---- Save Supplementary Figure: Q-Q plots ----
ggsave(
  filename = "output/figures/Supplementary_QQ_plots.tiff",
  plot     = qq_plot_figure,
  width    = 9,
  height   = 6,
  units    = "in",
  dpi      = 600,
  device   = "tiff",
  compression = "lzw"
)

## ---- Save Supplementary Figure: Residual vs Fitted ----
ggsave(
  filename = "output/figures/Supplementary_Residual_vs_Fitted.tiff",
  plot     = resid_fitted_figure,
  width    = 9,
  height   = 6,
  units    = "in",
  dpi      = 600,
  device   = "tiff",
  compression = "lzw"
)

######### YEAR-TO-YEAR TRAIT DISTRIBUTIONS ####################################################
## long data with duplicated "Across" group, six traits, abbreviated labels
plot_data_combined <- Architecture_data %>%
  select(Year, Genotype, all_of(architecture_traits)) %>%
  pivot_longer(all_of(architecture_traits), names_to = "Trait", values_to = "Value") %>%
  mutate(Year_group = as.character(Year)) %>%
  bind_rows(
    Architecture_data %>%
      select(Year, Genotype, all_of(architecture_traits)) %>%
      pivot_longer(all_of(architecture_traits), names_to = "Trait", values_to = "Value") %>%
      mutate(Year_group = "Across")
  ) %>%
  mutate(
    Year_group = factor(Year_group, levels = c("Across", "2018", "2019", "2020")),
    Trait      = factor(Trait, levels = architecture_traits,
                        labels = trait_abbreviations[architecture_traits])
  )

## Tukey-adjusted pairwise YEAR comparisons (genotype-aware), six traits
year_comparisons <- lapply(architecture_traits, function(trait_name) {
  m   <- lm(as.formula(paste0(trait_name, " ~ Year + Genotype")), data = Architecture_data)
  emm <- emmeans(m, ~ Year)
  pairs(emm, adjust = "tukey") %>% as.data.frame() %>% mutate(Trait = trait_name)
}) %>% bind_rows()

## brackets: SIGNIFICANT ONLY (drop ns entirely)
significance_brackets <- year_comparisons %>%
  mutate(
    group1 = sub("Year(\\d+).*", "\\1", contrast),
    group2 = sub(".*Year(\\d+)$", "\\1", contrast),
    stars  = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                 labels = c("***", "**", "*", "ns")),
    Trait  = factor(Trait, levels = architecture_traits,
                    labels = trait_abbreviations[architecture_traits])
  ) %>%
  filter(stars != "ns")   ## <-- significant only, no ns, blank where non-significant

## y-positions, stacked per panel
significance_brackets <- significance_brackets %>%
  left_join(
    plot_data_combined %>% group_by(Trait) %>% summarise(ymax = max(Value), .groups = "drop"),
    by = "Trait"
  ) %>%
  group_by(Trait) %>%
  mutate(y_position = ymax * (1.05 + 0.07 * (row_number() - 1))) %>%
  ungroup()

cols <- c("Across" = "#7FBF7F",
          "2018"   = "#C9A3E9",
          "2019"   = "#F58A87",
          "2020"   = "#9092EC")

architecture_violin_figure <- ggplot(plot_data_combined,
                                     aes(x = Year_group, y = Value, fill = Year_group)) +
  geom_violin(aes(fill = Year_group), trim = FALSE, scale = "width",
              colour = "grey67", linewidth = 1.0, show.legend = FALSE) +
  geom_boxplot(width = 0.13, fill = NA, colour = "white",
               linewidth = 1.0, outlier.colour = "black", outlier.size = 1.2) +
  stat_summary(fun = median, geom = "crossbar", width = 0.12,
               colour = "grey45", linewidth = 0.3, fatten = 0) +
  facet_wrap(~ Trait, scales = "free_y", ncol = 3, strip.position = "left") +
  stat_pvalue_manual(significance_brackets,
                     label = "stars",
                     xmin = "group1", xmax = "group2",
                     y.position = "y_position",
                     tip.length = 0.01,
                     inherit.aes = FALSE) +
  scale_fill_manual(values = cols) +
  labs(x = NULL, y = NULL, fill = "Year") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

architecture_violin_figure

ggsave(
  filename = "output/figures/Figure_trait_distribution_across_years.tiff",
  plot     = architecture_violin_figure,
  width    = 12,
  height   = 7,
  units    = "in",
  dpi      = 600,
  device   = "tiff",
  compression = "lzw")

cat("\n=== Descriptive statistics and ANOVA: tables and figures saved to output/ ===\n")
