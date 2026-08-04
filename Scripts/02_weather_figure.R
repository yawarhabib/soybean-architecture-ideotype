###############################################################################################
###               Multi-trait selection defines a plant-architecture ideotype               ###
###                           for soybean in the Russian Far East                           ###
###                                                                                         ###
###                                       ----------                                        ###
###                Script 02: Weather during the three growing seasons                      ###
###                                                                                         ###
###                            R script by Yawar Habib, Skoltech                            ###
###           Center for Bio and Medical Technologies, ORCID 0009-0008-1812-3910            ###
###############################################################################################

######### PURPOSE #############################################################################
## Figure 1. Daily precipitation, the 7-day rolling mean air temperature and the rolling
## minimum-maximum range for the three growing seasons (2018-2020) at the FERIA field
## station, with the growing season and the sowing-to-emergence period marked below each
## panel.

######### DATA ################################################################################
## data/rp5_daily.csv
## Daily observations for the Khabarovsk station, obtained from the Reliable Prognosis
## weather archive (rp5.ru). Columns: date, n_obs (observations per day), t_mean_c, t_min_c,
## t_max_c, precip_mm. The record covers 1 May to 30 September of each year.
##
## data/Architecture_63genotypes_2018_2020_clean.csv
## Used only to derive the end of the growing-season strip from the observed days to
## maturity (see the GROWING SEASON STRIPS section below).

######### OUTPUT ##############################################################################
## Figure 1  output/figures/Figure_climate.tiff

######### PREPARATION OF THE WORKING INTERFACE IN R ###########################################
### I. Set working directory
# Set the working directory to the root of this repository, i.e. the folder that contains
# the data/, Scripts/ and output/ folders.
# On RStudio: tab 'Session' -> Set Working Directory -> Choose Directory.

### II. R packages needed for the analysis
# dplyr, lubridate, zoo, ggplot2

### III. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
gc()

library(dplyr)
library(lubridate)
library(zoo)
library(ggplot2)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

######### WEATHER RECORD AND ROLLING MEANS ####################################################
## plot_doy maps every date onto a common year so the three panels share one x axis
df <- read.csv("data/rp5_daily.csv", header = TRUE) %>%
  mutate(date = as.Date(date), year = year(date)) %>%
  filter(year %in% c(2018, 2019, 2020)) %>%
  mutate(plot_doy = as.Date(paste0("2000-", format(date, "%m-%d")))) %>%
  group_by(year) %>% arrange(date) %>%
  mutate(
    t_mean_smooth = rollmean(t_mean_c, 7, fill = NA, align = "center"),
    t_min_smooth  = rollmean(t_min_c,  7, fill = NA, align = "center"),
    t_max_smooth  = rollmean(t_max_c,  7, fill = NA, align = "center")
  ) %>% ungroup()

precip_max <- max(df$precip_mm, na.rm = TRUE)
temp_max   <- max(df$t_max_c,   na.rm = TRUE)
scale_fac  <- precip_max / temp_max

######### GROWING SEASON STRIPS ###############################################################
## Sowing and full-emergence dates are the recorded field operations for each season.
## The end of the growing season is not a field operation but a trait: days to maturity was
## recorded from sowing, so the strip ends at sowing + the mean days to maturity of the 63
## accessions in that year. This is derived from the trait data rather than assumed.
sowing_dates <- data.frame(
  year        = c(2018, 2019, 2020),
  sowing      = as.Date(c("2018-05-17", "2019-05-27", "2020-05-13")),
  full_shoots = as.Date(c("2018-05-30", "2019-06-10", "2020-05-28"))
)

mean_days_to_maturity <- read.csv("data/Architecture_63genotypes_2018_2020_clean.csv",
                                  header = TRUE) %>%
  group_by(year = Year) %>%
  summarise(mean_DM = mean(Vegetation_period), .groups = "drop")

mean_days_to_maturity

## map the three seasons onto the common plotting year used by plot_doy
to_plot_year <- function(x) as.Date(paste0("2000-", format(x, "%m-%d")))

strips <- sowing_dates %>%
  left_join(mean_days_to_maturity, by = "year") %>%
  mutate(grow_end = sowing + round(mean_DM)) %>%
  mutate(across(c(sowing, full_shoots, grow_end), to_plot_year))

strips

col_precip    <- "#4197d8"   # periwinkle blue (your 2020 color) — precip bars
col_temp_line <- "#E5645F"   # deepened coral (your 2019 #F58A87 darkened) — temp line
col_temp_fill <- "#F58A87"   # coral (your 2019 color, lightened via alpha) — ribbon
col_grow      <- "#CDE3CD"   # pale green (your "Across" #7FBF7F lightened) — growing season
col_sowing    <- "#7FBF7F"   # green (your "Across" color) — sowing

strip_y <- -precip_max * 0.07
strip_h <- precip_max * 0.035

######### FIGURE 1 - WEATHER DURING THE THREE GROWING SEASONS #################################
fig1 <- ggplot(df, aes(x = plot_doy)) +
  geom_col(aes(y = precip_mm), fill = col_precip, colour = NA, width = 0.9) +
  geom_ribbon(aes(ymin = t_min_smooth * scale_fac, ymax = t_max_smooth * scale_fac),
              fill = col_temp_fill, alpha = 0.35) +
  geom_line(aes(y = t_mean_smooth * scale_fac), colour = col_temp_line, linewidth = 0.8) +
  # growing-season strip (sowing -> mean maturity)
  geom_rect(data = strips,
            aes(xmin = sowing, xmax = grow_end,
                ymin = strip_y - strip_h, ymax = strip_y + strip_h, fill = "Growing season"),
            inherit.aes = FALSE) +
  # sowing segment (sowing -> full shoots)
  geom_rect(data = strips,
            aes(xmin = sowing, xmax = full_shoots,
                ymin = strip_y - strip_h, ymax = strip_y + strip_h, fill = "Sowing"),
            inherit.aes = FALSE) +
  scale_fill_manual(name = NULL,
                    values = c("Growing season" = col_grow, "Sowing" = col_sowing)) +
  facet_wrap(~ year, nrow = 1) +
  scale_y_continuous(
    name = "Precipitation (mm)",
    expand = expansion(mult = c(0.03, 0.05)),
    sec.axis = sec_axis(~ . / scale_fac, name = "Temperature (°C)",
                        breaks = scales::pretty_breaks(5))
  ) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month", expand = expansion(add = 1)) +
  labs(x = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
    axis.title.y.left  = element_text(colour = col_precip, face = "bold", margin = margin(r = 8)),
    axis.text.y.left   = element_text(colour = col_precip),
    axis.ticks.y.left  = element_line(colour = col_precip),
    axis.title.y.right = element_text(colour = col_temp_line, face = "bold", margin = margin(l = 8)),
    axis.text.y.right  = element_text(colour = col_temp_line),
    axis.ticks.y.right = element_line(colour = col_temp_line),
    axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 12),
    panel.spacing      = unit(1.2, "lines"),
    plot.margin        = margin(12, 12, 12, 12),
    legend.position    = "bottom"
  )

fig1

######### SAVE OUTPUT #########################################################################
ggsave(
  filename = "output/figures/Figure_climate.tiff",
  plot = fig1, width = 12, height = 7, units = "in",
  dpi = 600, device = "tiff", compression = "lzw"
)

cat("\n=== Figure 1 saved to output/figures/Figure_climate.tiff ===\n")
