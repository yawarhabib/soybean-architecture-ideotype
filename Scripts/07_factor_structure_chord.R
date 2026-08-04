###############################################################################################
###               Multi-trait selection defines a plant-architecture ideotype               ###
###                           for soybean in the Russian Far East                           ###
###                                                                                         ###
###                                       ----------                                        ###
###                    Script 07: Trait-to-factor chord diagram                             ###
###                                                                                         ###
###                            R script by Yawar Habib, Skoltech                            ###
###           Center for Bio and Medical Technologies, ORCID 0009-0008-1812-3910            ###
###############################################################################################

######### PURPOSE #############################################################################
## Figure 7. The trait-to-factor loading structure of the two factor analyses underlying the
## selection indices, drawn as chord diagrams. Panel A is the MGIDI factor analysis, which
## retained two factors, and panel B the MTSI factor analysis, which retained three. Chord
## width is the absolute loading and chord colour its sign, so the panels show which traits
## define each factor and in which direction.

######### DATA ################################################################################
## output/tables/Table_MGIDI_factor_analysis.csv          (Trait, FA1, FA2, ...)
## output/tables/Supplementary_MTSI_factor_analysis.csv    (Trait, FA1, FA2, FA3, ...)
##
## Both are written by script 06, which must therefore be run first.

######### OUTPUT ##############################################################################
## Figure 7  output/figures/Figure_factor_structure_chord.tiff

######### PREPARATION OF THE WORKING INTERFACE IN R ###########################################
### I. Set working directory
# Set the working directory to the root of this repository, i.e. the folder that contains
# the data/, Scripts/ and output/ folders.
# On RStudio: tab 'Session' -> Set Working Directory -> Choose Directory.

### II. R packages needed for the analysis
# circlize

### III. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
gc()

library(circlize)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

## Panel A <- MGIDI factor analysis (two factors), panel B <- MTSI factor analysis (three).
## Both tables are produced by script 06.
mgidi <- read.csv("output/tables/Table_MGIDI_factor_analysis.csv",        check.names = FALSE)
mtsi  <- read.csv("output/tables/Supplementary_MTSI_factor_analysis.csv", check.names = FALSE)

mgidi
mtsi

######### SHARED DEFINITIONS ##################################################################
trait_order <- c("PH", "FPH", "NB", "PNMS", "PNLB", "DM")

col_pos <- "#5A8FAE"   # positive loading (blue)
col_neg <- "#C98BAA"   # negative loading (rose)

trait_grp_col <- c(PH = "#4C72B0", FPH = "#4C72B0", PNMS = "#4C72B0",  # main-axis / height
                   NB = "#55A868", PNLB = "#55A868",                    # branching
                   DM = "#DD8452")                                      # phenology
factor_col <- "#8073AC"  # slate-purple for FA sectors

######### CHORD PANEL FUNCTION ################################################################
## ---- one chord panel from a loadings df, using the given factor columns ----
draw_chord <- function(fa, fcols) {
  fa <- fa[match(trait_order, fa$Trait), ]
  edges <- do.call(rbind, lapply(fcols, function(f)
    data.frame(Trait = fa$Trait, Factor = f, load = fa[[f]])))
  edges$abs  <- abs(edges$load)
  edges$sign <- ifelse(edges$load >= 0, "pos", "neg")

  mat      <- matrix(0, length(trait_order), length(fcols),
                     dimnames = list(trait_order, fcols))
  link_col <- matrix(NA, length(trait_order), length(fcols), dimnames = dimnames(mat))
  for (i in seq_len(nrow(edges))) {
    mat[edges$Trait[i], edges$Factor[i]]      <- edges$abs[i]              # width = |loading|
    link_col[edges$Trait[i], edges$Factor[i]] <-
      adjustcolor(ifelse(edges$sign[i] == "pos", col_pos, col_neg), 0.7)   # colour = sign
  }
  grid_col <- c(trait_grp_col[trait_order],
                setNames(rep(factor_col, length(fcols)), fcols))

  n_tr <- length(trait_order)
  circos.clear()
  circos.par(gap.after = c(rep(3, n_tr - 1), 12, rep(3, length(fcols) - 1), 12),
             start.degree = 90, canvas.xlim = c(-1.18, 1.18), canvas.ylim = c(-1.18, 1.18),
             points.overflow.warning = FALSE)
  chordDiagram(mat, grid.col = grid_col, col = link_col,
               annotationTrack = "grid", preAllocateTracks = list(track.height = 0.10),
               transparency = 0.3, link.sort = TRUE, link.decreasing = TRUE)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.9, CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.95, font = 2)
  }, bg.border = NA)
}

######### FIGURE 7 - TRAIT TO FACTOR LOADING STRUCTURE ########################################
## ---- compose: two panels + corner tags + shared legend ----
tiff("output/figures/Figure_factor_structure_chord.tiff", width = 3200, height = 1750, res = 300,
     compression = "lzw")   # use png("...png", width=3200, height=1750, res=300) for a quick-view PNG
par(mfrow = c(1, 2), mar = c(1.5, 0.5, 1.5, 0.5), oma = c(2, 0, 0, 0))

draw_chord(mgidi, c("FA1", "FA2"))
mtext("A", side = 3, adj = 0.02, line = -0.5, cex = 1.5, font = 1)

draw_chord(mtsi, c("FA1", "FA2", "FA3"))
mtext("B", side = 3, adj = 0.02, line = -0.5, cex = 1.5, font = 1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlim = c(0, 1), ylim = c(0, 1))
legend(x = 0.5, y = 0.055, xjust = 0.5,
       legend = c("Positive loading", "Negative loading",
                  "Main-axis / height", "Branching", "Phenology (DM)", "Factor"),
       fill = c(adjustcolor(col_pos, 0.7), adjustcolor(col_neg, 0.7),
                "#4C72B0", "#55A868", "#DD8452", factor_col),
       border = NA, bty = "n", cex = 0.8, horiz = TRUE, text.width = 0.115)
dev.off()

cat("\n=== Figure 7 saved to output/figures/Figure_factor_structure_chord.tiff ===\n")
