# Load required packages
library(raster)
library(circlize)

# Set working directory and LULC classes
setwd("D:/Dissertation/Process/Classified_rasters/")
lulc_classes <- c("Forest","Rangeland","Cropland","Wetland","Water","Built-up")
lulc_colors <- setNames(
  c("#065106","#A6F040","#F3E537","#00CED1","#1555E4","#FF1000"),
  lulc_classes
)

# Input raster reading and validation
cat("Loading LULC rasters for 2010 and 2020...\n")
lulc_2010 <- raster("LULC2010.tif")
lulc_2020 <- raster("LULC2020.tif")

if (!compareRaster(lulc_2010, lulc_2020, extent = TRUE, rowcol = TRUE, crs = TRUE, res = TRUE)) {
  stop("ERROR: The rasters do not match! Please check their properties.")
} else {
  cat("Raster properties match.\n")
}

# Extract pixel values and build transition matrix
v1 <- getValues(lulc_2010)
v2 <- getValues(lulc_2020)

valid <- !is.na(v1) & !is.na(v2)
v1 <- v1[valid]
v2 <- v2[valid]

tm <- matrix(0, nrow = 6, ncol = 6, dimnames = list(lulc_classes, lulc_classes))

for (k in seq_along(v1)) {
  i <- v1[k]
  j <- v2[k]
  if (i >= 1 && i <= 6 && j >= 1 && j <= 6) {
    tm[i, j] <- tm[i, j] + 1
  }
}

write.csv(tm, "transition_matrix_2010_2020.csv", row.names = TRUE)

# Compute link start and end positions
total_pixels <- sum(tm)
row_totals <- rowSums(tm)
col_totals <- colSums(tm)

get_positions <- function(vals, total) {
  if (total == 0) return(list(starts = numeric(0), ends = numeric(0)))
  props <- vals / total
  starts <- cumsum(c(0, props))[1:length(vals)]
  ends <- starts + props
  list(starts = starts, ends = ends)
}

from_pos <- lapply(seq_len(6), function(i) get_positions(tm[i, ], row_totals[i]))
to_pos   <- lapply(seq_len(6), function(j) get_positions(tm[, j], col_totals[j]))

# Initialize output device and Circos layout
tiff("ChordDiagram_2010_2020.tiff",
     width  = 1122 / 300,
     height = 1122 / 300,
     units  = "in",
     res    = 300,
     compression = "lzw")

par(mar = rep(0, 4), oma = rep(0, 4), xpd = NA)

circos.clear()
circos.par(
  start.degree  = 90,
  gap.degree    = 3,
  track.margin  = c(0, 0),
  cell.padding  = c(0, 0, 0, 0),
  canvas.xlim   = c(-1.03, 1.03),
  canvas.ylim   = c(-1.03, 1.03)
)

circos.initialize(factors = lulc_classes, xlim = matrix(c(rep(0, 6), rep(1, 6)), ncol = 2))

# Draw sector tracks and class labels
circos.trackPlotRegion(
  track.index = 1,
  ylim        = c(0, 1),
  track.height = 0.18,
  bg.col      = lulc_colors,
  bg.border   = "black",
  panel.fun = function(x, y) {
    sec   <- get.cell.meta.data("sector.index")
    xlim  <- get.cell.meta.data("xlim")
    ylim  <- get.cell.meta.data("ylim")
    circos.text(
      x = mean(xlim), y = mean(ylim),
      labels = sec,
      facing = "bending.inside",
      niceFacing = TRUE,
      cex   = 1.0,
      col   = "black",
      font  = 2,
      family = "Arial"
    )
  }
)

circos.trackPlotRegion(
  track.index = 2,
  ylim        = c(0, 1),
  track.height = 0.05,
  bg.col      = NA,
  bg.border   = NA,
  panel.fun   = function(x, y) {}
)

# Draw transition links with scaled widths
min_lwd <- 0.5
max_lwd <- 15

for (i in seq_len(6)) {
  for (j in seq_len(6)) {
    val <- tm[i, j]
    if (val > 0) {
      rp <- from_pos[[i]]
      tp <- to_pos[[j]]
      x1 <- c(rp$starts[j], rp$ends[j])
      x2 <- c(tp$starts[i], tp$ends[i])
      
      prop <- val / total_pixels
      lwd_val <- min_lwd + prop * (max_lwd - min_lwd)
      
      circos.link(
        sector.index1 = lulc_classes[i],
        point1        = x1,
        sector.index2 = lulc_classes[j],
        point2        = x2,
        col           = adjustcolor(lulc_colors[lulc_classes[i]], alpha.f = 0.6),
        border        = NA,
        lwd           = lwd_val
      )
    }
  }
}

# Add percentage tick marks and labels
tick_positions <- seq(0, 1, length.out = 5)
label_positions <- tick_positions
label_texts <- c("0%", "25%", "50%", "75%")

op <- par(family = "Arial")

for (i in seq_along(lulc_classes)) {
  circos.axis(
    h             = "top",
    major.at      = tick_positions,
    labels        = rep("", length(tick_positions)),
    sector.index  = lulc_classes[i],
    track.index   = 1,
    col           = "black",
    lwd           = 0.6,
    major.tick.length = 0.3
  )
  
  for (j in seq_along(label_positions)) {
    x <- label_positions[j]
    y <- 1.4
    circos.text(
      x           = x,
      y           = y,
      labels      = label_texts[j],
      sector.index = lulc_classes[i],
      track.index  = 1,
      facing      = "outside",
      niceFacing  = TRUE,
      cex         = 0.6,
      col         = "black",
      font        = 2,
      family      = "Arial"
    )
  }
}

par(op)

# Save output and close device
dev.off()
cat("✔ Chord diagram saved as 'ChordDiagram_2010_2020.tiff' in the working directory.\n")