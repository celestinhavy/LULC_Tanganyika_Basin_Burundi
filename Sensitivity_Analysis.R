# OAT sensitivity analysis – one parameter at a time

input_dir  <- "D:/Dissertation/Input"
output_dir <- "D:/Dissertation/Output"

slope_thr <- 9
elev_thr <- 1372
mult <- c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)

library(terra)
library(ggplot2)
library(sf)

# Load reference raster and study area mask
lulc_t1 <- rast(file.path(input_dir, "LULC2020.tif"))
ref <- lulc_t1
study_area_sf <- st_read(file.path(input_dir, "Study_area.shp"), quiet = TRUE)
study_vect <- vect(study_area_sf)
study_mask <- rasterize(study_vect, ref, field = 1)
study_mask <- ifel(is.na(study_mask), NA, 1)

load_rast <- function(f) {
  r <- rast(f)
  r <- crop(r, ref)
  r <- resample(r, ref, method = "near")
  r <- mask(r, study_mask)
  return(r)
}

# Load remaining data and derive degraded areas
lulc_t2 <- rast(file.path(input_dir, "LULC2030.tif"))
dem <- rast(file.path(input_dir, "DEM.tif"))
slope <- rast(file.path(input_dir, "Slope.tif"))
pc <- load_rast(file.path(output_dir, "2030_Protection_Core.tif"))

N_classes <- c(1, 2, 4, 5)
A_classes <- c(3, 6)

degraded <- lulc_t2 %in% A_classes
degraded <- ifel(is.na(degraded), 0, degraded)
degraded <- ifel(degraded > 0, 1, 0)
degraded <- mask(degraded, study_mask)

conversion <- (lulc_t1 %in% N_classes) & (lulc_t2 %in% A_classes)
conversion <- ifel(is.na(conversion), 0, conversion)
conversion <- ifel(conversion > 0, 1, 0)
conversion <- mask(conversion, study_mask)

# Baseline buffers
HI_base <- load_rast(file.path(input_dir, "water_buffer_binary380m.tif"))
HI_base <- ifel(is.na(HI_base), 0, HI_base); HI_base <- ifel(HI_base > 0, 1, 0)
CI_base <- load_rast(file.path(input_dir, "cultural_buffer_binary1000m.tif"))
CI_base <- ifel(is.na(CI_base), 0, CI_base); CI_base <- ifel(CI_base > 0, 1, 0)
PA_base <- load_rast(file.path(input_dir, "protected_buffer_binary1000m.tif"))
PA_base <- ifel(is.na(PA_base), 0, PA_base); PA_base <- ifel(PA_base > 0, 1, 0)

# Evaluation function
eval_zones <- function(HI = HI_base, CI = CI_base, PA = PA_base,
                       s_thr = slope_thr, e_thr = elev_thr) {
  steep <- slope >= s_thr
  steep <- ifel(is.na(steep), 0, steep); steep <- mask(steep, study_mask)
  erosion <- (dem >= e_thr) & (slope >= s_thr)
  erosion <- ifel(is.na(erosion), 0, erosion); erosion <- mask(erosion, study_mask)
  priority <- erosion | (HI == 1) | (CI == 1) | (PA == 1) | steep
  priority <- ifel(is.na(priority), 0, priority); priority <- mask(priority, study_mask)
  rh <- priority & (degraded | conversion) & !pc
  rh <- ifel(is.na(rh), 0, rh); rh <- ifel(rh > 0, 1, 0); rh <- mask(rh, study_mask)
  su <- !pc & !rh
  su <- ifel(is.na(su), 0, su); su <- ifel(su > 0, 1, 0); su <- mask(su, study_mask)
  cell_area <- prod(res(ref)) / 1e6
  data.frame(
    PC = sum(values(pc), na.rm = TRUE) * cell_area,
    RH = sum(values(rh), na.rm = TRUE) * cell_area,
    SU = sum(values(su), na.rm = TRUE) * cell_area
  )
}

# Load buffer by multiplier
get_buffer <- function(type, m) {
  lab <- sprintf("%02d", round(m * 10))
  f <- file.path(input_dir, paste0(type, "_", lab, "x.tif"))
  r <- load_rast(f)
  r <- ifel(is.na(r), 0, r); r <- ifel(r > 0, 1, 0)
  return(r)
}

# Run OAT for each parameter (all using the same multiplier set)
res_all <- data.frame()

for (m in mult) {
  HI <- get_buffer("HI", m)
  tmp <- eval_zones(HI = HI)
  tmp$Param <- "Hydrological"
  tmp$X <- 380 * m
  res_all <- rbind(res_all, tmp)
}

for (m in mult) {
  CI <- get_buffer("CI", m)
  tmp <- eval_zones(CI = CI)
  tmp$Param <- "Cultural"
  tmp$X <- 1000 * m
  res_all <- rbind(res_all, tmp)
}

for (m in mult) {
  PA <- get_buffer("PA", m)
  tmp <- eval_zones(PA = PA)
  tmp$Param <- "Protected"
  tmp$X <- 1000 * m
  res_all <- rbind(res_all, tmp)
}

# Slope: use the same multipliers
slope_vals <- round(slope_thr * mult, 1)
for (s in slope_vals) {
  tmp <- eval_zones(s_thr = s)
  tmp$Param <- "Slope"
  tmp$X <- s
  res_all <- rbind(res_all, tmp)
}

# Elevation: use the same multipliers
elev_vals <- round(elev_thr * mult, 0)
for (e in elev_vals) {
  tmp <- eval_zones(e_thr = e)
  tmp$Param <- "Elevation"
  tmp$X <- e
  res_all <- rbind(res_all, tmp)
}

write.csv(res_all, file.path(output_dir, "OAT_Sensitivity.csv"), row.names = FALSE)

# Visualization (unchanged)
zone_colors <- c("PC" = "#2E86AB", "RH" = "#D95D39", "SU" = "#F0C27B")
zone_labels <- c("PC" = "PC", "RH" = "RH", "SU" = "SU")

make_plot <- function(df, pname, xlab) {
  d <- df[df$Param == pname, ]
  ggplot(d, aes(x = X)) +
    geom_line(aes(y = PC, color = "PC"), linewidth = 1) +
    geom_line(aes(y = RH, color = "RH"), linewidth = 1) +
    geom_line(aes(y = SU, color = "SU"), linewidth = 1) +
    scale_color_manual(values = zone_colors, labels = zone_labels) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                       expand = expansion(mult = c(0.05, 0.1))) +
    labs(x = xlab, y = expression("Area (km"^2*")")) +
    theme_classic() +
    theme(legend.position = c(0.95, 0.5),
          legend.justification = c(1, 0.5),
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
          legend.margin = margin(2, 2, 2, 2),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          plot.title = element_blank(),
          plot.margin = margin(5, 15, 5, 5))
}

panels <- list(
  list(data = res_all, param = "Hydrological", xlab = "Hydrological buffer (m)"),
  list(data = res_all, param = "Cultural",     xlab = "Cultural buffer (m)"),
  list(data = res_all, param = "Protected",    xlab = "Protected area buffer (m)"),
  list(data = res_all, param = "Slope",        xlab = "Slope threshold (°)"),
  list(data = res_all, param = "Elevation",    xlab = "Elevation threshold (m)")
)

for (i in seq_along(panels)) {
  p <- make_plot(panels[[i]]$data, panels[[i]]$param, panels[[i]]$xlab)
  fname <- file.path(output_dir, paste0("OAT_", panels[[i]]$param, ".png"))
  ggsave(fname, plot = p,
         width = 1122/300, height = 1122/300,
         dpi = 300, bg = "white")
}

cat("OAT sensitivity analysis completed. Outputs in", output_dir, "\n")