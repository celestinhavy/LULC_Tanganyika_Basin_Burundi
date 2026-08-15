# Spatial prioritization for sustainable freshwater basin management 
# based on projected land use and land cover (LULC).
# The analysis presented here applies the 2040 LULC projection scenario.
# ============================================================================
# USER-DEFINED PARAMETERS
input_dir  <- "D:/Dissertation/Input"
output_dir <- "D:/Dissertation/Output"

forest_class <- 1; rangeland_class <- 2; cropland_class <- 3
wetland_class <- 4; water_class <- 5; builtup_class <- 6

slope_threshold <- 9; water_buffer <- 380; cultural_buffer <- 1000; pa_buffer <- 1000

# LOAD PACKAGES
library(terra); library(segmented); library(sf); library(ggplot2)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# LOAD INPUT DATA
lulc_t1 <- rast(file.path(input_dir, "LULC2020.tif"))
lulc_t2 <- rast(file.path(input_dir, "LULC2040.tif"))
dem <- rast(file.path(input_dir, "DEM.tif"))
slope <- rast(file.path(input_dir, "Slope.tif"))
HI <- rast(file.path(input_dir, "water_buffer_binary380m.tif"))
CI <- rast(file.path(input_dir, "cultural_buffer_binary1000m.tif"))
PA <- rast(file.path(input_dir, "protected_buffer_binary1000m.tif"))
study_area_sf <- st_read(file.path(input_dir, "Study_area.shp"), quiet = TRUE)

# REFERENCE RASTER & MASK
ref <- lulc_t1
study_vect <- vect(study_area_sf)
study_mask <- rasterize(study_vect, ref, field = 1)
study_mask <- ifel(is.na(study_mask), NA, 1)

# HARMONIZE RASTERS
harmonize <- function(r) {
  r <- crop(r, ref); r <- resample(r, ref, method = "near"); r <- mask(r, study_mask); return(r)
}

lulc_t1 <- harmonize(lulc_t1); lulc_t2 <- harmonize(lulc_t2)
dem <- harmonize(dem); slope <- harmonize(slope)
HI <- harmonize(HI); CI <- harmonize(CI); PA <- harmonize(PA)

study_mask <- ifel(is.na(lulc_t1), NA, 1)

HI <- ifel(is.na(HI), 0, HI); HI <- ifel(HI > 0, 1, 0)
CI <- ifel(is.na(CI), 0, CI); CI <- ifel(CI > 0, 1, 0)
PA <- ifel(is.na(PA), 0, PA); PA <- ifel(PA > 0, 1, 0)

# ELEVATION THRESHOLD
natural <- lulc_t1 %in% c(forest_class, wetland_class, water_class)
e <- values(dem); n <- values(natural)
valid <- !is.na(e) & !is.na(n); e <- e[valid]; n <- n[valid]

bw <- 2 * IQR(e) / length(e)^(1/3)
n_bins <- max(20, ceiling(diff(range(e)) / bw))
breaks <- seq(min(e), max(e), length.out = n_bins + 1)
bins <- cut(e, breaks = breaks, include.lowest = TRUE)

prop_natural <- tapply(n, bins, mean, na.rm = TRUE)
mid_elev <- (head(breaks, -1) + tail(breaks, -1)) / 2
df <- data.frame(elev = mid_elev, prop = prop_natural); df <- na.omit(df)

lm_fit <- lm(prop ~ elev, data = df)
seg_fit <- try(segmented(lm_fit, seg.Z = ~elev), silent = TRUE)

if (inherits(seg_fit, "try-error")) {
  spline_fit <- smooth.spline(df$elev, df$prop, spar = 0.6)
  pred <- predict(spline_fit, df$elev)
  e_thr <- df$elev[which.max(diff(pred$y))]
} else {
  e_thr <- seg_fit$psi[2]
}

elevation_threshold <- round(e_thr, 0)

# PROTECTION CORE
N_classes <- c(forest_class, rangeland_class, wetland_class, water_class)
A_classes <- c(cropland_class, builtup_class)

pc <- (lulc_t1 == forest_class & lulc_t2 == forest_class) |
  (lulc_t1 == wetland_class & lulc_t2 == wetland_class) |
  (lulc_t1 == water_class & lulc_t2 == water_class)
pc <- ifel(is.na(pc), 0, pc); pc <- ifel(pc > 0, 1, 0); pc <- mask(pc, study_mask)
writeRaster(pc, file.path(output_dir, "2040_Protection_Core.tif"), overwrite = TRUE)

# DEGRADED AREAS & CONVERSION RISK
degraded <- lulc_t2 %in% A_classes
degraded <- ifel(is.na(degraded), 0, degraded); degraded <- ifel(degraded > 0, 1, 0); degraded <- mask(degraded, study_mask)

conversion <- (lulc_t1 %in% N_classes) & (lulc_t2 %in% A_classes)
conversion <- ifel(is.na(conversion), 0, conversion); conversion <- ifel(conversion > 0, 1, 0); conversion <- mask(conversion, study_mask)

# RESTORATION PRIORITY CONDITIONS
steep <- slope >= slope_threshold
steep <- ifel(is.na(steep), 0, steep); steep <- mask(steep, study_mask)

erosion <- (dem >= elevation_threshold) & (slope >= slope_threshold)
erosion <- ifel(is.na(erosion), 0, erosion); erosion <- mask(erosion, study_mask)

priority <- erosion | (HI == 1) | (CI == 1) | (PA == 1) | steep
priority <- ifel(is.na(priority), 0, priority); priority <- mask(priority, study_mask)

# RESTORATION HOTSPOTS
rh <- priority & (degraded | conversion) & !pc
rh <- ifel(is.na(rh), 0, rh); rh <- ifel(rh > 0, 1, 0); rh <- mask(rh, study_mask)
writeRaster(rh, file.path(output_dir, "2040_Restoration_Hotspots.tif"), overwrite = TRUE)

# SUSTAINABLE USE AREAS
su <- !pc & !rh
su <- ifel(is.na(su), 0, su); su <- ifel(su > 0, 1, 0); su <- mask(su, study_mask)
writeRaster(su, file.path(output_dir, "2040_Sustainable_Use.tif"), overwrite = TRUE)

# PRIORITY MAP
pm <- ifel(pc == 1, 1, ifel(rh == 1, 2, ifel(su == 1, 3, NA)))
pm <- mask(pm, study_mask)
levels(pm) <- data.frame(ID = 1:3, Zone = c("Protection Core", "Restoration Hotspot", "Sustainable Use"))
writeRaster(pm, file.path(output_dir, "2040_Priority_Map.tif"), overwrite = TRUE)

# AREA STATISTICS
calc_stats <- function(study_mask, pc, rh, su, ref, output_dir) {
  total_cells <- sum(values(study_mask), na.rm = TRUE)
  cell_area <- prod(res(ref)) / 1e6
  total_area <- total_cells * cell_area
  
  get_area <- function(r, name) {
    r <- mask(r, study_mask)
    cells <- sum(values(r), na.rm = TRUE)
    area <- cells * cell_area
    pct <- (area / total_area) * 100
    return(data.frame(Zone = name, Area_km2 = round(area, 2), Percentage = round(pct, 2)))
  }
  
  stats <- rbind(get_area(pc, "Protection Core"), 
                 get_area(rh, "Restoration Hotspot"), 
                 get_area(su, "Sustainable Use"))
  stats <- rbind(stats, data.frame(Zone = "TOTAL", Area_km2 = round(total_area, 2), Percentage = 100))
  
  write.csv(stats, file.path(output_dir, "2040_Area_Statistics.csv"), row.names = FALSE)
  
  cat("\n=== AREA STATISTICS 2040 ===\n")
  cat(sprintf("Study area: %.2f km²\n", total_area))
  cat("----------------------------------------\n")
  for (i in 1:(nrow(stats)-1)) {
    cat(sprintf("  %-20s: %9.2f km² (%5.2f %%)\n", stats$Zone[i], stats$Area_km2[i], stats$Percentage[i]))
  }
  cat("----------------------------------------\n")
  cat(sprintf("  %-20s: %9.2f km² (%5.2f %%)\n", stats$Zone[nrow(stats)], stats$Area_km2[nrow(stats)], stats$Percentage[nrow(stats)]))
  return(stats)
}
area_stats <- calc_stats(study_mask, pc, rh, su, ref, output_dir)

# BAR CHART VISUALIZATION
create_barchart <- function(area_stats, output_dir, dpi = 300) {
  stats_plot <- area_stats[1:3, ]
  stats_plot$Zone <- factor(stats_plot$Zone, 
                            levels = c("Protection Core", "Sustainable Use", "Restoration Hotspot"))
  
  p <- ggplot(stats_plot, aes(x = Zone, y = Area_km2, fill = Zone)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = paste0(round(Area_km2, 2), " km²\n(", Percentage, "%)")), 
              vjust = -0.25, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = c("#2E86AB", "#F0C27B", "#D95D39")) +
    scale_x_discrete(labels = c("PC", "SU", "RH")) +
    labs(x = "Priority Zone", y = "Area (km²)") +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.line = element_line(color = "black", linewidth = 0.5),
          plot.margin = margin(5, 5, 5, 5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  ggsave(filename = file.path(output_dir, "2040_Priority_Zone_Area.png"),
         plot = p, width = 1122/dpi, height = 1122/dpi, dpi = dpi, bg = "white")
  
  return(p)
}
create_barchart(area_stats, output_dir)
cat(" RUNNING COMPLETED\n")