# install and load packages
packages <- c("terra", "dplyr", "readr")
new_packages <- packages[!packages %in% rownames(installed.packages())]
if (length(new_packages) > 0) install.packages(new_packages, dependencies = TRUE)
library(terra)
library(dplyr)
library(readr)

# paths
input_dir <- "D:/Dissertation/Process/Classified_rasters"
output_dir <- "D:/Dissertation/Output/Tables"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# lulc classes and years
classes <- c("Forest", "Rangeland", "Cropland", "Wetland", "Water", "Built-up")
years <- c(2000, 2010, 2020, 2030, 2040)

# load lulc rasters
lulc <- lapply(years, function(y) rast(file.path(input_dir, paste0("LULC", y, ".tif"))))
names(lulc) <- as.character(years)

# pixel area
pixel_area <- prod(res(lulc[[1]])) / 1e6

# lulc area and percentage
area_table <- bind_rows(lapply(seq_along(lulc), function(i) {
  f <- as.data.frame(freq(lulc[[i]]))
  f <- f[!is.na(f$value) & f$value %in% 1:6, ]
  data.frame(year = years[i], class = classes[as.integer(f$value)],
             pixels = f$count, area_km2 = f$count * pixel_area)
})) %>%
  group_by(year) %>%
  mutate(percentage = area_km2 / sum(area_km2) * 100) %>%
  ungroup()

write_csv(area_table, file.path(output_dir, "LULC_area_percentage.csv"))

# temporal change
periods <- data.frame(t1 = c(2000, 2010, 2020, 2030),
                      t2 = c(2010, 2020, 2030, 2040))

change_table <- bind_rows(lapply(seq_len(nrow(periods)), function(i) {
  a <- area_table[area_table$year == periods$t1[i], c("class", "area_km2")]
  b <- area_table[area_table$year == periods$t2[i], c("class", "area_km2")]
  names(a)[2] <- "area_t1_km2"
  names(b)[2] <- "area_t2_km2"
  x <- merge(a, b, by = "class")
  x$t1 <- periods$t1[i]
  x$t2 <- periods$t2[i]
  x$change_km2 <- x$area_t2_km2 - x$area_t1_km2
  x$change_percent <- x$change_km2 / x$area_t1_km2 * 100
  x$annual_rate <- x$change_percent / (x$t2 - x$t1)
  x[c("t1", "t2", "class", "area_t1_km2", "area_t2_km2",
      "change_km2", "change_percent", "annual_rate")]
}))

write_csv(change_table, file.path(output_dir, "LULC_change.csv"))

# transition function
transition <- function(r1, r2, t1, t2) {
  a <- values(r1, mat = FALSE)
  b <- values(r2, mat = FALSE)
  keep <- !is.na(a) & !is.na(b) & a %in% 1:6 & b %in% 1:6
  a <- as.integer(a[keep])
  b <- as.integer(b[keep])
  tab <- table(factor(a, levels = 1:6), factor(b, levels = 1:6))
  x <- expand.grid(from = 1:6, to = 1:6)
  x$pixels <- as.vector(tab)
  x$area_km2 <- x$pixels * pixel_area
  x$from_class <- classes[x$from]
  x$to_class <- classes[x$to]
  x$t1 <- t1
  x$t2 <- t2
  x[c("t1", "t2", "from", "to", "from_class", "to_class",
      "pixels", "area_km2")]
}

# transition matrices
transition_table <- bind_rows(lapply(seq_len(nrow(periods)), function(i)
  transition(lulc[[as.character(periods$t1[i])]],
             lulc[[as.character(periods$t2[i])]],
             periods$t1[i], periods$t2[i])
))

write_csv(transition_table, file.path(output_dir, "LULC_transition_matrices.csv"))

# major transitions
major_transitions <- transition_table[
  transition_table$from != transition_table$to &
    transition_table$pixels > 0, ]

major_transitions <- major_transitions[
  order(major_transitions$t1, major_transitions$t2,
        -major_transitions$area_km2), ]

write_csv(major_transitions, file.path(output_dir, "LULC_major_transitions.csv"))

# projected 2020–2040 transitions
transition_2020_2040 <- transition(lulc[["2020"]], lulc[["2040"]], 2020, 2040)
transition_2020_2040 <- transition_2020_2040[
  transition_2020_2040$from != transition_2020_2040$to &
    transition_2020_2040$pixels > 0, ]
transition_2020_2040 <- transition_2020_2040[
  order(-transition_2020_2040$area_km2), ]

write_csv(transition_2020_2040,
          file.path(output_dir, "LULC_transitions_2020_2040.csv"))

cat("\nLULC analysis completed successfully.\n")
cat("Results saved to:", output_dir, "\n")