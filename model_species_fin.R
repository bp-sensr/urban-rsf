# Load packages

save.image("urban-rsf.RData")

library(tidyverse) #V2.0.0
library(sf) #V1.1-1
library(wildrtrax) #V1.5.0
library(tidyr) #V1.3.2
library(dplyr) #V1.2.1
library(terra) #V1.9-27
library(MuMIn) #V1.48.19
library(raster) #V3.6-32
library(sjPlot) #V2.9.0
library(PBSmapping) #V2.74.1
library(ggpubr) #V1.0.0
library(tidyterra)

# Authorize to WildTrax
Sys.setenv(WT_USERNAME = "", WT_PASSWORD = "")
wt_auth()

# WildTrax locations
locations_sf <- wt_get_sync("organization_locations", organization = 5466) |>
  distinct(location, longitude, latitude) |>
  drop_na() |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

locations_coords <- wt_get_sync("organization_locations", organization = 5466) |>
  distinct(location, latitude, longitude) |>
  drop_na()

locations_coords$X <- locations_coords$longitude
locations_coords$Y <- locations_coords$latitude
attributes(locations_coords)$projection <- "LL"
locations_coords1 <- convUL(locations_coords, km = TRUE)
locations_coords1$X <- locations_coords1$X * 1000
locations_coords1$Y <- locations_coords1$Y * 1000

# Coordinate system and 150 m buffers ----
# All spatial layers are transformed to the same projected CRS here.
locations_sf <- st_transform(locations_sf, crs = 32612)
locations_buff <- st_buffer(locations_sf, dist = 150)

# uPLVI landcover ----
uplvi_nat <- read_sf("/users/alexandremacphail/desktop/coe_geospatial/uPLVI/uplvi.shp") |>
  dplyr::select(id, primeclas1, landclas1, stype1, stype2, stype3, stypeper1, stypeper2, stypeper3, ecounit1, ecounit2, ecounit3, geometry)

uplvi_proj <- st_transform(uplvi_nat, crs = 32612)

intersected <- st_intersection(locations_buff, uplvi_proj)

land_class_area <- intersected |>
  mutate(area_m2 = as.numeric(st_area(geometry))) |>
  st_drop_geometry() |>
  group_by(location, primeclas1, landclas1, stype1) |>
  summarise(area_m2 = sum(area_m2)) |>
  ungroup()

buffer_area <- pi * 150^2

# Streetlights ----
light <- read_sf("/users/alexandremacphail/desktop/coe_geospatial/streetlight/streetlight.shp")
light_proj <- st_transform(light, crs = 32612)
intersected_light <- st_intersection(locations_buff, light_proj)

# Watt average was not in the data models, so added here
light_dens <- intersected_light |>
  st_drop_geometry() |>
  group_by(location, type) |>
  summarise(n = n_distinct(id),
            watt_average = mean(watts),
            sd_average = sd(watts)) |>
  ungroup()

# Noise ----
noise <- read_sf("/users/alexandremacphail/desktop/coe_geospatial/COE_2025_noise.shp")
noise_proj <- st_transform(noise, crs = 32612)
intersected_noise <- st_intersection(locations_buff, noise_proj)

noise_avg <- intersected_noise |>
  st_drop_geometry() |>
  group_by(location) |>
  summarise(noise_db = mean(noise)) |>
  ungroup()

# Connectivity ----
connec <- read_sf("/users/alexandremacphail/desktop/coe_geospatial/circuit2026.shp")
connec_proj <- st_transform(connec, crs = 32612)
intersected_connec <- st_intersection(locations_buff, connec_proj)

connec_avg <- intersected_connec |>
  st_drop_geometry() |>
  group_by(location) |>
  summarise(connectivity = mean(connectivi)) |>
  ungroup()

# River and creek distance ----
# River file is too large to load as a vector layer, so distances
# are extracted from rasters instead.
river_creek <- raster("/users/alexandremacphail/desktop/coe_geospatial/COE_Distance_to_RiverCreek_10by10.tif")

# Distance to North Saskatchewan River only
sask_river <- raster("/users/alexandremacphail/desktop/coe_geospatial/nsask_dist.tif")




###### Function and use

wt_rsf <- function(projects = c(1750, 4581, 4559, 3778, 4572, 4575), species) {

  # Load in used/unused data
  used <- map_dfr(projects, ~ wt_download_report(.x, "ARU", "main")) |>
    dplyr::select(location, latitude, longitude, recording_date_time, species_code) |>
    distinct() |>
    group_by(location) |>
    mutate(presence = case_when(species_code == species ~ 1, TRUE ~ 0)) |>
    ungroup() |>
    dplyr::select(location, latitude, longitude, presence) |>
    distinct()

  # Add in proportion cover
  oo <- land_class_area |>
    mutate(prop = area_m2 / buffer_area) |>
    dplyr::select(location, stype1, prop) |>
    distinct() |>
    left_join(used, relationship = "many-to-many") |>
    drop_na() |>
    relocate(presence, .after = location)

  # Adding in lightpost density metrics
  oo1 <- oo |>
    dplyr::select(location, presence, stype1, prop) |>
    distinct() |>
    left_join(light_dens, relationship = "many-to-many") |>
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

  # Add in noise information
  oo2 <- oo1 |>
    dplyr::select(location, presence, stype1, prop, n) |>
    distinct() |>
    left_join(noise_avg, relationship = "many-to-many")

  # Add in river distance information
  oo3 <- merge(oo2, locations_coords1, by = "location")
  oo3$river.dist <- raster::extract(sask_river, oo3[, 9:10])
  oo3$creek.dist <- raster::extract(river_creek, oo3[, 9:10])

  # Add in connectivity information
  oo_fin <- oo3 |>
    dplyr::select(location, presence, stype1, prop, n, noise_db, river.dist, creek.dist) |>
    distinct() |>
    left_join(connec_avg, relationship = "many-to-many")

  wide <- pivot_wider(oo_fin, names_from = stype1, values_from = prop)

  wide[is.na(wide)] <- 0

  # Prep RSF variables
  rsf_ready <- wide |>
    rename(AF  = `Aggregates and/or fill site`,
           AIH = `Transportation Surface`,
           AS  = `Acreage Subdivision`,
           AW  = `Anthropogenic Water Body`,
           BPC = `Building/Parking Complex`,
           CA  = `Annual Crops`,
           CDS = `Commercial/Industrial Development`,
           CP  = `Tame Pasture`,
           CPR = `Rough Pasture`,
           CS  = `Closed Shrub`,
           ECS = `Established Commercial / Industrial`,
           EMS = `Exposed Mineral Soil`,
           ERC = `Established Residental Community`,
           FS  = `Farmyard/Acreage`,
           FT  = `Forested`,
           GF  = `Grass Fen`,
           HG  = `Herbaceous Grass`,
           M   = `Marsh`,
           MG  = `Maintained Grass`,
           MS  = `Medial Shrub`,
           NG  = `Non maintained Grass/Shrubs`,
           NT  = `Nursery/Tree farm`,
           NW  = `Natural Water Body`,
           OG  = `Oil and/or Gas field`,
           RC  = `Recent Clearing`,
           RDS = `Residential Development`,
           SF  = `Shrub Fen`,
           TS  = `Treed Shelterbelt`,
           TT  = `Transplant Trees`)

  rsf_ready1 <- rsf_ready |>
    group_by(location) |>
    mutate(HIGH_DEV = AIH + BPC + RC, .keep = "unused") |>
    mutate(MID_DEV = ERC + ECS, .keep = "unused") |>
    mutate(MOD = TS + HG + CP + FS + TT + CA + NG + CPR + AS + NT, .keep = "unused") |>
    mutate(MOD_MG = MG, .keep = "unused") |>
    mutate(NAT_LAND = EMS, .keep = "unused") |>
    mutate(NAT_WATER = NW, .keep = "unused") |>
    mutate(NAW = FT, .keep = "unused") |>
    mutate(NNW = MS + CS, .keep = "unused") |>
    mutate(VAR_DEV = OG + CDS + RDS, .keep = "unused") |>
    mutate(WET = M + SF + GF, .keep = "unused") |>
    mutate(ANT = AW, .keep = "unused") |>
    mutate(FIL = AF, .keep = "unused")

  # Scale all variables
  scale_vars <- c("n", "HIGH_DEV", "MID_DEV", "MOD", "MOD_MG", "NAT_LAND", "NAT_WATER", "NAW", "NNW", "VAR_DEV", "WET", "ANT", "FIL", "noise_db", "connectivity", "river.dist", "creek.dist")
  rsf_ready1 <- rsf_ready1 |>
    ungroup() |>
    mutate(across(all_of(scale_vars), as.numeric))
  rsf_ready2 <- rsf_ready1 |>
    mutate(across(all_of(scale_vars), ~ as.numeric(scale(.x))))
  rsf_ready3 <- rsf_ready2 |>
    group_by(location) |>
    filter(presence == max(presence)) |>
    ungroup()

  # Run a logistic regression
  null <- glm(presence ~ -1, data = rsf_ready3, family = binomial(link = "logit"), na.action = na.pass)
  mod <- glm(presence ~ -1 + n + noise_db + connectivity + river.dist + creek.dist + HIGH_DEV + MOD + MOD_MG + NAT_LAND + NAT_WATER + NAW + NNW + VAR_DEV + WET + ANT + FIL, data = rsf_ready3, family = binomial(link = "logit"), na.action = na.pass)
  dredge <- dredge(mod, beta = "sd", evaluate = TRUE)

  # Cleanest model sans interactions. Currently modelled for each species manually.
  if (species == "BANS") {
    mod_V1 <- glm(presence ~ -1 + connectivity + noise_db + WET + river.dist, data = rsf_ready3, family = binomial(link = "logit"))
  } else if (species == "TRES") {
    mod_V1 <- glm(presence ~ -1 + FIL + HIGH_DEV + MOD + NNW + NAT_WATER + WET + river.dist, data = rsf_ready3, family = binomial(link = "logit"))
  } else if (species == "CAWA") {
    mod_V1 <- glm(presence ~ -1 + creek.dist + NAW + VAR_DEV, data = rsf_ready3, family = binomial(link = "logit"))
  } else if (species == "YEWA") {
    mod_V1 <- glm(presence ~ -1 + MOD + NAT_WATER + NNW + NAW + WET, data = rsf_ready3, family = binomial(link = "logit"))
  } else {
    stop("No more species.")
  }

  # Adding in interactions
  # NOTE: currently unused below (not plotted, predicted, or returned) — wire it in or remove
  mod_V2 <- glm(presence ~ -1 + MID_DEV + MOD_MG + NAT_LAND + NAW + NNW + river.dist + n + MOD * river.dist + NAT_LAND * river.dist + NAW * river.dist + NNW * river.dist, data = rsf_ready2, family = binomial(link = "logit"))

  # VIF check
  #car::vif(mod)
  #summary(mod)
  #summary(mod_V1)

  # Odds ratio plots
  odds_best <- plot_model(mod_V1, sort.est = TRUE,
                          title = paste(species, "- Best Model"),
                          vline.color = "transparent") +
    geom_hline(yintercept = 1, color = "black", linewidth = 1.5, linetype = "dashed") +
    font_size(title = 22, axis_title.x = 20, labels.x = 20, labels.y = 20)
  odds_full <- plot_model(mod, sort.est = TRUE,
                          title = paste(species, "- Full Model"),
                          vline.color = "transparent") +
    geom_hline(yintercept = 1, color = "black", linewidth = 1.5, linetype = "dashed") +
    font_size(title = 22, axis_title.x = 20, labels.x = 20, labels.y = 20)

  # Starting Map creation ----
  # base_map <- rast("/users/alexandremacphail/desktop/coe_geospatial/testh1.tif") # Get from Google Drive or from env
  names(base_map) <- c("ANT", "creek.dist", "noise_db", "connectivity", "FIL", "HIGH_DEV", "n", "MID_DEV", "MOD_MG", "MOD", "NAT_LAND", "NAT_WATER", "NAW", "NNW", "VAR_DEV", "WET", "river.dist")

  # Standardize raster using the same mean and SD used for model fitting
  base_map1 <- base_map
  for (v in scale_vars) {base_map1[[v]] <- (base_map[[v]] - mean(rsf_ready1[[v]], na.rm = TRUE)) / sd(rsf_ready1[[v]], na.rm = TRUE)}

  # Maps ----
  map.rsf  <- predict(base_map1, mod, type = "response")     # from FULL model
  map.rsf1 <- predict(base_map1, mod_V1, type = "response")  # from BEST/REDUCED model

  return(list(full = map.rsf,        # Full model (all predictors)
              reduced = map.rsf1,    # Best/reduced model (species-specific subset)
              dredge_odds = odds_best,
              full_odds = odds_full,
              model_full = mod,
              model_reduced = mod_V1,
              data = rsf_ready3))

}

bans_rsf <- wt_rsf(species = "BANS")
tres_rsf <- wt_rsf(species = "TRES")
cawa_rsf <- wt_rsf(species = "CAWA")
yewa_rsf <- wt_rsf(species = "YEWA")

# To plot RSF best model
#plot(bans_rsf[[1]])

