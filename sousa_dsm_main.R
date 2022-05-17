##%######################################################%##
#                                                          #
#'###     Contrasting aerial abundance estimates for     ####
#'###          two sympatric dolphin species at          ####
#'###          a regional scale using distance           ####
#'###       sampling and density surface modelling       ####
#                                                          #
##%######################################################%##

# Authors: Raudino et al. (2022)

# Required libraries ------------------------------------------------------

library(pacman)
pacman::p_load(raster,         # Geographic Data Analysis and Modeling 
               magrittr,       # Pipes
               sf,             # Simple features for r
               here,           # Simpler way to find files
               sp,             # Classes and Methods for Spatial Data
               Distance,       # Distance Sampling
               dsm,            # Density surface modelling
               ggpubr,         # 'ggplot2' Based Publication Ready Plots
               ggspatial,      # Spatial Data Framework for ggplot2
               patchwork,      # The Composer of Plots
               ggsn,           # North arrows and scale bars
               smoothr,        # Smoothing spatial features
               evaluate,       # Parsing and Evaluation Tools
               viridisLite,    # Colour palettes
               tidyverse)      # Tools and functions for data science

# Set tibble options
options(tibble.width = Inf) # All tibble columns shown
options(pillar.neg = FALSE) # No colouring negative numbers
options(pillar.subtle = TRUE)
options(pillar.sigfig = 5)

spatial.scale <- 5000
dolphin.species <- c("tursiops", "sousa")
survey.years <- list(2015, c(2016, 2017))
names(survey.years) <- c("2015", "2016_2017")

source("sousa_dsm_functions.R")
source("soap_check.R")

load("data/xy_area.RData")
load("data/sousa_mrds.RData") # MRDS analysis

set.seed(354870)

# GIS ---------------------------------------------------------------------

# Define geographic coordinate systems
crs.latlon <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
crs.utm <- sp::CRS("+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

#'------------------------------------------
# Survey transects 
#'------------------------------------------

# Import survey transects
transects <- suppressWarnings(purrr::map(
  .x = names(survey.years),
  .f = ~ {
    if (.x == "2015") {
      shp <- sf::st_read("gis/2015_transects.shp")
    } else {
      shp <- sf::st_read("gis/2016_transects.shp")
    }
    sf::st_crs(shp) <- crs.utm
    shp %>% dplyr::filter(!st_is_empty(.))
  }
) %>%
  purrr::set_names(x = ., nm = names(survey.years)))

# Assign route IDs
transects$`2015`$route <- paste0("T2015-", 1:nrow(transects$`2015`))
transects$`2016_2017`$route <- paste0("T201617-", 1:nrow(transects$`2016_2017`))

transects$`2015` <- dplyr::select(transects$`2015`, X, Y, Length, geometry, route)
transects$`2016_2017` <- dplyr::select(transects$`2016_2017`, X, Y, Length, geometry, route)

#'------------------------------------------
# Study area
#'------------------------------------------

study_area <- purrr::map(
  .x = survey.years,
  .f = ~ sf::st_read(paste0("gis/pilbara_", paste0(.x, collapse = "_"), ".shp")) %>%
    sf::st_transform(x = ., crs = crs.utm)
) %>%
  purrr::set_names(x = ., nm = names(survey.years))

# Add a small buffer (will be useful to avoid edge effects later)
sbuffer <- purrr::map(.x = study_area, .f = ~sf::st_buffer(x = .x, dist = 5000))

#'------------------------------------------
# WA coastline
#'------------------------------------------

# Import coastline / landmass shapefile, project to UTM and crop
wa <- sf::st_read("gis/Australian_mainland.shp")
wa <- sf::st_transform(x = wa, crs = crs.utm)
wa <- sf::st_crop(x = wa, y = raster::extent(setNames(c(117756.7, 838920.7, 7444076, 7873106), c("xmin", "xmax", "ymin", "ymax"))))

# Create spatial object to define the corresponding area of ocean
wa_extent <- purrr::map(.x = sbuffer, .f = ~sf::st_as_sfc(sf::st_bbox(.x)))
wa_ocean <- purrr::map(.x = wa_extent, .f = ~sf::st_difference(st_union(.x), st_union(wa)))

# Quick example plot
par(mfrow = c(2, 1))
purrr::walk(
  .x = names(survey.years),
  .f = ~ {
    plot(sf::st_geometry(wa_ocean[[.x]]), col = "lightblue")
    plot(sf::st_geometry(study_area[[.x]]), add = TRUE)
  }
)
par(mfrow = c(1, 1))

# Refined study area
study_area_coastline <- purrr::map2(
  .x = wa_ocean,
  .y = study_area,
  .f = ~ sf::st_intersection(x = .x, y = .y)
)

#'------------------------------------------
# Geographic coordinates
#'------------------------------------------

# Create rasters of x and y
xy_area <- purrr::map(.x = raster_area, .f = ~raster::as.data.frame(x = .x, xy = TRUE))

x_area <- purrr::map(.x = names(survey.years),
                     .f = ~{xy_area[[.x]] %>%
                         dplyr::select(x, y) %>% 
                         dplyr::mutate(layer = x) %>% 
                         raster::rasterFromXYZ(.) %>% 
                         raster::crop(x = ., y = raster_area[[.x]]) %>% 
                         raster::mask(x = ., mask = raster_area[[.x]])}) %>% 
  purrr::set_names(x = ., nm = names(survey.years))

y_area <- purrr::map(.x = names(survey.years),
                     .f = ~{xy_area[[.x]] %>%
                         dplyr::select(x, y) %>% 
                         dplyr::mutate(layer = y) %>% 
                         raster::rasterFromXYZ(.) %>% 
                         raster::crop(x = ., y = raster_area[[.x]]) %>% 
                         raster::mask(x = ., mask = raster_area[[.x]])}) %>% 
  purrr::set_names(x = ., nm = names(survey.years))

#'------------------------------------------
# Combine rasters into a stack
#'------------------------------------------

env <- purrr::map(
  .x = names(survey.years),
  .f = ~ raster::stack(x_area[[.x]], y_area[[.x]])
) %>% purrr::set_names(x = ., nm = names(survey.years))

env <- purrr::map(.x = env, .f = ~ {
  names(.x) <- c("x", "y")
  .x})

#'------------------------------------------
# Hex grid
#'------------------------------------------

# Create hexagon grid for predictions
hexgrid <- purrr::map(.x = study_area,
                      .f = ~{
                        sf::st_make_grid(x = .x,
                                         cellsize = spatial.scale, # in m
                                         what = "polygons",
                                         crs = sf::st_crs(wa),
                                         square = FALSE) %>% 
                          sf::st_sf(ID = 1:length(lengths(.)), .)}) # Add index

# Quick plots
par(mfrow = c(2, 1))
purrr::walk(.x = names(survey.years), .f = ~{
  plot(raster_area[[.x]])
  plot(st_geometry(hexgrid[[.x]]), add = TRUE)
  plot(st_geometry(wa), add = TRUE, col = "grey")
  plot(st_geometry(study_area[[.x]]), add = TRUE, col = "transparent")
})
par(mfrow = c(1, 1))

# Extract mean values of covariates within each hexagon
hex_env <- purrr::map2(
  .x = hexgrid,
  .y = env,
  .f = ~ raster::extract(
    x = .y,
    y = sf::as_Spatial(.x),
    fun = mean, na.rm = TRUE, df = TRUE
  )
)

hexgrid <- purrr::map2(
  .x = hexgrid,
  .y = hex_env,
  .f = ~ {
    dplyr::left_join(x = .x, y = .y, by = "ID") %>%
      dplyr::filter(is.na(x) == FALSE)
  }
)

# Crop grid
area_ocean <- purrr::map(.x = study_area, .f = ~sf::st_difference(st_union(.x), st_union(wa)))
hexgrid <- purrr::map2(.x = hexgrid, .y = area_ocean, .f = ~sf::st_intersection(x = .x, y = .y))

# Calculate the area of each hexagon (in m)
hexgrid <- purrr::map(.x = hexgrid, .f = ~{
  .x$area <- as.numeric(sf::st_area(.x))
  .x})

# Example maps
purrr::walk(.x = names(survey.years), .f = ~{
  
  p1 <- ggplot(data = hexgrid[[.x]]) + 
    geom_sf(aes(fill = area), col = "transparent")
  
  gg.add <- function(pp){
    res <- pp +
      geom_sf(data = wa) +
      scale_fill_gradientn(colours = pals::viridis(1000)) +
      geom_sf(data = study_area[[.x]], fill = "transparent") +
      geom_sf(data = wa) +
      coord_sf(xlim = sf::st_bbox(study_area[[.x]])[c(1,3)], 
               ylim = sf::st_bbox(study_area[[.x]])[c(2,4)], expand = FALSE) +
      theme_minimal()
    return(res)}
  
  print(gg.add(p1))
  
})

# Data -------------------------------------------------------------

pred.df <- purrr::map(.x = hexgrid, .f = ~{
    tmp <- sf::st_centroid(.x) %>% sf::st_coordinates(.) %>% data.frame(.) %>%
      cbind(., .x$area)
    names(tmp) <- c("x", "y", "area")
    tmp})

#'------------------------------------------
# Segment data
#'------------------------------------------
# Import and compile the segment data
seg.df <- purrr::map_df(
  .x = unlist(survey.years),
  .f = ~ {
    readr::read_csv(paste0("data/", .x, "_datafiles/", .x, "_segdata", spatial.scale, ".csv")) %>%
      dplyr::mutate(
        year = .x,
        Transect.Label = paste0("T", year, "-", Transect.Label),
        Sample.Label = paste0("S", year, "-", Sample.Label)
      )
  }
)

seg.df <- seg.df %>% 
  dplyr::select(-latitude, -longitude, -dplyr::starts_with("depth")) %>% 
  sf::st_as_sf(x = ., coords = c("x", "y"), remove = FALSE) 

seg.df <- split(x = seg.df, f = seg.df$year)
seg.df$`2016_2017` <- rbind(seg.df$`2016`, seg.df$`2017`)
seg.df$`2016` <- seg.df$`2017` <- NULL

seg.df <- purrr::map(.x = names(survey.years),
                     .f = ~{seg.df[[.x]] %>%
                         sf::st_set_crs(x = ., value = crs.utm)}) %>% 
  purrr::set_names(x = ., nm = names(survey.years))

# Retrieve route IDs
route_id <- purrr::map2(.x = seg.df, 
                        .y = transects, 
                        .f = ~.y$route[sf::st_nearest_feature(.x, .y)])

seg.df <- purrr::map2(.x = seg.df,
                      .y = route_id,
                      .f = ~ { .x$route <- .y
                      .x})

# Check that routes have been correctly allocated
purrr::map(.x = names(seg.df),
           .f = ~{
             ggplot() + 
               geom_sf(data = transects[[.x]]["route"], aes(col = as.factor(route))) +
               geom_sf(data = seg.df[[.x]][, "route"], aes(col = as.factor(route)), shape = 1) +
               theme(legend.position = "none")})

# Convert back to standard data.frame
seg.df <- purrr::map(.x = seg.df, .f = ~sf::st_drop_geometry(.x))

#'------------------------------------------
# Distance data
#'------------------------------------------
# Import and compile the distance sampling data
dist.files <- purrr::map(
  .x = dolphin.species,
  .f = ~ {
    lapply(X = survey.years, FUN = function(x) 
      here::here(paste0("data/", x, "_datafiles/", x, "_distdata_", .x, spatial.scale, "_truncated.csv")))}) %>% purrr::set_names(x = ., nm = dolphin.species)

dist.df <- 
  lapply(X = dolphin.species, FUN = function(u){
    purrr::map_df(
      .x = 1:length(unlist(dist.files[[u]])),
      .f = ~ {
        if (file.exists(unlist(dist.files[[u]])[[.x]])) {
          readr::read_csv(unlist(dist.files[[u]])[[.x]]) %>%
            janitor::clean_names(.) %>%
            dplyr::select(object, size, distance, effort, detected,
                          beaufort, x, y, turbidity) %>%
            dplyr::mutate(year = unlist(survey.years)[.x],
                          object = paste0("obj", year, "-", object)) %>% 
            dplyr::rename(Effort = effort) %>% 
            dplyr::mutate(dist_band = cut_width(distance, width = 100, boundary = 0))
        }})
  }) %>% purrr::set_names(., nm = dolphin.species)

#'------------------------------------------
# Observation data
#'------------------------------------------
# Import and compile the observation data

obs.files <- purrr::map(
  .x = dolphin.species,
  .f = ~ {
    lapply(X = survey.years, FUN = function(x) here::here(paste0("data/", x, "_datafiles/", x, "_obsdata_", .x, spatial.scale, "_truncated.csv")))}) %>% purrr::set_names(x = ., nm = dolphin.species)

obs.df <- 
  lapply(X = dolphin.species, FUN = function(u){
    purrr::map_df(
      .x = 1:length(unlist(obs.files[[u]])),
      .f = ~ {
        if (file.exists(unlist(obs.files[[u]])[[.x]])) {
          readr::read_csv(unlist(obs.files[[u]])[[.x]]) %>%
            janitor::clean_names(.) %>%
            dplyr::rename(Effort = effort, Sample.Label = sample_label) %>% 
            dplyr::mutate(year = unlist(survey.years)[.x],
                          object = paste0("obj", year, "-", object)) %>% 
            dplyr::mutate(Sample.Label = paste0("S", year, "-", Sample.Label))
        }})
  }) %>% purrr::set_names(., nm = dolphin.species)

# Dolphin sightings
sightings <- purrr::map(.x = setNames(dolphin.species, dolphin.species),
                        .f = ~ dplyr::left_join(x = obs.df[[.x]], 
                                                y = dist.df[[.x]][, c("x", "y", "object")], by = "object"))

# Exploratory data analysis -------------------------------------------------------------

# Observation distance
purrr::walk(.x = dolphin.species, 
            .f = ~{
              p <- ggplot(dist.df[[.x]], aes(distance)) + 
                geom_histogram(col = "black", binwidth = 30) +
                facet_grid(rows = vars(year), scales = "free_y") + ggtitle(.x)
              print(p)})

# Group size
purrr::walk(.x = dolphin.species, 
            .f = ~{
              p <- ggplot(dist.df[[.x]], aes(size)) + 
                geom_histogram(col = "black", binwidth = 1) +
                facet_grid(rows = vars(year), scales = "free_y") + ggtitle(.x)
              print(p)})

par(mfrow = c(2,2))

# Group size as a function of distance from the trackline and Beaufort sea state
purrr::walk(.x = dolphin.species, 
            .f = ~{
              
              # plots of distance vs. cluster size
              plot(dist.df[[.x]]$distance, dist.df[[.x]]$size, main = .x, xlab = "Distance (m)",
                   ylab = "Group size", pch = 19, cex = 0.5, col = gray(0.7))
              
              # lm fit
              l.dat <- data.frame(distance = seq(0, 500, len = 1000))
              lo <- lm(size~distance, data = dist.df[[.x]])
              lines(l.dat$distance, as.vector(predict(lo, l.dat)))
              
              plot(dist.df[[.x]]$distance, dist.df[[.x]]$beaufort, main = .x, xlab = "Distance (m)",
                   ylab = "Beaufort sea state", pch = 19, cex = 0.5, col = gray(0.7))
              
            })

purrr::map(.x = dolphin.species, 
           .f = ~{
             p <- ggplot(data = dist.df[[.x]]) + 
               geom_boxplot(aes(x = dist_band, y = size), fill = "grey", notch = FALSE) +
               ggtitle(.x) + xlab("Distance band") + ylab("Group size") +
               theme_minimal(); p})


# Detection function models -------------------------------------------------------------

trunc.dist <- list(tursiops = 520, sousa = 450)

detfc.models <- c("HN", # Half-Normal
                  "HR", # Hazard-Rate
                  "HN_Beaufort", # HN with Beaufort
                  "HR_Beaufort", # HR with Beaufort
                  "HN_size", # HN with group size
                  "HR_size", # HR with group size
                  "Unif_cosine") # Uniform with cosine adjustment

detfc <- purrr::map(.x = dolphin.species,
                    .f = ~{
                      
                      dat <- data.frame(dist.df[[.x]])
                      
                      m1 <- Distance::ds(data = dat, 
                                         truncation = trunc.dist[[.x]],
                                         key = "hn", adjustment = NULL)
                      m2 <- Distance::ds(data = dat, 
                                         truncation = trunc.dist[[.x]],
                                         key = "hr", adjustment = NULL)
                      m3 <- Distance::ds(data = dat, 
                                         truncation = trunc.dist[[.x]],
                                         formula = ~as.factor(beaufort),
                                         key = "hn",
                                         adjustment = NULL)
                      m4 <- Distance::ds(data = dat, 
                                         truncation = trunc.dist[[.x]],
                                         formula = ~as.factor(beaufort),
                                         key = "hr", adjustment = NULL)
                      m5 <- Distance::ds(data = dat, 
                                         truncation = trunc.dist[[.x]],
                                         formula = ~size,
                                         key = "hn", adjustment = NULL)
                      m6 <- Distance::ds(data = dat, 
                                         truncation = trunc.dist[[.x]],
                                         formula = ~size,
                                         key = "hr", adjustment = NULL)
                      
                      ll <- list(m1, m2, m3, m4, m5, m6)
                      names(ll) <- detfc.models[1:6]
                      ll
                    }) %>% purrr::set_names(., dolphin.species)

# Model comparison tables
detfc.comp <- lapply(
  X = dolphin.species,
  FUN = function(p) {
    purrr::map_df(.x = detfc[[p]], .f = ~ AIC(.x)) %>%
      dplyr::mutate(model = detfc.models[1:6]) %>%
      dplyr::mutate(delta_aic = AIC - min(AIC)) %>%
      dplyr::arrange(delta_aic)
  }
) %>%
  purrr::set_names(., dolphin.species)

CvM <- lapply(
  X = dolphin.species,
  FUN = function(p) {
    purrr::map_df(.x = detfc[[p]], .f = ~ Distance::gof_ds(.x)$dsgof$CvM)})

detfc.comp <- purrr::map2(.x = detfc.comp, .y = CvM, 
                          .f = ~dplyr::select(.y, p) %>% 
                            dplyr::bind_cols(.x, .))
row.names(detfc.comp[[1]]) <- NULL
row.names(detfc.comp[[2]]) <- NULL

# Identify best models according to the AIC
detfc.best <- lapply(X = dolphin.species,
                     FUN = function(p){
                       detfc.comp[[p]][1, ]$model}) %>%
  purrr::set_names(., dolphin.species)

# Model checks - AIC is not always a cookbook
plot_detfc(model.rank = c(1, 1))
plot_detfc(model.rank = c(1, 4))

# detfc.best$sousa <- "HR_size"

# Soap film smoother ---------------------------------------------------------------

# Convert spatial objects to data.frame
survey_area <- purrr::map(.x = study_area_coastline, 
                          .f = ~{
                            s_a <- sf::as_Spatial(.x)
                            s_a <- data.frame(s_a@polygons[[1]]@Polygons[[1]]@coords)
                            names(s_a) <- c("x", "y")
                            s_a})

# Create soap film grid - dsm models not fitting if more than 30 knots
soap.knots <- purrr::map(.x = survey_area, .f = ~ dsm::make.soapgrid(bnd = .x, n.grid = c(20, 15)))

# Identify segments inside boundaries
seg.onoff <- purrr::map2(
  .x = seg.df,
  .y = survey_area,
  .f = ~ {
    x <- .x$x
    y <- .x$y
    mgcv::inSide(x = x, y = y, bnd = as.list(.y))
  }
)

# Identify prediction cells inside boundaries
pred.onoff <- purrr::map2(
  .x = pred.df,
  .y = survey_area,
  .f = ~ {
    x <- .x$x
    y <- .x$y
    mgcv::inSide(x = x, y = y, bnd = as.list(.y))
  }
)

# Identify knots inside boundaries
knots.onoff <- purrr::map2(
  .x = soap.knots,
  .y = survey_area,
  .f = ~ {
    x <- .x$x
    y <- .x$y
    mgcv::inSide(x = x, y = y, bnd = as.list(.y))
  }
)

soap.knots <- purrr::map2(.x = soap.knots, .y = knots.onoff, .f = ~.x[.y, ])
seg.df.soap <- purrr::map2(.x = seg.df, .y = seg.onoff, .f = ~.x[.y, ])
pred.df.soap <- purrr::map2(.x = pred.df, .y = pred.onoff, .f = ~.x[.y, ])

# Soap knots need to be stored as separate objects for DSM fitting
knots.2015 <- soap.knots$`2015`
knots.2016_2017 <- soap.knots$`2016_2017`

# Add factor variables for route ID and year
seg.df.soap <- purrr::map(.x = seg.df.soap, .f = ~{
  .x$route_f <- factor(.x$route)
  .x$year_f <- factor(.x$year)
  .x})

# Manually remove knots that are on or outside boundary
text.counter <- NULL

# Fix knots for 2015
while (is.null(text.counter)){
  
  test.dsm <- evaluate::evaluate("dsm(formula = abundance.est ~
                  s(x, y, bs = \"so\", xt = list(bnd = list(survey_area$`2015`))) +
                  s(route_f, bs = \"re\"),
                ddf.obj = detfc[[\"tursiops\"]][[detfc.best[[\"tursiops\"]]]],
                segment.data = data.frame(seg.df.soap$`2015`),
                knots = knots.2015,
                observation.data = data.frame(obs.df$tursiops[obs.df$tursiops$year == 2015, ]),
                method = \"REML\",
                family = tw())")
  
  if(stringr::str_sub(as.character(test.dsm[[2]]), 1, 21) == "Error in crunch.knots"){
    
    error.msg <- as.character(test.dsm[[2]])
    error.msg <- stringr::str_sub(error.msg, gregexpr(":", error.msg)[[1]][1] + 2, nchar(error.msg))
    problem.knot <- readr::parse_number(error.msg)
    knots.2015 <- knots.2015[-problem.knot, ]
    cat("Removing knot:", problem.knot, "\n")
    text.counter <- NULL
    
  } else {
    
    text.counter <- "Stop"
  }
  
} # End while loop

# Fix knots for 2016-2017
text.counter <- NULL

while (is.null(text.counter)){
  
  test.dsm <- evaluate::evaluate("dsm(formula = abundance.est ~
                  s(x, y, bs = \"so\", xt = list(bnd = list(survey_area$`2016_2017`))) +
                  s(route_f, bs = \"re\"),
                ddf.obj = detfc[[\"sousa\"]][[detfc.best[[\"sousa\"]]]],
                segment.data = data.frame(seg.df.soap$`2016_2017`),
                knots = knots.2016_2017,
                observation.data = data.frame(obs.df$sousa[obs.df$sousa$year == 2016, ]),
                method = \"REML\",
                family = tw())")
  
  if(stringr::str_sub(as.character(test.dsm[[2]]), 1, 21) == "Error in crunch.knots"){
    
    error.msg <- as.character(test.dsm[[2]])
    error.msg <- stringr::str_sub(error.msg, gregexpr(":", error.msg)[[1]][1] + 2, nchar(error.msg))
    problem.knot <- readr::parse_number(error.msg)
    knots.2016_2017 <- knots.2016_2017[-problem.knot, ]
    cat("Removing knot:", problem.knot, "\n")
    text.counter <- NULL
    
  } else {
    text.counter <- "Stop"
  }
}

# test.dsm.soap <- dsm(formula = abundance.est ~
#                   s(x, y, bs = "so", k = 30, xt = list(bnd = list(survey_area$`2016_2017`))),
#                 ddf.obj = detfc[["tursiops"]][[detfc.best[["tursiops"]]]],
#                 segment.data = data.frame(seg.df.soap$`2016_2017`),
#                 knots = knots.2016_2017,
#                 observation.data = data.frame(obs.df$tursiops[obs.df$tursiops$year == 2016, ]),
#                 method = "REML",
#                 family = tw())

# test.dsm <- dsm(formula = abundance.est ~
#                   s(x, y, bs = "tp", k = 30), 
#                 ddf.obj = detfc[["tursiops"]][[detfc.best[["tursiops"]]]],
#                 segment.data = data.frame(seg.df.soap$`2016_2017`),
#                 # knots = knots.2016_2017,
#                 observation.data = data.frame(obs.df$tursiops[obs.df$tursiops$year == 2016, ]),
#                 method = "REML",
#                 family = tw())

# preddata.varprop <- split(pred.df.soap$`2016_2017`, 1:nrow(pred.df.soap$`2016_2017`))
# testvar <- dsm.var.gam(test.dsm.soap, pred.data=preddata.varprop, off.set=pred.df.soap$`2016_2017`$area)
# testvar
# 
# # Prediction data.frame
# newdat <- pred.df.soap$`2016_2017`
# preddf <- hexgrid$`2016_2017`[pred.onoff$`2016_2017`,] # As hexgrid is an sf object
# preddf$Nhat <- predict(object = test.dsm.soap, 
#                        newdata = newdat, 
#                        off.set = newdat$area)
# ggplot(data = preddf) + 
#   geom_sf(fill = "white", col = "transparent") +
#   geom_sf(aes(fill = Nhat), col = "transparent") +
#   geom_point(data = sightings$tursiops[sightings$tursiops$year == 2016,], aes(x,y,col=factor(year), size = size)) +
#   geom_sf(data = study_area_coastline$`2016_2017`,
#           fill = "transparent", col = "black", size = 0.25) +
#   scale_fill_gradientn(colours = pals::viridis(1000)) + 
#   theme(legend.text = element_text(colour = "black", size = 14),
#         legend.title = element_blank())

# Checks the soap film
par(mfrow = c(1, 2))
purrr::map2(
  .x = survey_area,
  .y = list(knots.2015, knots.2016_2017),
  .f = ~ soap_check(bnd = list(list(x = .x$x, y = .x$y)), knots = .y, tol = 1e-6))

hexgridpred <- purrr::map2(.x = hexgrid, .y = pred.onoff, .f = ~.x[.y,])


# Spatial models -----------------------------------------------------------

# Define model formulae - use abundance.est here when the detection function model
# contains a observation-level covariate (i.e., group size). The response variable is therefore no
# longer the raw counts per segment but the estimated abundance per segment 
# (using a Horvitz-Thompson-like estimator). Use count otherwise.

model.formulae <- list(
  " ~ s(x, y, bs = \"so\", xt = list(bnd = list(survey_area$`2015`))) +
  s(route_f, bs = \"re\")",
  " ~ year_f + s(x, y, m = 2, bs = \"so\", xt = list(bnd = list(survey_area$`2016_2017`))) +
  s(x, y, by = year_f, bs = \"so\", xt = list(bnd = list(survey_area$`2016_2017`)), m = 1) + s(route_f, bs = \"re\")"
) %>% purrr::set_names(x = ., nm = names(survey.years))


bias.correction <- mrds.df %>% 
  dplyr::group_by(species) %>%
  dplyr::summarise(p0 = mean(p0)) %>%
  dplyr::rename(perc = p0) %>%
  # From Brown et al. and Durden et al. respectively
  dplyr::mutate(avail = c(0.48, mean(c(0.64, 0.7784, 0.8)))) %>% 
  dplyr::mutate(bias = avail * perc)

# Forcada et al. 2004: 0.7784
# Lauriano et al. 2014: 0.8

# Fit the DSMs
dsm.models <- purrr::map(
  .x = setNames(dolphin.species, dolphin.species),
  .f = ~ {
    tmp <- lapply(
      X = names(survey.years),
      FUN = function(yr) {
        species.name <- .x
        mformulae <- as.formula(paste0(ifelse(grepl(pattern = "_", x = detfc.best[[species.name]]),
                                              "abundance.est", "count"), model.formulae[[yr]]))
        if (class(mformulae) == "formula") mformulae <- list(xy = mformulae)
        
        if (species.name == "sousa" & yr == "2015") {
          list()
        } else {
          purrr::map(
            .x = 1:length(mformulae),
            .f = ~ {
              cat(paste0(
                "Fitting DSM for ", species.name,
                " (", yr, ") \n"
              ))
              observ <- data.frame(obs.df[[species.name]])
              if(yr == "2015"){
                dsm::dsm(
                  formula = mformulae[[.x]],
                  ddf.obj = detfc[[species.name]][[detfc.best[[species.name]]]],
                  segment.data = data.frame(seg.df.soap[[yr]]),
                  knots = knots.2015,
                  observation.data = observ[observ$distance <= trunc.dist[[species.name]], ],
                  method = "REML", 
                  availability = bias.correction[bias.correction$species == species.name, ]$bias,
                  family = tw())
              } else {
                dsm::dsm(
                  formula = mformulae[[.x]],
                  ddf.obj = detfc[[species.name]][[detfc.best[[species.name]]]],
                  segment.data = data.frame(seg.df.soap[[yr]]),
                  knots = knots.2016_2017,
                  observation.data = observ[observ$distance <= trunc.dist[[species.name]], ],
                  availability = bias.correction[bias.correction$species == species.name, ]$bias,
                  method = "REML", 
                  family = tw())
              }}
          ) %>% purrr::set_names(., names(mformulae))
        }
      }
    )
    names(tmp) <- names(survey.years)
    tmp
  }
)

# Extract AIC values
# dsm.aic <- purrr::map_depth(.x = dsm.models, .depth = 3, .f = ~AIC(.x)) %>% 
# purrr::map_depth(.x = ., .depth = 2, .f = ~unlist(.x))

# Models that minimise the AIC
# best.dsm <- purrr::map_depth(.x = dsm.aic, .depth = 2, .f = ~names(which.min(.x)))
# best.dsmNames <- unlist(best.dsm) %>% purrr::map2(.x = ., .y = names(.), .f = ~paste0(.y, ".", .x))
# best.dsm <- purrr::cross2(.x = dolphin.species, .y = names(survey.years)) %>% 
#   purrr::map(.x = ., .f = ~{
#     whichmod <- best.dsm[[.[[1]]]][[.[[2]]]]
#     if(!is.null(whichmod)) dsm.models[[.[[1]]]][[.[[2]]]][[whichmod]]}) %>% 
#   purrr::compact() %>% 
#   purrr::set_names(x = ., nm = names(best.dsmNames))

# Perform model checks
par(mfrow = c(1, 3))
purrr::map(.x = dolphin.species,
           .f = ~{
             lapply(X = names(survey.years), 
                    FUN = function(m){
                      mobject <- dsm.models[[.x]][[m]]
                      if(length(mobject) > 0) {dsm_check(mobject[[1]], rq = FALSE)
                        plot(1, type = 'n', axes = FALSE, ann = FALSE)
                        legend(1, 1, legend = paste0(.x, " - ", m))}
                    })})

# Auto-correlation
par(mfrow = c(1, 3))
purrr::map(.x = dolphin.species,
           .f = ~{
             lapply(X = names(survey.years), 
                    FUN = function(m){
                      mobject <- dsm.models[[.x]][[m]]
                      if(length(mobject) > 0) {
                        dsm.cor(mobject[[1]], max.lag = 10, 
                                Segment.Label = "Sample.Label", main = paste0(.x, " - ", m))}
                    })})

par(mfrow = c(1,1))

# Estimates of abundance
N_est <- purrr::cross2(.x = dolphin.species, .y = survey.years) %>% 
  purrr::map(.x = ., .f = ~
               lapply(X = .[[2]], FUN = function(d){
                 ab.est <- get_N(species = .[[1]], yr = d)
                 if(!is.null(ab.est)) names(ab.est) <- paste0(.[[1]], " - ", d)
                 ab.est
               })) %>% unlist(.) %>% 
  tibble::as_tibble(., rownames = "species_year") %>% 
  tidyr::separate(col = species_year, into = c("species", "year"), sep = " - ", convert = TRUE) %>% 
  dplyr::rename(Nhat = value)

# Variance estimation -----------------------------------------------------

# Note that for models where there are covariates at the individual level we cannot calculate the variance via the variance propagation method (dsm.var.prop) of Williams et al (2011). Instead we can use a GAM uncertainty estimation and combine it with the detection function uncertainty via the delta method (dsm.var.gam), which simply sums the squared coefficients of variation to get a total coefficient of variation (and therefore assumes that the detection process and spatial process are independent). There are no restrictions on the form of the detection function when using dsm.var.gam.

# In a nutshell:

# dsm.var.gam can be used in any situation BUT you're ignoring the covariance between the detection function covariates and the spatial model (e.g., that Beaufort varies spatially). That's fine for the case where there's
# no covariates in the detection function.

# dsm.var.prop takes care of the case where the detection function covariates have a covariance with the spatial model (spatial beaufort etc.) but that covariate can only change at the level of the segment (no animal/group level covariates).

# When using estimated abundance as the response, you can only use dsm.var.gam.

# Here, we implement the analytic approach (sometimes called the delta method) described in 
# see eqn 6 of https://link.springer.com/article/10.1007%2Fs13253-021-00438-2#Sec3).
# The other option would be to implement a simulation-based approach (the next section
# in that paper), which can be easier to code-up 
# (see https://github.com/densitymodelling/nefsc_fin_mrds_dsm/blob/master/dsm_fin_analysis.Rmd#L163)

CV_est <- purrr::cross2(.x = dolphin.species, .y = survey.years) %>% 
  purrr::map(.x = ., .f = ~
               lapply(X = .[[2]], FUN = function(d){
                 cv.est <- get_variance(species = .[[1]], yr = d, preddata = pred.df.soap) 
                 if(!is.null(cv.est)) do.call(rbind, cv.est) %>% 
                   tibble::as_tibble(., rownames = "type") %>% 
                   dplyr::rename(CV = V1) %>% 
                   dplyr::mutate(year = d, species = .x[[1]])
               })) %>% purrr::flatten() %>% 
  do.call(rbind, .)

N_est <- dplyr::left_join(x = N_est, y = CV_est %>% dplyr::filter(type == "total"), by = c("year", "species")) %>% dplyr::select(-type)

# Delta method approximate asymptotic CI
N_est <- N_est %>% dplyr::mutate(lower.95 = Nhat / exp(1.96*sqrt(log(1 + CV^2))),
                                 upper.95 = Nhat * exp(1.96*sqrt(log(1 + CV^2))))

# Finally, add %Deviance explained
N_est$Dev <- purrr::map(.x = dsm.models, .f = ~{
  
  # 2015
  if(purrr::is_empty(.x[[1]])){
    s1 <- NULL
  } else {  
    s1 <- summary(.x[[1]][[1]]) 
  }
  s2 <- summary(.x[[2]][[1]]) # 2016_2017
  c(ifelse(is.null(s1), NA, s1$dev.expl), s2$dev.expl, s2$dev.expl)
  
}) %>% do.call(c, .) %>% .[!is.na(.)]

# write.csv(N_est, file = "out/Abundance_estimates.csv", row.names = FALSE)

# Maps ----------

# Set up list object to store predictions
prediction.list <- list(species = dolphin.species, year = unlist(survey.years)) %>% purrr::cross_df()

# Generate model predictions and estimates of variance
prediction.maps <- purrr::map(
  .x = dolphin.species,
  .f = ~ {
    df <- prediction.list %>% dplyr::filter(species == .x)
    dfl <- purrr::map(
      .x = 1:nrow(df),
      .f = ~ get_map(
        species = df$species[.x],
        yr = df$year[.x],
        preddata = pred.df.soap,
        hex.grid = hexgridpred,
        plot.sightings = TRUE
      )
    ) %>%
      purrr::set_names(., nm = df$year)
  }
) %>%
  purrr::set_names(., nm = dolphin.species)

# Percentage of IUCN range ----------------------------------------------------------------

IUCN_sousa <- sf::st_read("gis/IUCN_range/sousa_range.shp")
(sf::st_area(x = study_area$`2015`) / 1000000) / (sf::st_area(x = IUCN_sousa) / 1000000)
(sf::st_area(x = study_area$`2016_2017`) / 1000000) / (sf::st_area(x = IUCN_sousa) / 1000000)
