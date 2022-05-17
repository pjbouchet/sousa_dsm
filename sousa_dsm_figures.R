##%######################################################%##
#                                                          #
#'###     Contrasting aerial abundance estimates for     ####
#'###          two sympatric dolphin species at          ####
#'###          a regional scale using distance           ####
#'###       sampling and density surface modelling       ####
#                                                          #
##%######################################################%##

# Authors: Raudino et al. (2022)

# Figure 1 ----------------------------------------------------------------

n.insets <- 6

# Convert sightings df to simple feature for plotting
sightings.sf <- purrr::map2(.x = sightings, 
                            .y = dolphin.species,
                            .f = ~dplyr::mutate(.x, species = .y)) %>% 
  do.call(rbind, .) %>% 
  sp::SpatialPointsDataFrame(coords = cbind(.$x, .$y), 
                             data = ., proj4string = crs.utm) %>% 
  sf::st_as_sf(.) %>% 
  dplyr::mutate(species = ifelse(species == "sousa", "Sousa saluhensis", "Tursiops aduncus"))

# Import transects common to both years
transects.overlap <- sf::read_sf("gis/2015_2016_transects.shp")

# Insets
insets <- purrr::map(.x = LETTERS[1:n.insets], 
                     .f = ~{
                       sf::read_sf(paste0("gis/inset_", .x, ".shp")) %>% 
                         sf::st_transform(x = ., crs = crs.utm)})
names(insets) <- LETTERS[1:n.insets]

inset.coords <- sf::read_sf("gis/inset_coordinates.shp") %>% 
  sf::st_transform(x = ., crs = crs.utm) %>% 
  sf::as_Spatial()
inset.coords <- cbind(inset.coords$inset_ID, coordinates(inset.coords))

# This function outputs all the relevant maps to the "out" folder
# These are then compiled in Affinity Designer (external programme)
make_fig1()

# Figure 2 ----------------------------------------------------------------

d1 <- gg.detfc("sousa", "HR_size", group.size = c(1, 5, 20))
d2 <- gg.detfc("tursiops", "HR_size", group.size = c(1, 5, 20))
fig2 <- (d1 + d2) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold", size = 30))
ggsave(filename = "out/fig2.pdf", plot = fig2, height = 5, width = 11, dpi = 450)

# Figure 3 ----------------------------------------------------------------

# Save prediction maps to disk
purrr::walk(
  .x = 1:nrow(prediction.list),
  .f = ~ {
    mn <- c("Nhat", "CVhat")
    for (j in mn) {
      fn <- paste0(
        "out/map_", prediction.list$species[.x], "_",
        prediction.list$year[.x], "_", j, ".pdf"
      )
      ln <- paste0(
        "out/legend_", prediction.list$species[.x], "_",
        prediction.list$year[.x], "_", j, ".pdf"
      )
      fm <- prediction.maps[[prediction.list$species[.x]]][[as.character(prediction.list$year[.x])]]$plot[[j]]
      lg <- prediction.maps[[prediction.list$species[.x]]][[as.character(prediction.list$year[.x])]]$lgd[[j]]
      if (!is.null(fm)) ggsave(filename = fn, plot = fm, height = 5, width = 11, dpi = 450)
      if (!is.null(lg)) ggsave(filename = ln, plot = lg, dpi = 450)
    }
  }
)
