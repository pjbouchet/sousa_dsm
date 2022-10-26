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

# Import rivers
rivers <- sf::st_read("gis/Four_rivers.shp") %>% 
  sf::st_transform(x = ., crs = crs.utm)

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

d1 <- gg.detfc("sousa", "HR_size", group.size = c(1, 5, 10, 20))
d2 <- gg.detfc("tursiops", "HR_size", group.size = c(1, 5, 10, 20))
fig2 <- (d1 + d2) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold", size = 30))

ggsave(filename = "out/Figure2.pdf", plot = fig2, height = 5, width = 11, dpi = 450)

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


# Export the predictions as raster files
# Create empty raster with desired resolution

purrr::map(.x = dolphin.species,
           .f = ~{
             
      if(.x == "tursiops"){
      r <- sf::st_drop_geometry(prediction.maps[[.x]]$`2015`$plot$Nhat$data[, c("x", "y", "Nhat")])
      ras <- raster::raster(raster::extent(hexgrid$`2015`), res = spatial.scale)
      raster::crs(ras) <- crs.utm
      outrast <- raster::rasterize(r[,1:2], y = ras, field = r$Nhat, fun = mean_ras)
      raster::writeRaster(x = outrast, filename = "out/tursiops_Nhat_2015", format = "GTiff")
      }
             
      ras <- raster::raster(raster::extent(hexgrid$`2016_2017`), res = spatial.scale)
      raster::crs(ras) <- crs.utm
      r1 <- sf::st_drop_geometry(prediction.maps[[.x]]$`2016`$plot$Nhat$data[, c("x", "y", "Nhat")])
      r2 <- sf::st_drop_geometry(prediction.maps[[.x]]$`2017`$plot$Nhat$data[, c("x", "y", "Nhat")])
      
      outrast1 <- raster::rasterize(r1[,1:2], y = ras, field = r1$Nhat, fun = mean_ras)
      outrast2 <- raster::rasterize(r2[,1:2], y = ras, field = r2$Nhat, fun = mean_ras)
      
    raster::writeRaster(x = outrast1, filename = paste0("out/", .x, "_Nhat_2016"), format = "GTiff")
    raster::writeRaster(x = outrast2, filename = paste0("out/", .x, "_Nhat_2017"), format = "GTiff") })

purrr::map(.x = dolphin.species,
                 .f = ~{
                   
      if(.x == "tursiops"){
      r <- sf::st_drop_geometry(prediction.maps[[.x]]$`2015`$plot$CVhat$data[, c("x", "y", "CVhat")])
      ras <- raster::raster(raster::extent(hexgrid$`2015`), res = spatial.scale)
      raster::crs(ras) <- crs.utm
      outrast <- raster::rasterize(r[,1:2], y = ras, field = r$CVhat, fun = mean_ras)
      raster::writeRaster(x = outrast, filename = "out/tursiops_CV_2015", format = "GTiff")
      }
                   
      ras <- raster::raster(raster::extent(hexgrid$`2016_2017`), res = spatial.scale)
      raster::crs(ras) <- crs.utm
      r1 <- sf::st_drop_geometry(prediction.maps[[.x]]$`2016`$plot$CVhat$data[, c("x", "y", "CVhat")])
      r2 <- sf::st_drop_geometry(prediction.maps[[.x]]$`2017`$plot$CVhat$data[, c("x", "y", "CVhat")])
                   
      outrast1 <- raster::rasterize(r1[,1:2], y = ras, field = r1$CVhat, fun = mean_ras)
      outrast2 <- raster::rasterize(r2[,1:2], y = ras, field = r2$CVhat, fun = mean_ras)
                   
    raster::writeRaster(x = outrast1, filename = paste0("out/", .x, "_CV_2016"), format = "GTiff")
    raster::writeRaster(x = outrast2, filename = paste0("out/", .x, "_CV_2017"), format = "GTiff") })



# Results ----------------------------------------------------------------

# Total No. sightings by species and year with and without truncation

purrr::map(.x = 1:2, .f = ~dist.df[[.x]] |> 
             dplyr::mutate(species = dolphin.species[.x])) %>%
  do.call(rbind, .) |> 
  dplyr::mutate(grpsize = ifelse(size <=4, "1–4", "5+")) |> 
  dplyr::group_by(species, year, grpsize) |> 
  dplyr::summarise(n()) |> 
  dplyr::arrange(desc(species))

purrr::map(.x = 1:2, .f = ~dist.df[[.x]] |> 
             dplyr::mutate(species = dolphin.species[.x]) |> 
             dplyr::filter(distance<= trunc.dist[[.x]])) %>%
  do.call(rbind, .) |> 
  dplyr::mutate(grpsize = ifelse(size <=4, "1–4", "5+")) |> 
  dplyr::group_by(species, year, grpsize) |> 
  dplyr::summarise(n()) |> 
  dplyr::arrange(desc(species))

# Effort

seg.df %>% do.call(rbind, .) |> 
  dplyr::group_by(year) |> 
  dplyr::summarise(n(), tot = sum(Effort))

purrr::map(.x = 1:2, .f = ~obs.df[[.x]] |> dplyr::mutate(species = dolphin.species[.x])) %>%
  do.call(rbind, .) |> 
  dplyr::group_by(species, year) |> 
  dplyr::summarise(n = length(unique(Sample.Label))) |> 
  dplyr::arrange(desc(species))

p.tursiops <- detfc$tursiops$HN_size$ddf$fitted |> 
  tibble::enframe() |> 
  dplyr::rename(object = name, p = value)

p.sousa <- detfc$sousa$HN_size$ddf$fitted |> 
  tibble::enframe() |> 
  dplyr::rename(object = name, p = value)

obs.df_trunc <- purrr::map2(.x = obs.df, .y = trunc.dist, .f = ~dplyr::filter(.x, distance <= .y))

purrr::map2(.x = obs.df_trunc, .y = list(p.tursiops, p.sousa),
                      .f = ~dplyr::left_join(.x, .y, by = "object")) %>%
  purrr::map(.x = ., .f = ~dplyr::mutate(.x, N_hat = size / p)) %>%
  purrr::map2(.x = ., .y = dolphin.species, .f = ~dplyr::mutate(.x, species = .y)) %>%
  do.call(rbind, .) |> 
  dplyr::group_by(species, year) |> 
  dplyr::summarise(N_max = max(N_hat)) |> 
  dplyr::arrange(desc(species))
  
purrr::map(.x = unlist(survey.years), .f = ~{
  tmp <- prediction.maps$tursiops[[as.character(.x)]]$plot$Nhat$data$Nhat/(prediction.maps$tursiops[[as.character(.x)]]$plot$Nhat$data$area/1000000)
  mean(tmp)
  # sd(tmp)
})

purrr::map(.x = unlist(survey.years), .f = ~{
  tmp <- prediction.maps$sousa[[as.character(.x)]]$plot$Nhat$data$Nhat/(prediction.maps$tursiops[[as.character(.x)]]$plot$Nhat$data$area/1000000)
  mean(tmp)
  sd(tmp)
})
