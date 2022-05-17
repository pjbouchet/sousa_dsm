##%######################################################%##
#                                                          #
#'###     Contrasting aerial abundance estimates for     ####
#'###          two sympatric dolphin species at          ####
#'###          a regional scale using distance           ####
#'###       sampling and density surface modelling       ####
#                                                          #
##%######################################################%##

# Authors: Raudino et al. (2022)

##%######################################################%##
#                                                          #
#'###           Wrapper function to get basic            ####
#'###         plots of detection function & GOF          ####
#                                                          #
##%######################################################%##

#' @param model.rank Ranks of the models to be plotted for tursiops and sousa respectively

plot_detfc <- function(model.rank = c(1, 1)){
  
  par(mfrow = c(2, 2))
  purrr::walk(.x = 1:length(dolphin.species), 
              .f = ~{
                model.name <- detfc.comp[[.x]]$model[model.rank[.x]]
                species.name <- dolphin.species[.x]
                plot(detfc[[.x]][[model.name]], 
                     showpoints = FALSE, pl.den = 0, lwd = 2,
                     main = paste0(species.name, " (", model.name, ")"))
                ddf.gof(detfc[[.x]][[model.name]]$ddf)
              }) 
}

##%######################################################%##
#                                                          #
#'###             Function to plot detection             ####
#'###         function for different group sizes         ####
#                                                          #
##%######################################################%##

#' @param species Species whose detection function is to be plotted
#' @param ds.model Detection function model. Can be set to "best" to grab the model with lowest AIC.
#' @param group.size Group sizes for which to plot separate detection functions. Defaults to NULL if size is not included in the model. 
#' @param col.by.species # Use manuscript colours for sousa and tursiops
#' @param col.p # Colour of the average detection curve
#' @param col.hist # Fill colour for histogram bars

gg.detfc <- function(species,
                     ds.model,
                     group.size = NULL,
                     col.by.species = TRUE,
                     col.p = "black",
                     col.hist = "grey70"){
  
  # Hazard rate key function
  key.fct.hz <- function(distance, key.scale, key.shape) {
    return(1 - exp(-(distance / key.scale)^(-key.shape)))
  }
  
  # Half normal key function
  key.fct.hn <- function (distance, key.scale) 
  {exp(-((distance/(sqrt(2) * key.scale))^2))}
  
  # Scaling factor needed for correct plotting (from mrds:::scalevalue)
  scale.value <- function(key.scale, z) {
    exp(as.matrix(z) %*% key.scale)
  }
  
  # Define required parameters
  breaks <- NULL
  vname <- "distance"
  point <- TRUE
  pdf <- FALSE
  
  # Extract model and data  
  if(ds.model == "best") ds.model <- detfc.best[[species]]
  model <- detfc[[species]][[ds.model]]
  hazrate <- grepl(pattern = "HR", x = ds.model)
  if(!grepl(pattern = "size", x = ds.model)) group.size <- NULL
  
  dist.data <- dist.df[[species]]
  ddfobj <- model$ddf$ds$aux$ddfobj
  dat <- model$ddf$data
  
  # Distance bounds
  width <- model$ddf$meta.data$width
  left <- model$ddf$meta.data$left
  
  # Vector of distances used to re-create the detection function
  distances <- seq(0, width, length.out = 100)

  if(!is.null(group.size)){
    if(hazrate){
     
      # Shape parameter
      key.shape <- scale.value(key.scale = ddfobj$shape$parameters, z = matrix(1, nrow = 100, 1))
     
      # Scale parameter
      mat.hr <- purrr::map(.x = group.size,
                           .f = ~{
                             mat.tmp <- matrix(1, nrow = 100, ncol = 2)
                             # mat.tmp[, 2] <- log(.x) # when log(size) included in detection function
                             mat.tmp[, 2] <- .x
                             scale.value(key.scale = ddfobj$scale$parameters, z = mat.tmp)})
    } else{
      mat.hr <- purrr::map(.x = group.size,
                           .f = ~{
                             mat.tmp <- matrix(1, nrow = 100, ncol = 2)
                             mat.tmp[, 2] <- .x
                             scale.value(key.scale = ddfobj$scale$parameters, z = mat.tmp)})
    }
  }

  # Detection probability values and abundance estimate
  selected <- rep(TRUE, nrow(ddfobj$xmat))
  xmat <- ddfobj$xmat[selected, ]
  
  # Detection probability for each fitted value & Nhat estimate
  if (length(model$ddf$fitted) == 1) {
    pdot <- rep(model$ddf$fitted, sum(as.numeric(selected)))
  } else {
    pdot <- model$ddf$fitted[selected]
    Nhat <- sum(1 / pdot)
  }
  
  # Histogram breaks
  n <- length(xmat$distance)
  nc <- round(sqrt(n), 0)
  binned <- FALSE
  if (is.null(model$ds$aux$int.range)) {
    int.range <- c(0, trunc.dist[[species]])
  } else {
    int.range <- model$ds$aux$int.range
  }
  if (is.matrix(int.range)) {
    max.range <- as.vector(int.range[1, ])
    int.range <- int.range[2:dim(int.range)[1], ]
    range.varies <- TRUE
  } else {
    max.range <- int.range
    normalize <- FALSE
    range.varies <- FALSE
  }
  breaks <- c(max(0, (max.range[1])), max.range[1] + ((max.range[2] - max.range[1])/nc) * (1:nc))

  if (breaks[1] > left) {
    breaks <- c(left, breaks)
    nc <- nc + 1
  }
  
  nc <- length(breaks) - 1
  lower <- min(breaks)
  upper <- max(breaks)
  dat <- dat[selected, ]
  keep <- dat[, vname] >= lower & dat[, vname] <= upper
  hist.obj <- hist(dat[, vname][keep], breaks = breaks, plot = FALSE)
  ymax <- max(hist.obj$counts)
  
  # Histogram values
  expected.counts <- (breaks[2:(nc + 1)] - breaks[1:nc]) * (Nhat/breaks[nc + 1])
  if (!(pdf & point)) {
    hist.obj$density <- hist.obj$counts/expected.counts
    hist.obj$density[expected.counts == 0] <- 0
  }
  
  point_vals <- mrds:::detfct(xmat$distance, ddfobj, select = selected, width = width)
  if (is.null(ylim)) ylim <- c(0, max(hist.obj$density, max(point_vals)))

  hist.data <- data.frame(xmin = breaks[1:length(breaks)-1],
                          xmax = breaks[2:length(breaks)],
                          ymin = 0,
                          ymax = hist.obj$density)

  # Compute average detection probability function
  finebr <- seq(0, width, length.out = 101)
  xgrid <- NULL
  linevalues <- NULL
  newdat <- xmat
  for (i in 1:(length(finebr) - 1)) {
    x <- (finebr[i] + finebr[i + 1])/2
    xgrid <- c(xgrid, x)
    newdat$distance <- rep(x, nrow(newdat))
    detfct.values <- detfct(newdat$distance, ddfobj, select = selected, width = width)
    linevalues <- c(linevalues, sum(detfct.values/pdot)/sum(1/pdot))
  }
  y.avg <- linevalues
  
  # Detection function for each group size
  if(!is.null(group.size)){
    if(hazrate){
      y.val <- purrr::map(.x = 1:length(group.size), 
                          .f = ~{
                            key.fct.hz(distances, key.scale = mat.hr[[.x]], key.shape = key.shape) %>% 
                              as.data.frame(.) %>% 
                              dplyr::mutate(Size = group.size[.x]) %>% 
                              dplyr::rename(value = V1)})
    } else {
      y.val <- purrr::map(.x = 1:length(group.size), 
                          .f = ~{
                            key.fct.hn(distances, key.scale = mat.hr[[.x]]) %>% 
                              as.data.frame(.) %>% 
                              dplyr::mutate(Size = group.size[.x]) %>% 
                              dplyr::rename(value = V1)}) 
    }
    
    y.val.df <- do.call(rbind, y.val) %>% 
      dplyr::mutate(distce = rep(distances, length(group.size))) %>% 
      as_tibble(.)
  } else {
    y.val.df <- data.frame(distce = distances, value = y.avg) 
  }
  
  if (!col.by.species) {
  if (!is.null(group.size) & length(group.size) < 5) stop("Must specify 5 or more size classes when col.by.species = FALSE")
  line.cols <- viridisLite::viridis(n = length(group.size), begin = 0.2)
} else {
  if (species == "sousa") species.colour <- "#E7B800" else species.colour <- "#65ABB7"
  line.cols <- purrr::map(
    .x = seq(0.3, 1, length.out = length(group.size)),
    .f = ~ hexa2hex(input.colour = species.colour, opacity = .x)
  ) %>%
    do.call(c, .)
}

gPlot <- ggplot(data = hist.data) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = col.hist, color = "black") +
  {
    if (!is.null(group.size)) {
      geom_line(
        data = y.val.df, size = 0.7,
        aes(x = distce, y = value, group = factor(Size), colour = factor(Size))
      )
    }
  } +
  {
    if (is.null(group.size)) {
      geom_line(
        data = y.val.df, size = 0.7, aes(x = distce, y = value)
      )
    }
  } +
  geom_line(
    data = data.frame(distce = distances, value = y.avg),
    aes(x = distce, y = value), colour = col.p, size = 0.7
  ) +
  xlab("Distance (m)") +
  ylab("Detection probability") +
  labs(colour = "Group size") +
  scale_colour_manual(values = line.cols) +
  theme(
    plot.margin = unit(c(0.5, 1.5, 0.5, 1), "cm"),
    axis.text.x = element_text(size = 11, colour = "black"),
    axis.text.y = element_text(size = 11, angle = 90, hjust = 0.5, vjust = 0.5, colour = "black"),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 14),
    axis.title.x = element_text(margin = margin(t = 20, r = 00, b = 0, l = 0), size = 14),
    legend.text = element_text(size = 11),
    legend.position = "top"
  )

return(gPlot)
  
}

##%######################################################%##
#                                                          #
#'###             Wrapper function to create             ####
#'###            diagnostics plots for a DSM             ####
#                                                          #
##%######################################################%##

#' @param model Input DSM
#' @param rq Logical, defaults to TRUE. Whether to create plots using randomized quantile residuals.
#' @param gg.plot Logical, use ggplot for graphics if TRUE.

dsm_check <- function(model, rq = TRUE, gg.plot = FALSE){
  
  print(summary(model))
  par(mfrow = c(3,2))
  if(rq) dsm::rqgam.check(model) else gam.check(model)
  
  environment(model$family$variance)$p <- model$family$getTheta(TRUE)
  qres <- statmod::qres.tweedie(model)
  
  if(gg.plot){
    ggplot(aes(x, y),
           data = seg.df.soap) + 
      geom_point(aes(color=qres))+
      facet_wrap(~as.factor(year), nrow = 3, ncol = 1)+
      scale_color_viridis_c()+
      theme_bw()}
  
}

##%######################################################%##
#                                                          #
#'###      Function to retrieve abundance estimates      ####
#                                                          #
##%######################################################%##

#' @param species Species of interest
#' @param yr Survey year

get_N <- function(species, yr){
  
  which.year <- purrr::map_lgl(.x = survey.years, .f = ~any(.x == yr))
  
  if(sum(which.year) == 0){
    
    res <- NULL 
    
  } else {
    
    if(length(dsm.models[[species]][[names(survey.years)[which.year]]]) > 0){
      
      route_id_tmp <- seg.df[[names(survey.years)[which.year]]] %>% 
        dplyr::filter(year == yr) %>% 
        dplyr::pull(route) %>% 
        unique() %>% 
        sample(., size = 1)
      
      newdat <- pred.df[[names(survey.years)[which.year]]] %>% 
        dplyr::mutate(year_f = as.character(yr), route_f = route_id_tmp) %>% 
        data.frame(.)
      
      res <- predict(object = dsm.models[[species]][[names(survey.years)[which.year]]][[1]], 
                     newdata = newdat, off.set = newdat$area,
                     exclude = "s(route_f)") %>% sum(na.rm = TRUE) 
      
    } else { res <- NULL }
  }
  return(res)
}

##%######################################################%##
#                                                          #
#'###         Function to map model predictions          ####
#                                                          #
##%######################################################%##

#' @param species Species name
#' @param yr Survey year
#' @param preddata Prediction grid object
#' @param hex.grid Hexagonal grid object
#' @param plot.sightings Logical. Whether to overlay sightings on the map.
#' @param col.sightings Colour to use to represent sightings.
#' @param show.gridlines Logical. Whether to show grid lines.
#' @param show.labels Logical. Whether to show coordinate labels.

get_map <- function(species, 
                    yr, 
                    preddata,
                    hex.grid,
                    plot.sightings = TRUE, 
                    col.sightings = "black",
                    show.gridlines = FALSE,
                    show.labels = FALSE){
  
  which.year <- purrr::map_lgl(.x = survey.years, .f = ~any(.x == yr))
  
  if(length(dsm.models[[species]][[names(survey.years)[which.year]]]) > 0){
  
    # Extract DSM model
    dsm.model <- dsm.models[[species]][[names(survey.years)[which.year]]][[1]]
    
    # Random route
    route_id_tmp <- seg.df[[names(survey.years)[which.year]]] %>% 
      dplyr::filter(year == yr) %>% 
      dplyr::pull(route) %>% 
      unique() %>% 
      sample(., size = 1)
    
    # Prediction data.frame
    newdat <- preddata[[names(survey.years)[which.year]]] %>% 
      dplyr::mutate(year_f = as.character(yr), route_f = route_id_tmp) %>% 
      data.frame(.)
  
    preddf <- hex.grid[[names(survey.years)[which.year]]] # As hexgrid is an sf object
    preddf$Nhat <- predict(object = dsm.model, 
                        newdata = newdat, 
                        off.set = newdat$area,
                        exclude = "s(route_f)")
  
    CVhat <- get_variance(species = species, yr = yr, grid = TRUE, preddata = preddata)
    preddf$CVhat <- as.numeric(CVhat$total)
    preddf$transp <- 1 - scales::rescale(preddf$CVhat, c(0.1, 1))
    
    preddf <- preddf %>% 
      dplyr::mutate(CVhat_bin = cut(CVhat, breaks = c(0,0.2,0.4,0.6,0.8,1))) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(CVhat_bin = ifelse(CVhat > 1, ">1", as.character(CVhat_bin))) %>%
      dplyr::ungroup()
    
    if(plot.sightings) sight <- sightings[[species]] %>% 
      dplyr::filter(year == yr)
    
    # Generate map
    gPlot <- list()
    
    # gPlot$biscale <- ggplot(data = preddf) + 
    #   geom_sf(fill = "white", col = "transparent") +
    #   geom_sf(aes(fill = Nhat, alpha = transp), col = "transparent") +
    #   geom_sf(data = study_area_coastline[[names(survey.years)[which.year]]],
    #           fill = "transparent", col = "black", size = 0.25) +
    #   scale_fill_gradientn(colours = pals::viridis(1000)) +
    #   theme(legend.text = element_text(colour = "black", size = 14),
    #         legend.title = element_blank())
    
    gPlot$Nhat <- ggplot(data = preddf) + 
      geom_sf(fill = "white", col = "transparent") +
      geom_sf(aes(fill = Nhat), col = "transparent") +
      geom_sf(data = study_area_coastline[[names(survey.years)[which.year]]],
              fill = "transparent", col = "black", size = 0.25) +
      scale_fill_gradientn(colours = pals::viridis(1000)) + 
      theme(legend.text = element_text(colour = "black", size = 14),
            legend.title = element_blank())
    
    gPlot$CVhat <- ggplot(data = preddf) + 
      geom_sf(fill = "white", col = "transparent") +
      geom_sf(aes(fill = as.factor(CVhat_bin)), col = "transparent") +
      geom_sf(data = study_area_coastline[[names(survey.years)[which.year]]],
              fill = "transparent", col = "black", size = 0.25) +
      # scale_fill_gradientn(colours = pals::brewer.ylorbr(1000)) +
      scale_fill_manual(values = pals::brewer.ylorbr(6))
      theme(legend.text = element_text(colour = "black", size = 14),
            legend.title = element_blank())
    
    glegends <- purrr::map(.x = gPlot, .f = ~ggpubr::get_legend(.x) %>% as_ggplot(.))
    
    gPlot <- purrr::map(.x = gPlot, 
                        .f = ~{
         .x +
      {if(plot.sightings) geom_point(data = sight, aes(x = x, y = y, size = size), colour = col.sightings, 
                                     shape = 1)} +
      coord_sf(xlim = c(sf::st_bbox(study_area[[names(survey.years)[which.year]]])[1],
                        sf::st_bbox(study_area[[names(survey.years)[which.year]]])[3]),
               ylim = c(sf::st_bbox(study_area[[names(survey.years)[which.year]]])[2],
                        sf::st_bbox(study_area[[names(survey.years)[which.year]]])[4]),
               expand = FALSE) +
      
      scale_size_continuous(limits = c(1, 35), name = "Group size",
                            breaks = c(1, 2, 5, 10, 15, 20, 30), range = c(1, 10),
                            guide = guide_legend(override.aes = list(linetype = rep("blank", 7), 
                                                                     shape = rep(1, 7), fill = "black"))) +
      theme_minimal() + 
      xlab("") + ylab("") +
      theme(legend.position = "none",
            axis.text.x = element_text(colour = "black", size = 12),
            axis.text.y = element_text(colour = "black", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
            panel.grid.major = element_line(size = 0, linetype = 'solid', colour = "#93bac4")) + 
                            
     {if(show.gridlines) theme(panel.grid.major = element_line(size = 0.1, 
                                                               linetype = 'solid', colour = "#93bac4"))} +
                            {if(!show.labels) theme(axis.text.x = element_blank(),
                             axis.text.y = element_blank())}
                        })
            
  } else {
    gPlot <- NULL
    glegends <- NULL
  } # End if
  return(list(plot = gPlot, lgd = glegends)) 
}

##%######################################################%##
#                                                          #
####        Function to get variance estimates           ####
#                                                          #
##%######################################################%##

#' @param species Species
#' @param yr Survey year
#' @param grid Logical. If TRUE, computes CV for each cell of a prediction grid.
#' @param preddata Prediction grid.

get_variance <- function(species, yr, grid = FALSE, preddata){
  
  which.year <- purrr::map_lgl(.x = survey.years, .f = ~any(.x == yr))
  
  if(length(dsm.models[[species]][[names(survey.years)[which.year]]]) > 0){
    
    # CV of detectability component
    detfc.summary <- summary(detfc[[species]][[detfc.best[[species]]]])
    cv.detfc <- as.numeric(detfc.summary$ds$average.p.se / detfc.summary$ds$average.p)
    
    # Extract DSM model
    dsm.model <- dsm.models[[species]][[names(survey.years)[which.year]]][[1]]
    
    # Random route
    route_id_tmp <- seg.df[[names(survey.years)[which.year]]] %>% 
      dplyr::filter(year == yr) %>% 
      dplyr::pull(route) %>% 
      unique() %>% 
      sample(., size = 1)
    
    # Prediction data.frame
    newdat <- preddata[[names(survey.years)[which.year]]] %>% 
      dplyr::mutate(year_f = as.character(yr), route_f = route_id_tmp) %>% 
      data.frame(.)
    
    if(!grid){
    
    # CV for the spatial component
    cv.gam <- get_cv_gam(dsm.object = dsm.model, pred.data = newdat)
    overall.cv <- sqrt(cv.detfc^2 + cv.gam^2)
    
    } else {
      
    newgrid <- split(newdat, 1:nrow(newdat))
    pb <- progress::progress_bar$new(format = " Computing [:bar] :percent eta: :eta ",
                                     total = length(newgrid), clear = FALSE, width = 80)
    cv.gam <- purrr::map_dbl(.x = newgrid, .f = ~{
      pb$tick(); get_cv_gam(dsm.object = dsm.model, pred.data = .x)}) 
    overall.cv <- purrr::map_dbl(.x = cv.gam, .f = ~sqrt(cv.detfc^2 + .x^2))
    }

    
    res <- list(det = cv.detfc, gam = cv.gam, total = overall.cv)
    
    
  } else {
    
    res <-  NULL
    
  }
  return(res)
}

##%######################################################%##
#                                                          #
####        Function to compute the variance of         ####
####           the spatial component of a DSM           ####
#                                                          #
##%######################################################%##

#' @param dsm.object DSM model.
#' @param pred.data Prediction data.frame.

get_cv_gam <- function(dsm.object, pred.data){
  
  # Inverse link functions
  tmfn <- dsm.object$family$linkinv
  dtmfn <- function(eta){ifelse(is.na(eta), NA, grad(tmfn, ifelse(is.na(eta), 0, eta)))}
  
  # from mvbutils
  "%**%"<-function(x, y){
    dimnames(x) <- NULL
    dimnames(y) <- NULL
    
    if(length(dim(x)) == 2 && length(dim(y)) == 2 && dim(x)[2] ==
       1 && dim(y)[1] == 1){
      return(c(x) %o% c(y))}
    
    if((!is.null(dim(x)) && any(dim(x) == 1))){
      dim(x) <- NULL}
    
    if((!is.null(dim(y)) && any(dim(y) == 1))){
      dim(y) <- NULL}
    
    if(is.null(dim(x)) && is.null(dim(y))){
      if(length(x) == length(y)){
        x <- x %*% y
      }else{
        if ((length(x) != 1) && (length(y) != 1)){
          stop(paste("lengths of x (",length(x),") and y (",
                     length(y),") are incompatible",sep=""))
        }else{
          x <- x * y
        }
      }
    }else{
      x <- x %*% y
    }
    
    if((!is.null(dim(x)) && any(dim(x) == 1))){
      dim(x) <- NULL
    }
    return(x)
  }
  
  cft <- coef(dsm.object) # Extract model coefficients
  cft[which(startsWith(x = names(cft), prefix = "s(route_f)"))] <- 0 # Zero out random effects
  off.set <- pred.data$area # Extract offset
  
  # Convert to list objects
  if (is.data.frame(pred.data) & is.vector(off.set)) {
    pred.data <- list(pred.data)
    off.set <- list(off.set)
  }
  
  # Set up prediction matrices
  preddo <- list(length(pred.data))
  dpred.db <- matrix(0, length(pred.data), length(cft))
  
  for (ipg in seq_along(pred.data)) {
    
    # Initialise offset vector
    pred.data[[ipg]]$off.set <- rep(0, nrow(pred.data[[ipg]]))
    
    # Make predictions - lpmatrix yields the values of the linear predictor (minus any offset) 
    # when postmultiplied by the parameter vector
    lpmat <- predict(dsm.object, newdata = pred.data[[ipg]], type = "lpmatrix", exclude = "s(route_f)")
    lppred <- lpmat %**% cft # Linear predictor
    
    if (length(off.set[[ipg]]) == 1) {
      this.off.set <- rep(off.set[[ipg]], nrow(pred.data[[ipg]]))
    } else {
      this.off.set <- off.set[[ipg]]
    }
    # Apply inverse functions
    preddo[[ipg]] <- this.off.set %**% tmfn(lppred)
    dpred.db[ipg, ] <- this.off.set %**% (dtmfn(lppred) * lpmat)
  }
  
  vpred <- dpred.db %**% tcrossprod(vcov(dsm.object), dpred.db)
  cv.gam <- sqrt(vpred) / sum(unlist(preddo))
  return(cv.gam)

}

##%######################################################%##
#                                                          #
#'###           Function to generate Figure 1            ####
#                                                          #
##%######################################################%##

make_fig1 <- function(){
  
  # Cities
  cities <- rbind(tibble::tibble(city = "Exmouth", lat = -21.9313, lon = 114.1237)) %>% 
    rbind(., tibble::tibble(city = "Onslow", lat = -21.6627, lon = 115.1070)) %>% 
    rbind(., tibble::tibble(city = "Port Hedland", lat = -20.3110, lon = 118.6095)) %>% 
    rbind(., tibble::tibble(city = "Karratha", lat = -20.7337, lon = 116.8447)) %>% 
    dplyr::arrange(city) %>% 
    sp::SpatialPointsDataFrame(coords = cbind(.$lon, .$lat), data = ., proj4string = crs.latlon) %>% 
    sp::spTransform(x = ., CRSobj = crs.utm) %>% 
    as_tibble() %>% 
    dplyr::mutate(coords.x1.label = coords.x1, coords.x2.label = coords.x2)
  
  cities[cities$city == "Exmouth", ]$coords.x1.label <- cities[cities$city == "Exmouth", ]$coords.x1 + 25000
  cities[cities$city == "Exmouth", ]$coords.x2.label <- cities[cities$city == "Exmouth", ]$coords.x2 + 4500
  cities[cities$city == "Onslow", ]$coords.x1.label <- cities[cities$city == "Onslow", ]$coords.x1 + 16000
  cities[cities$city == "Onslow", ]$coords.x2.label <- cities[cities$city == "Onslow", ]$coords.x2 - 13000
  cities[cities$city == "Karratha", ]$coords.x1.label <- cities[cities$city == "Karratha", ]$coords.x1 - 10500
  cities[cities$city == "Karratha", ]$coords.x2.label <- cities[cities$city == "Karratha", ]$coords.x2 - 19000
  cities[cities$city == "Port Hedland", ]$coords.x1.label <- 
    cities[cities$city == "Port Hedland", ]$coords.x1 + 28000
  cities[cities$city == "Port Hedland", ]$coords.x2.label <- 
    cities[cities$city == "Port Hedland", ]$coords.x2 - 13000
  
  # Combine all transects and add year
  all.transects <- purrr::map2(.x = transects, 
                               .y = names(transects), 
                               .f = ~dplyr::mutate(.x, year = .y)) %>% 
    do.call(rbind, .)
  
  all.transects$year <- ifelse(all.transects$year == "2016_2017", "2016 & 2017", all.transects$year)
  
  fig1 <- ggplot() +
    geom_sf(data = all.transects, aes(colour = year), fill = "transparent", lwd = 0.4, show.legend = "line") +
    geom_sf(data = transects.overlap, linetype = "dashed") +
    geom_sf(data = wa, col = "#232323", fill = "#e1d4b2", lwd = 0.25) +
    geom_point(data = cities, aes(x = coords.x1, y = coords.x2)) + 
    geom_text(data = cities, aes(x = coords.x1.label, y = coords.x2.label, label = city), size = 4)
  
  for(layer in 1:length(insets)){
    fig1 <- fig1 + 
      geom_sf(data = insets[[layer]], fill = "transparent", col = "black") + 
      annotate("text", x = as.numeric(inset.coords[layer, 2]) + 4000,
               y = as.numeric(inset.coords[layer, 3]) + 10000, 
               label = LETTERS[layer], fontface = 2, size = 6)}
  
  fig1 <- fig1 +
    
    annotate("text", x = 245268, y = 7763498, label = "Indian Ocean", 
             fontface = 'italic', colour = "#3d7899", angle = 25, size = 5) +
    annotate("text", x = 214461, y = 7646577, label = "North West \n Cape", 
             fontface = 'italic', colour = "black", angle = 0, size = 3.5) + 
    annotate("text", x = 366590, y = 7807484, label = "Montebello \n Islands", 
             fontface = 'italic', colour = "black", angle = 0, size = 3.5) + 
    annotate("text", x = 477127, y = 7807484, label = "Dampier \n Archipelago", 
             fontface = 'italic', colour = "black", angle = 0, size = 3.5) +
    
    geom_segment(aes(x = 197328, y = 7576897, xend = 210461, yend = 7626577)) + 
    geom_segment(aes(x = 346564, y = 7740123, xend = 360590, yend = 7787484)) + 
    geom_segment(aes(x = 456844, y = 7727229, xend = 474127, yend = 7787484)) + 
    
    scale_colour_manual(values = c("black", "white"),  name = "Surveys",
                        guide = guide_legend(override.aes = list(linetype = c("solid", "solid"), 
                                                                 shape = c(NA, NA)), size = 2)) +
    scale_fill_manual(values = c("#E7B800", "#65ABB7"), name = "Species",
                      guide = guide_legend(override.aes = list(linetype = c("blank", "blank"), 
                                                               shape = c(21, 21), size = 5))) + 
    scale_size_continuous(limits = c(1, 35), name = "Group size",
                          breaks = c(1, 2, 5, 10, 15, 20, 30), range = c(1, 10),
                          guide = guide_legend(override.aes = list(linetype = rep("blank", 7), 
                                                                   shape = rep(21, 7), fill = "black"))) +
    coord_sf(expand = TRUE, 
             xlim = c(sf::st_bbox(all.transects)[1],
                      sf::st_bbox(all.transects)[3] - 25000),
             ylim = c(sf::st_bbox(all.transects)[2] - 4000,
                      sf::st_bbox(all.transects)[4])) +
    theme_minimal() +
    xlab("") + ylab("") +
    theme(legend.position = "none", #c(0.82, 0.35)
          legend.box = "horizontal",
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          axis.text.x = element_text(colour = "black", size = 12),
          axis.text.y = element_text(colour = "black", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
          panel.background = element_rect(fill = "#bad8e0", colour = "#bad8e0", size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "#93bac4"))
  
  fig1.legend <- fig1 +
    geom_sf(data = sightings.sf, aes(size = size, fill = species), pch = 21, show.legend = "point") +
    coord_sf(expand = TRUE, 
             xlim = c(sf::st_bbox(all.transects)[1],
                      sf::st_bbox(all.transects)[3] - 25000),
             ylim = c(sf::st_bbox(all.transects)[2] - 4000,
                      sf::st_bbox(all.transects)[4])) +
    theme_minimal() +
    xlab("") + ylab("") +
    theme(legend.position = c(0.5, 0.35),
          legend.box = "horizontal",
          panel.background = element_rect(fill = "#bad8e0", colour = "#bad8e0", size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "#93bac4"))
  
  fig1 <- fig1 + ggsn::north(data = all.transects, location = "topleft", anchor = unlist(list(x = 180691, y = 7825095)), symbol = 12, scale = 0.05) + 
    annotate("text", label = "N", x = 196700, y = 7833600, size = 4.5, fontface = "bold") +
    ggsn::scalebar(all.transects, dist = 50, dist_unit = "km", transform = FALSE, model = "WGS84",
                   border.size = 0.35, st.dist = 0.03, st.size = 4, 
                   anchor = unlist(list(x = 556700, y = 7515000)))
  
  fig1.legend <- ggpubr::get_legend(fig1.legend)
  fig1.legend <- as_ggplot(fig1.legend)
  
  inset.pos <- purrr::map(.x = LETTERS[1:n.insets], 
                          .f = ~{
    if(.x == "A") res <- c(sf::st_bbox(insets[[.x]])[1] + 2500, sf::st_bbox(insets[[.x]])[4] - 2500)
    if(.x == "B") res <- c(sf::st_bbox(insets[[.x]])[1] + 1500, sf::st_bbox(insets[[.x]])[4] - 2000)
    if(.x == "C") res <- c(sf::st_bbox(insets[[.x]])[1] + 2000, sf::st_bbox(insets[[.x]])[4] - 2500)
    if(.x == "D") res <- c(sf::st_bbox(insets[[.x]])[1] + 500, sf::st_bbox(insets[[.x]])[4] - 2500)
    if(.x == "E") res <- c(sf::st_bbox(insets[[.x]])[1] + 0, sf::st_bbox(insets[[.x]])[4] - 2500)
    if(.x == "F") res <- c(sf::st_bbox(insets[[.x]])[1] - 2500, sf::st_bbox(insets[[.x]])[4] - 4500)
    res
  }) %>% purrr::set_names(x = ., nm = LETTERS[1:n.insets])

  
  inset.plots <- purrr::map2(.x = insets,
                             .y = LETTERS[1:n.insets],
                             .f = ~{
                              
                              p <- ggplot() +
                                geom_sf(data = all.transects, aes(colour = year), 
                                        fill = "transparent", lwd = 0.4, show.legend = FALSE) +
                                geom_sf(data = transects.overlap, linetype = "dashed") +
                                geom_sf(data = sightings.sf, 
                                        aes(size = size, fill = species), pch = 21, show.legend = FALSE) +
                                geom_sf(data = wa, col = "#232323", fill = "#e1d4b2", lwd = 0.25) +
                                scale_fill_manual(values = c("#E7B800", "#65ABB7")) +
                                scale_color_manual(values = c("black", "white")) +
                                coord_sf(expand = TRUE, 
                                         xlim = c(sf::st_bbox(.x)[1], sf::st_bbox(.x)[3]),
                                         ylim = c(sf::st_bbox(.x)[2], sf::st_bbox(.x)[4])) +
                                theme_minimal() + xlab("") + ylab("") +
                                theme(panel.background = element_rect(fill = "#bad8e0", 
                                                                      colour = "#bad8e0", size = 0.5, linetype = "solid"),
                                      panel.grid.major = element_line(size = 0, linetype = 'solid', 
                                                                      colour = "#93bac4"),
                                      axis.title.x = element_blank(),
                                      axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.title.y = element_blank(),
                                      axis.text.y = element_blank(),
                                      axis.ticks.y = element_blank(),
                                      panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
                                scale_size_continuous(limits = c(1, 35), 
                                                      breaks = c(1, 2, 5, 10, 15, 20, 30), 
                                                      range = c(1, 10),
                                                      guide = guide_legend(override.aes = 
                                                                             list(shape = 16, fill = "black"))) +
                                ggspatial::annotation_scale(location = "br", text_cex = 0.85) +
                                annotate("text", label = .y, x = inset.pos[[.y]][1], 
                                         y = inset.pos[[.y]][2], size = 8, fontface = "bold")
                              
                              p})  
  names(inset.plots) <- LETTERS[1:n.insets]
  

  
  ggsave(filename = "out/figure1.pdf", plot = fig1, dpi = 450, height = 7, width = 10)
  ggsave(filename = "out/figure1_legend.pdf", plot = fig1.legend, dpi = 450)
  
  plot.width <- plot.height <- list()
  plot.width$A <- 2.5
  plot.height$A <- 2.5
  plot.width$B <- 3
  plot.height$B <- 3
  plot.width$C <- 2.5
  plot.height$C <- 2.5
  plot.width$D <- 4
  plot.height$D <- 3
  plot.width$E <- 4
  plot.height$E <- 4
  plot.width[["F"]] <- 5
  plot.height[["F"]] <- 4
  
  # ggsave(filename = "out/figure1_inset_A.pdf", plot = inset.plots$A, dpi = 450, width = 4, height = 4)
  for(pl in LETTERS[1:n.insets]) ggsave(filename = paste0("out/figure1_inset_", pl, ".pdf"), dpi = 450,
                                 plot = inset.plots[[pl]],
                                 width = plot.width[[pl]],
                                 height = plot.height[[pl]])
}

##%######################################################%##
#                                                          #
#'### Function to simulate transparency using HEX colours   ####
#                                                          #
##%######################################################%##

#' @param input.colour Input colour
#' @param opacity Transparency level
#' @param bg.colour Background colour

hexa2hex <- function(input.colour, 
                     opacity, 
                     bg.colour = "white"){
  
  #'---------------------------------------------
  # PARAMETERS
  #'---------------------------------------------
  #' @param input.colour Initial colour.
  #' @param opacity Desired level of transparency (number between 0 and 1).
  #' @param bg.colour Colour of the background. Defaults to 'white'.
  #'---------------------------------------------
  
  # White background
  
  bg <- grDevices::col2rgb(bg.colour, alpha = FALSE)
  
  # Convert input colour to RGB
  
  rgbcol <- grDevices::col2rgb(input.colour, alpha = FALSE)
  
  # Calculate red, green, blue values corresponding to input colour at chosen transparency level
  
  rc <- (1 - opacity) * bg[1,] + opacity * rgbcol[1,]
  gc <- (1 - opacity) * bg[2,] + opacity * rgbcol[2,]
  bc <- (1 - opacity) * bg[3,] + opacity * rgbcol[3,]
  
  # Convert back to hex
  
  rgb2hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)
  return(rgb2hex(r = rc, g = gc, b = bc))
}