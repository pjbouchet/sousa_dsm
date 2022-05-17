##%######################################################%##
#                                                          #
#'###     Contrasting aerial abundance estimates for     ####
#'###          two sympatric dolphin species at          ####
#'###          a regional scale using distance           ####
#'###       sampling and density surface modelling       ####
#                                                          #
##%######################################################%##

# Authors: Raudino et al. (2022)

library(tidyverse)
library(patchwork)

angles <- read.csv("data/aerial horizontal angles_2015_2016.csv")

aircraft.speed <- 100 * 0.514 # m/sec

# Horizontal angles - following Forcada et al. 2004
# Bottlenose dolphin abundance in the NW Mediterranean: addressing heterogeneity in distribution
alpha.backward <- quantile(angles$horizontal.angle, 0.90)
alpha.forward <- 180 - alpha.backward

# Define function to calculate time window (in sec)
timew <- function(x, forward, backward, speed){
  theta.forward <- abs(90-forward)
  theta.backward <- backward - 90
  x * (tan(theta.forward*(pi / 180)) + tan(theta.backward*(pi / 180))) / speed
} 

calc_rt <- function(species, t.window, empirical = TRUE){
  
  dmax <- round(max(dist.df[[species]]$distance), -2) + 100
  
  # Vector of right-truncation distances (m0)
  rt <- seq(100, dmax, 10)
  
  # Calculate mean time window for each right-truncation distance
  timewindows <- sapply(X = rt, FUN = function(a) {
  if(empirical) xx <- dist.df[[species]][dist.df[[species]]$distance < a, ]$distance else
    xx <- seq(100, dmax, 1)[seq(100, dmax, 1) < a]
    mean(timew(
      x = xx,
      forward = alpha.forward,
      backward = alpha.backward,
      speed = aircraft.speed
    ))
  })
  
  rt.gam <- mgcv::gam(timewindows~s(rt))
  opt.fun <- function(x){abs(t.window-predict(rt.gam, data.frame(rt=x)))}
  opt.value <- optimize(f = opt.fun, interval = c(100, dmax))$minimum
  
  p1 <- ggplot(data = data.frame(t = predict(rt.gam, data.frame(rt)), rt = rt), aes(x = rt, y = t)) +
    geom_hline(yintercept = t.window, col = "grey", linetype = "dashed") +
    geom_vline(xintercept = opt.value, col = "grey", linetype = "dashed") +
    geom_point(pch = 16, size = 2.5) +
    geom_line() +
    xlab("Right-truncation distance") +
    ylab("Time window") + 
    theme(axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 12, colour = "black"),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 00, l = 0)),
          plot.margin = margin(t = 1.5, r = 1, b = 1, l = 1, "cm"),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 12))
  
  # Percentage of observations discarded after right-truncation
  exclude.perc <- sum(dist.df[[species]]$distance > opt.value)/nrow(dist.df[[species]])
  
  p2 <- ggplot(data = data.frame(x = dist.df[[species]]$distance), aes(x)) +
    geom_histogram(colour = 'black', fill = "lightgrey", breaks = hist(dist.df[[species]]$distance, breaks = 10)$breaks) +
    theme(axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 12, colour = "black"),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 00, l = 0)),
          plot.margin = margin(t = 1.5, r = 1, b = 1, l = 1, "cm"),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 12)) +
    geom_vline(xintercept = opt.value, colour = "steelblue") +
    xlab("Distance")
  
  print(p1 / p2)
  
  return(list(rt = opt.value, perc = exclude.perc))
  
}

trunc.sousa <- calc_rt("sousa", 7, TRUE)
# calc_rt("tursiops", 7.5)

##%######################################################%##
#                                                          #
####     Sousa: 445 m (4.72%) Tursiops: 520 m      ####
#                                                          #
##%######################################################%##

