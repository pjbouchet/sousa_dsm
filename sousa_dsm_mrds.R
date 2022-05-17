##%######################################################%##
#                                                          #
#'###     Contrasting aerial abundance estimates for     ####
#'###          two sympatric dolphin species at          ####
#'###          a regional scale using distance           ####
#'###       sampling and density surface modelling       ####
#                                                          #
##%######################################################%##

# Authors: Raudino et al. (2022)

# Initialisation ----------------------------------------------------------

library(mrds)
library(tidyverse)

trunc.dist <- list(tursiops = 520, sousa = 445) # From sousa_dsm_truncation.R

mrds.df <- data.frame(
  species = c(rep("tursiops", 3), rep("sousa", 2)), 
  year = c(2015:2017, 2016:2017),
  p0 = NA,
  p0.SE = NA,
  FI = NA,
  FI.cov = NA,
  PI = NA,
  PI.cov = NA,
  final.model = NA)

init.list <- setNames(vector(mode = "list", length = 3), c(2015:2017))
mrds.gof.list <- mrds.model.list <- 
  list(tursiops = list(FI = init.list, PI = init.list),
       sousa = list(FI = init.list, PI = init.list))

for(i in 1:nrow(mrds.df)){
  
  plot.title <- paste0(stringr::str_to_sentence(mrds.df$species[i]), " (", mrds.df$year[i], ")")
  cat(plot.title, "--------------------\n")
  
  # Import the dolphin data
  detections <- read.csv(paste0("data/mrds/", mrds.df$year[i], "_", mrds.df$species[i], "_MRDS.csv"))
  
  # Investigate the structure of the dataset
  # str(detections)
  
  # Format variables
  detections$observer <- as.factor(detections$observer)
  detections$object <- as.numeric(detections$object)
  detections$detected <- as.numeric(detections$detected)
  detections$size <- as.numeric(detections$size)
  detections$turbidity <- as.factor(detections$turbidity)
  detections$beaufort <- as.factor(detections$beaufort)
  detections$distance <- detections$distance - 100
  detections <- na.omit(detections)

  pdf(file = paste0("data/mrds/", mrds.df$species[i], "_", mrds.df$year[i], "_plots.pdf"))
  par(mfrow = c(2, 2))
  
  # Full independence (FI) ------------------------------------------------------
  
  cat("\nFitting full independence models ...\n")
  
  # IO configuration with full independence model as a starting point
  fi.mr.dist <- mrds.model.list[[mrds.df$species[i]]]$FI[[as.character(mrds.df$year[i])]][["nocovariate"]] <- 
    mrds::ddf(method ='io.fi', 
                          mrmodel = ~glm(link = 'logit', formula = ~distance),
                          data = detections,
                          meta.data = list(width = trunc.dist[[mrds.df$species[i]]]))
  
  # summary(fi.mr.dist)
  
  # Produce goodness of fit statistics and a qq plot
  mrds.gof.list[[mrds.df$species[i]]]$FI[[as.character(mrds.df$year[i])]][["~distance"]] <- 
    mrds::ddf.gof(fi.mr.dist, main = paste0(plot.title, " - FI model"), pch = 16)
  
  # Set up list with possible models containing covariates
  mr.formula <- c("~distance","~distance+size","~distance+turbidity","~distance+beaufort",
                  "~distance+size+turbidity","~distance+size+beaufort","~distance+turbidity+beaufort",
                  "~distance+size+turbidity+beaufort")
  
  num.mr.models <- length(mr.formula)
  
  # Create data.frame to store results
  fi.results <- data.frame(MRmodel = mr.formula, AIC = rep(NA, num.mr.models))
  
  # Loop through all MR models
  for (mod in 1:num.mr.models) {
    print(mr.formula[mod])
    tryCatch(fi.model  <- mrds::ddf(method = 'io.fi', 
                           mrmodel = ~glm(link='logit', formula = as.formula(mr.formula[mod])),
                           data = detections, 
                           meta.data = list(width = trunc.dist[[mrds.df$species[i]]])),
             error = function(e) NA)
    fi.results$AIC[mod] <- summary(fi.model)$aic
  }
  
  # Calculate delta AIC
  fi.results$deltaAIC <- fi.results$AIC - min(fi.results$AIC, na.rm = TRUE)
  fi.results <- fi.results[order(fi.results$deltaAIC), ]
  # knitr::kable(fi.results, digits = 2)
  
  # Fit chosen model
  fi.mr.dist.cov <- mrds.model.list[[mrds.df$species[i]]]$FI[[as.character(mrds.df$year[i])]][["covariates"]] <- 
    mrds::ddf(
    method = "io.fi",
    mrmodel = ~ glm(link = "logit", formula = as.formula(fi.results$MRmodel[1])),
    data = detections,
    meta.data = list(width = trunc.dist[[mrds.df$species[i]]]))
  
  mrds.gof.list[[mrds.df$species[i]]]$FI[[as.character(mrds.df$year[i])]][[fi.results$MRmodel[1]]] <- ddf.gof(fi.mr.dist.cov,
    main = paste0(plot.title, " - FI model\n", fi.results$MRmodel[1]), pch = 16)
  
  mrds.df$FI[i] <- AIC(fi.mr.dist)$AIC
  mrds.df$FI.cov[i] <- AIC(fi.mr.dist.cov)$AIC
  
  # if(mrds.df$FI[i] < mrds.df$FI.cov[i]){
  #   mrds.df$FI.chsq[i] <- round(gof.results.FI$chisquare$pooled.chi$p, 5)
  #   mrds.df$FI.CvM[i] <- round(gof.results.FI$dsgof$CvM$p, 5)
  #   mrds.df$FI.AIC[i] <- AIC(fi.mr.dist)$AIC
  # } else {
  #   mrds.df$FI.chsq[i] <- round(gof.results.FIcov$chisquare$pooled.chi$p, 5)
  #   mrds.df$FI.CvM[i] <- round(gof.results.FIcov$dsgof$CvM$p, 5)
  #   mrds.df$FI.AIC[i] <- AIC(fi.mr.dist.cov)$AIC
  # }
  
  # Point independence (PI) ------------------------------------------------------
  
  cat("\nFitting point independence models ...\n")
  
  # With the PI assumption, the overall detection probability function is obtained 
  # by combining the shape of the detection function estimated from the DS model and 
  # the intercept parameter obtained from the MR model.
  
  # Fit trial configuration with point independence model
  pi.mr.dist <- mrds.model.list[[mrds.df$species[i]]]$PI[[as.character(mrds.df$year[i])]][["nocovariate"]] <- ddf(method = 'io', 
                    mrmodel = ~glm(link = 'logit', formula = as.formula(fi.results$MRmodel[1])),
                    dsmodel = ~cds(key = 'hn'), 
                    data = detections,
                    meta.data = list(width = trunc.dist[[mrds.df$species[i]]]))
  
  mrds.gof.list[[mrds.df$species[i]]]$PI[[as.character(mrds.df$year[i])]][["~distance"]] <- 
    ddf.gof(pi.mr.dist, main = paste0(plot.title, " - PI model"), pch = 16)
  
  # Point independence model that include covariates in DS model
  # Use selected MR model, iterate across DS models
  ds.formula <- c("~size","~turbidity","~beaufort","~size+turbidity","~size+beaufort","~turbidity+beaufort", "~size+turbidity+beaufort")
  num.ds.models <- length(ds.formula)
  
  # Create dataframe to store results
  pi.results <- data.frame(DSmodel = ds.formula, AIC = rep(NA, num.ds.models))
  
  # Loop through ds models - use selected MR model from earlier
  for (mod in 1:num.ds.models) {
    print(ds.formula[mod])
    tryCatch(pi.model <- mrds::ddf(
      method = "io", 
      mrmodel = ~ glm(link = "logit", formula = as.formula(fi.results$MRmodel[1])),
      dsmodel = ~ mcds(key = "hn", formula = as.formula(ds.formula[mod])),
      data = detections, 
      meta.data = list(width = trunc.dist[[mrds.df$species[[i]]]])),
    error = function(e) NA)
    pi.results$AIC[mod] <- summary(pi.model)$AIC
  }
  
  # Calculate delta AIC
  pi.results$deltaAIC <- pi.results$AIC - min(pi.results$AIC, na.rm = TRUE)
  pi.results <- pi.results[order(pi.results$deltaAIC), ]
  # knitr::kable(pi.results, digits = 2)
  
  # Fit chosen model
  pi.mr.dist.cov <- mrds.model.list[[mrds.df$species[i]]]$PI[[as.character(mrds.df$year[i])]][["covariates"]] <- 
    mrds::ddf(
      method = "io",
      mrmodel = ~ glm(link = "logit", formula = as.formula(fi.results$MRmodel[1])),
      dsmodel = ~ mcds(key = "hn", formula = as.formula(pi.results$DSmodel[1])),
      data = detections,
      meta.data = list(width = trunc.dist[[mrds.df$species[i]]]))
  
  mrds.gof.list[[mrds.df$species[i]]]$PI[[as.character(mrds.df$year[i])]][[pi.results$DSmodel[1]]] <- 
   ddf.gof(pi.mr.dist.cov,
    main = paste0(plot.title, " - PI model\n", pi.results$DSmodel[1]), pch = 16)
  
  mrds.df$PI[i] <- AIC(pi.mr.dist)$AIC 
  mrds.df$PI.cov[i] <- AIC(pi.mr.dist.cov)$AIC
  
  if(min(c(mrds.df$FI[i], mrds.df$FI.cov[i])) > min(c(mrds.df$PI[i], mrds.df$PI.cov[i])) |
     (min(c(mrds.df$FI[i], mrds.df$FI.cov[i])) < min(c(mrds.df$PI[i], mrds.df$PI.cov[i]))) &
     (min(c(mrds.df$FI[i], mrds.df$FI.cov[i])) - min(c(mrds.df$PI[i], mrds.df$PI.cov[i])) <= 2))
     {
    mrds.df$final.model[i] <- "Point"
    if(mrds.df$PI.cov[i] < mrds.df$PI[i]) mysum <- summary(pi.mr.dist.cov) else mysum <- summary(pi.mr.dist)
    mrds.df$p0[i] <- round(mysum$mr.summary$average.p0, 3)
    mrds.df$p0.SE[i] <- round(mysum$mr.summary$average.p0.se, 3)
    } else {
    mrds.df$final.model[i] <- "Full"
    if(mrds.df$FI.cov[i] < mrds.df$FI[i]) mysum <- summary(fi.mr.dist.cov) else mysum <- summary(fi.mr.dist)
    mrds.df$p0[i] <- round(mysum$average.p0, 3)
    mrds.df$p0.SE[i] <- round(mysum$average.p0.se, 3)
    }

  dev.off()

  if(i == 3){
  mrds.df$p0[3] <- round(summary(pi.mr.dist)$mr.summary$average.p0, 3)
  mrds.df$p0.SE[3] <- round(summary(pi.mr.dist)$mr.summary$average.p0.se, 3)
  }
  
  cat("\n")
  
} # End for i

save(mrds.df, mrds.gof.list, mrds.model.list, file = "data/sousa_mrds.RData")