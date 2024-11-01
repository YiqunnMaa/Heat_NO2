### Key R code for the study
### Spatial heterogeneity in synergistic effects of extreme heat and NO2 exposures
### on cardiorespiratory hospitalizations in California
library(dplyr)
library(data.table)
library(lubridate)
library(gnm)
library(msm)
library(tictoc)


## 1. state-wide case-crossover analysis
# RERI function
reri_interaction <- function(b, v, exposure1, exposure2) {
  reri <- exp(b[1] + b[2] + b[3]) - exp(b[1]) - exp(b[2]) + 1
  
  k1 <- exp(b[1] + b[2] + b[3]) - exp(b[1])
  k2 <- exp(b[1] + b[2] + b[3]) - exp(b[2])
  k3 <- exp(b[1] + b[2] + b[3])
  vreri <- v[1, 1] * k1 * k1 + 
    v[2, 2] * k2 * k2 + 
    v[3, 3] * k3 * k3 + 
    2 * v[1, 2] * k1 * k2 +
    2 * v[1, 3] * k1 * k3 +
    2 * v[2, 3] * k2 * k3
  se_reri <- sqrt(vreri)
  
  reri_ci95_l <- reri - 1.96 * se_reri
  reri_ci95_u <- reri + 1.96 * se_reri
  
  est <- exp(b)
  se <- sqrt(c(v[1, 1], v[2, 2], v[3, 3]))
  est_l <- exp(b - 1.96 * se)
  est_u <- exp(b + 1.96 * se)
  
  k3_se <- deltamethod(~exp(x1+x2+x3), b, v)
  k3_l <- k3 - 1.96 * k3_se
  k3_u <- k3 + 1.96 * k3_se
  
  temp <- data.frame(
    est=c(est, reri, k3),
    ll=c(est_l, reri_ci95_l, k3_l),
    ul=c(est_u, reri_ci95_u, k3_u),
    se=c(se, se_reri, k3_se))
  rownames(temp) <- c(exposure1, exposure2, "interaction_multiplicative", "interaction_reri", "joint_effet")
  return(temp)
}

# run the model
model_interaction <- gnm(csdresp ~ EH.this*NO2.this,
                         data = data, na.action = na.exclude,
                         family = quasipoisson, eliminate = stratum)

# extract the results
model_interaction.coef <- model_interaction$coefficients[c("EH.this","NO2.this","EH.this:NO2.this")]
model_interaction.vcov <- vcov(model_interaction, dispersion = NULL, with.eliminate = TRUE)

table.RERI <- reri_interaction(b=model_interaction.coef, v=model_interaction.vcov, "extreme_heat only", "NO2 exposure only")


## 2. within-community matched design and BHM
# define functions for the monthly weighting method
month.control <- function(exposed, control, nbuffer) { 
  if (length(exposed) > 0) {
    out <- lapply(exposed, function(e_day) {
      e_buffer <- e_day + (-nbuffer:nbuffer) ## potential days for control (buffer day window from exposure)
      c_pool <- control[control %in% e_buffer]  ## potential control days after excluding days close to other exposure
      return(c(e_day, c_pool))
    })
    names(out) <- exposed
  } else {
    out <- numeric()
  }
  return(out)
}

month.wt.analysis <- function(exposure, event, c.list, outcome.dt) { 
  if (length(c.list) > 0) {
    out <- numeric()
    for (i in 1:length(c.list)) {
      if (length(c.list[[i]]) > 1) {
        e_day <- c.list[[i]][1]
        c_pool <- data.table(date=c.list[[i]], exposed = c(1, rep(0, length(c.list[[i]])-1)))
        c_pool$wt <- 1/round(abs(as.numeric(difftime(c_pool$date, e_day, units = "days"))), digits = 0)
        c_pool$wt[1] <- 0
        c_pool <- merge(c_pool, outcome.dt, by="date", all.x=TRUE)
        
        rr <- c_pool[exposed==1, eval(as.name(event))]/c_pool[, sum(eval(as.name(event))*wt)/sum(wt)]
        rr <- ifelse(is.infinite(rr), NA, rr)
        out <- c(out, rr)
      }
    }
    
    if (length(out)==0) {
      temp <- simpleError(paste("No available control day"))
    } else {
      temp <- c(mean(out, na.rm=TRUE), sum(!is.na(out)))
      if (is.na(temp[1])) {
        temp <- simpleError(paste("No non-infinite RR"))
        temp$call <- NULL
      }
    }
    
    
  } else {
    temp <- simpleError(paste("No available exposed day"))
    temp$call <- NULL
  }
  return(temp)
}

# run the analysis for each exposure
for (i in 1:nrow(bar)) {
  data <- data.all[data.all$ZCTA==bar$ZCTA[i],]
  buffer_days <- data$date[data$EH.this==1 | data$NO2.this==1]
  buffer_days <- unique(c(buffer_days-3, buffer_days-2, buffer_days-1, buffer_days, buffer_days + 1, buffer_days + 2, buffer_days + 3))
  day00_bobb <- data$date[!(data$date %in% buffer_days)]
  bar$bobb_control_pool[i] <- length(day00_bobb)
  
  for (exposure_ in c("EH1NO21", "EH1NO20", "EH0NO21")) {
    # identify exposed and control days
    days <- data$date[data$EH.this==as.numeric(substring(exposure_, 3, 3)) & data$NO2.this==as.numeric(substring(exposure_, 7, 7))]
    bar[i, exposure_] <- length(days)
    baz2 <- month.control(days, day00_bobb, 30)
    
    # weighted analysis
    temp4 <- month.wt.analysis(exposure_, event_, baz2, data)
    if (!inherits(temp4, what = "condition")) {
      bar[i, paste0(c("month_wt_rr_", "month_wt_ngrp_"), exposure_)] <- temp4
    } else {
      fail <- rbind(fail, data.frame(ZCTA=bar$ZCTA[i], event = event_, exposure = exposure_, ngrps = length(baz2), method = "month_wt", error=as.character(temp4)))
      bar[i, paste0("month_wt_ngrp_", exposure_)] <- length(baz2)
    }
  }
}

# BHM
bef.sp <- spLM(est ~ 1, data = bayesDF, coords = as.matrix(bayesDF[,.(FinalLon, FinalLat)]),
               starting = list("phi" = 2.5, "sigma.sq" = 0.4, "tau.sq" = 0.17), 
               tuning = list("phi" = 0.125, "sigma.sq" = 0.02, "tau.sq" = 0.0085),
               priors = list("phi.Unif" = c(0.001, 6), 
                             "sigma.sq.IG" = c(2, 1/0.4), 
                             "tau.sq.IG" = c(2, 1/0.17)),
               cov.model = "spherical", n.samples = n.samples, verbose = TRUE, n.report=2000)


burn.in <- floor(0.75*n.samples)
bef.sp <- spRecover(bef.sp, start = burn.in, n.report=1000)


## 3. meta regression
for (i in 1:nrow(result.table)) {
  var.i <- result.table[i,"var"]
  
  # meta-regression of reri after bayesian spatial pooling
  m2 <- metagen(TE = RERI_pooled, seTE = RERI_pooled_SD, studlab = ZCTA, data = data.i)
  g2 <- metareg(m2, as.formula(paste0(" ~ ", var.i)))
  result.table[i, c("meta_coef", "meta_SE")] <- c(g2$beta[2], g2$se[2])
}

