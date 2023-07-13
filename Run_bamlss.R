rm(list = ls())
if (!require(BayesX)) {install.packages("BayesX")}
if (!require(parallel)) {install.packages("parallel")}
if (!require(bamlss)) {install.packages("bamlss")}

# ## Source function for data generation
source("./Assistance_functions.R") 

## ############################################################################################# ##
## Spatial predictor in etas --------------------------------------------------------------------
## ############################################################################################# ##

## Data generation  -----------------------------------------------------------------------------#
set.seed(1234561)

dat <- sim(n= 200, geo_pos = "etas",
           geo_file = file.path("./westerngermany.bnd"))

dat <- merge(dat$df, dat$df.id, all = T)

dat$id2 <- factor(dat$id) # convert id into factor
dat$regions1 <- factor(dat$regions) # convert regions into factor

# Ggenerate penalty matrix for estimation of spatial effect
map <- read.bnd(file.path("./westerngermany.bnd")) 
K <- bnd2gra(map)
dat$regions1 <- factor(dat$regions) # re-assign factor to have original n = 200 levels
IDcol <- IDrow <- colnames(K) %in% levels(dat$regions1)

K <- K[IDrow, IDcol]

## Estimation -----------------------------------------------------------------------------------#
# set formula for estimation
f <- list(
  Surv2(T_surv, delta, obs = y) ~ s(T_surv,k=20),
  gamma ~ Xs1 + s(Xs2, bs="ps", k=20) +
    te(regions1, bs = "mrf", xt = list(penalty = K)),
  mu ~ Xl1 + s(Xl2, bs="ps", k=20) + Xls1 + s(Xls2, bs="ps", k=20) + Xls3 + T_long +
    s(id2, bs="re") + s(id2, T_long,bs="re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

# estimate posterior modes
set.seed(321)
m_mode <- tryCatch({bamlss(f, data = dat, family = "jm",
                           timevar = "T_long", idvar = "id2", maxit = 5000,
                           verbose = F)
  },
  error = function(e) e)
save(m_mode, file = file.path("./etas_bamlss_mode.Rdata"))

# estimate posterior means with posterior modes as starting values
set.seed(321)
m <- tryCatch({bamlss(f, data = dat, family = "jm",
                      timevar = "T_long", idvar = "id2", optimizer = FALSE,
                      jm.start = coef(m_mode), n.iter = 12000, burnin = 2000, thin = 10,
                      verbose = F)
  },
  error = function(e) e)

save(m, file = file.path("./etas_bamlss_means.Rdata"))

## ############################################################################################# ##
## Spatial predictor in etals --------------------------------------------------------------------
## ############################################################################################# ##
rm(list = grep("sim", ls(), invert = T, value = T))
## Data generation  -----------------------------------------------------------------------------#
set.seed(1234561)

dat <- sim(n= 200, geo_pos = "etals",
           geo_file = file.path("./westerngermany.bnd"))

dat <- merge(dat$df, dat$df.id, all = T)

dat$id2 <- factor(dat$id) # convert id into factor
dat$regions1 <- factor(dat$regions) # convert regions into factor

# Ggenerate penalty matrix for estimation of spatial effect
map <- read.bnd(file.path("./westerngermany.bnd")) 
K <- bnd2gra(map)
dat$regions1 <- factor(dat$regions) # re-assign factor to have original n = 200 levels
IDcol <- IDrow <- colnames(K) %in% levels(dat$regions1)

K <- K[IDrow, IDcol]

## Estimation -----------------------------------------------------------------------------------#
# set formula for estimation
f <- list(
  Surv2(T_surv, delta, obs = y) ~ s(T_surv,k=20),
  gamma ~ Xs1 + s(Xs2, bs="ps",k=20),
  mu ~ Xl1 + s(Xl2, bs="ps", k=20) + Xls1 + s(Xls2, bs="ps", k=20) + Xls3 + T_long +
    te(regions, bs = "mrf", xt = list(penalty = K)) +
    s(id2, bs="re") + s(id2, T_long,bs="re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
  )

# estimate posterior modes
set.seed(321)
m_mode <- tryCatch({bamlss(f, data = dat, family = "jm",
                           timevar = "T_long", idvar = "id2", maxit = 5000,
                           verbose = F)
},
error = function(e) e)
save(m_mode, file = file.path("./etals_bamlss_mode.Rdata.Rdata"))

# estimate posterior means with posterior modes as starting values
set.seed(321)
m <- tryCatch({bamlss(f, data = dat, family = "jm",
                      timevar = "T_long", idvar = "id2", optimizer = FALSE,
                      jm.start = coef(m_mode), n.iter = 12000, burnin = 2000, thin = 10,
                      verbose = F)
},
error = function(e) e)

save(m, file = file.path("./etals_bamlss_mode.Rdata.Rdata"))

## ############################################################################################# ##
## Spatial predictor in etal --------------------------------------------------------------------
## ############################################################################################# ##
rm(list = grep("sim", ls(), invert = T, value = T))
## Data generation  -----------------------------------------------------------------------------#
set.seed(1234561)

dat <- sim(n= 200, geo_pos = "etal",
           geo_file = file.path("./westerngermany.bnd"))

dat <- merge(dat$df, dat$df.id, all = T)

dat$id2 <- factor(dat$id) # convert id into factor
dat$regions1 <- factor(dat$regions) # convert regions into factor

# Ggenerate penalty matrix for estimation of spatial effect
map <- read.bnd(file.path("./westerngermany.bnd")) 
K <- bnd2gra(map)
dat$regions1 <- factor(dat$regions) # re-assign factor to have original n = 200 levels
IDcol <- IDrow <- colnames(K) %in% levels(dat$regions1)

K <- K[IDrow, IDcol]

## Estimation -----------------------------------------------------------------------------------#
# set formula for estimation
f <- list(
  Surv2(T_surv, delta, obs = y) ~ s(T_surv,k=20),
  gamma ~ Xs1 + s(Xs2, bs="ps",k=20),
  mu ~ Xl1 + s(Xl2, bs="ps", k=20) + Xls1 + s(Xls2, bs="ps", k=20) + Xls3 + T_long +
    te(regions1, bs = "mrf", xt = list(penalty = K)) +
    s(id2, bs="re") + s(id2, T_long,bs="re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
  )

# estimate posterior modes
set.seed(321)
m_mode <- tryCatch({bamlss(f, data = dat, family = "jm",
                           timevar = "T_long", idvar = "id2", maxit = 5000,
                           verbose = F)
},
error = function(e) e)
save(m_mode, file = file.path("./etal_bamlss_mode.Rdata.Rdata"))

# estimate posterior means with posterior modes as starting values
set.seed(321)
m <- tryCatch({bamlss(f, data = dat, family = "jm",
                      timevar = "T_long", idvar = "id2", optimizer = FALSE,
                      jm.start = coef(m_mode), n.iter = 12000, burnin = 2000, thin = 10,
                      verbose = F)
},
error = function(e) e)

save(m, file = file.path("./etal_bamlss_mode.Rdata.Rdata"))