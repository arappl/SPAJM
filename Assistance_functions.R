## Simulation function ##################################################
#' Function to simulate a data structure suitable for joint modelling
#'
#' @param n number of individuals
#' @param geo_pos position of spatial effect in predictors, takes values c("none", "etal", "etals", "etas"). "" means "no spatial effect in model"
#' @param geo_file path to .bnd file, from with 

sim <- function(n, geo_pos, geo_file) {
  # initialize values
  geo_l <- geo_ls <- geo_s <- 0
  
  # n = 200
  n_i = 6
  rho <- 0
  alpha = -0.3
  
  int <- 0
  betal <- -0.5 #c(-0.5, 0.7, 1.3, 0.3, 0.5)
  betat <-  0.4
  betals <- 0.9 #c(0.9, 0.3, -1, 0.2, -0.4)
  betas <-  0.1
  sigma2 <- 0.5
  
  switch(geo_pos,
         "etal" = {geo_l <- 1},
         "etals" = {geo_ls <- 1},
         "etas" = {geo_s <- 1})
  
  # Start sim
  ### generate id vector
  id = rep(1:n, each = n_i)
  first = rep(c(1, rep(0, n_i-1)), n)
  
  ### generate covariate matrices, standardize columns to mean = 0, sd = 1
  ### note, that Xls contains duplicated values of the same measurement, Xls_un not
  if(all(betal == 0) == TRUE){Xl = matrix(0, n*n_i, 1)}else{
    pl = length(betal)
    Xl <- runif(n * n_i * pl, -1, 1)
    # generate correlation matrix
    S <- rho^as.matrix(dist(1:pl))
    # generate design matrix
    Xl <- matrix(Xl, ncol= pl, nrow = n*n_i) %*% chol(S)
    colnames(Xl) <- paste0("Xl", 1:pl)
    Xl <- scale(Xl)
  }
  # for saving in output
  fl1 <- data.frame(Reduce(cbind, lapply(1:pl, function(x) Xl[, x]*betal[x])))
  names(fl1) <- paste0("fl", 1:ncol(fl1), "_lin")
  
  if(all(betas == 0) == TRUE){Xs = matrix(0, n, 1)}else{
    ps = length(betas)
    Xs <- runif(n * ps, -1, 1)
    # generate correlation matrix
    S <- rho^as.matrix(dist(1:ps))
    # generate design matrix
    Xs <- matrix(Xs, ncol= ps, nrow = n) %*% chol(S)
    colnames(Xs) <- paste0("Xs", 1:ps)
    Xs <- scale(Xs)
  }
  fs1 <- data.frame(Reduce(cbind, lapply(1:ps, function(x) Xs[, x]*betas[x])))
  names(fs1) <- paste0("fs", 1:ncol(fs1), "_lin")
  
  if(all(betals == 0) == TRUE){Xls = matrix(0, n*n_i, 1); Xls_un = matrix(0, n, 1)}else{
    pls = length(betals)
    Xls <- rep(runif(n * pls, -1, 1), each = n_i)
    # generate correlation matrix
    S <- rho^as.matrix(dist(1:pls))
    # generate design matrix
    Xls <- matrix(Xls, ncol= pls, nrow = n*n_i) %*% chol(S)
    colnames(Xls) <- paste0("Xls", 1:pls)
    Xls <- scale(Xls)
    Xls_un = Xls[first==1,, drop = FALSE]
  }
  fls1 <- data.frame(Reduce(cbind, lapply(1:pls, function(x) Xls[, x]*betals[x])))
  names(fls1) <- paste0("fls", 1:ncol(fls1), "_lin")
  
  # generate vars and non-linear effects
  Xl2 <- matrix(0, n*n_i, 1)
  Xl2[,1] <- scale(runif(n*n_i, -1, 1))
  fl2 <- function(e) {0.5*e + 15* dnorm(2*(e-0.2)) - dnorm(e + 0.4)}
  fl2 <- fl2(Xl2[,1])
  
  Xs2 <- matrix(0, n, 1)
  Xs2[,1] <- scale(runif(n, -pi, pi))
  fs2 <- 0.5*sin(Xs2[,1])
  
  Xls2 <- matrix(0, n*n_i, 1)
  Xls2[,1] = rep(scale(runif(n, -pi, pi)), each = n_i)
  fls2 <- -0.5*sin(Xls2[,1])
  
  Xls3 <- scale(runif(n*n_i, -1, 1))
  fls3 <- -0.5*Xls3
  
  # generate geo effect
  # map <- read.bnd(paste0(geo_path, "westerngermany.bnd"))
  map <- BayesX:::read.bnd(geo_file)
  centroids <- BayesX:::get.centroids(map)
  centroids <- centroids[!duplicated(rownames(centroids)), ]
  # standardize centroids == spatial variable s
  centroids_xsd <- (centroids[, "centroid_x"]-mean(centroids[,"centroid_x"]))/sd(centroids[,"centroid_x"])
  centroids_ysd <- (centroids[,"centroid_y"]-mean(centroids[,"centroid_y"]))/sd(centroids[,"centroid_y"])
  # geo effect
  f_geo0 <- sin(centroids_xsd) * cos(0.5 * centroids_ysd)
  f_geo0 <- f_geo0 - mean(f_geo0)
  
  n_geo <- nrow(centroids)
  
  if(n >= n_geo) {
    regions <- unlist(lapply(1:floor(n/n_geo), function(z) sample(rownames(centroids), size = n_geo, replace = F)))
    regions <- c(regions, sample(rownames(centroids), size = n %% n_geo, replace=F))
  } else {
    regions <- sample(rownames(centroids), size = n, replace=F)
  }
  
  f_geo <- f_geo0[match(regions, rownames(centroids))]
  
  
  ### generate time points
  day_in_year = sample(1:365, n*n_i, replace = TRUE)
  time = rep(seq(0,(n_i-1)*365, 365), n)
  time = time + day_in_year
  for(i in 1:n){time[id==i] = time[id==i] - min(time[id==i])}
  T_long = time/(n_i*365)
  
  ### generate random effects and compute predictor vectors
  gamma0 = rnorm(n, 0, sqrt(2))
  gamma1 = rnorm(n, 0, sqrt(2))
  etal = int + Xl%*%betal + fl2 + geo_l*rep(f_geo, each = n_i)
  etas = Xs%*%betas + fs2 + geo_s*f_geo
  etals = Xls%*%betals + fls2 + fls3 + geo_ls*rep(f_geo, each = n_i) + rep(gamma0, each=n_i) + T_long*(rep(gamma1, each=n_i) + betat)
  
  ### simulate longitudinal outcome
  y = rnorm(n*n_i, etal + etals, sqrt(sigma2))
  # y <- scale(y, scale = F)
  
  if (!geo_pos %in% c("etas", "none")) {
    df <- data.frame(id, y, Xl, Xl2, Xls, Xls2, Xls3, T_long, fl1, fl2, fls1, fls2, fls3, "ft" = betat*T_long,
                     "regions" = rep(regions, each = n_i), "f_geo" = rep(f_geo, each = n_i),
                     "raneff" = rep(gamma0, each = n_i), "ranslopes" = rep(gamma1, each = n_i),
                     etals, etal)
  } else {
    df <- data.frame(id, y, Xl, Xl2, Xls, Xls2, Xls3, fl1, fl2, fls1, fls2, fls3, T_long,
                     "raneff" = rep(gamma0, each = n_i), "ranslopes" = rep(gamma1, each = n_i),
                     etals, etal)
  }
  
  scale = 0.4; shape = 1.5
  lambda = function(t){scale*shape*t^(shape - 1)}
  
  # etals.h <- Xls%*%betals
  etals.h <- etals - (rep(gamma0, each=n_i) + T_long*(rep(gamma1, each=n_i) + betat)) # subtract time, since the hazard part takes time directly into account
  
  t.dom = seq(0, 15, 0.001)
  
  # Define repetition numbers of etals for haz.dom
  etals_reps <- tapply(findInterval(T_long, t.dom, left.open = T), id, function(x) diff(sort(c(x, length(t.dom)))))
  # ## explanation
  # T_long_in_t.dom <- findInterval(T_long, t.dom, left.open = T) # positions of individual T_longs in t.dom
  # # ## diff(sort(c(x, length(t.dom))) The differences between interval bounds (positions and the upper limit of t.dom [length(t.dom)])
  # # ## amount to the number of times each value of etals has to be repeated
  # etals_reps <- tapply(T_long_in_t.dom, id, function(x) diff(sort(c(x, length(t.dom)))))
  
  haz.dom = lambda(t.dom)
  haz.dom = matrix(haz.dom, n, length(haz.dom), byrow = TRUE)
  
  etals.dom <- rep(etals.h, unlist(etals_reps))
  etals.dom <-  matrix(etals.dom, n, length(t.dom), byrow = TRUE)
  
  haz.dom = haz.dom * exp(as.numeric(etas) + alpha*(etals.dom + gamma0 + (betat + gamma1) %*% t(t.dom)))
  
  # haz.dom = matrix(haz.dom, n, length(haz.dom), byrow = TRUE)
  # haz.dom = cumsum(haz.dom*.001)
  
  haz.dom = t(apply(haz.dom, 1, function(x) cumsum(x*.001)))
  
  St = exp(-haz.dom)
  
  # T_surv = c()
  # for(j in 1:n){
  #   u = runif(1)
  #   T_surv[j] = t.dom[sum(St[j,] > u)]
  # }
  
  T_surv = c()
  for(j in 1:n){
    u = runif(1)
    if(sum(St[j,] > u) == 0) {
      T_surv[j] <- 0
    } else {
      T_surv[j] = t.dom[sum(St[j,] > u)]
    }
  }
  if(any(T_surv == 0)) {T_surv[T_surv == 0] <- 0.000001}
  
  delta = rep(1, n)
  censoring = 1
  delta[T_surv>censoring] = 0
  T_surv[delta==0] = censoring
  
  # ## but 50% of the censored had a censoring time before 1
  n_cens <- floor(0.5*sum(delta == 0))
  ind <- sample(1:length(T_surv[delta == 0]), n_cens)
  T_surv[delta == 0][ind] <- runif(n_cens)
  
  lambda_surv <- unlist(lapply(1:n, function(x) {haz.dom[x, sum(t.dom <= T_surv[x])]}))
  
  if (geo_pos == "etas") {
    df.id <- data.frame("id" = unique(id), T_surv, delta, Xs, Xs2, "lambda0" = lambda(T_surv), fs1, fs2, etas, lambda_surv,
                        "regions" = regions, "f_geo" = f_geo)
  } else {
    df.id <- data.frame("id" = unique(id), T_surv, delta, Xs, Xs2, "lambda0" = lambda(T_surv), fs1, fs2, etas, lambda_surv)
  }
  
  time_zero <- which(T_long > rep(T_surv, each = n_i))
  
  df <- df[-time_zero, ]
  
  names(df) <- gsub("\\.", "", names(df))
  names(df.id) <- gsub("\\.", "", names(df.id))
  
  input <- list("n" = n, "n_i" = n_i, "rho" = rho, 
                "int" = int, "betal" = betal, "betat" = betat, "betals" = betals, "betas" = betas, "alpha" = alpha, 
                "fl" = "0.5*x + 15* dnorm(2*(x-0.2)) - dnorm(x + 0.4)", 
                "fs" = "0.5*sin(x)", 
                "fls" = "-0.5*sin(x)",
                "fls(t)_lin" = -0.5,
                "ft" = "betat * T_long",
                "f_geo" = list("function" = "sin(centroids_x) * cos(0.5 * centroids_y)", "position" = geo_pos, "map" = geo_file),
                "lambda0" = c("scale" = 0.4, "shape" = 1.5, "function" = "scale*shape*t^(shape - 1)"),
                "sigma2" = sigma2, "B" = matrix(c(2,0,0,2), 2, 2))
  
  if (geo_pos == "none") {map_save <- NULL} else {map_save <- map}
  return(list("input" = input, "df" = df, "df.id" = df.id, "map" = map_save))
}
