rm(list = ls())
if (!require(pammtools)) {install.packages("pammtools")}
if (!require(BayesX)) {install.packages("BayesX")}
if (!require(parallel)) {install.packages("parallel")}

## #################### !!! ################### ##
## Set path to your local BayesX installation
BXpath <- ""
## ############################################ ##

# ## Source function for data generation
source("./Assistance_functions.R") 

## ############################################################################################# ##
## Spatial predictor in etas --------------------------------------------------------------------
## ############################################################################################# ##
## Data generation  -----------------------------------------------------------------------------#
set.seed(1234561)

dat <- sim(n= 200, geo_pos = "etas",
           geo_file = file.path("./westerngermany.bnd"))

# ## Generate augmented data set for use in BayesX ############################ ##
# rename datasets for BayesX compatible variable names
names(dat$df) <- gsub("\\.", "", names(dat$df))
names(dat$df.id) <- gsub("\\.", "", names(dat$df.id))

# augment data set
ped_df <- as_ped(
  data = list(dat$df.id, dat$df),
  formula = as.formula(paste0("Surv(T_surv, delta) ~ . + concurrent(", paste(names(dat$df)[-grep("id|T_long", names(dat$df))], collapse = ", ") , ", tz_var = 'T_long')")),
  id = "id",
  cut = sort(c(unique(dat$df$T_long), unique(dat$df.id$T_surv[dat$df.id$delta == 1])))
)


names(ped_df) <- gsub("offset", "off", names(ped_df))
ped_df <- ped_df[, -grep("interval", names(ped_df))]

dummy_df <- dat$df
names(dummy_df)[grep("T_long", names(dummy_df))] <- "tend"
dummy_df$tstart <- unlist(tapply(dummy_df$tend, INDEX = dummy_df$id, function(x) c(0, x[-length(x)])))
NAM <- names(ped_df)[!names(ped_df) %in% names(dummy_df)]
MAT <- matrix(0, nrow = nrow(dummy_df), ncol = length(NAM))
colnames(MAT) <- NAM
dummy_df <- cbind(dummy_df, MAT)
dummy_df$regions <- ped_df$regions[match(dummy_df$id, ped_df$id)]

# create weights to be used for correct estimation in BayesX
dummy_df$wy <- 1
ped_df$wy <- 0

ped_df <- rbind(ped_df, dummy_df)
ped_df$wh <- as.numeric(ped_df$wy == 0)

rm(NAM, MAT)

# save data set in ".raw"-format for BayesX use
write.table(ped_df, paste0("./data_etas.raw"), col.names=TRUE, row.names=FALSE, sep=" ", quote=FALSE)

## Batch Creation --------------------------------------------------------------------------------#
path <- "./"
dir <- "Results_etas"
if(!dir.exists(file.path(path, dir))) {
  dir.create(file.path(path, dir))
} 

options(scipen = 999)
iter = 70000
burnin = 10000
step = 60

jmProg <- paste0(path, "batch_etas.txt")
write(NULL, jmProg)
write(paste("% usefile",jmProg),jmProg)
write(" ",jmProg,append=T)
write("delimiter = ;",jmProg,append=T)
write(" ",jmProg,append=T)

write("dataset d;",jmProg,append=T)
write(paste0("d.infile using ", path, "data_etas.raw;"),jmProg,append=T)
write("map g;",jmProg,append=T)
write(paste0("g.infile using ", path, "westerngermany.bnd;"),jmProg,append=T)

write("mcmcreg m;" ,jmProg,append=T)
write(paste0("m.outfile =  ", path, dir, "/r ;") ,jmProg,append=T)
write(paste0("logopen, replace using ", path, dir, "/log.txt ;") ,jmProg,append=T)

write(paste0("m.hregress y = Xls1 + Xls3 + Xls2(pspline) + tend + id(random) + tend*id(random), family=gaussian_shared equationtype=mu_shared iterations=",iter," burnin=",burnin," step=",step," IWLSlineff using d;"), jmProg, append=T)
write(paste0("m.hregress y = const + Xl1 + Xl2(pspline) weight wy, family=gaussian equationtype=mu using d;"), jmProg, append=T)
write(paste0("m.hregress ped_status = const + tend(pspline) + off(offset) + Xs1 + Xs2(pspline) + regions(spatial, map = g) weight wh,
  family=poisson equationtype=lambda setseed = 123 using d;"), jmProg, append=T)

write("m.getsample; \n", jmProg, append=T)
write("drop d m;",jmProg,append=T)
write(" ",jmProg,append=T)
write("delimiter = return;",jmProg,append=T)

## Run BayesX --------------------------------------------------------------------------------#
batch <- list.files(pattern ="batch_etas", recursive = T, full.names = T)

times <- system.time(
  system(
    paste(BXpath, batch, sep = " "),
    wait = TRUE)
)
save(times, file = file.path("./BX_run_times_etas.Rdata"))

## ############################################################################################# ##
## Spatial predictor in etals --------------------------------------------------------------------
## ############################################################################################# ##
rm(list = grep("sim|BXpath", ls(), value = T, invert = T))
## Data generation  -----------------------------------------------------------------------------#
set.seed(1234561)

dat <- sim(n= 200, geo_pos = "etals",
           geo_file = file.path("./westerngermany.bnd"))

# ## Generate augmented data set for use in BayesX ############################ ##
# rename datasets for BayesX compatible variable names
names(dat$df) <- gsub("\\.", "", names(dat$df))
names(dat$df.id) <- gsub("\\.", "", names(dat$df.id))

# augment data set
ped_df <- as_ped(
  data = list(dat$df.id, dat$df),
  formula = as.formula(paste0("Surv(T_surv, delta) ~ . + concurrent(", paste(names(dat$df)[-grep("id|T_long", names(dat$df))], collapse = ", ") , ", tz_var = 'T_long')")),
  id = "id",
  cut = sort(c(unique(dat$df$T_long), unique(dat$df.id$T_surv[dat$df.id$delta == 1])))
)

names(ped_df) <- gsub("offset", "off", names(ped_df)) # rename variable "offset" to avoid confusion with "offset" operator in BayesX
ped_df <- ped_df[, -grep("interval", names(ped_df))] # eliminate non-quantitative variables

dummy_df <- dat$df
names(dummy_df)[grep("T_long", names(dummy_df))] <- "tend"
dummy_df$tstart <- unlist(tapply(dummy_df$tend, INDEX = dummy_df$id, function(x) c(0, x[-length(x)])))
NAM <- names(ped_df)[!names(ped_df) %in% names(dummy_df)]
MAT <- matrix(0, nrow = nrow(dummy_df), ncol = length(NAM))
colnames(MAT) <- NAM
dummy_df <- cbind(dummy_df, MAT)

# create weights to be used for correct estimation in BayesX
dummy_df$wy <- 1
ped_df$wy <- 0

ped_df <- rbind(ped_df, dummy_df)
ped_df$wh <- as.numeric(ped_df$wy == 0)

rm(NAM, MAT)

# save data set in ".raw"-format for BayesX use
write.table(ped_df, paste0("./data_etals.raw"), col.names=TRUE, row.names=FALSE, sep=" ", quote=FALSE)

## Batch Creation --------------------------------------------------------------------------------#
path <- "./"
dir <- "Results_etals"
if(!dir.exists(file.path(path, dir))) {
  dir.create(file.path(path, dir))
} 

options(scipen = 999)
iter = 70000
burnin = 10000
step = 60

jmProg <- paste0(path, "batch_etals.txt")
write(NULL, jmProg)
write(paste("% usefile",jmProg),jmProg)
write(" ",jmProg,append=T)
write("delimiter = ;",jmProg,append=T)
write(" ",jmProg,append=T)

write("dataset d;",jmProg,append=T)
write(paste0("d.infile using ", path, "data_etals.raw;"),jmProg,append=T)
write("map g;",jmProg,append=T)
write(paste0("g.infile using ", path, "westerngermany.bnd;"),jmProg,append=T)

write("mcmcreg m;" ,jmProg,append=T)
write(paste0("m.outfile =  ", path, dir, "/r ;") ,jmProg,append=T)
write(paste0("logopen, replace using ", path, dir, "/log.txt ;") ,jmProg,append=T)

write(paste0("m.hregress y = Xls1 + Xls3 + Xls2(pspline) + tend + regions(spatial, map = g) + id(random) + tend*id(random), family=gaussian_shared equationtype=mu_shared iterations=",iter," burnin=",burnin," step=",step," IWLSlineff using d;"), jmProg, append=T)
write(paste0("m.hregress y = const + Xl1 + Xl2(pspline) weight wy, family=gaussian equationtype=mu using d;"), jmProg, append=T)
write(paste0("m.hregress ped_status = const + tend(pspline) + off(offset) + Xs1 + Xs2(pspline) weight wh,
family=poisson equationtype=lambda setseed = 321 using d;"), jmProg, append=T)

write("m.getsample; \n", jmProg, append=T)
write("drop d m;",jmProg,append=T)
write(" ",jmProg,append=T)
write("delimiter = return;",jmProg,append=T)

## Run BayesX --------------------------------------------------------------------------------#
batch <- list.files(pattern ="batch_etals", recursive = T, full.names = T)

times <- system.time(
  system(
    paste(BXpath, batch, sep = " "),
    wait = TRUE)
)
save(times, file = file.path("./BX_run_times_etals.Rdata"))

## ############################################################################################# ##
## Spatial predictor in etal --------------------------------------------------------------------
## ############################################################################################# ##
rm(list = grep("sim|BXpath", ls(), value = T, invert = T))
## Data generation  -----------------------------------------------------------------------------#
set.seed(1234561)

dat <- sim(n= 200, geo_pos = "etal",
           geo_file = file.path("./westerngermany.bnd"))

# ## Generate augmented data set for use in BayesX ############################ ##
# rename datasets for BayesX compatible variable names
names(dat$df) <- gsub("\\.", "", names(dat$df))
names(dat$df.id) <- gsub("\\.", "", names(dat$df.id))

# augment data set
ped_df <- as_ped(
  data = list(dat$df.id, dat$df),
  formula = as.formula(paste0("Surv(T_surv, delta) ~ . + concurrent(", paste(names(dat$df)[-grep("id|T_long", names(dat$df))], collapse = ", ") , ", tz_var = 'T_long')")),
  id = "id",
  cut = sort(c(unique(dat$df$T_long), unique(dat$df.id$T_surv[dat$df.id$delta == 1])))
)

names(ped_df) <- gsub("offset", "off", names(ped_df)) # rename variable "offset" to avoid confusion with "offset" operator in BayesX
ped_df <- ped_df[, -grep("interval", names(ped_df))] # eliminate non-quantitative variables

dummy_df <- dat$df
names(dummy_df)[grep("T_long", names(dummy_df))] <- "tend"
dummy_df$tstart <- unlist(tapply(dummy_df$tend, INDEX = dummy_df$id, function(x) c(0, x[-length(x)])))
NAM <- names(ped_df)[!names(ped_df) %in% names(dummy_df)]
MAT <- matrix(0, nrow = nrow(dummy_df), ncol = length(NAM))
colnames(MAT) <- NAM
dummy_df <- cbind(dummy_df, MAT)

# create weights to be used for correct estimation in BayesX
dummy_df$wy <- 1
ped_df$wy <- 0

ped_df <- rbind(ped_df, dummy_df)
ped_df$wh <- as.numeric(ped_df$wy == 0)

rm(NAM, MAT)

# save data set in ".raw"-format for BayesX use
write.table(ped_df, paste0("./data_etal.raw"), col.names=TRUE, row.names=FALSE, sep=" ", quote=FALSE)

## Batch Creation --------------------------------------------------------------------------------#
path <- "./"
dir <- "Results_etal"
if(!dir.exists(file.path(path, dir))) {
  dir.create(file.path(path, dir))
} 

options(scipen = 999)
iter = 70000
burnin = 10000
step = 60

jmProg <- paste0(path, "batch_etal.txt")
write(NULL, jmProg)
write(paste("% usefile",jmProg),jmProg)
write(" ",jmProg,append=T)
write("delimiter = ;",jmProg,append=T)
write(" ",jmProg,append=T)

write("dataset d;",jmProg,append=T)
write(paste0("d.infile using ", path, "data_etal.raw;"),jmProg,append=T)
write("map g;",jmProg,append=T)
write(paste0("g.infile using ", path, "westerngermany.bnd;"),jmProg,append=T)

write("mcmcreg m;" ,jmProg,append=T)
write(paste0("m.outfile =  ", path, dir, "/r ;") ,jmProg,append=T)
write(paste0("logopen, replace using ", path, dir, "/log.txt ;") ,jmProg,append=T)

write(paste0("m.hregress y = Xls1 + Xls3 + Xls2(pspline) + tend + id(random) + tend*id(random), family=gaussian_shared equationtype=mu_shared iterations=",iter," burnin=",burnin," step=",step," IWLSlineff using d;"), jmProg, append=T)
write(paste0("m.hregress y = const + Xl1 + Xl2(pspline) + regions(spatial, map = g) weight wy, family=gaussian equationtype=mu using d;"), jmProg, append=T)
write(paste0("m.hregress ped_status = const + tend(pspline) + off(offset) + Xs1 + Xs2(pspline) weight wh,
family=poisson equationtype=lambda setseed = 321 using d;"), jmProg, append=T)

write("m.getsample; \n", jmProg, append=T)
write("drop d m;",jmProg,append=T)
write(" ",jmProg,append=T)
write("delimiter = return;",jmProg,append=T)

## Run BayesX --------------------------------------------------------------------------------#
batch <- list.files(pattern ="batch_etal", recursive = T, full.names = T)

times <- system.time(
  system(
    paste(BXpath, batch, sep = " "),
    wait = TRUE)
)
save(times, file = file.path("./BX_run_times_etal.Rdata"))