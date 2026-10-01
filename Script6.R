adm_location <- "d:/data/ms/Lorentz-Workshop/book/data/output/columns_adm.csv" # adjust if you have a different output location
data_adm <- read.csv(adm_location) # read age-depth model data
key <- "mid"
t <-  data_adm[,which(grepl("Myr",names(data_adm)))] # time step
h <- data_adm[,which(grepl(key,names(data_adm)))] # heights at the location with key in the name
# define age-depth model
adm <- admtools::tp_to_adm(t = t,
                          h = h,
                          T_unit = "Myr",
                          L_unit = "m")
plot(adm,
     lwd_destr = 1, # line width destructive interval (hiatus)
     lwd_acc = 2, # line width accumulative interval
     col_destr = NULL,
     col_acc = "black")
admtools::T_axis_lab()
admtools::L_axis_lab()
title(main = "Age-depth model")
seed <- 42
set.seed(seed)
lambda <- 2 # origination rate [Myr^-1]
mu <- 0.3 # extinction rate [Myr^-1]
psi <- 70 # fossil sampling rate [fossils per Myr per lineage]
# timescale of observation - depends on duration of your forward simulation
t_max <- admtools::max_time(adm)
tree <- TreeSim::sim.bd.age(age = t_max,
                           numbsim = 1,
                           lambda = lambda,
                           mu = mu,
                           frac = 1,
                           mrca = FALSE,
                           complete = TRUE)[[1]]
fossils <- FossilSim::sim.fossils.poisson(rate = psi,
                                         tree = tree,
                                         root.edge = FALSE)
plot(fossils, 
     tree = tree,
     main = "FBD tree")
t_min <- admtools::min_time(adm)
h_max <- admtools::max_height(adm)
h_min <- admtools::min_height(adm)
wd_location <- "d:/data/ms/Lorentz-Workshop/book/data/output/columns_wd.csv" # adjust if you have a different output location
data_wd <- read.csv(wd_location) # read age-depth model data
# select water depth column that contains the key
wd <- data_wd[,which(grepl(key,names(data_wd)))]
plot(x = t,
     y = wd,
     type = "l",
     xlab = "Time [Myr]",
     ylab = "Water depth [m]",
     main = "Water depth")
wd_opt <- 5 # pref. water depth 5 m
wd_tolerance <- 10 # water depth tolerance 10 m (sd)
cutoff_val <- 0  # does not live above water 
niche_def <- StratPal::snd_niche(opt = wd_opt, 
                                tol = wd_tolerance, 
                                cutoff_val = cutoff_val) 
wd_vals <- seq(from = min(wd),
              to = max(wd),
              length.out = 100)
plot(x = wd_vals,
     y = niche_def(wd_vals),
     type = "l",
     xlab = "Water depth [m]",
     ylab = "Sampling probability",
     main = "Niche definition")
# define function that specifies water depth as a function of time
gradient <- approxfun(x = t,
                     y = wd)
plot(x = t,
     y = gradient(t),
     type = "l",
     xlab = "Time [Myr]",
     ylab = "Water depth [m]",
     main = "Water depth as gradient")
plot(x = t,
     y = t |> gradient() |> niche_def(),
     type = "l",
     xlab = "Time [Myr]",
     ylab = "Sampling probability",
     main = "Sampling probability of taxa, time domain")
plot(x = t |> gradient() |> niche_def(),
     y = t |> admtools::time_to_strat(adm),
     type = "l",
     xlab = "Sampling probability",
     ylab = "Stratigrapic position [m]",
     main = "Sampling probability of taxa, stratigraphic domain")
# baseline scenarion without ecological effects
StratPal::p3_var_rate(x = c(0, t_max), # fossil abundance doubles over the time interval of observation 
                      y = c(1,2),
                      from = 0,
                      to = 1,
                      f_max = 3,
                      n = 1000) |>  # sample 1000 fossils to reduce noise
  hist(breaks = seq(from = t_min, to = t_max, length.out = 30),
       main = "Fossil abundance \nBaseline without ecological effect",
       xlab = "Time",
       ylab = "Fossil occurrences")
StratPal::p3_var_rate(x = c(0,1),
                      y = c(1,2),
                      from = 0,
                      to = 1,
                      f_max = 3,
                      n = 1000) |>
  StratPal::apply_niche(niche_def = niche_def,
                        gc = gradient) |>
  hist(breaks = seq(t_min, t_max, length.out = 30),
       main = "Fossil abundance with ecological effect",
       xlab = "Time",
       ylab = "Fossil occurrences")
StratPal::p3_var_rate(x = c(0,1),
                      y = c(1,2),
                      from = t_min,
                      to = t_max,
                      f_max = 3,
                      n = 1000) |>
  StratPal::apply_niche(niche_def = niche_def,
                        gc = gradient) |>
  admtools::time_to_strat(adm) |>
  hist(breaks = seq(from = h_min, to = h_max, length.out = 30),
       main = "Fossil abundance with ecological effects\nin the stratigraphic domain",
       xlab = "Stratigraphic position [m]",
       ylab = "Fossil occurrences")
# As baseline
plot(fossils, 
     tree = tree,
     main = "Full fossil record without ecological effects",
     rho = 0)
fossils_after_ecol <- fossils |>
  admtools::rev_dir(ref = t_max) |>
  StratPal::apply_niche(niche_def = niche_def, gc = gradient) |>
  admtools::rev_dir(ref = t_max)
plot(fossils_after_ecol, 
     tree = tree,
     main = "FBD tree after ecological effects",
     rho = 0)
# defines taphonomic filter that removes all fossils coinciding with gaps
f <- StratPal::strat_filter(adm)
fossils_filtered <- fossils |>
  admtools::rev_dir(ref = t_max) |>
  StratPal::apply_taphonomy(identity, f) |>
  admtools::rev_dir(ref = t_max)
plot(fossils_filtered,
     tree = tree,
     main = "Fossils during gaps removed")
adm_true <- adm # "truth"
adm_assumed <- admtools::tp_to_adm(t = c(t_min, t_max),
                                  h = c(h_min, h_max),
                                  T_unit = "Myr",
                                  L_unit = "m")
plot(adm_true)
admtools::T_axis_lab()
admtools::L_axis_lab()
title(main = "True age-depth model")
plot(adm_assumed)
admtools::T_axis_lab()
admtools::L_axis_lab()
title(main = "Assumed age-depth model")
true_age <- 0.75 # pick an arbitrary age
assumed_age <- true_age |>
  admtools::time_to_strat(adm_true) |>
  admtools::strat_to_time(adm_assumed)
assumed_age
StratPal::p3(rate = 1,
             from = t_min,
             to = t_max,
             n = 1000) |>
  admtools::time_to_strat(adm_true, destructive = FALSE) |>
  admtools::strat_to_time(adm_assumed) |>
  hist(breaks = seq(from = t_min, to = t_max, length.out = 30),
       xlab = "Time [Myr]",
       main = "Inferred extinction dymanics",
       freq = FALSE,
       ylab = "Last occcurrences")
# inferred time
inferred_time <- t |> 
  admtools::time_to_strat(adm_true) |>
  admtools::strat_to_time(adm_assumed) 
# plot error
plot(x = t,
     y = t - inferred_time, # real time minus inferred time
     type = "l",
     xlab = "Time [Myr]",
     ylab = "Error [Myr]")
# no error as reference
lines(x = c(t_min, t_max),
      y = c(0,0),
      lty = 3)
fossils_misassigned_age <- fossils |>
  admtools::rev_dir(ref = t_max) |>
  admtools::time_to_strat(adm_true) |>
  admtools::strat_to_time(adm_assumed) |>
  admtools::rev_dir(ref = t_max)
plot(fossils_misassigned_age,
     tree = tree,
     main = "Tree with wrong reconstructed fossil ages")

fossils_misassigned_age <- fossils |>
  admtools::rev_dir(ref = t_max) |>
  admtools::time_to_strat(adm_true) |>
  admtools::strat_to_time(adm_assumed) |>
  admtools::rev_dir(ref = t_max)
plot(fossils_misassigned_age,
     tree = tree,
     main = "Tree with wrong reconstructed fossil ages")
