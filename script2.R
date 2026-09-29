adm_name = "data/output/columns_adm.csv" # adjust if you have a different output location
data_adm = read.csv(adm_name) # read age-depth model data
key = "proximal"
t =  data_adm[,which(grepl("Myr",names(data_adm)))] # time step
h = data_adm[,which(grepl(key,names(data_adm)))] # heights at the location with key in the nameadm = admtools::tp_to_adm(t = t,
adm = admtools::tp_to_adm(t = t,
                          h = h,
                          T_unit = "Myr",
                          L_unit = "m")
adm
summary(adm)
plot(adm,
     lwd_destr = 1, # line width destructive interval (hiatus)
     lwd_acc = 2, # line width accumulative interval
     col_destr = "blue",
     col_acc = "black")
admtools::T_axis_lab() # add time axis label
admtools::L_axis_lab() # add length axis label (stratigraphic domain)
title(main = paste("Age-depth model at ", key, "location"))
# no of hiatuses
admtools::get_hiat_no(adm)
# total duration
admtools::get_total_duration(adm)
# total thickness of sediment accumulated
admtools::get_total_thickness(adm)
# stratigraphic (in)completeness
# note this is incompleteness on the timescale of the time steps of your forward simulation
# here on the scale of your forward simulation
admtools::get_incompleteness(adm)
admtools::get_completeness(adm)
# constants for convenience
t_min = admtools::min_time(adm)
t_max = admtools::max_time(adm)
h_min = admtools::min_height(adm)
h_max = admtools::max_height(adm)
t_step = 0.01 # 10 kyr resolution
h_step = 1 # meter resolution
# plotting values
t_vals = seq(t_min, t_max, by = t_step)
h_vals = seq(h_min, h_max, by = h_step)
# Time domain plot
plot(x = t_vals,
     y = admtools::sed_rate_t(adm, t_vals),
     type = "l",
     xlab = "Time [Myr]",
     ylab = "Sediment accumulation rate [m/Myr]",
     main = "Sediment accumulation (Time domain)")
# Stratigraphic domain plot
plot(x = admtools::sed_rate_l(adm, h_vals),
     y = h_vals,
     type = "l",
     xlab = "Sediment accumulation [m/Myr]",
     ylab = "Stratigraphic position [m]",
     main = "Sediment accumulation (Stratigraphic domain)")
plot(x = admtools::condensation(adm, h_vals),
     y = h_vals,
     type = "l",
     xlab = "Condensation [Myr/m]",
     ylab = "Stratigraphic position [m]",
     main = "Sediment condensation (Stratigraphic domain)")
# times
t_example = c(0.2, 0.5)
# transform into stratigraphic positions
h_example = admtools::time_to_strat(t_example,                           adm)
# display
h_example
h_example = admtools::time_to_strat(t_example,
                                    adm,
                                    destructive = FALSE)
# display
h_example
# example stratigraphic positions
h_example = c(4,30)
t_example = admtools::strat_to_time(obj = h_example,
                                    x = adm)
t_example
# set seed for reproducibility
seed = 42
set.seed(seed)
# small constant
eps = 0.01
max_rate = 1000
# fossil occurrences peak 10-fold beteen 0.5 and 0.55 Myr 
rate = approxfun(x = c(t_min, 0.45, 0.45 + eps, 0.53, 0.53 + eps, t_max),
                 y = 100 * c(1, 1, 10, 10, 1, 1),
                 rule = 2)
plot(x = t,
     y = rate(t),
     type = "l",
     xlab = "Time [Myr]",
     ylab = "Abundance [Fossils/Myr]",
     main = "Fossil abundance",
     ylim = c(0, max_rate))
sim_fossils = StratPal::p3_var_rate(x = rate,
                                    from = t_min,
                                    to = t_max,
                                    f_max = max_rate)
hist(sim_fossils,
     breaks = seq(t_min, t_max, length.out = 30),
     xlab = "Time [Myr]",
     main = " Simulated fossil abundance in the time domain")
sim_fossils_strat = sim_fossils |> 
  admtools::time_to_strat(adm)
hist(sim_fossils_strat,
     breaks = seq(h_min, h_max, length.out = 30),
     xlab = "Stratigraphic position [m]",
     main = "Fossil abundance in the stratigraphic domain")
StratPal::p3(rate = 10,
             from = h_min,
             to = h_max) |>
  admtools::strat_to_time(x = adm) |>
  hist(breaks = seq(t_min, t_max, length.out = 30),
       xlab = "Time [Myr]",
       ylab = "Fossil abundance",
       main = "Reconstructed fossil abundance in time domain")
n_lfo = 1000 # no of last occurrences
f_lo = approxfun(x = c(t_min, t_max),
                 y = c(1,2), # double over timescale of observation
                 rule = 2)
plot(x = t,
     y = f_lo(t),
     type = "l",
     xlab = "Time [Myr]",
     ylab = "Extinction rate [-]",
     main = "Extinction rate",
     yaxt = "n",
     ylim = c(0,max(f_lo(t))))
last_occ = StratPal::p3_var_rate(x = f_lo,
                                 from = t_min,
                                 to = t_max,
                                 n = n_lfo, # condition on the required number of last occurrences
                                 f_max = 2)
hist(last_occ,
     breaks = seq(t_min, t_max, length.out = 30),
     xlab = "Time [Myr]",
     ylab = "Last occurrences [#/Myr]",
     main = "Last occurrences in the time domain")
last_occ |>
  admtools::time_to_strat(x = adm,
                          destructiv = FALSE) |>
  hist(breaks = seq(h_min, h_max, length.out = 30),
       xlab = "Stratigraphic position [m]",
       main = "Last occurrences in the stratigraphic domain")
t_res = 0.01 # temporal resolution
t_sampled = seq(from = t_min, # times where the lineage is sampled
                to = t_max,
                by = t_res)
# model parameters
sigma = 1 # variance
mu = 3 # directionality
lineage_t = StratPal::random_walk(t = t_sampled,
                                  sigma = sigma,
                                  mu = mu)
plot(lineage_t,
     xlab = "Time [Myr]",
     ylab = "Trait value",
     type = "l",
     main = "Trait evolution in the time domain")
lineage_t |>
  admtools::time_to_strat(x = adm) |>
  plot(type = "l",
       xlab = "Trait value",
       ylab = "Stratigraphic position [m]",
       main = "Trait evolution in the stratigraphic domain")
key1 = "proximal"
key2 = "midplatform"
key3 = "slope"
h1 = data_adm[,which(grepl(key1,names(data_adm)))] 
h2 = data_adm[,which(grepl(key2,names(data_adm)))] 
h3 = data_adm[,which(grepl(key3,names(data_adm)))] 
adm_proximal = admtools::tp_to_adm(t = t,
                                   h = h1,
                                   T_unit = "Myr",
                                   L_unit = "m")
adm_mid = admtools::tp_to_adm(t = t,
                              h = h2,
                              T_unit = "Myr",
                              L_unit = "m")
adm_slope = admtools::tp_to_adm(t = t,
                                h = h3,
                                T_unit = "Myr",
                                L_unit = "m")
ddc = admtools::adm_to_ddc(adm1 = adm_proximal,
                           adm2 = adm_slope)
plot(ddc,
     type = "l",
     xlab = "Stratigraphic position proximal [m]",
     ylab = "Stratigraphic postion on slope [m]",
     main = "Depth-depth curve for correlation")