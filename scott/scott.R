# Replicate Excel spreadsheet developed by Scott Meyer
# Meyer_SPR_Kodiak black rf_March_2018.xlsx 

# load ----
source("R/helper.R")

# data ----

# observed commercial and sport fishery age data from the spreadsheet
# the plus group for Scott's analysis was 50
plus_group <- 50

# data can be evaluated for commercial (comm), or sport or both
# select one of the these combinations

fisheries <- c("comm")
# fisheries <- c("sport")
# fisheries <- c("comm", "sport")

# data are weighted by sample size by age/year
# then set as a proportion that sums to 1

read_csv(here::here("scott/brf_age_prop.csv")) %>%
  filter(year>2006, fishery %in% fisheries) %>%
  mutate(age = ifelse(age>plus_group, plus_group, age)) %>% 
  uncount(count) %>% 
  mutate(tot_n = n()) %>% 
  group_by(year) %>% 
  mutate(annual_n = n()) %>% 
  group_by(age, year) %>% 
  mutate(age_n = n()) %>% 
  group_by(year, age) %>% 
  summarise_each(list(mean)) %>%
  group_by(age) %>% 
  summarise(prop = sum(age_n) / mean(tot_n)) %>% 
  left_join(data.frame(age = 0:plus_group), .) %>% 
  mutate(prop = replace_na(prop, 0)) -> dat

# inputs from spreadsheet

# estimate length/weight relationship
# a * length^b

lw_int = -10.37118
lw_slope = 2.84756

# maturity at age from Kodiak study 
# A50 - 10.47 years 

mat_int = -7.521637
mat_slope = 0.717806

# length at age 
# 2 parameter von Bertalanffy
# Worton and Rosenkrantz 2003

Linf =	54.6525
kappa =	0.178
Lzero =	7

# generate data for input to model based upon above parameters

tibble(age = 0:plus_group) %>% 
  mutate(length = Linf -(Linf - Lzero) * exp(-kappa * age),
         length = ifelse(length<0, 0.01, length),
         weight = exp(lw_int + lw_slope * log(length)),
         mature = 1 / (1 + exp(mat_int + mat_slope * age)) * exp(mat_int + mat_slope * age)) -> ins  

# TMB model ----
setwd(here::here("scott"))

# must first compile the C++ code (only need done once)
compile("scott.cpp")
# load the compiled code
dyn.load(dynlib("scott"))

# data input for the C++ model

data = list(ages = dat$age, 
            paaC = dat$prop,
            mu_M = 0.123,    # Prior on M
            sd_M = 0.05,     # Prior on M sd
            mu_F = 0.1,      # Prior on F
            sd_F = 0.05,     # Prior on F sd
            mu_mu = 8.0,    # Prior on age at fully selected
            sd_mu = 0.05,  # Prior on age at fully selected sd
            mu_ups = 1.4,   # Prior on selectivity slope
            sd_ups = 0.05, # prior on selectivity slope sd
            laa = ins$length,
            maa = ins$mature,
            waa = ins$weight)

# starting values for parameters
# note - parameters in log space
params = list(logM = log(0.123),	
              logF = log(0.10),
              logmu = log(8.03),		    
              logupsilon = log(1.43),
              logsigR = 0.02)

# map is used to "fix" parameters 
# this does not fix any parameters
# map = list()

# this map would fix M at whatever starting value was chosen
map = list(logM = factor(NA))


# parameter bounds
# L <- c(logM = log(0.04),
#        logF = log(.02),
#        logmu = log(7),
#        logupsilon = log(.03),
#        logsigR = log(0.0001))
# 
# U <- c( logM = log(0.3),
#         logF = log(.4),
#         logmu = log(12),
#         logupsilon = log(4),
#         logsigR = log(10))

L <- c(
       logF = log(.02),
       logmu = log(7),
       logupsilon = log(.03),
       logsigR = log(0.0001))

U <- c(
       logF = log(.4),
       logmu = log(12),
       logupsilon = log(4),
       logsigR = log(10))

# build model
model <- MakeADFun(data = data, 
                   parameters = params, 
                   DLL="scott", 
                   map = map)


# optimize the model
fit <- nlminb(model$par, 
              model$fn, 
              model$gr,
              lower = L,
              upper = U)

best <- model$env$last.par.best
rep <- sdreport(model)
(ss = summary(rep, select = "all", p.value = T))


propC <- model$report()$propC
(M <- model$report()$M)
(F <- model$report()$F)
Fa <- model$report()$Fa
saC <- model$report()$saC
(mu <- model$report()$mu)
(upsilon <- model$report()$upsilon)
Ca <- model$report()$Ca 
Ua <- model$report()$Va
Na <- model$report()$Na
(spr <- model$report()$spr)

tibble(age = 0:plus_group) %>% 
  mutate(length = Linf -(Linf - Lzero) * exp(-kappa * age),
         length = ifelse(length<0, 0.01, length),
         weight = exp(lw_int + lw_slope * log(length)),
         mature = 1 / (1 + exp(mat_int + mat_slope * age)) * exp(mat_int + mat_slope * age),
         unfished = Ua * weight * mature,
         fished = Na * weight * mature,
         Ca = Ca * weight,
         expN = Na * saC,
         tot_bio = Na * weight,
         exp_bio = expN * weight,
         Fa = Fa,
         Sa = saC) -> report

report %>% 
  dplyr::select(age, Sa, mature) %>% 
  pivot_longer(-age, names_to = "measure", "value") %>% 
  ggplot(aes(age, value, color = measure)) + 
  geom_line() +
  scale_color_viridis_d(name = "", end = 0.75) +
  theme(legend.position = c(0.8, 0.2))


report %>% 
  ggplot(aes(age, propC)) + 
  geom_point() +
  geom_bar(aes(y = dat$prop), stat = "identity", alpha = 0.3) 

