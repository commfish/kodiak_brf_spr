# Black rockfish SPR and selectivity-at-age
# ben.williams@noaa.gov, jane.sullivan1@alaska.gov
# Last updated May 1, 2020

# Developed for Westward Region
# Contact: carrie.worton@alaska.gov

# Load ----
source(here::here("R/helper.r"))

# Data ----

# combine sport and comm age data
plus_group <- 30

# data can be evaluated for commercial (comm), or sport or both
# select one of the these combinations

# Female only data for commercial fisheries in Afognak, Eastside, Northeast, Southeast
# mimics Excel spreadsheet from Scott Meyer

read_csv(here::here("data/docksideAWLSM.csv"), guess = 50000) %>% 
  rename_all(tolower) %>% # rename all to lowercase 
  # select() example: new name = old name
  dplyr::select(date = sample_date, Area = mgmt_area, Section = section,
                gear = gear_code, sex, length, maturity, weight, age = final_age) %>% 
  # mutate adds new column, lubridate changes times and dates to different formats
  mutate(date = lubridate::mdy(date), # mdy = month, day, year
         year = lubridate::year(date),
         Year = factor(year), # capitilized objects are factors
         Age = factor(age),
         # case_when() similar to ifelse()
         Sex = case_when(sex==1 ~ "male", 
                         sex==2 ~ "female"),
         weight = ifelse(weight==0, NA, weight), # cleaning false 0s, making them NAs
         age = ifelse(age>plus_group, plus_group, age)) %>% 
  # Filtering is like subsetting
  filter(age>0, # only want ages greater than 0
         !is.na(Sex), # !is.na = "cannot equal NA", removing NAs sex
         # %in% similar to == except for a list or a vector of options
         Section %in% c("Afognak", "Eastside", "Northeast", "Southeast"), 
         # Only keeping females!
         sex==2) %>% 
  # select retains specific columns
  dplyr::select(Section, sex, length, weight, age, year, Year, Age, Sex) -> brf

# estimate length/weight relationship
# weight = a * length ^ b

# lm() = fits a linear model
# $ extracts coefficients or parameter estimates from the linear model
# unname() = cleans up the output
lw <- unname(lm(log(weight) ~ log(length), data = brf)$coef)
lw_int <- (lw[1])
lw_slope <- (lw[2])

data.frame(length = 1:70) %>% 
  # adding a new column for predicted weights - transformation out of log-space
  mutate(weight = exp(lw_int + lw_slope * log(length)))  %>% 
  ggplot(aes(x = length, y = weight)) + 
  geom_line() + 
  geom_point(data = brf, alpha = 0.2) # alpha makes overlapping point transluscent

# Logistic parameters for maturity at age from Kodiak study (note that these
# values may be updated in 2021): Worton, C., and D. Urban. 2005. Life history
# parameters of black rockfish in the Kodiak Area. Pages 43-53. [In]: Nearshore
# marine research in Alaska (IV): Final comprehensive progress report. NOAA
# Cooperative Agreement NA16FN2808. Alaska Department of Fish and Game, Division
# of Commercial Fisheries, Juneau.
mat_int = -7.521637
mat_slope = 0.717806

data.frame(age = 1:plus_group) %>% 
  mutate(maturity = exp(mat_int + mat_slope * age) / (1 + exp(mat_int + mat_slope * age)))  %>% 
  ggplot(aes(age, maturity)) + 
  geom_line() 

# von b - these parameter estimates were created using the vonb.R script
vonb <- read_csv(here::here("output/vonb.csv")) # reading in the results from vonb.R
Linf <- vonb$value[vonb$param=="f_Linf"] # subseting the individual parameter estimates
kappa <- vonb$value[vonb$param=="f_kappa"] # this looks high
t0 <- vonb$value[vonb$param=="f_t0"] # may want to revisit, usually t0 are negative

# Visualizes the vonB fit to the data
brf %>% 
  mutate(fit = Linf * (1 - exp(-kappa * (age - t0)))) %>% 
  ggplot(aes(age, length)) + 
  geom_point() + 
  geom_line(aes(y = fit)) +
  expand_limits(y = 0, x = 0)

# Store all model inputs into a single object that spans the ages used in the model
tibble(age = 0:plus_group) %>% # tibble() creates a data.frame()
  mutate(length = Linf * (1 - exp(-kappa * (age - t0))), # growth 
         length = ifelse(length<0, 0.01, length), # making any negative lengths that were artifacts of the vonB model = small values
         weight = exp(lw_int + lw_slope * log(length)), # weight-length relationship
         # maturity
         mature = exp(mat_int + mat_slope * age) / (1 + exp(mat_int + mat_slope * age))) -> ins  

# potential values for natural mortality (M) - you can find these and similar M
# estimators using the barefoot ecologist tool:
# http://barefootecologist.com.au/shiny_m.html
4.899 * seq(30,50,5) ^-0.916
4.118 * kappa^ 0.73 * Linf^-0.33

# Age composition data ----

# Weight catch by area and year - the goal is to make the age composition as
# representative of your catch as possible

# To illustrate this concept, I created a fake dataset of catch in earch area
# and year. You'll want to use actual harvest data instead.
fake_catch <- data.frame(expand.grid(Section = unique(brf$Section),
            year = unique(brf$year))) 
fake_catch <- fake_catch %>% 
  mutate(catch = runif(min = 20, max = 300, n = nrow(fake_catch)))
write_csv(fake_catch, path = "data/kodiak_catch_example.csv")

# Create a look-up tables of weighting factors for each area-year combo that is
# catch normalized (scaled between 0 and 1) 
areayear_weights <- fake_catch %>% 
  mutate(weight = scales::rescale(catch)) %>% 
  select(-catch)

# Age composition data
select_dat <- brf %>% 
  count(year, Section, age) %>%
  # join to area-year look-up table 
  left_join(areayear_weights, by = c("year", "Section")) %>% 
  # adjust sample size by multiplying by area-year weighting factor
  mutate(adj_n = n * weight,
         # total adjusted sample size
         total_adj_n = sum(adj_n)) %>% 
  # within each age, add up all the samples and divide by the total sample size
  # to get proportions-at-age
  group_by(age) %>%
  summarise(prop = round(sum(adj_n) / unique(total_adj_n), 5)) %>% 
  ungroup() %>% 
  # join back to full dataframe of all ages for input to the model
  full_join(data.frame(age = 0:plus_group)) %>% 
  mutate(prop = replace_na(prop, 0)) %>% 
  arrange(age)

# check that comp sums to 1
sum(select_dat$prop)

# plot of combined comp
ggplot(select_dat, aes(x = age, y = prop)) +
  geom_bar(stat = "identity",
           position = "dodge") 

# TMB model ----
setwd(here::here("TMB"))

compile("target_spr.cpp")
dyn.load(dynlib("target_spr"))

data = list(ages = select_dat$age,
            paaC = select_dat$prop,
            # Normal prior on estimated natural mortality Normal ~ (mu, sd)
            mu_M = 0.18, # 
            sd_M = 0.05,
            # Normal prior on estimated fishing mortality Normal ~ (mu, sd)
            mu_F = 0.1,
            sd_F = 0.05,
            # Normal prior on Fspr (fishing mortality rate at target SPR)
            mu_Fspr = 0.1,
            sd_Fspr = 0.05,
            # Priors for logistic selectivity parameters mu (a50) and upsilson
            # (rate/scale), Normal ~ (mu, sd)
            mu_mu = 8,
            sd_mu = 0.25,
            mu_ups = 1.5,
            sd_ups = 0.05,
            laa = ins$length, # this model doesn't use length-at-age
            maa = ins$mature, #
            waa = ins$weight,
            target_spr = 0.50)


params = list(logM = log(0.183),	
              logF = log(0.08),
              logmu = log(8),		    
              logupsilon = log(1.3),
              logsigR = 0.02,
              logFspr = log(0.07))

map = list()


# upper and lower parameter bounds

L <- c(logM = log(0.02),
       logF = log(.02), 
       logmu = log(5), 
       logupsilon = log(.07), 
       logsigR = log(0.0001),
       logFspr = log(0.02))

U <- c(logM = log(0.4),
       logF = log(.4),
       logmu = log(12), 
       logupsilon = log(5), 
       logsigR = log(10),
       logFspr = log(0.4))

# build model
model <- MakeADFun(data = data, 
                   parameters = params, 
                   DLL="target_spr", 
                   map = map)


# optimize the model
fit <- nlminb(model$par, 
              model$fn, 
              model$gr,
              lower = L,
              upper = U)

# Output of parameter estimates and maximum gradient component (mgc), which
# should be < 0.001, preferrably < 0.00001. If mgc is high or the rep says
# "Hessian not positive definite" DO NOT USE MODEL OUTPUT. It means the model
# has not converged
rep <- sdreport(model)
rep 

# Summary table of estimated parameters and SEs (all are in log space)
(ss = summary(rep, select = "all", p.value = TRUE))

# use model$report() to extract different estimated or derived quantities from
# the model
propC <- model$report()$propC # catch age comps
(M <- model$report()$M) # natural mortality estimate
(F <- model$report()$F) # fishing mortality estimate
(Fa <- model$report()$Fa) # fishing mortality at age (fully selected)
(Fspr <- model$report()$Fspr) # estimated fishing mortality at target SPR 
saC <- model$report()$saC # selectivity-at-age
(mu <- model$report()$mu) # a50: age at 50% selectivity
(upsilon <- model$report()$upsilon) # selectivity slope
(Ca <- model$report()$Ca) # catch-at-age in numbers
Ua <- model$report()$Va # unfished numbers-at-age
Na <- model$report()$Na # fished numbers-at-age (at current F)
(spr <- model$report()$spr) # spr estimate 
(catch_target_spr <- model$report()$catch_target_spr)

(priors <- model$report()$priors) # likelihood component related to priors
(comp_nll <- model$report()$comp_nll) # likelihood component related to fit of age comps
(spr_nll <- model$report()$spr_nll) # likelihood component related to spr penalty
(tot_nll <- model$report()$tot_nll) # likelihood component related to spr penalty

data.frame(saf = saC,
           M = M,
           F = F) %>%
#   write_csv(here::here("data/select.csv"))

# Report file
tibble(age = 0:plus_group) %>% 
  mutate(length = Linf * (1 - exp(-kappa * (age - t0))),
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
         selectivity = saC) -> report

report %>% 
  dplyr::select(age, selectivity, mature) %>% 
  pivot_longer(-age, names_to = "measure", "value") %>% 
  ggplot(aes(age, value, color = measure)) + 
  geom_line() +
  scale_color_viridis_d(name = "", end = 0.75) +
  theme(legend.position = c(0.8, 0.2))

sum(report$Ca)

report %>% 
  ggplot(aes(age, propC)) + 
  geom_point() +
  geom_bar(aes(y =select_dat$prop), stat = "identity", alpha = 0.3) 

report %>% 
  ggplot(aes(age, unfished)) + 
  geom_line() +
  geom_line(aes(y = fished), col = 4)

report %>% 
  summarise(spr = sum(fished) / sum(unfished))

