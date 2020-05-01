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
         # If you only wanted to look at certain years:
         filter(year > 2012) %>% 
         # Only keeping females!
         sex==2) %>% 
  # select retains specific columns
  dplyr::select(Section, sex, length, weight, age, year, Year, Age, Sex) -> brf

# catch by district used to weight age comps

# whole weight in lb
catch <- read_csv("data/district_year_catch.csv", guess_max = 100000) %>% 
  rename(Section = district)

# Allomtry ----

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

# Maturity ----

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

# Growth ----

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

# Create a look-up tables of weighting factors for each area-year combo that is
# catch normalized (scaled between 0 and 1) 

areayear_weights <- catch %>% 
  # subsets catch in years when we have age data
  filter(year %in% unique(brf$year)) %>% 
  mutate(weight = scales::rescale(catch)) %>%
  # - removes column catch
  select(-catch)

# Age composition data
select_dat <- brf %>% 
  # Creates a table of counts of ages by year, Section, and age
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
  mutate(prop = replace_na(prop, 0)) %>% # replace NAs in prop column with 0
  arrange(age) # sorts by age

# check that comp sums to 1
sum(select_dat$prop)

# plot of combined comp
ggplot(select_dat, aes(x = age, y = prop)) +
  geom_bar(stat = "identity",
           position = "dodge") 

# TMB (Template Model Builder) model ----
setwd(here::here("TMB"))

compile("target_spr.cpp") # if model compiles, "0" will print to the console
dyn.load(dynlib("target_spr")) # links compiled cpp to R

# Prepping data for TMB, which accepts only lists for data, parameters, and
# upper/lower bounds
data = list(ages = select_dat$age, # ages
            paaC = select_dat$prop, # proportions at age, "catch comps"
            
            # B Williams selected these priors. Please contact him for more information.
            
            # Normal prior on estimated natural mortality Normal ~ (mu, sd)
            mu_M = 0.18, 
            sd_M = 0.05,
            # Normal prior on estimated fishing mortality Normal ~ (mu, sd)
            mu_F = 0.1,
            sd_F = 0.05,
            # Normal prior on Fspr (estimated fishing mortality rate at target SPR)
            mu_Fspr = 0.1,
            sd_Fspr = 0.05,
            # Priors for logistic selectivity parameters mu (a50, age at 50%
            # selectivity) and upsilson (rate or scale), Normal ~ (mu, sd)
            mu_mu = 8,
            sd_mu = 0.25,
            mu_ups = 1.5,
            sd_ups = 0.05,
            
            # Additional data inputs 
            laa = ins$length, # length-at-age, this model doesn't use length-at-age
            maa = ins$mature, # maturity-at-age
            waa = ins$weight, # weight-at-age (kg)
            
            # Target SPR - management decision. Using 0.5 as an example, but
            # it's subject to change
            target_spr = 0.50)

# Parameters that are estimated in the model and their associated starting
# values. They are all in log space - log space is used for parameters that are
# also going to be positive. Helpful for model stabilization. Jittering or
# making subtle changes to these starting values can help you test the stability
# of this model. If small changes in starting values lead to large changes in
# results, you've likely found a local minimum (= false convergence)
params = list(logM = log(0.183),	# M 
              logF = log(0.08),   # Fishing mortality 
              logmu = log(8),		  # mu = a50 selectivity
              logupsilon = log(1.3), # upsilon = rate that asymptotic selectivity is reached
              # standard deviation in the age compositions, which Ben has
              # assumed are lognormally distributed. Potential area for future
              # development to move towards a multinomial likelihood.
              logsigR = 0.02,     
              # Fspr is the estimated fishing mortality rate at your target SPR
              # (e.g. SPR = 0.5)
              logFspr = log(0.07))

map = list() # map is used to fix parameters - TMB documentation for map

# upper and lower parameter bounds. These lists need to be in the same order as
# the parameter list. When you're examining parameter estimates, make sure they
# aren't hitting bounds. If a parameter converges on a bound, it is a sign that
# the model has not convergences.
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

# build model (AD = autodifferentiation, makes estimation of nonlinear models
# super quick)
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

# Results ----

# Output of parameter estimates and maximum gradient component (mgc), which
# should be < 0.001, preferrably < 0.00001. If mgc is high or the rep says
# "Hessian not positive definite" or parameter estimates are hitting a bound or
# standard errors are NaNs, DO NOT USE MODEL OUTPUT. It means the model has not
# converged!

# In nonlinear models with lots of parameters, there are multiple levels of
# convergence: (1)  Data to fit well (2)  Stay within acceptable parameters
# space defined by upper and lower bounds (3)  Gradient has to do with how
# informative the data are to the estimates. If gradient is steep (low), it
# means the estimates are well-informed
rep <- sdreport(model)
rep 

# Summary table of estimated parameters and SEs (all are in log space)
(ss = summary(rep, select = "all", p.value = TRUE))

# use model$report() to extract different estimated or derived quantities from
# the model
propC <- model$report()$propC # catch age comps
(M <- model$report()$M) # natural mortality estimate
(F <- model$report()$F) # current fishing mortality estimate
(Fa <- model$report()$Fa) # fishing mortality at age (F * selectivity-at-age)
(Fspr <- model$report()$Fspr) # estimated fishing mortality at target SPR 
(saC <- model$report()$saC) # selectivity-at-age
(mu <- model$report()$mu) # a50: age at 50% selectivity
(upsilon <- model$report()$upsilon) # selectivity slope
(Ca <- model$report()$Ca) # catch-at-age in numbers
(Ua <- model$report()$Va) # unfished numbers-at-age (theoretical)
(Na <- model$report()$Na) # fished numbers-at-age (at current F)
(spr <- model$report()$spr) # current spr estimate 

# Likelihood compenents
(priors <- model$report()$priors) # likelihood component related to priors
# your comp_nll should be the largest absolute likelihood component
(comp_nll <- model$report()$comp_nll) # likelihood component related to fit of age comps 
(spr_nll <- model$report()$spr_nll) # likelihood component related to spr penalty
(tot_nll <- model$report()$tot_nll) # likelihood component related to spr penalty

# data.frame(saf = saC,
#            M = M,
#            F = F) %>%
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

# Comparison of selectivity and maturity
report %>% 
  dplyr::select(age, selectivity, mature) %>% 
  pivot_longer(-age, names_to = "measure", "value") %>% 
  ggplot(aes(age, value, color = measure)) + 
  geom_line() +
  scale_color_viridis_d(name = "", end = 0.75) +
  theme(legend.position = c(0.8, 0.2))

# Relative catch under current estimated F compared to estimated F under SPR
# target. If SPR target is an OFL, you would want to be way below the target. If
# SPR target is a true management target (e.g. ABC/ACL), then you might want to
# increase or decrease based on the percent different between these values.
sum(report$Ca) # current catch
(catch_target_spr <- model$report()$catch_target_spr) # catch under target spr

# Fit to the age comps: a couple issues to look for include over or
# underestimation over a several sequential ages
report %>% 
  ggplot(aes(age, propC)) + 
  geom_point() +
  geom_bar(aes(y =select_dat$prop), stat = "identity", alpha = 0.3) 

# Plot of unfished (black) relative to fished (blue)
report %>% 
  ggplot(aes(age, unfished)) + 
  geom_line() +
  geom_line(aes(y = fished), col = 4)

