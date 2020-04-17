# black rockfish selectivity at age
# ben.williams@alaska.gov

# load ----
source(here::here("R/helper.r"))

# data ----

# note if the lengths are not rounded there are occasionaly more precise estimates
# these will not allow the model to converge
# combine sport and comm age data
plus_group <- 30

# data can be evaluated for commercial (comm), or sport or both
# select one of the these combinations

# Female only data for commercial fisheries in Afognak, Eastside, Northeast, Southeast
# mimics Excel spreadsheet from Scott Meyer

read_csv(here::here("data/docksideAWLSM.csv"), guess = 50000) %>% 
  rename_all(tolower) %>% 
  dplyr::select(date = sample_date, Area = mgmt_area, Section = section,
                gear = gear_code, sex, length, maturity, weight, age = final_age) %>% 
  mutate(date = lubridate::mdy(date),
         year = lubridate::year(date),
         Year = factor(year),
         Age = factor(age),
         Sex = case_when(sex==1 ~ "male",
                         sex==2 ~ "female"),
         weight = ifelse(weight==0, NA, weight),
         age = ifelse(age>plus_group, plus_group, age)) %>% 
  filter(age>0, 
         !is.na(Sex),
         Section %in% c("Afognak", "Eastside", "Northeast", "Southeast"), sex==2) %>% 
  dplyr::select(Section, sex, length, weight, age, year, Year, Age, Sex) -> brf


# estimate length/weight relationship
# a * length^b

lw <- unname(lm(log(weight) ~ log(length), data = brf)$coef)
lw_int = (lw[1])
lw_slope = (lw[2])

m_lw_int = (lw[1])
m_lw_slope = (lw[2])


data.frame(length = 1:70) %>% 
  mutate(weight = exp(lw_int + lw_slope * log(length)))  %>% 
  ggplot(aes(length, weight)) + 
  geom_line() + 
  geom_point(data = brf, alpha = 0.2)

# maturity at age from Kodiak study 
# A50 - 10.47 years may not apply...
mat_int = -7.521637
mat_slope = 0.717806

data.frame(age = 1:plus_group) %>% 
  mutate(maturity = exp(mat_int + mat_slope * age) / (1 + exp(mat_int + mat_slope * age)))  %>% 
  ggplot(aes(age, maturity)) + 
  geom_line() 

# von b
vonb <- read_csv(here::here("output/vonb.csv"))
Linf <- vonb$value[vonb$param=="f_Linf"]
kappa <- vonb$value[vonb$param=="f_kappa"]
t0 <- vonb$value[vonb$param=="f_t0"]

brf %>% 
  mutate(fit = Linf * (1 - exp(-kappa * (age - t0)))) %>% 
  ggplot(aes(age, length)) + 
  geom_point() + 
  geom_line(aes(y = fit)) +
  expand_limits(y = 0, x = 0)

tibble(age = 0:plus_group) %>% 
  mutate(length = Linf * (1 - exp(-kappa * (age - t0))),
         length = ifelse(length<0, 0.01, length),
         weight = exp(lw_int + lw_slope * log(length)),
         mature = 1 / (1 + exp(mat_int + mat_slope * age)) * exp(mat_int + mat_slope * age)) -> ins  

# potential M
4.899 * seq(30,50,5) ^-0.916
4.118 * kappa^ 0.73 * Linf^-0.33

# clean data ----
brf %>% 
  dplyr::select(age, year) %>% 
  mutate(total_n = n()) %>% 
  group_by(year) %>% 
  mutate(annual_n = n()) %>% 
  group_by(age, year) %>% 
  mutate(n = n()) %>% 
  ungroup %>% 
  mutate(prop = range01(n / annual_n) * annual_n) %>% 
  group_by(age) %>% 
  summarise(prop = sum(prop) / mean(total_n)) %>% 
  mutate(prop = range01(prop)) %>% 
  left_join(data.frame(age = 0:plus_group), .) %>% 
  mutate(prop = replace_na(prop, 0)) -> select_dat


# TMB model ----
setwd(here::here("TMB"))

compile("target_spr.cpp")
dyn.load(dynlib("target_spr"))

data = list(ages = select_dat$age,
            paaC = select_dat$prop,
            mu_M = 0.18,
            sd_M = 0.05,
            mu_F = 0.1,
            sd_F = 0.05,
            mu_Fspr = 0.1,
            sd_Fspr = 0.05,
            mu_mu = 8,
            sd_mu = 0.25,
            mu_ups = 1.5,
            sd_ups = .05,
            laa = ins$length,
            maa = ins$mature,
            waa = ins$weight,
            target_spr = 0.50)


params = list(logM = log(0.183),	
              logF = log(0.08),
              logmu = log(8),		    
              logupsilon = log(1.3),
              logsigR = 0.02,
              logFspr = log(0.07))

map = list()


# parameter bounds

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

rep <- sdreport(model)
(ss = summary(rep, select = "all", p.value = T))


propC <- model$report()$propC
(M <- model$report()$M)
(F <- model$report()$F)
Fa <- model$report()$Fa
(Fspr <- model$report()$Fspr)
saC <- model$report()$saC
(mu <- model$report()$mu)
(upsilon <- model$report()$upsilon)
Ca <- model$report()$Ca 
Ua <- model$report()$Va
Na <- model$report()$Na
(spr <- model$report()$spr)
(spr_nll <- model$report()$spr_nll)
# data.frame(saf = saC, 
#            M = M, 
#            F = F) %>% 
#   write_csv(here::here("data/select.csv"))


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