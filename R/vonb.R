# von Bertalanffy growth rate
# log space via TMB
# set at age plus group of 30

# load ----
source(here::here("R/helper.R"))

# dat ----

read_csv(here::here("data/docksideAWLSM.csv"), guess = 50000) %>% 
  rename_all(tolower) %>% 
  dplyr::select(date = sample_date, Area = mgmt_area, Section = section,
                gear = gear_code, sex, length, maturity, weight, age = final_age) %>% 
  mutate(date = lubridate::mdy(date),
         year = lubridate::year(date),
         Year = factor(year),
         Age = factor(age),
         Sex = case_when(sex==1 ~ "male",
                         sex==2 ~ "female")) %>% 
  filter(length>20, age>0, 
         Section %in% c("Afognak", "Eastside", "Northeast", "Southeast"), sex%in%1:2) -> brf 

filter(brf, sex == 2) -> fem
filter(brf, sex == 1) -> mal

# TMB model ----
setwd(here::here("tmb"))

compile("vonb.cpp")
dyn.load(dynlib("vonb"))


# female model ----
data = list(age = fem$age,
            length = fem$length)

# starting parameters
params <- list(logLinf = log(55), logkappa = log(0.15), t0 = -2, logSigma = 0.001)
ll <- list(logLinf = log(25), logkappa = log(0.05), t0 = -5, logSigma = log(0.0001))
ul <- list(logLinf = log(75), logkappa = log(0.45), t0 = 5, logSigma = log(10))

# build model
model <- MakeADFun(data = data, 
                   parameters = params, 
                   DLL="vonb")

# optimize the model
fit <- nlminb(model$par, 
              model$fn, 
              model$gr,
              lower = ll,
              upper = ul)

best <- model$env$last.par.best
rep <- sdreport(model)
rep

best
f_summary <- summary(rep, select = "all", p.value = T)

# female output ------------------
f_fit <- model$report()$fit
(f_Linf <- model$report()$Linf)
f_Linf_se = (f_summary[1,2])
(f_kappa <- model$report()$kapp)
f_kappa_se = (f_summary[2,2])
(f_t0 <- model$report()$t0)
f_t0_se = (f_summary[3,2])
(f_Sigma <- model$report()$Sigma)

# female figs -------------------
data.frame(age = fem$age, length = fem$length, fit = f_fit) %>% 
  ggplot(aes(age, length)) + 
  geom_point() + 
  geom_line(aes(y = fit))+
  expand_limits(x = 0, y = 0)


# male model ----
data = list(age = mal$age,
            length = mal$length)

# build model
model <- MakeADFun(data = data, 
                   parameters = params, 
                   DLL="vonb")

# optimize the model
fit <- nlminb(model$par, 
              model$fn, 
              model$gr,
              lower = ll,
              upper = ul)

rep <- sdreport(model)
m_summary <- summary(rep, select = "all", p.value = T)

# male output ------------------
m_fit <- model$report()$fit
(m_Linf <- model$report()$Linf)
m_Linf_se = (m_summary[1,2])
(m_kappa <- model$report()$kapp)
m_kappa_se = (m_summary[2,2])
(m_t0 <- model$report()$t0)
m_t0_se = (m_summary[3,2])
(m_Sigma <- model$report()$Sigma)

# male figs -------------------
data.frame(age = mal$age, length = mal$length, fit = m_fit) %>% 
  ggplot(aes(age, length)) + 
  geom_point() + 
  geom_line(aes(y = fit))+
  expand_limits(x = 0, y = 0)


# output ----------
data.frame(param = c("f_Linf", "f_Linf_se", "f_kappa", "f_kappa_se", "f_t0", "f_t0_se",
                     "m_Linf", "m_Linf_se", "m_kappa", "m_kappa_se", "m_t0", "m_t0_se"),
           value = c(f_Linf, f_Linf_se, f_kappa, f_kappa_se, f_t0, f_t0_se,
                     m_Linf, m_Linf_se, m_kappa, m_kappa_se, m_t0, m_t0_se)) %>% 
  write_csv(here::here("output/vonb.csv"))
