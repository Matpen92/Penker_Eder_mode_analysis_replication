# header ---------------------------
#
# Author: Matpen
#
#
# Date Created: 2022-28-09
#
#
# Script Name: mode-effect analysis
#
#
# Script Description: apply bayesian methods (CFA lin.-reg.) to detect 
# mixed-mode effect on two constructs of the ISSP-envionmental module
#
#
# inputs: .sav file --> ssö data 
# 
# 
# outputs two .doc files including regression tabels and two .png files 
#  showing posterior densities
# 
#
#
# Notes:
#
#either use renv::restore() to restore packages inlcuded in lockfile or
#or if error --> manually install list of packages using code below
#
## Packages  -------------------------

renv::restore()


packages <- c("tidyverse", "rstan", "rstanarm", "haven",
              "lavaan", "blavaan", "patchwork",
              "bayesplot", "sjPlot", "janitor", 
              "tidybayes", "table1", "sjmisc")


lapply(packages, install.packages , character.only = TRUE)

lapply(packages, library, character.only = TRUE)

## Packages  -------------------------

library(dplyr)
library(rstan)
library(brms)
library(haven)
library(lavaan)
library(janitor)
library(rstanarm)
library(blavaan)
library(bayesplot)
library(blavaan)
library(sjPlot)
library(tidybayes)
library(table1)
library(sjmisc)
library(flextable)
library(DiagrammeR)

# data --------------------------------------------------------------------

issp <-
  read_sav("data/ISSP_env2021_cleaned.sav") %>% 
  clean_names(.)

# table one ---------------------------------------------------------------

issp$SurveyMode <- factor(issp$quelle, labels = c("CATI", "CAPI"))


issp <-
  issp %>% 
  mutate(d27_rec = case_when(
    d27 == 1 ~ 1,
    d27 == 2 ~ 2,
    d27 == 3 ~ 3,
    d27 == 4 ~ 4,
    d27 == 5 ~ 5,
    d27 >=  6 ~ 6,
    is.na(d27) ~ 7))#


issp <- 
  issp %>% 
  mutate(d3_rec = case_when(
    d3 <= 5 ~ 1,
    d3 > 5 & d3 < 9 ~ 2,
    d3 >= 9 ~ 3 ) )


issp$d27_rec <-  factor(issp$d27_rec, 
                        labels = c("ÖVP", "SPÖ", "FPÖ" ,
                                   "Green Party", "NEOS", 
                                   "Other Party", "No Answer"))
issp$d3_rec<-  factor(issp$d3_rec, 
                      labels = c("< Secondary Education ", "Secondary Education",#
                                 "Tertiary Education" ))

label(issp$age) <- "Age"
label(issp$d1)  <- "Gender"
label(issp$d27_rec) <- "Vote last Election"
label(issp$b15) <- "Political Interest (high-low)"
label(issp$d35) <- "Personal Net-Income"
label(issp$d28) <- "Left-Right Self-Assesment (left-right)"
label(issp$d3_rec) <- "Highest completed degree of education"


t1 <- issp %>% 
  table1(~ to_factor(d1) + age + to_factor(d3_rec)  + d35 + b15 + to_factor(d27_rec)+ d28  | SurveyMode, data =.)

t1flex(t1) %>% 
  flextable::save_as_docx(path="output/tables/table1_paper.docx")

# mean table --------------------------------------------------------------

invertItem <- function(items, maxvalue = 5, minvalue = 1){
  if(!is.numeric(maxvalue)) stop("maxvalue must be numeric!")
  if(!is.numeric(minvalue)) stop("minvalue must be numeric!")
  if(minvalue >= maxvalue) warning("Warning in invertItem: minvalue >= maxvalue")
  return((maxvalue+minvalue) - items)
}

issp <- 
  issp %>% 
  mutate_at(c("c11_1", "c11_2", "c11_3"),invertItem)


issp %>% 
  select(c5_1:c5_4,c11_1:c11_3, quelle) %>% 
  table1(~.  | to_factor(quelle),
         data = .)  %>% 
  table1::t1flex(.) %>% 
  width(.,j = 1, 3,  unit = "cm") %>% 
  flextable::save_as_docx(path="output/tables/mean_table_paper.docx")


# Bayesian sem insti --------------------------------------------------------

mod1 <- ' 

trust =~  c5_1 + c5_2 + c5_3  + c5_4

'

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


null.model <- c(paste0("c5_", 1:4, " ~~ c5_", 1:4), paste0("c5_", 1:4, " ~ 1"))

fit00 <-
  bsem(null.model,
       n.chains = 4,
       burnin = 1500,
       sample = 1500,
       target = "stan",
       seed = 123,
       data= issp)

fit1 <-
  bsem(mod1,
       n.chains = 4,
       burnin = 1500,
       sample = 1500,
       target = "stan",
       seed = 123,
       group = "SurveyMode",
       data= issp)

fit2 <-
  bsem(mod1,
       n.chains = 4,
       burnin = 1500,
       sample = 1500,
       target = "stan",
       seed = 123,
       group = "SurveyMode",
       group.equal = c("loadings"),
       data= issp)

fit3 <- 
  bsem(mod1,
       n.chains = 4,
       burnin = 1500,
       sample = 1500,
       target = "stan",
       seed = 123,
       group = "SurveyMode",
       group.equal = c("intercepts", "loadings"),
       data= issp)


blavFitIndices(fit1, baseline.model = fit00, pD = "dic") %>% 
  summary(.,  central.tendency = c("mean","median"), prob = .95)

blavFitIndices(fit2, baseline.model = fit00, pD = "dic") %>% 
  summary(.,  central.tendency = c("mean","median"), prob = .95)

blavFitIndices(fit3, baseline.model = fit00, pD = "dic") %>% 
  summary(.,  central.tendency = c("mean","median"), prob = .95)

fitMeasures(fit1)
fitMeasures(fit2)
fitMeasures(fit3)

# bayesian sem will -------------------------------------------------------

mod1 <- ' 

willing =~ c11_1 + c11_2 + c11_3 

'


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

null.model <- c(paste0("c11_", 1:3, " ~~ c11_", 1:3), paste0("c11_", 1:3, " ~ 1"))


fit0 <-
  bsem(null.model,
       n.chains = 4,
       burnin = 1500,
       sample = 1500,
       target = "stan",
       seed = 123,
       data= issp)

fita <-
  bsem(mod1,
       n.chains = 4,
       burnin = 1500,
       sample = 1500,
       target = "stan",
       seed = 123,
       group = "SurveyMode",
       data= issp)

fitb <-
  bsem(mod1,
       n.chains = 4,
       burnin = 1500,
       sample = 1500,
       target = "stan",
       seed = 123,
       group = "SurveyMode",
       group.equal ="loadings",
       data= issp)

fitc <- 
  bsem(mod1,
       n.chains = 4,
       target = "stan",
       seed = 123,
       burnin = 1500,
       sample = 1500,
       group = "SurveyMode",
       group.equal = c("intercepts", "loadings"),
       data= issp)

blavFitIndices(fita, baseline.model = fit0, pD = "dic") %>% 
  summary(.,  central.tendency = c("mean","median"), prob = .95)

blavFitIndices(fitb, baseline.model = fit0, pD = "dic") %>% 
  summary(.,  central.tendency = c("mean","median"), prob = .95)

blavFitIndices(fitc, baseline.model = fit0, pD = "dic") %>% 
  summary(.,  central.tendency = c("mean","median"), prob = .95)

fitMeasures(fita)
fitMeasures(fitb)
fitMeasures(fitc)

# observed values test willigness ----------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
# 
# invertItem <- function(items, maxvalue = 5, minvalue = 1){
#   if(!is.numeric(maxvalue)) stop("maxvalue must be numeric!")
#   if(!is.numeric(minvalue)) stop("minvalue must be numeric!")
#   if(minvalue >= maxvalue) warning("Warning in invertItem: minvalue >= maxvalue")
#   return((maxvalue+minvalue) - items)
# }


issp <- 
  issp %>% 
  rowwise() %>% 
  mutate(will = mean(c(c11_1 , c11_2 , c11_3), na.rm = T),
         trust = mean( c(c5_1 , c5_2 , c5_3 , c5_4) , na.rm = T),
         d27_rec = case_when(d27 == 1 ~ 1,
                             d27 == 2 ~ 2,
                             d27 == 3 ~ 3,
                             d27 == 4 ~ 4,
                             d27 == 5 ~ 5,
                             d27 >  5 ~ 6,
                             is.na(d27) ~ 7))

###

mod1 <- 
  issp %>% 
  stan_glm(trust ~  quelle,
           family = gaussian(link = "identity"),
           seed = 12345,
           prior = normal(0, 2),
           chains = 4,
           data = .) 


mod2 <- 
  issp %>% 
  stan_glm(will ~  quelle,
           family = gaussian(link = "identity"),
           seed = 12345,
           prior = normal(0, 2),
           chains = 4,
           data = .) 


mod3 <- 
  issp %>% 
  stan_glm( trust ~  quelle +  as.factor(d27_rec)   + d28  ,
            family = gaussian(link = "identity"),
            seed = 12345,
            prior = normal(0, 2),
            chains = 4,
            data = .) 


mod4 <-
  issp %>% 
  stan_glm( will ~  quelle +  as.factor(d27_rec)   + d28   ,
            family = gaussian(link = "identity"),
            seed = 12345,
            prior = normal(0, 2),
            chains = 4,
            data = .) 

# Outputs -----------------------------------------------------------------

tab_model(mod3, mod1,  show.intercept = F, 
          string.pred = "Coefficients",
          pred.labels = c( " Survey mode CAPI (CATI=ref.)"),
          dv.labels = c("Institutional Trust", "Willigness to Sacrifice"),
          file = "output/tables/only_mode_effect.doc")


tab_model(mod1, mod3, mod2, mod4, 
          show.intercept = F, 
          collapse.ci = T,
          string.pred = "Coefficients" ,
          pred.labels = c( "Survey mode CAPI (CATI=ref.)", 
                           "SPÖ", "FPÖ", "Green Party",
                           "NEOS", "Other Party", "No Answer",
                           "Left-Right Self-Assesment (left-right)"),
          file = "output/tables/will_mode_and_controls.doc")

tab_model(mod1, mod3, mod2, mod4, 
          show.intercept = F, 
          collapse.ci = T,
          string.pred = "Coefficients" ,
          pred.labels = c( "Survey mode CAPI (CATI=ref.)", 
                           "SPÖ", "FPÖ", "Green Party",
                           "NEOS", "Other Party", "No Answer",
                           "Left-Right Self-Assesment (left-right)"))


p1 <-
  plot(mod3,  pars = "quelle" ,  plotfun = "hist") +
  xlab("Institutional Trust") +
  theme(axis.title.x  = element_text(family = "Comic Sans MS") )+
  theme_blank()


p2 <- 
  plot(mod1, pars = "quelle" , plotfun = "hist") + 
  xlab("Willigness to Sacrifice") +
  theme(axis.title.x  = element_text(family = "Comic Sans MS")) + 
  theme_blank()

library(patchwork)

p_expo <- 
  p1 + p2 +
  plot_annotation(title = "Posterior Distributions of Mode Effects (CATI = ref.)",
                  theme = theme(plot.title = element_text(size = 14,
                                                          family = "Comic Sans MS")))

ggsave(filename = "Posteriors_modes.png",
       plot = p_expo,
       path = "output/plots/")


# Tidy draws and overlapping densities ------------------------------------

df1 <- 
  mod2 %>% 
  spread_draws(quelle) %>% 
  select(quelle) %>% 
  mutate(marker = rep(1, 4000))

df2 <- 
  mod4 %>% 
  spread_draws(quelle) %>% 
  select(quelle) %>% 
  mutate(marker = rep(2, 4000))

df3 <- full_join(df1, df2)

cols <- c("1" = "#E69F00", "2" = "#0072B2")

p3 <- 
  df3 %>% 
  ggplot(aes(x= quelle, fill = as.factor(marker))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols, 
                    labels= c("Unadjusted Mode Differences", "Adjusted Mode Differences")) +
  theme_blank() +
  theme(legend.title = element_blank()) + 
  labs(x = "")  






