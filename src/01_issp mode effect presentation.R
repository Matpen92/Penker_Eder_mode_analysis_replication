# header ---------------------------
#
# Author: Matpen
#
#
# Date Created: 2022-05-30
#
#
# Script Name:
#
#
# Script Description:
#
#
# inputs
# 
# 
# outputs
#
#
#
# Notes:
#   
#
## Packages  -------------------------

packages <- c("tidyverse", "rstan", "rstanarm", "haven",
              "lavaan", "blavaan", "patrchwork",
              "bayesplot", "sjPlot", "janitor", 
              "tidybayes", "table1", "sjmisc")


lapply(packages, install.packages , character.only = TRUE)

lapply(packages, library, character.only = TRUE)


# data --------------------------------------------------------------------

issp <-
  read_sav("data/ISSP_env2021_cleaned.sav") %>% 
  clean_names(.)

# data preparation one --------------------------------------------------------


issp$SurveyMode <- factor(issp$quelle, labels = c("CATI", "CAPI"))

issp$d27_rec <-  factor(issp$d27_rec, 
                        labels = c("ÖVP", "SPÖ", "FPÖ" ,
                                   "Green Party", "NEOS", 
                                   "Other Party", "No Answer"))

invertItem <- function(items, maxvalue = 5, minvalue = 1){
  if(!is.numeric(maxvalue)) stop("maxvalue must be numeric!")
  if(!is.numeric(minvalue)) stop("minvalue must be numeric!")
  if(minvalue >= maxvalue) warning("Warning in invertItem: minvalue >= maxvalue")
  return((maxvalue+minvalue) - items)
}

issp <- 
  issp %>% 
  rowwise() %>% 
  mutate(will = mean( c11_1 , c11_2 , c11_3, na.rm = T),
         trust = mean( c5_1 , c5_2 , c5_3 , c5_4 , na.rm = T),
         d27_rec = case_when(d27 == 1 ~ 1,
                             d27 == 2 ~ 2,
                             d27 == 3 ~ 3,
                             d27 == 4 ~ 4,
                             d27 == 5 ~ 5,
                             d27 >  5 ~ 6,
                             is.na(d27) ~ 7),
         will_rec = invertItem(will ))


label(issp$age) <- "Age"
label(issp$d1)  <- "Gender"
label(issp$d27_rec) <- "Vote last Election"
label(issp$b15) <- "Political Interest (high-low)"
label(issp$d35) <- "Personal Net-Income"
label(issp$d28) <- "Left-Right Self-Assesment (left-right)"



# desciptive table --------------------------------------------------------


issp %>% 
  table1(~ d1 + age +  + d35 + b15 + to_factor(d27_rec)+ d28  | SurveyMode, data =.)

  
# Bayesian sem insti --------------------------------------------------------

mod1 <- ' 

trust =~ c5_1 + c5_2 + c5_3  + c5_4

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


mod1 <- 
  issp %>% 
    stan_glm( will_rec ~  quelle,
              family = gaussian(link = "identity"),
              seed = 12345,
              prior = normal(0, 2),
              chains = 4,
              data = .) 

mod2 <- 
  issp %>% 
  stan_glm( will_rec ~  quelle +  as.factor(d27_rec) + b15  + d28  ,
            family = gaussian(link = "identity"),
            seed = 12345,
            prior = normal(0, 2),
            chains = 4,
            data = .) 


mod3 <-
  issp %>% 
  stan_glm( trust ~  quelle  ,
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


tab_model(mod1, mod2, show.intercept = F, 
          string.pred = "Coefficients" ,
          pred.labels = c( "Survey mode CAPI (CATI=ref.)", 
                          "SPÖ", "FPÖ", "Green Party",
                          "NEOS", "Other Party", "No Answer",
                          "Political Interes (high - low)",
                          "Left-Right Self-Assesment (left-right"),
          file = "output/tables/will_mode_and_controls.doc")


p1 <- 
  plot(mod3,  pars = "quelle" ,  plotfun = "hist") +
  xlab("Institutional Trust") +
  theme_blank()
 
                 
p2 <- 
  plot(mod1, pars = "quelle" , plotfun = "hist") + 
  xlab("Willigness to Sacrifice") +
  theme_blank()


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
  mod1 %>% 
  spread_draws(quelle) %>% 
  select(quelle) %>% 
  mutate(marker = rep(1, 4000))

df2 <- 
  mod2 %>% 
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
  labs(x = "")  +
  ggtitle("Unadjusted vs. Adjusted CATI vs. CAPI Mode Differences: Willigness to Sacrifice")

p3

# export p3 using RSTUDIO plot-pane, other export function ignore theme
