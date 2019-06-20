### load required libraries ----
library("lme4")
library("MuMIn")
library("tidyverse")

### print package version ----
packageVersion("lme4") # 1.1.21
packageVersion("MuMIn") # 1.43.6
packageVersion("tidyverse") # 1.2.1

#################################################################################################################
# function to test overdispersion (https://rdrr.io/github/markushuff/PsychHelperFunctions/src/R/overdisp_fun.R) #
#################################################################################################################
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
#################################################################################################################


### read data ----
data_flu <- read.csv("../avian_flu_data/flu_data_working.csv") %>% mutate(is_flu = 1)
data <- read_csv("../pca_and_predictions/data_pc1.csv")
df <- left_join(data, data_flu, by = c("iba_label", "species_name")) %>%
  mutate(is_flu = replace_na(is_flu, 0))

### explore data ----
flu_incidence <- inner_join(
  inner_join(
    df %>% select(iba_label, species_name),
    data_flu %>% select(iba_label) %>% distinct,
    by = "iba_label") %>%
    group_by(species_name) %>%
    summarise(n = n()),
  data_flu %>% filter(is_flu == 1) %>%
    group_by(species_name) %>%
    summarise(n_flu = n()),
  by = "species_name") %>%
  arrange(desc(n_flu))

# save output to file (we will need it to select the top IBAs from our predictions)
write_csv(flu_incidence, "flu_incidence.csv")

### fit a generalized linear mixed-effects model using the first PC ----
df_flu <- inner_join(
  inner_join(
    df,
    data_flu %>% filter(is_flu == 1) %>% select(iba_label) %>% distinct,
    by = "iba_label"),
  data_flu %>% filter(is_flu == 1) %>% select(species_name) %>% distinct,
  by = "species_name")
# write_csv(df_flu, "df_flu.csv")

flu_glmer <- glmer(is_flu ~ 1
                   + pc1 # fixed factor (first PC)
                   + (1|species_name) # random factor
                   + (1|iba_label), # random factor
                   data = df_flu,
                   family="binomial",
                   glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(flu_glmer)

##############################################################
#Random effects:
#  Groups       Name        Variance Std.Dev.
#iba_label    (Intercept) 0.3908   0.6251  
#species_name (Intercept) 0.3010   0.5486  
#Number of obs: 441, groups:  iba_label, 35; species_name, 22
##############################################################
#Fixed effects:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -1.81017    0.23566  -7.681 1.57e-14 ***
#  pc1         -0.16220    0.07403  -2.191   0.0285 * 
##############################################################

### test overdispersion of the response variable for this binomial distribution ----
overdisp_fun(flu_glmer)
# overdispersion of the distribution that we have used: ratio = 0.802; p-value = 0.999 (there is no overdispersion)

# compare our full glmer model with one without the fixed-effects factor
anova(update(flu_glmer, .~. -pc1), flu_glmer, test = "LR")
# chi-squared = 4.57; df = 1; p-value = 0.032

### plots ----
plot(log.betweenness ~ log.degree_in, data = df_flu)
plot(log.strength_in ~ log.degree_in, data = df_flu)
plot(log.betweenness ~ log.strength_in, data = df_flu)
# positive relationships among the three predictors

ggplot(data = df_flu, aes(x = pc1, y = is_flu)) +
  geom_jitter(height = 0.05, width = 0) +
  geom_smooth(method = "glm", method.args=list(family = binomial))
