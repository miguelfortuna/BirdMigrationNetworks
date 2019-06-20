library("car")
library("tidyverse")

# read data ----
df <- read_csv("sampling_effort.csv")

# linear model ----
m <- lm(birds ~ ringers * ibas, data = df)

summary(m)
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  -3059.8579  6074.1224  -0.504    0.619  
#ringers         17.9356    15.8305   1.133    0.268  
#ibas           160.2727    80.3233   1.995    0.057 .
#ringers:ibas    -0.1016     0.1194  -0.851    0.403  

Anova(m, type = 2)
#Anova Table (Type II tests)
#Response: birds
#Sum Sq Df F value   Pr(>F)   
#ringers       147128746  1  0.6089 0.442542   
#ibas         2039566211  1  8.4404 0.007573 **
#ringers:ibas  174925857  1  0.7239 0.402944   
#Residuals    6041115274 25  
