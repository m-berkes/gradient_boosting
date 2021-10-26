########### Libraries #######################################

# Check.Packages
# Function to simplify package install
# Requires 'packages' as array of strings using 'c'

check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# check dependencies and install if needed

packages<-c("ggplot2",  "plyr", "dplyr", "tidyverse", "car", "DescTools", "plotrix", "lsr", "e1071", "patchwork", "RColorBrewer",
            "emmeans", "ggpubr", "fastmatch", "BayesianFirstAid", "extrafont", "MatchIt", "FactoMineR", "modelr", "infer", "interactions",
            "factoextra", "twang", "gbm", "survey", "lattice", "cobalt", "xtable", "modelsummary")
check.packages(packages)

# run this after installing new fonts on system
font_import()
loadfonts()

## needed for extrafonts if saving graph to PDF
## Needed only on Windows - run once per R session
## Adjust the path to match your installation of Ghostscript - so obviously you also need to download Ghostscript
Sys.setenv(R_GSCMD = "D:/Program Files/gs9.26/bin/gswin32c.exe")


################# Set up #####################

## load file
biomarkers <- read_csv(
  file = "boosting.csv",         #change this location
  col_names=TRUE, na = "."
  )

names(biomarkers)[names(biomarkers) == "PTGENDER"] <- "SEX"
biomarkers$SEX <- factor(biomarkers$SEX)
biomarkers$LANGUSE <- factor(biomarkers$LANGUSE, labels = c("ML", "BL"))
biomarkers %>%
  distinct(DIAGNOSIS) %>%
  mutate(DIAGNOSIS = fct_relevel(DIAGNOSIS, c("CN", "MCI", "AD"))) %>%
  arrange(DIAGNOSIS)
biomarkers$DIAGNOSIS <- factor(biomarkers$DIAGNOSIS)
biomarkers$ETHCAT <- factor(biomarkers$ETHCAT)
biomarkers$RACCAT <- factor(biomarkers$RACCAT)
biomarkers$PTAU[ biomarkers$PTAU == "<8" ] <- 8
biomarkers$PTAU[ biomarkers$PTAU == ">120" ] <- 120
biomarkers$TAU[ biomarkers$TAU == "<80" ] <- 80
biomarkers$TAU[ biomarkers$TAU == ">1300" ] <- 1300
biomarkers$PTAU <- as.numeric(biomarkers$PTAU)
biomarkers$TAU <- as.numeric(biomarkers$TAU)

#adds in columns for presence of 1 or 2 APOE4 alleles
biomarkers <- biomarkers %>%
  mutate(APOE = if_else(biomarkers$APOEGEN1 == '4' | biomarkers$APOEGEN2 == '4', '1', '0'))
biomarkers <- biomarkers %>%
  mutate(APOE2 = if_else(biomarkers$APOEGEN1 == '4' & biomarkers$APOEGEN2 == '4', '1', '0'))
biomarkers$APOE <- as.factor(biomarkers$APOE)
biomarkers$APOE2 <- as.factor(biomarkers$APOE2)

# assign new numbers based on CogProfile: 0 is healthy, 1 is MCI or AD
biomarkers$CogProfileNum <- as.numeric(biomarkers$DXCURRENT)
biomarkers$CogProfileNum[ biomarkers$CogProfileNum == 7 ] <- 1
biomarkers$CogProfileNum[ biomarkers$CogProfileNum == 1 ] <- 0
biomarkers$CogProfileNum[ biomarkers$CogProfileNum >1 ] <- 1

# remove Hispanic/Latino and unknown ethnicity monolinguals
biomarkers <- biomarkers[!(biomarkers$LANGUSE == 'ML' & biomarkers$ETHCAT != '2' ), ]

# subset of only monolinguals
biomono <- subset(biomarkers, biomarkers$LANGUSE == "ML")
biomono <- subset(biomono, select = c("PID","SITEID","SEX","AGE","EDU","LANGUSE","ETHCAT", "RACCAT", "DXCURRENT","DIAGNOSIS", "ABETA40","ABETA42",
                                      "ABRATIO","PTAU","TAU","APOEGEN1","APOEGEN2","APOE", "APOE2", "CDR_GLOBAL","CDR_SUM","MMSE",
                                      "ADAS11","ADAS13","FA","AD","RD","LIMMTOTAL","LDELTOTAL","DSPANFOR", "DSPANFLTH",
                                      "DSPANBAC", "DSPANBLTH","CATANIMSC", "CATVEGESC", "TRAASCOR", "TRABSCOR",
                                      "BNTTOTAL", "ANARTERR", "MINTTOTAL","CogProfileNum"))

# subset of only monolinguals with DTI data
biomonoDTI <- subset(biomono, !is.na(FA))



################## Base descriptives monolinguals ################


ddply(biomono, .(DIAGNOSIS), summarise, 
      N = length(DIAGNOSIS),
      Fem = length(which(SEX == '2')),
      Cog = mean(CogProfileNum, na.rm = TRUE),
      Age_mean = mean(AGE),
      Age_sd = sd(AGE),
      FA_mean = mean(FA, na.rm = TRUE),
      FA_sd = sd(FA, na.rm = TRUE),
      AD_mean = mean(AD, na.rm = TRUE), 
      AD_sd = sd(AD, na.rm = TRUE),
      RD_mean = mean(RD, na.rm = TRUE),
      RD_sd = sd(RD, na.rm = TRUE),
      Edu_mean = mean(EDU, na.rm = TRUE),
      Edu_sd = sd(EDU, na.rm = TRUE),
      MMSE_mean = mean(MMSE, na.rm = TRUE),
      MMSE_sd = sd(MMSE, na.rm = TRUE),
      AB42 = mean(ABETA42, na.rm = TRUE),
      AB42_sd = sd(ABETA42, na.rm = TRUE),
      AB40 = mean(ABETA40, na.rm = TRUE),
      AB40_sd = sd(ABETA40, na.rm = TRUE),
      ABRatio = mean(ABRATIO, na.rm = TRUE),
      PTau = mean(PTAU, na.rm = TRUE),
      PTau_sd = sd(PTAU, na.rm = TRUE),
      Tau = mean(TAU, na.rm = TRUE),
      Tau_sd = sd(TAU, na.rm = TRUE),
      APOE = length(which(APOE==1)),
      APOEboth = length(which(APOE2==1))
      )



############### Monolingual GBM Analyses ########################

# PCA(X, scale.unit = TRUE, ncp = 5, graph = TRUE), requires FactoMiner package.
# Do not run PCA with missing values, it will impute values and assign to subjects - has no meaning.

for.monopca <- subset(biomonoDTI, select = c("FA", "AD", "RD"))
res.monopca <- PCA(for.monopca, graph = TRUE)
scoresmonoPCA <- res.monopca$ind$coord[,1]
biomonoDTI <- cbind(biomonoDTI, scoresmonoPCA)
biomonofull <-left_join(biomono, biomonoDTI, by = NULL, all.x = TRUE)

ddply(biomonofull, .(DIAGNOSIS), summarise, 
      PCA = mean(scoresmonoPCA, na.rm = TRUE)
      )



## propensity scores on multiple diagnosis groups - TWANG PACKAGE REQUIRED
set.seed(1)

# set to optimal iterations based on previous runs and the summary(mnps.mono$psList) command
mnps.mono <- mnps(DIAGNOSIS ~ SEX + AGE + EDU + RACCAT,
                  data = biomonofull,
                  estimand = "ATE",
                  verbose = FALSE,
                  stop.method = c("es.mean", "ks.max"),
                  n.trees = 6450)

# list of PS scores
mnps.mono$psList$CN$ps$es.mean.ATE
mnps.mono$psList$MCI$ps$es.mean.ATE
mnps.mono$psList$AD$ps$es.mean.ATE
get.weights(mnps.mono, stop.method = "es.mean", estimand = "ATE") # propensity score weights

# graphical assessments of balance
m1p3 <- love.plot(bal.tab(mnps.mono), which.treat = NULL, stars = "raw", size = 5, alpha = 0.5) #requires cobalt package
m1p1 <- plot(mnps.mono, plots = 1)
m1p2 <- plot(mnps.mono, plots = 2, subset = "es.mean")
plot(mnps.mono$psList$CN, plots = 3)#,  pairwiseMax = FALSE, figureRows = 3)
plot(mnps.mono$psList$MCI, plots = 3)#,  pairwiseMax = FALSE, figureRows = 3)
plot(mnps.mono$psList$AD, plots = 3)#,  pairwiseMax = FALSE, figureRows = 3)
plot(mnps.mono$psList$CN, plots = 4)
plot(mnps.mono$psList$MCI, plots = 4)
plot(mnps.mono$psList$AD, plots = 4)
plot(mnps.mono$psList$CN, plots = 6)
plot(mnps.mono$psList$MCI, plots = 6)
plot(mnps.mono$psList$AD, plots = 6)

# tabular assessments of balance
GBM1.bal.tab.CN <- bal.table(mnps.mono$psList$CN, digits = 3)
GBM1.bal.tab.MCI <- bal.table(mnps.mono$psList$MCI, digits = 3)
GBM1.bal.tab.AD <- bal.table(mnps.mono$psList$AD, digits = 3)
summary(mnps.mono$psList$CN) #3198 iterations
summary(mnps.mono$psList$MCI) #5138 iterations
summary(mnps.mono$psList$AD) #6405 iterations

write.table(GBM1.bal.tab.CN, file = "GBM1_bal_tab_CN.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)
write.table(GBM1.bal.tab.MCI, file = "GBM1_bal_tab_MCI.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)
write.table(GBM1.bal.tab.AD, file = "GBM1_bal_tab_AD.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)


# relative influence of predictors
summary(mnps.mono$data) # overall info
summary(mnps.mono$psList$CN$gbm.obj, n.trees=mnps.mono$psList$CN$n.trees, plot=TRUE)
summary(mnps.mono$psList$MCI$gbm.obj, n.trees=mnps.mono$psList$MCI$n.trees, plot=TRUE)
summary(mnps.mono$psList$AD$gbm.obj, n.trees=mnps.mono$psList$AD$n.trees, plot=TRUE)

# get weights
biomonofull$w <- get.weights(mnps.mono, stop.method = "es.mean")
design.mnps <- svydesign(ids=~1, weights=~w, data = biomonofull)

# examining brain PC scores based on diagnosis after re-weighting
glm1.uw <- lm(scoresmonoPCA ~ DIAGNOSIS, data = biomonofull) # unweighted comparisons
contrasts(biomonofull$DIAGNOSIS) # how diagnoses were coded
glm1 <- svyglm(scoresmonoPCA ~ DIAGNOSIS, design = design.mnps) # weighted comparisons
summary(glm1.uw)
summary(glm1) # there is no F statistic for svyglm object as the test is not fit by maximum likelihood
regTermTest(glm1, ~DIAGNOSIS) # this provides a Wald test of the hypothesis that all coefficients associated with a particular regression term are zero

# use the following if needed
# -sum(coef(glm1.uw)[-1]) # unweighted AD estimate
# -sum(coef(glm1)[-1]) # AD estimate, which is a linear combination of other two estimates

sqrt(c(-1,-1) %*% summary(glm1)$cov.scaled[-1,-1] %*% c(-1,-1)) # standard error of this estimate, i.e., AD estimate
anova(glm1.uw, glm1)

# some relationships of interest
summary(lm(MMSE ~ DIAGNOSIS, data= biomonofull))
summary(svyglm(MMSE ~ DIAGNOSIS, design = design.mnps))
regTermTest(svyglm(MMSE ~ DIAGNOSIS, design = design.mnps), ~DIAGNOSIS)


## GBM analysis looking at demographic and biomarker variables
set.seed(1)
mono.gbm2 <- subset(biomonofull, select = c("DIAGNOSIS", "SEX","AGE","EDU","ETHCAT", "RACCAT", "ABETA40","ABETA42",
                                            "ABRATIO","PTAU","TAU","APOE", "APOE2","MMSE",
                                            "ADAS11","ADAS13","scoresmonoPCA","LIMMTOTAL","LDELTOTAL","DSPANFOR", "DSPANFLTH",
                                            "DSPANBAC", "DSPANBLTH","CATANIMSC", "CATVEGESC", "TRAASCOR", "TRABSCOR",
                                            "BNTTOTAL", "ANARTERR", "MINTTOTAL"))

# using optimal iterations
mnps.mono.gbm2 <- mnps(DIAGNOSIS ~ SEX+AGE+EDU+RACCAT+ABETA42+ABRATIO+PTAU+TAU+APOE,
                       data = mono.gbm2,
                       estimand = "ATE",
                       verbose = FALSE,
                       stop.method = c("es.mean"),
                       n.trees = 1100)

# graphical assessments of balance
m2p3 <- love.plot(bal.tab(mnps.mono.gbm2), which.treat = NULL, stars = "raw", size = 5, alpha = 0.5)
m2p1 <- plot(mnps.mono.gbm2, plots = 1)
m2p2 <- plot(mnps.mono.gbm2, plots = 2, subset = "es.mean")
plot(mnps.mono.gbm2$psList$CN, plots = 3,  pairwiseMax = FALSE, figureRows = 3)
plot(mnps.mono.gbm2$psList$MCI, plots = 3,  pairwiseMax = FALSE, figureRows = 3)
plot(mnps.mono.gbm2$psList$AD, plots = 3,  pairwiseMax = FALSE, figureRows = 3)
plot(mnps.mono.gbm2$psList$CN, plots = 4)
plot(mnps.mono.gbm2$psList$MCI, plots = 4)
plot(mnps.mono.gbm2$psList$AD, plots = 4)
plot(mnps.mono.gbm2$psList$CN, plots = 6)
plot(mnps.mono.gbm2$psList$MCI, plots = 6)
plot(mnps.mono.gbm2$psList$AD, plots = 6)

# tabular assessments of balance
bal.table(mnps.mono.gbm2, digits = 3)

GBM2.bal.tab.CN <- bal.table(mnps.mono.gbm2$psList$CN, digits = 3)
write.csv(bal.table(mnps.mono.gbm2$psList$CN, digits = 3), "gbm2_CN.csv")
GBM2.bal.tab.MCI <- bal.table(mnps.mono.gbm2$psList$MCI, digits = 3)
write.csv(bal.table(mnps.mono.gbm2$psList$MCI, digits = 3), "gbm2_MCI.csv")
GBM2.bal.tab.AD <- bal.table(mnps.mono.gbm2$psList$AD, digits = 3)
write.csv(bal.table(mnps.mono.gbm2$psList$AD, digits = 3), "gbm2_AD.csv")

write.table(GBM2.bal.tab.CN, file = "GBM2_bal_tab_CN.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)
write.table(GBM2.bal.tab.MCI, file = "GBM2_bal_tab_MCI.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)
write.table(GBM2.bal.tab.AD, file = "GBM2_bal_tab_AD.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)


summary(mnps.mono.gbm2$psList$CN) #1046 iterations
summary(mnps.mono.gbm2$psList$MCI) #346 iterations
summary(mnps.mono.gbm2$psList$AD) #248 iterations

# relative influence of predictors
summary(mnps.mono.gbm2$data) # overall info
summary(mnps.mono.gbm2$psList$CN$gbm.obj, n.trees=mnps.mono.gbm2$psList$CN$n.trees, plot=TRUE)
summary(mnps.mono.gbm2$psList$MCI$gbm.obj, n.trees=mnps.mono.gbm2$psList$MCI$n.trees, plot=TRUE)
summary(mnps.mono.gbm2$psList$AD$gbm.obj, n.trees=mnps.mono.gbm2$psList$AD$n.trees, plot=TRUE)

# get weights
mono.gbm2$w <- get.weights(mnps.mono.gbm2, stop.method = "es.mean")
design.mnps.gbm2 <- svydesign(ids=~1, weights=~w, data = mono.gbm2)
# design.mnps.gbm2 <- update( design.mnps.gbm2 , DIAGNOSIS = relevel(DIAGNOSIS,"AD") ) # relevels the diagnosis factor to CN as baseline variable. Use if needed

# examining brain PC scores based on diagnosis after re-weighting
glm2.uw <- lm(scoresmonoPCA ~ DIAGNOSIS, data = mono.gbm2) # unweighted comparisons
glm2 <- svyglm(scoresmonoPCA ~ DIAGNOSIS, design = design.mnps.gbm2) # weighted comparisons
summary(glm2.uw)
summary(glm2)
regTermTest(glm2, ~DIAGNOSIS)

# as with model 1, use if needed
# -sum(coef(glm2)[-1]) #AD estimate, which is a linear combination of other two estimates
# sqrt(c(-1,-1) %*% summary(glm2)$cov.scaled[-1,-1] %*% c(-1,-1)) # standard error of this estimate

# some relationships of interest follow
summary(lm(MMSE ~ DIAGNOSIS, data = mono.gbm2))
summary(svyglm(MMSE ~ DIAGNOSIS, design = design.mnps.gbm2))
regTermTest(svyglm(MMSE ~ DIAGNOSIS, design = design.mnps.gbm2), ~DIAGNOSIS)

summary(svyglm(ABETA42 ~ APOE, design = design.mnps.gbm2)) #relation between APOE and AB42
summary(svyglm(ABETA42 ~ APOE + DIAGNOSIS, design = design.mnps.gbm2)) #relation between APOE and DIAGNOSIS on AB42
regTermTest(svyglm(ABETA42 ~ APOE, design = design.mnps.gbm2), ~APOE)

summary(svyglm(PTAU ~ APOE, design = design.mnps.gbm2)) #relation between APOE and TAU
summary(svyglm(TAU ~ APOE + DIAGNOSIS, design = design.mnps.gbm2)) #relation between APOE and DIAGNOSIS on TAU
regTermTest(svyglm(PTAU ~ APOE, design = design.mnps.gbm2), ~APOE)

ddply(mono.gbm2, .(APOE), summarise, 
      tau_mean = mean(TAU), na.rm = TRUE,
      ptau_mean = mean(PTAU), na.rm = TRUE,
      AB42_mean = mean(ABETA42), na.rm = TRUE)





############## Both language group set up #################

# all subjects, less variables
bio.all <- subset(biomarkers, select = c("PID","SEX","AGE","EDU","LANGUSE","ETHCAT", "RACCAT","DIAGNOSIS", "ABETA40","ABETA42",
                                      "ABRATIO","PTAU","TAU","APOE", "APOE2","MMSE",
                                      "ADAS11","ADAS13","FA","AD","RD","LIMMTOTAL","LDELTOTAL","DSPANFOR", "DSPANFLTH",
                                      "DSPANBAC", "DSPANBLTH","CATANIMSC", "CATVEGESC", "TRAASCOR", "TRABSCOR",
                                      "BNTTOTAL", "ANARTERR", "MINTTOTAL","CogProfileNum"))

# subset of only those with DTI data
bothDTI <- subset(biomarkers, !is.na(FA))
# select only necessary data
bothDTI <- subset(bothDTI, select = c("PID","SEX","AGE","EDU","LANGUSE","ETHCAT","RACCAT","DIAGNOSIS","ABETA40","ABETA42",
                                    "ABRATIO","PTAU","TAU","APOE","APOE2","MMSE","ADAS11","ADAS13",
                                    "FA","AD","RD","LIMMTOTAL","LDELTOTAL","DSPANFOR", "DSPANFLTH",
                                    "DSPANBAC", "DSPANBLTH","CATANIMSC", "CATVEGESC", "TRAASCOR", "TRABSCOR",
                                    "BNTTOTAL", "ANARTERR", "MINTTOTAL","CogProfileNum"))





################## Base descriptives ################

ddply(bio.all, .(LANGUSE), summarise, 
      N = length(LANGUSE),
      Fem = length(which(SEX == '2')),
      Cog = mean(CogProfileNum, na.rm = TRUE),
      CN = length(which(DIAGNOSIS == "CN")),
      MCI = length(which(DIAGNOSIS == "MCI")),
      Alz = length(which(DIAGNOSIS == "AD")),
      Age_mean = mean(AGE),
      Age_sd = sd(AGE),
      FA_mean = mean(FA, na.rm = TRUE),
      FA_sd = sd(FA, na.rm = TRUE),
      AD_mean = mean(AD, na.rm = TRUE), 
      AD_sd = sd(AD, na.rm = TRUE),
      RD_mean = mean(RD, na.rm = TRUE),
      RD_sd = sd(RD, na.rm = TRUE),
      Edu_mean = mean(EDU, na.rm = TRUE),
      Edu_sd = sd(EDU, na.rm = TRUE),
      MMSE_mean = mean(MMSE, na.rm = TRUE),
      MMSE_sd = sd(MMSE, na.rm = TRUE),
      AB42 = mean(ABETA42, na.rm = TRUE),
      AB42_sd = sd(ABETA42, na.rm = TRUE),
      AB42_N = length(which(ABETA42 != 'NA')),
      AB40 = mean(ABETA40, na.rm = TRUE),
      AB40_sd = sd(ABETA40, na.rm = TRUE),
      AB40_N = length(which(ABETA40 != 'NA')),
      ABRatio = mean(ABRATIO, na.rm = TRUE),
      PTau = mean(PTAU, na.rm = TRUE),
      PTau_sd = sd(PTAU, na.rm = TRUE),
      PTau_N = length(which(PTAU != 'NA')),
      Tau = mean(TAU, na.rm = TRUE),
      Tau_sd = sd(TAU, na.rm = TRUE),
      Tau_N = length(which(TAU != 'NA')),
      APOE = length(which(APOE=="1"))
      )

chisq_test(bio.all, SEX ~ LANGUSE)
summary(lm(AGE ~ LANGUSE, data = bio.all))
t.test(AGE ~ LANGUSE, data = bio.all)
summary(lm(EDU ~ LANGUSE, data = bio.all))
t.test(EDU ~ LANGUSE, data = bio.all)
summary(lm(MMSE ~ LANGUSE, data = bio.all))


ddply(bothDTI, .(LANGUSE), summarise, 
      N = length(LANGUSE),
      Fem = length(which(SEX == '2')),
      Cog = mean(CogProfileNum, na.rm = TRUE),
      CN = length(which(DIAGNOSIS == "CN")),
      MCI = length(which(DIAGNOSIS == "MCI")),
      Alz = length(which(DIAGNOSIS == "AD")),
      Age_mean = mean(AGE),
      Age_sd = sd(AGE),
      FA_mean = mean(FA, na.rm = TRUE),
      FA_sd = sd(FA, na.rm = TRUE),
      AD_mean = mean(AD, na.rm = TRUE), 
      AD_sd = sd(AD, na.rm = TRUE),
      RD_mean = mean(RD, na.rm = TRUE),
      RD_sd = sd(RD, na.rm = TRUE),
      Edu_mean = mean(EDU, na.rm = TRUE),
      Edu_sd = sd(EDU, na.rm = TRUE),
      MMSE_mean = mean(MMSE, na.rm = TRUE),
      MMSE_sd = sd(MMSE, na.rm = TRUE),
      AB42 = mean(ABETA42, na.rm = TRUE),
      AB42_sd = sd(ABETA42, na.rm = TRUE),
      AB40 = mean(ABETA40, na.rm = TRUE),
      AB40_sd = sd(ABETA40, na.rm = TRUE),
      ABRatio = mean(ABRATIO, na.rm = TRUE),
      PTau = mean(PTAU, na.rm = TRUE),
      PTau_sd = sd(PTAU, na.rm = TRUE),
      Tau = mean(TAU, na.rm = TRUE),
      Tau_sd = sd(TAU, na.rm = TRUE),
      APOE = length(which(APOE==1)),
      APOEboth = length(which(APOE2==1))
      )





############### DTI PCA for both groups ########################

# PCA(X, scale.unit = TRUE, ncp = 5, graph = TRUE), requires FactoMiner package.
# Do not run PCA with missing values, it will impute values and assign to subjects - has no meaning.

for.pca <- subset(bothDTI, select = c("FA", "AD", "RD"))
res.pca <- PCA(for.pca, graph = TRUE)
scoresPCA <- res.pca$ind$coord[,1]
bothDTI <- cbind(bothDTI, scoresPCA)
biofull <-left_join(bio.all, bothDTI, by = NULL, all.x = TRUE) # all subjects plus PC scores for those with DTI
biofull <- subset(biofull, !is.na(DIAGNOSIS)) #eliminates bilinguals without diagnoses

ddply(biofull, .(LANGUSE, DIAGNOSIS), summarise,
      N = length(which(scoresPCA != "NA")),
      PCA = mean(scoresPCA, na.rm = TRUE),
      PCA_sd = sd(scoresPCA, na.rm = TRUE)
      )



################ Analyses for full dataset #########################

reg.all2 <- lm(scoresPCA ~ LANGUSE + DIAGNOSIS + EDU + AGE + SEX + APOE + APOE2 + ABETA42 + ABRATIO + PTAU + TAU,
               data = biofull)
summary(reg.all2)
EtaSq(reg.all2)

# to compare regression 1 against 2, the number of subjects needs to be the same.
# Run regression 2 first, so that regression 1 uses the right amount of subjects.
reg.all1 <- lm(scoresPCA ~ LANGUSE + DIAGNOSIS + EDU + AGE + SEX, data = na.omit(biofull[ , all.vars(formula(reg.all2))])) 
summary(reg.all1)
EtaSq(reg.all1)

# comparison
anova(reg.all1, reg.all2)



# propensity scores - TWANG PACKAGE REQUIRED
set.seed(1)
biofull$LANGUSE <- as.integer(biofull$LANGUSE) # ps function required integer for the target variable
biofull$LANGUSE <- (biofull$LANGUSE-1)


# using more optimal iterations
ps.all <- ps(LANGUSE ~ SEX  + AGE + EDU + ETHCAT + RACCAT + ABETA42 + ABRATIO + PTAU + TAU +
               APOE,
             data = as.data.frame(biofull), #setting as dataframe seems necessary for ps, but wasn't needed for mnps
             estimand = "ATE",
             verbose = FALSE,
             stop.method = c("es.mean"),
             n.trees = 1000)

# list of PS scores
ps.all$ps
get.weights(ps.all, stop.method = "es.mean", estimand = "ATE") # propensity score weights

# graphical assessments of balance
m3p3 <- love.plot(bal.tab(ps.all), stars = "raw", size = 5, alpha = 0.5)
m3p1 <- plot(ps.all, plots = 1)
m3p2 <- plot(ps.all, plots = 2)
plot(ps.all, plots = 3,  pairwiseMax = FALSE, figureRows = 3)
plot(ps.all, plots = 4)
plot(ps.all, plots = 5)
plot(ps.all, plots = 6)

# tabular assessments of balance
GBM3.bal.tab <- bal.table(ps.all, digits = 3)
write.table(GBM3.bal.tab, file = "GBM3_bal_tab.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

summary(ps.all) #99 iterations

# relative influence of predictors
summary(ps.all$data) # overall info
summary(ps.all$gbm.obj, n.trees=ps.all$desc$es.mean.ATE$n.trees, plot=TRUE)
influplot3 <- summary(ps.all$gbm.obj, n.trees=ps.all$desc$es.mean.ATE$n.trees, plot=TRUE)

# get weights
biofull$w <- get.weights(ps.all, stop.method = "es.mean")
design.ps <- svydesign(ids=~1, weights=~w, data=biofull)

# examining brain PC scores based on diagnosis after re-weighting
glm3.uw <- lm(scoresPCA ~ LANGUSE * DIAGNOSIS, data = biofull)
alias(glm3.uw)
summary(glm3.uw)
glm3 <- svyglm(scoresPCA ~ LANGUSE * DIAGNOSIS, design = design.ps)
alias(glm3)
summary(glm3)
regTermTest(glm3, ~LANGUSE * DIAGNOSIS)

glm3.treat <- svyglm(scoresPCA ~ LANGUSE * DIAGNOSIS, design = design.ps, contrast=list(DIAGNOSIS=contr.sum)) # estimate the causal effect of each treatment relative to the average potential outcome of all the treatments
summary(glm3.treat)
-sum(coef(glm3.treat)[-1]) # MCI estimate, which is a linear combination of other two estimates
sqrt(c(-1,-1) %*% summary(glm4.treat)$cov.scaled[-1,-1] %*% c(-1,-1)) # standard error of this estimate

# design.ps <- update( design.ps , DIAGNOSIS = relevel(DIAGNOSIS,"CN") ) # relevel baseline value to CN if needed

# some relationships of interest
summary(svyglm(MMSE ~ LANGUSE * DIAGNOSIS, design = design.ps))
regTermTest(svyglm(MMSE ~ LANGUSE * DIAGNOSIS, design = design.ps), ~LANGUSE * DIAGNOSIS)

summary(svyglm(ADAS13 ~ LANGUSE * DIAGNOSIS, design = design.ps))
regTermTest(svyglm(ADAS13 ~ LANGUSE*DIAGNOSIS, design = design.ps), ~LANGUSE * DIAGNOSIS)

summary(svyglm(LIMMTOTAL ~ LANGUSE * DIAGNOSIS, design = design.ps))
summary(svyglm(LDELTOTAL ~ LANGUSE * DIAGNOSIS, design = design.ps))
summary(svyglm(PTAU ~ LANGUSE * DIAGNOSIS, design = design.ps))
summary(svyglm(TAU ~ LANGUSE * DIAGNOSIS, design = design.ps))
summary(svyglm(ABETA42 ~ LANGUSE * DIAGNOSIS, design = design.ps))


# glm3 descriptives
ddply(biofull, .(LANGUSE), summarise, 
     CN = length(which(DIAGNOSIS=='CN')),
     MCI = length(which(DIAGNOSIS=='MCI')),
     Alz = length(which(DIAGNOSIS=='AD')),
     PC = mean(scoresPCA, na.rm = TRUE),
     AB42 = mean(ABETA42, na.rm = TRUE),
     Ptau = mean(PTAU, na.rm = TRUE),
     Ttau = mean(TAU, na.rm = TRUE),
     MMSE_m = mean(MMSE, na.rm = TRUE),
     MMSE_sd = sd(MMSE, na.rm = TRUE),
     LM_Total = mean(LIMMTOTAL, na.rm = TRUE),
     LM_Delay = mean(LDELTOTAL, na.rm = TRUE),
     ANART = mean(ANARTERR, na.rm = TRUE)
      )

ddply(biofull, .(DIAGNOSIS), summarise, 
      N = length(DIAGNOSIS),
      PC = mean(scoresPCA, na.rm = TRUE),
      AB42 = mean(ABETA42, na.rm = TRUE),
      Ptau = mean(PTAU, na.rm = TRUE),
      Ttau = mean(TAU, na.rm = TRUE),
      MMSE_m = mean(MMSE, na.rm = TRUE),
      MMSE_sd = sd(MMSE, na.rm = TRUE),
      LM_Total = mean(LIMMTOTAL, na.rm = TRUE),
      LM_Delay = mean(LDELTOTAL, na.rm = TRUE),
      ANART = mean(ANARTERR, na.rm = TRUE),
      ADAS_13 = mean(ADAS13, na.rm = TRUE)
      )

ddply(subset(biofull, LANGUSE=='1'), .(DIAGNOSIS), summarise, 
      N = length(DIAGNOSIS),
      Fem = length(which(SEX == '2')),
      Cog = mean(CogProfileNum, na.rm = TRUE),
      Age_mean = mean(AGE),
      Age_sd = sd(AGE),
      Edu_mean = mean(EDU, na.rm = TRUE),
      Edu_sd = sd(EDU, na.rm = TRUE),
      MMSE_m = mean(MMSE, na.rm = TRUE),
      MMSE_sd = sd(MMSE, na.rm = TRUE),
      PC = mean(scoresPCA, na.rm = TRUE),
      PC_n = length(which(scoresPCA!='NA')),
      FA_mean = mean(FA, na.rm = TRUE),
      FA_sd = sd(FA, na.rm = TRUE),
      AD_mean = mean(AD, na.rm = TRUE), 
      AD_sd = sd(AD, na.rm = TRUE),
      RD_mean = mean(RD, na.rm = TRUE),
      RD_sd = sd(RD, na.rm = TRUE),
      AB42 = mean(ABETA42, na.rm = TRUE),
      AB42_sd = sd(ABETA42, na.rm = TRUE),
      AB42_N = length(which(ABETA42 != 'NA')),
      AB40 = mean(ABETA40, na.rm = TRUE),
      AB40_sd = sd(ABETA40, na.rm = TRUE),
      AB40_N = length(which(ABETA40 != 'NA')),
      ABRatio = mean(ABRATIO, na.rm = TRUE),
      PTau = mean(PTAU, na.rm = TRUE),
      PTau_sd = sd(PTAU, na.rm = TRUE),
      PTau_N = length(which(PTAU != 'NA')),
      Tau = mean(TAU, na.rm = TRUE),
      Tau_sd = sd(TAU, na.rm = TRUE),
      Tau_N = length(which(TAU != 'NA')),
      APOE_N = length(which(APOE==1)),
      APOE_Tot = length(which(APOE != 'NA'))
      )

ddply(subset(biofull, LANGUSE=='0'), .(DIAGNOSIS), summarise, 
      N = length(DIAGNOSIS),
      Fem = length(which(SEX == '2')),
      Cog = mean(CogProfileNum, na.rm = TRUE),
      Age_mean = mean(AGE),
      Age_sd = sd(AGE),
      Edu_mean = mean(EDU, na.rm = TRUE),
      Edu_sd = sd(EDU, na.rm = TRUE),
      MMSE_m = mean(MMSE, na.rm = TRUE),
      MMSE_sd = sd(MMSE, na.rm = TRUE),
      PC = mean(scoresPCA, na.rm = TRUE),
      PC_n = length(which(scoresPCA!='NA')),
      FA_mean = mean(FA, na.rm = TRUE),
      FA_sd = sd(FA, na.rm = TRUE),
      AD_mean = mean(AD, na.rm = TRUE), 
      AD_sd = sd(AD, na.rm = TRUE),
      RD_mean = mean(RD, na.rm = TRUE),
      RD_sd = sd(RD, na.rm = TRUE),
      AB42 = mean(ABETA42, na.rm = TRUE),
      AB42_sd = sd(ABETA42, na.rm = TRUE),
      AB42_N = length(which(ABETA42 != 'NA')),
      AB40 = mean(ABETA40, na.rm = TRUE),
      AB40_sd = sd(ABETA40, na.rm = TRUE),
      AB40_N = length(which(ABETA40 != 'NA')),
      ABRatio = mean(ABRATIO, na.rm = TRUE),
      PTau = mean(PTAU, na.rm = TRUE),
      PTau_sd = sd(PTAU, na.rm = TRUE),
      PTau_N = length(which(PTAU != 'NA')),
      Tau = mean(TAU, na.rm = TRUE),
      Tau_sd = sd(TAU, na.rm = TRUE),
      Tau_N = length(which(TAU != 'NA')),
      APOE_N = length(which(APOE==1)),
      APOE_Tot = length(which(APOE != 'NA'))
      )





########### Plot editing ###################

m1p3 + theme(text = element_text(family = "Segoe UI"),
             legend.title = element_text(size = 12, colour = "#08306B", face = 2),
             legend.text = element_text(size = 12, face = 2),
             axis.title.x = element_text(colour = "#08306B", size = 14, face = 2),
             axis.title.y = element_text(colour = "#08306B", size = 14, face = 2),
             axis.text = element_text(colour="darkslategray", size = 12, face = 2),
             strip.text.x = element_text(size = 12, colour = "darkslateblue", face = "bold"),
             strip.text.y = element_text(size = 12, colour = "darkslateblue", face = "bold", angle = 180)
              )


m2p3 + theme(text = element_text(family = "Segoe UI"),
             legend.title = element_text(size = 12, colour = "#08306B", face = 2),
             legend.text = element_text(size = 12, face = 2),
             axis.title.x = element_text(colour = "#08306B", size = 14, face = 2),
             axis.title.y = element_text(colour = "#08306B", size = 14, face = 2),
             axis.text = element_text(colour="darkslategray", size = 12, face = 2),
             strip.text.x = element_text(size = 12, colour = "darkslateblue", face = "bold"),
             strip.text.y = element_text(size = 12, colour = "darkslateblue", face = "bold", angle = 180)
              )


m3p3 + theme(text = element_text(family = "Segoe UI"),
             legend.title = element_text(size = 12, colour = "#08306B", face = 2),
             legend.text = element_text(size = 12, face = 2),
             axis.title.x = element_text(colour = "#08306B", size = 14, face = 2),
             axis.title.y = element_text(colour = "#08306B", size = 14, face = 2),
             axis.text = element_text(colour="darkslategray", size = 12, face = 2),
             strip.text.x = element_text(size = 12, colour = "darkslateblue", face = "bold"),
             strip.text.y = element_text(size = 12, colour = "darkslateblue", face = "bold", angle = 180)
              )

  
  
  
  
  
  
  
