rm(list=ls())

set.seed(58)

library(dplyr)
library(metafor)
library(mice)
library(knitr)
library(here)
library(ggplot2)
library(raster)
library(RColorBrewer)


test <- readRDS(here("data/effectsheet.rds"))

#Generate a count for which effect is which
names(test)[1] <- "Study"
test <- test %>% group_by(Study) %>% mutate(Effect = row_number(Study))


test$logodds <- NA
test$r <- NA
test$d <- NA
test$b <- NA
test$other <- NA
test$logodds.se <- NA
test$r.se <- NA
test$d.se <- NA
test$b.se <- NA
test$other.se <- NA
test$original.type <- test$Type

#pubyear = year since 1996
test$pubyear <- test$pubyear - 1996


#need to get these effects log-transformed before doing anything.
for (o in 1:nrow(test)) {
  if (test$Type[o] == "OR" | test$Type[o] == "aOR") {
    test$Type[o] <- "logodds"
    test$Value[o] <- log(test$Value[o])
  }
}


#Now reverse everything that needs to be reverse-coded
test$Value <- test$Value * (1 - 2*test$REV1) * (1 - 2*test$REV2)


#Now we will put every effect size and SE in its own field
log.list <- c("logodds")
original.log.list <- c("log-odds ratio (logistic regression coefficient)", "logistic regression B", 
                       "logistic regression coefficient", "logistic regression", "logistic regression coefficient (PH model)")
r.list <- c("r", "latent variable correlations", "partial correlation", "z")
d.list <- c("d", "Cohens d**")
b.list <- c("beta", "b", "standardized regression coefficient", "standardized beta")


for (e in 1:nrow(test)) {
  if ((test$Type[e] %in% c(log.list, original.log.list)) == TRUE) {
    test$logodds[e] <- test$Value[e]
    test$logodds.se[e] <- test$SE[e]
  } else if ((test$Type[e] %in% r.list) == TRUE) {
    test$r[e] <- test$Value[e]
    test$r.se[e] <- test$SE[e]
  } else if ((test$Type[e] %in% d.list) == TRUE) {
    test$d[e] <- test$Value[e]
    test$d.se[e] <- test$SE[e]
  } else if ((test$Type[e] %in% b.list) == TRUE) {
    test$b[e] <- test$Value[e]
    test$b.se[e] <- test$SE[e]
  } else {
    test$other[e] <- test$Value[e]
    test$other.se[e] <- test$SE[e]
  }
}



#imputation function
impute.p <- function(old.p) {
  new.p <- old.p
  if (is.na(old.p) == FALSE) {
    if (old.p == "ns" | old.p == "p >.05" | old.p == "p > .05") {
      new.p <- runif(1, .05, 1)
    }
    if (old.p == "p < .05" | old.p == ".05" | old.p == "0.05") {
      new.p <- runif(1, .01, .05)
    }
    if (old.p == "p < .01" | old.p == ".01" | old.p == "0.01") {
      new.p <- runif(1, .001, .01)
    }
    if (old.p == "p < .001" | old.p == ".001" | old.p == "0.001" | old.p == "0.0" | old.p == "0") {
      new.p <- runif(1, .0001, .001)
    }
  }
  as.numeric(new.p)
}


#Now let's get the SE's for odds ratios which entails either:
#1. transforming the delta-method SE given by the authors, typically calculated as OR*se(logodds)
#(so to get se(logodds) we do se(OR)/OR)
#2. or, more commonly, using the CI to calculate the SE
for (r in 1:nrow(test)) {
  if (((test$Type[r] %in% log.list) == TRUE) & (is.na(test$SE[r]) == FALSE)) {
    test$SE[r]  <- test$SE[r]/test$original.value[r]
  }
  if (((test$Type[r] %in% log.list) == TRUE) & (is.na(test$SE[r]) == TRUE)) {
    test$SE[r]  <- (log(as.numeric(test$CI_upper[r])) - log(as.numeric(test$CI_lower[r])))/ (2*1.96)
  }
}


#Finally, we will get the SE's for the correlations
for (p in 1:nrow(test)) {
  if((test$Type[p] %in% r.list) == TRUE) {
    test$r.se[p] <- sqrt((1-test$r[p]^2)/(test$N[p] - 2))
  }
}

test$r.v <- test$r.se^2


#Now we create all the levels of a variable for the peer connectedness and substance use codes.

test$PCMeasure <- NA
test$PCMeasure[test$SP == 1] <- "Sociometric Popularity"
test$PCMeasure[test$SL == 1] <- "Sociometric Likeability"
test$PCMeasure[test$SM == 1 & test$FR == 1 & test$ISO != 1] <- "aaSociometric Friendship"
test$PCMeasure[test$ISO == 1] <- "Sociometric Isolation"
test$PCMeasure[test$SM == 0 & test$FR == 1 & test$ISO != 1] <- "aSelf-Reported Number of Friends"
test$PCMeasure[test$LON == 1] <- "Loneliness"
test$PCMeasure[test$TIME == 1] <- "Time Spent with Friends"
test$PCMeasure[test$SUP == 1] <- "Support from Friends"
test$PCMeasure[test$GENSUP == 1] <- "Support from All Peers"
test$PCMeasure[test$SRP == 1] <- "Self-Rated Popularity"
test$PCMeasure[test$COM == 1] <- "Social Competency"


test$Substance <- NA
test$SUMeasure <- NA
test$Timeframe <- NA
test$Substance[test$ALC == 1] <- "Alcohol Use"
test$Substance[test$BINGE == 1] <- "Binge Drinking or Drunkenness"
test$Substance[test$MJ == 1] <- "Marijuana Use"
test$Substance[test$TOB == 1] <- "Nicotine Use"
test$Substance[test$OD == 1] <- "Other Drug Use"
test$Substance[test$COMP == 1] <- "Compound Measure of Drug Use"
test$SUMeasure[test$Q == 1] <- "Quantity of Use"
test$SUMeasure[test$F == 1] <- "Frequency of Use"
test$SUMeasure[test$ANY == 1] <- "Any Use"
test$SUMeasure[test$PROB == 1] <- "Problem"
test$SUMeasure[test$Q == 0 & test$F == 0 & test$ANY == 0 & test$PROB == 0] <- "aaNo Measure Specified"
test$Timeframe[test$REC == 1] <- "Recent Use"
test$Timeframe[test$LT == 1] <- "Lifetime Use"
test$Timeframe[test$LT == 0 & test$REC == 0] <- "aNo Frequency Specified"


#We will generate 5 imputations, mostly to impute data for standard errors and ages. 
#The rationale for this is described in the Supplemental Materials. 
#But the first thing we'll do is just generate the 5 datasets.


dat.list <- vector(mode = "list", length = 5)

for (rep1 in 1:5) {

  dat.list[[rep1]] <- test
  
  #Now we will impute SE's on the basis of imputed p. values
#Note that we are doing this for all values, including log-odds as well as r and d
#but we had to log-transform OR's first because the SE's only make sense for the log-odds
dat.list[[rep1]]$sig.impute <- NA
for (s in 1:nrow(dat.list[[rep1]])) {
  if(nchar(dat.list[[rep1]]$Sig.[s]) > 0) {
    dat.list[[rep1]]$sig.impute[s] <- impute.p(dat.list[[rep1]]$Sig.[s])
    dat.list[[rep1]]$SE[s] <- abs(dat.list[[rep1]]$Value[s])/qnorm(1 - dat.list[[rep1]]$sig.impute[s]/2)
  }
}

#Note that we have some SE's that come in as 0
#We will fix them to some arbitrarily small value instead
for (v in 1:nrow(dat.list[[rep1]])) {
  if (is.na(dat.list[[rep1]]$SE[v]) == FALSE) {
    if (dat.list[[rep1]]$SE[v] == 0) {
      dat.list[[rep1]]$SE[v] <- .00001
    }
    if (dat.list[[rep1]]$SE[v] > abs(dat.list[[rep1]]$Value[v]*10)) {
      dat.list[[rep1]]$SE[v] <- abs(dat.list[[rep1]]$Value[v]*10)
    }
  }
}



#Now impute age
for (a in 1:nrow(dat.list[[rep1]])) {
  if (is.na(dat.list[[rep1]]$Age[a]) == TRUE) {
    if(is.na(dat.list[[rep1]]$Age_lower[a]) == FALSE) {
      dat.list[[rep1]]$Age[a] <- runif(1,dat.list[[rep1]]$Age_lower[a],dat.list[[rep1]]$Age_upper[a])
    } else {
      dat.list[[rep1]]$Age[a] <- runif(1, 10, 19)
    }
  }
}


#Now repeat the step that puts every different type of effect in its own column  
for (e in 1:nrow(dat.list[[rep1]])) {
  if ((dat.list[[rep1]]$Type[e] %in% c(log.list, original.log.list)) == TRUE) {
    dat.list[[rep1]]$logodds[e] <- dat.list[[rep1]]$Value[e]
    dat.list[[rep1]]$logodds.se[e] <- dat.list[[rep1]]$SE[e]
  } else if ((dat.list[[rep1]]$Type[e] %in% r.list) == TRUE) {
    dat.list[[rep1]]$r[e] <- dat.list[[rep1]]$Value[e]
    dat.list[[rep1]]$r.se[e] <- dat.list[[rep1]]$SE[e]
  } else if ((dat.list[[rep1]]$Type[e] %in% d.list) == TRUE) {
    dat.list[[rep1]]$d[e] <- dat.list[[rep1]]$Value[e]
    dat.list[[rep1]]$d.se[e] <- dat.list[[rep1]]$SE[e]
  } else if ((dat.list[[rep1]]$Type[e] %in% b.list) == TRUE) {
    dat.list[[rep1]]$b[e] <- dat.list[[rep1]]$Value[e]
    dat.list[[rep1]]$b.se[e] <- dat.list[[rep1]]$SE[e]
  } else {
    dat.list[[rep1]]$other[e] <- dat.list[[rep1]]$Value[e]
    dat.list[[rep1]]$other.se[e] <- dat.list[[rep1]]$SE[e]
  }
}
dat.list[[rep1]]$logodds.v <- dat.list[[rep1]]$logodds.se^2

}


#Now it's time to actually run the models! These are the models for zero-order effects, 
#the ones shown in Table 1 (with all zero-order effects) and Table 2 (with just correlation coefficients). 
#Note that this involves converting among all the different types of zero-order effects. 
#This entails an ugly run of code where we first convert Cohen's $d$ to correlation coefficient $r$. 
#We then convert the log-odds, which we got from contingency tables, to $d$ and then to $r$. 


#Make the objects which will hold all the imputed datasets
zero.order <- vector(mode = "list", length = 5)
just.r <- vector(mode = "list", length = 5)

#Make the objects which will hold all the model results
zo.model <- vector(mode = "list", length = 5)
r.model <- vector(mode = "list", length = 5)

for (rep2 in 1:5) {
  
  zero.order[[rep2]] <- subset(dat.list[[rep2]], dat.list[[rep2]]$Number == 0)
  
  #A bit hackish, but we're going to give converted values of r and d the "converted" label to distinguish them from things that were in that scale originally
  #So we will initialize values of "r.converted" etc. to the original value of r (or r.v, or d, or d.v) if it already existed; this will be blank for cases in which it didn't exist
  
  zero.order[[rep2]]$r.converted <- zero.order[[rep2]]$r
  zero.order[[rep2]]$r.converted.v <- zero.order[[rep2]]$r.v
  zero.order[[rep2]]$d.converted <- zero.order[[rep2]]$d
  zero.order[[rep2]]$d.converted.v <- zero.order[[rep2]]$d.se^2
  
  
  for (z in 1:nrow(zero.order[[rep2]])) {
    #convert d to r
    if (zero.order[[rep2]]$Type[z] == "d" | zero.order[[rep2]]$Type[z] == "Cohens d**") {
      zero.order[[rep2]]$r.converted[z] <- zero.order[[rep2]]$d.converted[z]/sqrt(4 + zero.order[[rep2]]$d.converted[z]^2)
      zero.order[[rep2]]$r.converted.v[z] <- (16*zero.order[[rep2]]$d.converted.v[z])/((4 + zero.order[[rep2]]$d.converted[z]^2)^3)
    }
    #convert logodds to d in cases in which there are no covariates (i.e., zero-order = true)
    #then convert those d's to r's
    #note that the last 2 lines of this for loop replicate what we did above; just want to make sure it gets converted for everyone
    if (zero.order[[rep2]]$Type[z] == "logodds" | zero.order[[rep2]]$Type[z] == "logistic regression coefficient") {
      zero.order[[rep2]]$d.converted[z] <- (sqrt(3)/pi)*zero.order[[rep2]]$Value[z]
      zero.order[[rep2]]$d.converted.v[z] <- (zero.order[[rep2]]$SE[z]^2)*(3/(pi^2))
      zero.order[[rep2]]$r.converted[z] <- zero.order[[rep2]]$d.converted[z]/sqrt(4 + zero.order[[rep2]]$d.converted[z]^2)
      zero.order[[rep2]]$r.converted.v[z] <- (16*zero.order[[rep2]]$d.converted.v[z])/((4 + zero.order[[rep2]]$d.converted[z]^2)^3)
    }
  }
  
  zero.order[[rep2]] <-  escalc(measure = "ZCOR", ri = r.converted, ni = N, data = zero.order[[rep2]])
  
  
  #There weren't enough cases with "Support from All Peers" in any of the zero-order measures to allow the model to be estimated
  #Then, within the ones that were in correlation metric originally, 
  zero.order[[rep2]] <- zero.order[[rep2]][which(zero.order[[rep2]]$PCMeasure != "Support from All Peers"),]
  #and there weren't enough cases with "Sociometric Isolation" when you consider only the effect sizes that were in correlation metric originally
  just.r[[rep2]] <- zero.order[[rep2]][which(zero.order[[rep2]]$Type == "r"),]
  just.r[[rep2]] <- just.r[[rep2]][which(just.r[[rep2]]$PCMeasure != "Sociometric Isolation"),]
  
  
  ## MODEL ESTIMATION SECTION ##
  
  #Model 1: All zero-order effects, even those which were converted, included
  zo.model[[rep2]] <- rma.mv(yi, 
                             vi, 
                             random = ~ 1 | Study/Effect, 
                             data = zero.order[[rep2]],
                             mods = ~ PCMeasure + 
                               Substance + Timeframe + SUMeasure + 
                               Age + Gender + Europe + AAO + AfricaME + LatinAmerica + Longitudinal + PC_outcome + 
                               pubyear + AddHealth)
  
  
  r.model[[rep2]] <- rma.mv(yi, 
                            vi, 
                            random = ~ 1 | Study/Effect, 
                            data = just.r[[rep2]],
                            mods = ~ PCMeasure + 
                              Substance + Timeframe + SUMeasure + 
                              Age + Gender + Europe + AAO + AfricaME + LatinAmerica + Longitudinal + PC_outcome + 
                              pubyear + AddHealth)
  
  
}

#We now run the models for partial effects. 
#This includes the models in Table 3 (for logistic regression) and Table 4 (for linear regression). 
#Two things merit attention here, from the perspective of recreating our decision-making. 
#First, note that for both logistic regression and linear regression, we eliminate certain categories. 
#This is because there were few enough effects that we couldn't actually run the models (e.g., only one or two papers report logistic regression with certain measures of peer connectedness). 
#Second, note that for log-odds we only used effects with variances greater than .0001 or less than 10. 
#Third, as noted in the paper, the Cochran test of moderators in the aggregate (i.e., $Q_M$) is not significant for the model for linear regressions, even though the specific effect of sociometric popularity is. 
#(This will become clear in the next chunk of code.) 
#Thus, we ran a pared-down model which did not include any of the control variables, as none of these were significant. 
#As will be shown in the next chunk, the test for moderators is highly significant when the pared-down model is run. 


#Make the objects that will hold all the imputed datasets
partial <- vector(mode = "list", length = 5)
logreg.partial <- vector(mode = "list", length = 5)
linreg.partial <- vector(mode = "list", length = 5)
linreg.subset <- vector(mode = "list", length = 5)

#Make the objects that will hold all the model results
linreg.model <- vector(mode = "list", length = 5)
logreg.model <- vector(mode = "list", length = 5)
linreg.pared.down <- vector(mode = "list", length = 5)

for (rep3 in 1:5){
  partial[[rep3]] <- subset(dat.list[[rep3]], dat.list[[rep3]]$Number > 0)
  
  
  #Further split the partial dataset into one for logistic regression and one for linear regression
  logreg.partial[[rep3]] <- subset(partial[[rep3]], partial[[rep3]]$Type %in% c(original.log.list, log.list))
  linreg.partial[[rep3]] <- subset(partial[[rep3]], partial[[rep3]]$Type == "b" | partial[[rep3]]$Type == "beta" | partial[[rep3]]$Type == "standardized beta" | partial[[rep3]]$Type == 
                                     "standardized regression coefficient")
  linreg.partial[[rep3]]$t.reg <- linreg.partial[[rep3]]$Value/linreg.partial[[rep3]]$SE
  
  
  linreg.partial[[rep3]] <- escalc(measure = "ZPCOR", ti = t.reg, ni = N, mi = Number, data = linreg.partial[[rep3]])
  linreg.partial[[rep3]] <- linreg.partial[[rep3]][which(abs(linreg.partial[[rep3]]$t.reg) < 10) , ]
  
  
  ##Linear regression model
  which.linreg <- which(linreg.partial[[rep3]]$PCMeasure %in% 
                          c("aaSociometric Friendship", "Sociometric Likeability", 
                            "Sociometric Popularity", "Support from Friends", "Time Spent with Friends", "aSelf-Reported Number of Friends", "Self-Reported Number of Friends"))
  linreg.subset[[rep3]] <- linreg.partial[[rep3]][which.linreg,]
  
  linreg.model[[rep3]] <- rma.mv(yi, 
                                 vi, 
                                 random = ~ 1 | Study/Effect, 
                                 data = linreg.subset[[rep3]],
                                 mods = ~ PCMeasure +
                                   Substance + Timeframe + SUMeasure + 
                                   Age + Gender + Europe + AAO + AfricaME + LatinAmerica + Longitudinal + pubyear + AddHealth +
                                   Number + PC_outcome + Including_peer_use + Including_other_peer_connectedness)
  
  linreg.pared.down[[rep3]] <- rma.mv(yi, 
                                      vi, 
                                      random = ~ 1 | Study/Effect, 
                                      data = linreg.subset[[rep3]],
                                      mods = ~ PCMeasure +
                                        Substance + Timeframe + SUMeasure + Including_peer_use)
  
  
  ##Logistic regression model
  for (q in 1:nrow(logreg.partial[[rep3]])) {
    if (is.na(logreg.partial[[rep3]]$logodds.v[q]) == FALSE) {
      if(logreg.partial[[rep3]]$logodds.v[q] < .0001 | 
         logreg.partial[[rep3]]$logodds.v[q] > 10)
      {logreg.partial[[rep3]]$logodds.v[q] <- NA}
    }
  }
  
  #Logistic regression model
  logreg.model[[rep3]] <- rma.mv(logodds, 
                                 logodds.v, 
                                 random = ~ 1 | Study/Effect, 
                                 data = logreg.partial[[rep3]],
                                 mods = ~ PCMeasure +
                                   Substance + Timeframe + SUMeasure + 
                                   Age + Gender + Europe + AAO + AfricaME + LatinAmerica + Longitudinal + pubyear + AddHealth + 
                                   Number + PC_outcome + Including_peer_use + Including_other_peer_connectedness)
  
  
  
  
}


for (rep4 in 1:5) {
  #### Now write out each imputed dataset
  zero.order.out <- zero.order[[rep4]][,c("Study", "Effect", "r.converted", "r.converted.v", "PCMeasure", "Substance", "Timeframe", "SUMeasure", "Age", "Gender", "Europe", "AAO", "AfricaME", "LatinAmerica", "Longitudinal", "pubyear", "Number", "PC_outcome", "Including_peer_use", "Including_other_peer_connectedness", "AddHealth")]
  names(zero.order.out)[c(3,4)] <- c("yi", "vi")
  write.table(zero.order.out, paste0("data/zeroorder_table_imp",rep4,".csv"), sep = ",", row.names = FALSE, col.names = TRUE)
  
  just.r.out <- just.r[[rep4]][,c("Study", "Effect", "r.converted", "r.converted.v", "PCMeasure", "Substance", "Timeframe", "SUMeasure", "Age", "Gender", "Europe", "AAO", "AfricaME", "LatinAmerica", "Longitudinal", "pubyear", "Number", "PC_outcome", "Including_peer_use", "Including_other_peer_connectedness", "AddHealth")]
  names(just.r.out)[c(3,4)] <- c("yi", "vi")
  write.table(just.r.out, paste0("data/justr_table_imp",rep4,".csv"), sep = ",", row.names = FALSE, col.names = TRUE)
  
  linreg.subset.out <- linreg.subset[[rep4]][,c("Study", "Effect", "yi", "vi", "PCMeasure", "Substance", "Timeframe", "SUMeasure", "Age", "Gender", "Europe", "AAO", "AfricaME", "LatinAmerica", "Longitudinal", "pubyear", "Number", "PC_outcome", "Including_peer_use", "Including_other_peer_connectedness", "AddHealth")]
  write.table(linreg.subset.out, paste0("data/linreg_table_imp",rep4,".csv"), sep = ",", row.names = FALSE, col.names = TRUE)
  
  
  logreg.subset.out <- logreg.partial[[rep4]][,c("Study", "Effect", "logodds", "logodds.v", "PCMeasure", "Substance", "Timeframe", "SUMeasure", "Age", "Gender", "Europe", "AAO", "AfricaME", "LatinAmerica", "Longitudinal", "pubyear", "Number", "PC_outcome", "Including_peer_use", "Including_other_peer_connectedness", "AddHealth")]
  names(logreg.subset.out)[c(3,4)] <- c("yi", "vi")
  write.table(logreg.subset.out, paste0("data/logreg_table_imp",rep4,".csv"), sep = ",", row.names = FALSE, col.names = TRUE)
}


#Get glanced information from the pool function
#Get the number of studies and effects in each model
zo.model[[1]]
linreg.model[[1]]
logreg.model[[1]]


#Get the number of studies reporting each of the following types of effects (i.e., k in Tables 1-4)

table(zero.order[[1]]$PCMeasure)
table(zero.order[[1]]$Substance)
table(zero.order[[1]]$SUMeasure)
table(zero.order[[1]]$Timeframe)

table(just.r[[1]]$PCMeasure)
table(just.r[[1]]$Substance)
table(just.r[[1]]$SUMeasure)
table(just.r[[1]]$Timeframe)

table(linreg.subset[[1]]$PCMeasure)
table(linreg.subset[[1]]$Substance)
table(linreg.subset[[1]]$SUMeasure)
table(linreg.subset[[1]]$Timeframe)

table(logreg.partial[[1]]$PCMeasure)
table(logreg.partial[[1]]$Substance)
table(logreg.partial[[1]]$SUMeasure)
table(logreg.partial[[1]]$Timeframe)

#Now we pool the results
pooled.zo <- pool(zo.model)
pooled.zo$glanced
round(colMeans(pooled.zo$glanced),3)
kable(summary(pooled.zo))
write.table(summary(pooled.zo), file = "data/pooled_zo_1007.csv", sep = ",", row.names = FALSE, col.names = TRUE)

pooled.r <- pool(r.model)
pooled.r$glanced
round(colMeans(pooled.r$glanced),3)
kable(summary(pooled.r))
write.table(summary(pooled.r), file = "data/pooled_r_1007.csv", sep = ",", row.names = FALSE, col.names = TRUE)

pooled.linreg <- pool(linreg.model)
pooled.linreg$glanced
round(colMeans(pooled.linreg$glanced),3)
kable(summary(pooled.linreg))
write.table(summary(pooled.linreg), file = "data/pooled_linreg_1007.csv", sep = ",", row.names = FALSE, col.names = TRUE)


pooled.linreg.pared.down <- pool(linreg.pared.down)
pooled.linreg.pared.down$glanced
round(colMeans(pooled.linreg.pared.down$glanced),3)
kable(summary(pooled.linreg.pared.down))
write.table(summary(pooled.linreg.pared.down), file = "data/pooled_linreg_pareddown_1007.csv", sep = ",", row.names = FALSE, col.names = TRUE)


pooled.logreg <- pool(logreg.model)
pooled.logreg$glanced
round(colMeans(pooled.logreg$glanced),3)
kable(summary(pooled.logreg))
write.table(summary(pooled.logreg), file = "data/pooled_logreg_1007.csv", sep = ",", row.names = FALSE, col.names = TRUE)
