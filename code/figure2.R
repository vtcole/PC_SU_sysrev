rm(list=ls())

library(dplyr)
library(metafor)
library(here)
library(ggplot2)
library(raster)
library(RColorBrewer)

##################

## Data management -- same steps as the regression analyses

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

#################3

test.justzo <- subset(test, test$Number == 0)

test.justzo$r.converted <- test.justzo$r
test.justzo$r.converted.v <- test.justzo$r.v
test.justzo$d.converted <- test.justzo$d
test.justzo$d.converted.v <- test.justzo$d.se^2

for (z in 1:nrow(test.justzo)) {
  if (test.justzo$Type[z] == "d" | test.justzo$Type[z] == "Cohens d**") {
    test.justzo$r.converted[z] <- test.justzo$d.converted[z]/sqrt(4 + test.justzo$d.converted[z]^2)
    test.justzo$r.converted.v[z] <- (16*test.justzo$d.converted.v[z])/((4 + test.justzo$d.converted[z]^2)^3)
  }
  if (test.justzo$Type[z] == "logodds" | test.justzo$Type[z] == "logistic regression coefficient") {
    test.justzo$d.converted[z] <- (sqrt(3)/pi)*test.justzo$Value[z]
    test.justzo$d.converted.v[z] <- (test.justzo$SE[z]^2)*(3/(pi^2))
    test.justzo$r.converted[z] <- test.justzo$d.converted[z]/sqrt(4 + test.justzo$d.converted[z]^2)
    test.justzo$r.converted.v[z] <- (16*test.justzo$d.converted.v[z])/((4 + test.justzo$d.converted[z]^2)^3)
  }
}



just.effects <- test.justzo[,c("Study", "PCMeasure", "SUMeasure", "Substance", "Timeframe", "r.converted", "r.converted.v", "N")]


#These are the functions for calculating the mean effects for each study's combination of substance use and peer connectedness 
#(i.e., each circle within a given subplot), and assigning them positions to graph them.

#This is just a helper function to get the vertical and horizontal positions for the plot. This basically positions each of the means we'll be calculating into a 3x4 matrix for each smaller subplot.

make.12 <- function(x) {
  x$horiz <- NA
  x$vert <- NA
  len <- nrow(x)
  for (row in 1:len) {
    x$horiz[row] <- 3-((12-row)%%4)
    x$vert[row] <- (row-x$horiz[row] - 1)/4
  }
  x
}


#Somewhat hackish function for laying out the different mean effects within each study. This calculates the mean value of all of a study's effects corresponding to a given combination of substance use and peer connection. It also assigns the vertical and horizontal position to each of these means, for plotting purposes.

mat.fn <- function(PCvar, which.su, SUvar, PC.title.var, SU.title.var) {
  the.mat <- subset(just.effects, just.effects$PCMeasure == PCvar & just.effects[,c(which.su)] == SUvar)
  the.effects <- aggregate(the.mat[,c("r.converted", "N")], by = list(the.mat$Study), mean)
  the.effects <- make.12(the.effects)
  the.effects$horiz <- the.effects$horiz/3 + .1
  the.effects$vert <- the.effects$vert/3 + .1
  the.effects$N[the.effects$N == max(the.effects$N)]
  the.effects$PC <- PCvar
  the.effects$SU <- SUvar
  the.effects$PC.title <- PC.title.var
  the.effects$SU.title <- SU.title.var
  the.effects
}


#Need an empty plot for any combination of substance use and peer connectedness taht doesn't have any effects.
empty <- data.frame()
empty.plot <- ggplot(empty) + geom_point() + xlim(0, 10) + ylim(0, 100) + 
  theme_void() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))


#Now we will use metafor to find the aggregate effects for every different level of the peer connectedness and substance use variables. Note that we use the original correlation metric here -- i.e., unlike in the regression models we don't convert them to z-scores. 
#That's so that we can get the actual number in correlation metric.


#This is how we find the correlation for each value of peer connectedness- or substance use- related variables in the aggregate -- i.e., the correlations that are listed in the margins.

#Sociometric friendship
sf.model <- rma.mv(r.converted, 
                   r.converted.v, 
                   random = ~ 1 | Study/Effect, 
                   data = (test.justzo[which(test.justzo$PCMeasure == "aaSociometric Friendship"),]))
sf.model.results <- paste("Sociometric Friendship", "r = .0229*", "14 studies", "86 effects", sep = "\n")

#Sociometric isolation
iso.model <- rma.mv(r.converted, 
                    r.converted.v, 
                    random = ~ 1 | Study/Effect, 
                    data = (test.justzo[which(test.justzo$PCMeasure == "Sociometric Isolation"),]))
summary(iso.model)
iso.model.results <- paste("Sociometric Isolation", "r = -.0690", "6 studies", "33 effects", sep = "\n")

#Sociometric likeability
sl.model <- rma.mv(r.converted, 
                   r.converted.v, 
                   random = ~ 1 | Study/Effect, 
                   data = (test.justzo[which(test.justzo$PCMeasure == "Sociometric Likeability"),]))
summary(sl.model)
sl.model.results <- paste("Sociometric Likeability", "r = .0403", "20 studies", "121 effects", sep = "\n")

#Sociometric popularity
sp.model <- rma.mv(r.converted, 
                   r.converted.v, 
                   random = ~ 1 | Study/Effect, 
                   data = (test.justzo[which(test.justzo$PCMeasure == "Sociometric Popularity"),]))
summary(sp.model)
sp.model.results <- paste("Sociometric Popularity", "r = .2022***", "16 studies", "86 effects", sep = "\n")

#Support from friends
sup.model <- rma.mv(r.converted, 
                    r.converted.v, 
                    random = ~ 1 | Study/Effect, 
                    data = (test.justzo[which(test.justzo$PCMeasure == "Support from Friends"),]))
summary(sup.model)
sup.model.results <- paste("Support from Friends", "r = .0209", "12 studies", "90 effects", sep = "\n")

#Time spent with friends
time.model <- rma.mv(r.converted, 
                     r.converted.v, 
                     random = ~ 1 | Study/Effect, 
                     data = (test.justzo[which(test.justzo$PCMeasure == "Time Spent with Friends"),]))
summary(time.model)
time.model.results <- paste("Time with Friends", "r = .1063", "10 studies", "89 effects", sep = "\n")

#Alcohol use
alc.model <- rma.mv(r.converted, 
                    r.converted.v,
                    random = ~ 1 | Study/Effect,
                    data = (test.justzo[which(test.justzo$Substance == "Alcohol Use"),]))
alc.model.results <- paste("Alcohol Use", "r = .0815***", "30 studies", "199 effects", sep = "\n")

#Binge drinking
binge.model <- rma.mv(r.converted, 
                      r.converted.v,
                      random = ~ 1 | Study/Effect,
                      data = (test.justzo[which(test.justzo$Substance == "Binge Drinking or Drunkenness"),]))
binge.model.results <- paste("Binge / Drunkenness", "r = .1045***", "11 studies", "49 effects", sep = "\n")

#Marijuana use
mj.model <- rma.mv(r.converted, 
                   r.converted.v,
                   random = ~ 1 | Study/Effect,
                   data = (test.justzo[which(test.justzo$Substance == "Marijuana Use"),]))
mj.model.results <- paste("Marijuana Use", "r = .0723*", "13 studies", "87 effects", sep = "\n")

#Nicotine use
nic.model <- rma.mv(r.converted, 
                    r.converted.v,
                    random = ~ 1 | Study/Effect,
                    data = (test.justzo[which(test.justzo$Substance == "Nicotine Use"),]))
nic.model.results <- paste("Nicotine Use", "r = .0606*", "27 studies", "158 effects", sep = "\n")

#Compound measure of substance use
com.model <- rma.mv(r.converted, 
                    r.converted.v,
                    random = ~ 1 | Study/Effect,
                    data = (test.justzo[which(test.justzo$Substance == "Compound Measure of Drug Use"),]))
com.model.results <- paste("Composite Measures", "r = .0398", "17 studies", "83 effects", sep = "\n")

#Measure = any substance use
any.model <- rma.mv(r.converted, 
                    r.converted.v,
                    random = ~ 1 | Study/Effect,
                    data = (test.justzo[which(test.justzo$SUMeasure == "Any Use"),]))
any.model.results <- paste("Any Use", "r = .0589", "16 studies", "191 effects", sep = "\n")

#Measure = frequency of substance use
freq.model <- rma.mv(r.converted, 
                     r.converted.v,
                     random = ~ 1 | Study/Effect,
                     data = (test.justzo[which(test.justzo$SUMeasure == "Frequency of Use"),]))
freq.model.results <- paste("Frequency of Use", "r = .0911***", "39 studies", "304 effects", sep = "\n")

#Measure = problems with substance use
prob.model <- rma.mv(r.converted, 
                     r.converted.v,
                     random = ~ 1 | Study/Effect,
                     data = (test.justzo[which(test.justzo$SUMeasure == "Problem"),]))
prob.model.results <- paste("Problem Use", "r = .0478", "5 studies", "48 effects", sep = "\n")

#Measure = quantity of substance use
q.model <- rma.mv(r.converted, 
                  r.converted.v,
                  random = ~ 1 | Study/Effect,
                  data = (test.justzo[which(test.justzo$SUMeasure == "Quantity of Use"),]))
q.model.results <- paste("Quantity Use", "r = .0419", "6 studies", "14 effects", sep = "\n")

#Timeframe = lifetime substance use
lt.model <- rma.mv(r.converted, 
                   r.converted.v,
                   random = ~ 1 | Study/Effect,
                   data = (test.justzo[which(test.justzo$Timeframe == "Lifetime Use"),]))
lt.model.results <- paste("Lifetime Use", "r = .0813**", "10 studies", "78 effects", sep = "\n")


#Timeframe = recent substance use
rec.model <- rma.mv(r.converted, 
                    r.converted.v,
                    random = ~ 1 | Study/Effect,
                    data = (test.justzo[which(test.justzo$Timeframe == "Recent Use"),]))
rec.model.results <- paste("Recent Use", "r = .0778***", "46 studies", "434 effects", sep = "\n")


#Apply the matrix function we used above to get the means and positions for the plots.

SP.ALC <- mat.fn("Sociometric Popularity", "Substance", "Alcohol Use", sp.model.results, alc.model.results)
SP.BINGE <- mat.fn("Sociometric Popularity", "Substance", "Binge Drinking or Drunkenness", sp.model.results, binge.model.results)
SP.COM <- mat.fn("Sociometric Popularity", "Substance", "Compound Measure of Drug Use", sp.model.results, com.model.results)
SP.MJ <- mat.fn("Sociometric Popularity", "Substance", "Marijuana Use", sp.model.results, mj.model.results)
SP.NIC <- mat.fn("Sociometric Popularity", "Substance", "Nicotine Use", sp.model.results, nic.model.results)
SP.ANY <- mat.fn("Sociometric Popularity", "SUMeasure", "Any Use", sp.model.results, any.model.results)
SP.FREQ <- mat.fn("Sociometric Popularity", "SUMeasure", "Frequency of Use", sp.model.results, freq.model.results)
SP.PROB <- mat.fn("Sociometric Popularity", "SUMeasure", "Problem", sp.model.results, prob.model.results)
SP.Q <- mat.fn("Sociometric Popularity", "SUMeasure", "Quantity of Use", sp.model.results, q.model.results)
SP.LT <- mat.fn("Sociometric Popularity", "Timeframe", "Lifetime Use", sp.model.results, lt.model.results)
SP.REC <- mat.fn("Sociometric Popularity", "Timeframe", "Recent Use", sp.model.results, rec.model.results)


SL.ALC <- mat.fn("Sociometric Likeability", "Substance", "Alcohol Use", sl.model.results, alc.model.results)
SL.BINGE <- mat.fn("Sociometric Likeability", "Substance", "Binge Drinking or Drunkenness", sl.model.results, binge.model.results)
SL.COM <- mat.fn("Sociometric Likeability", "Substance", "Compound Measure of Drug Use", sl.model.results, com.model.results)
SL.MJ <- mat.fn("Sociometric Likeability", "Substance", "Marijuana Use", sl.model.results, mj.model.results)
SL.NIC <- mat.fn("Sociometric Likeability", "Substance", "Nicotine Use", sl.model.results, nic.model.results)
SL.ANY <- mat.fn("Sociometric Likeability", "SUMeasure", "Any Use", sl.model.results, any.model.results)
SL.FREQ <- mat.fn("Sociometric Likeability", "SUMeasure", "Frequency of Use", sl.model.results, freq.model.results)
SL.PROB <- mat.fn("Sociometric Likeability", "SUMeasure", "Problem", sl.model.results, prob.model.results)
SL.Q <- mat.fn("Sociometric Likeability", "SUMeasure", "Quantity of Use", sl.model.results, q.model.results)
SL.LT <- mat.fn("Sociometric Likeability", "Timeframe", "Lifetime Use", sl.model.results, lt.model.results)
SL.REC <- mat.fn("Sociometric Likeability", "Timeframe", "Recent Use", sl.model.results, rec.model.results)


SF.ALC <- mat.fn("aaSociometric Friendship", "Substance", "Alcohol Use", sf.model.results, alc.model.results)
SF.BINGE <- mat.fn("aaSociometric Friendship", "Substance", "Binge Drinking or Drunkenness", sf.model.results, binge.model.results)
SF.COM <- mat.fn("aaSociometric Friendship", "Substance", "Compound Measure of Drug Use", sf.model.results, com.model.results)
SF.MJ <- mat.fn("aaSociometric Friendship", "Substance", "Marijuana Use", sf.model.results, mj.model.results)
SF.NIC <- mat.fn("aaSociometric Friendship", "Substance", "Nicotine Use", sf.model.results, nic.model.results)
SF.ANY <- mat.fn("aaSociometric Friendship", "SUMeasure", "Any Use", sf.model.results, any.model.results)
SF.FREQ <- mat.fn("aaSociometric Friendship", "SUMeasure", "Frequency of Use", sf.model.results, freq.model.results)
SF.PROB <- mat.fn("aaSociometric Friendship", "SUMeasure", "Problem", sf.model.results, prob.model.results)
SF.Q <- mat.fn("aaSociometric Friendship", "SUMeasure", "Quantity of Use", sf.model.results, q.model.results)
SF.LT <- mat.fn("aaSociometric Friendship", "Timeframe", "Lifetime Use", sf.model.results, lt.model.results)
SF.REC <- mat.fn("aaSociometric Friendship", "Timeframe", "Recent Use", sf.model.results, rec.model.results)

#Important note! Because of the data management we did, we need to re-reverse isolation. 
#That is, for the models, we reverse scored it so we're re-reversing it here. Why did we do this?
#Mostly just to maintain consistency across different parts of the code. 
#We didn't want to alter anything in the data management section of one and not another. 
ISO.ALC <- mat.fn("Sociometric Isolation", "Substance", "Alcohol Use", iso.model.results, alc.model.results)
ISO.ALC$r.converted <- ISO.ALC$r.converted*-1
ISO.BINGE <- empty.plot
ISO.COM <- mat.fn("Sociometric Isolation", "Substance", "Compound Measure of Drug Use", iso.model.results, com.model.results)
ISO.COM$r.converted <- ISO.COM$r.converted*-1
ISO.MJ <- mat.fn("Sociometric Isolation", "Substance", "Marijuana Use", iso.model.results, mj.model.results)
ISO.MJ$r.converted <- ISO.MJ$r.converted*-1
ISO.NIC <- mat.fn("Sociometric Isolation", "Substance", "Nicotine Use", iso.model.results, nic.model.results)
ISO.NIC$r.converted <- ISO.NIC$r.converted*-1
ISO.ANY <- mat.fn("Sociometric Isolation", "SUMeasure", "Any Use", iso.model.results, any.model.results)
ISO.ANY$r.converted <- ISO.ANY$r.converted*-1
ISO.FREQ <- mat.fn("Sociometric Isolation", "SUMeasure", "Frequency of Use", iso.model.results, freq.model.results)
ISO.FREQ$r.converted <- ISO.FREQ$r.converted*-1
ISO.PROB <- mat.fn("Sociometric Isolation", "SUMeasure", "Problem", iso.model.results, prob.model.results)
ISO.PROB$r.converted <- ISO.PROB$r.converted*-1
ISO.Q <- mat.fn("Sociometric Isolation", "SUMeasure", "Quantity of Use", iso.model.results, q.model.results)
ISO.Q$r.converted <- ISO.Q$r.converted*-1
ISO.LT <- mat.fn("Sociometric Isolation", "Timeframe", "Lifetime Use", iso.model.results, lt.model.results)
ISO.LT$r.converted <- ISO.LT$r.converted*-1
ISO.REC <- mat.fn("Sociometric Isolation", "Timeframe", "Recent Use", iso.model.results, rec.model.results)
ISO.REC$r.converted <- ISO.REC$r.converted*-1

SUP.ALC <- mat.fn("Support from Friends", "Substance", "Alcohol Use", sup.model.results, alc.model.results)
SUP.BINGE <- mat.fn("Support from Friends", "Substance", "Binge Drinking or Drunkenness", sup.model.results, binge.model.results)
SUP.COM <- mat.fn("Support from Friends", "Substance", "Compound Measure of Drug Use", sup.model.results, com.model.results)
SUP.MJ <- mat.fn("Support from Friends", "Substance", "Marijuana Use", sup.model.results, mj.model.results)
SUP.NIC <- mat.fn("Support from Friends", "Substance", "Nicotine Use", sup.model.results, nic.model.results)
SUP.ANY <- mat.fn("Support from Friends", "SUMeasure", "Any Use", sup.model.results, any.model.results)
SUP.FREQ <- mat.fn("Support from Friends", "SUMeasure", "Frequency of Use", sup.model.results, freq.model.results)
SUP.PROB <- empty.plot
SUP.Q <- empty.plot
SUP.LT <- mat.fn("Support from Friends", "Timeframe", "Lifetime Use", sup.model.results, lt.model.results)
SUP.REC <- mat.fn("Support from Friends", "Timeframe", "Recent Use", sup.model.results, rec.model.results)

TIME.ALC <- mat.fn("Time Spent with Friends", "Substance", "Alcohol Use", time.model.results, alc.model.results)
TIME.BINGE <- mat.fn("Time Spent with Friends", "Substance", "Binge Drinking or Drunkenness", time.model.results, binge.model.results)
TIME.COM <- mat.fn("Time Spent with Friends", "Substance", "Compound Measure of Drug Use", time.model.results, com.model.results)
TIME.MJ <- mat.fn("Time Spent with Friends", "Substance", "Marijuana Use", time.model.results, mj.model.results)
TIME.NIC <- mat.fn("Time Spent with Friends", "Substance", "Nicotine Use", time.model.results, nic.model.results)
TIME.ANY <- mat.fn("Time Spent with Friends", "SUMeasure", "Any Use", time.model.results, any.model.results)
TIME.FREQ <- mat.fn("Time Spent with Friends", "SUMeasure", "Frequency of Use", time.model.results, freq.model.results)
TIME.PROB <- empty.plot
TIME.Q <- empty.plot
TIME.LT <- mat.fn("Time Spent with Friends", "Timeframe", "Lifetime Use", time.model.results, lt.model.results)
TIME.REC <- mat.fn("Time Spent with Friends", "Timeframe", "Recent Use", time.model.results, rec.model.results)

mat.data <- rbind(SF.ALC, SF.BINGE, SF.COM, SF.MJ, SF.NIC,
                  SP.ALC, SP.BINGE, SP.COM, SP.MJ, SP.NIC,
                  SL.ALC, SL.BINGE, SL.COM, SL.MJ, SL.NIC,
                  ISO.ALC, ISO.BINGE, ISO.COM, ISO.MJ, ISO.NIC,
                  SUP.ALC, SUP.BINGE, SUP.COM, SUP.MJ, SUP.NIC,
                  TIME.ALC, TIME.BINGE, TIME.COM, TIME.MJ, TIME.NIC)

mat.data.2 <- rbind(SF.ANY, SF.FREQ, SF.PROB, SF.Q, SF.LT, SF.REC,
                    SP.ANY, SP.FREQ, SP.PROB, SP.Q, SP.LT, SP.REC,
                    SL.ANY, SL.FREQ, SL.PROB, SL.Q, SL.LT, SL.REC,
                    ISO.ANY, ISO.FREQ, ISO.PROB, ISO.Q, ISO.LT, ISO.REC,
                    SUP.ANY, SUP.FREQ, SUP.PROB, SUP.Q, SUP.LT, SUP.REC,
                    TIME.ANY, TIME.FREQ, TIME.PROB, TIME.Q, TIME.LT, TIME.REC)


names(mat.data)[2] <- "Correlation"
names(mat.data.2)[2] <- "Correlation"


#Now we use ggplot2 to make our plots. This includes Figure 2 and the plot on the Supplemental Material. 


#Figure 2
ggplot(mat.data, aes(x = horiz, y = vert, fill = Correlation)) + geom_point(aes(size = N), pch = 21)+ 
  scale_size(trans = "log10") + xlim(c(-.01, 1.25)) + ylim(c(-.01, .85)) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-.5, .5)) +
  theme_void() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        legend.position = "left") +
  facet_grid(PC.title ~ SU.title)


#Supplemental plot
ggplot(mat.data.2, aes(x = horiz, y = vert, fill = Correlation)) + geom_point(aes(size = N), pch = 21)+ 
  scale_size(trans = "log10") + xlim(c(-.01, 1.25)) + ylim(c(-.01, .85)) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-.5, .5)) +
  theme_void() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        legend.position = "left") +
  facet_grid(PC.title ~ SU.title)

