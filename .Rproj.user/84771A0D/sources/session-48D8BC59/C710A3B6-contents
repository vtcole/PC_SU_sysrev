rm(list=ls())

library(ggplot2)

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

#Get study-level dataset so that you can get the publication year of each study
studylevel <- test[,c("Study", "pubyear", "Longitudinal", "Europe", "LatinAmerica", "AAO", "AfricaME")]
studylevel$pubyear <- studylevel$pubyear + 1996

#Get just the number of studies per year
yeardata <- unique(studylevel[,c("Study", "pubyear")])
yeardata <- as.data.frame(table(yeardata$pubyear))
names(yeardata) <- c("Year", "Frequency")

#Plot the number of studies per year
ggplot(data=yeardata, aes(x = Year, y = Frequency, group = 1)) + 
  geom_line() + 
  geom_point() + 
  theme_minimal()

#Tabulate study design features
table(test$Longitudinal)
table(test$PC_outcome)

#Tabulate age and gender
mean(test$Age, na.rm = TRUE)
table(is.na(test$Age))

table(test$Gender) #87 effects were on male participants exclusively; 90 were on female participants exclusively
genderdata <- subset(test, test$Gender != 0 & test$Gender != 100)
hist(genderdata$Gender, main = "Frequency of gender (coded as percent male) for each effect", xlab = "Percent male")

hist(test$Age, xlab = "Age in Years", main = "Frequency of average age for each effect")
