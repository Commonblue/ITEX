library(tidyverse);library(metafor)

#prepare sept
##############################
plots<-read.table("Metadata_general.csv", header=T, sep=";", dec='.')
hyf<-read.table("Hyphal_length_analyses_ITEX.csv", header=T, sep=";", dec=',')
names(hyf)[names(hyf) == "Hyphal_length_uns"] <- "Unsep"
names(hyf)[names(hyf) == "Hyphal_length_sept"] <- "Sep"
names(hyf)[names(hyf) == "Total_hyphal_length"] <- "Total"

sephyf <- hyf[, c("Label", "Sep", "Treatment")]
plots <- plots[, c("Label", "Site", "Vegetationtype")]


df<-merge(plots,sephyf, by="Label")
df <- df[, !(names(df) == "Label")]
names(df)[names(df) == "Treatment"] <- "tt"



meanC <- aggregate(Sep ~ Site, data = df[df$tt == "Control", ], FUN = mean , na.rm = TRUE)
names(meanC)[2] <- "meanC"
stdevC <- aggregate(Sep ~ Site, data = df[df$tt == "Control", ], FUN= sd , na.rm = TRUE)
names(stdevC)[2] <- "stdevC"
nC <- aggregate(Sep ~ Site, data = df[df$tt == "Control", ], FUN=length)
names(nC)[2] <- "nC"

meanOTC <- aggregate(Sep ~ Site, data = df[df$tt == "OTC", ], FUN = mean , na.rm = TRUE)
names(meanOTC)[2] <- "meanOTC"
stdevOTC <- aggregate(Sep ~ Site, data = df[df$tt == "OTC", ], FUN= sd , na.rm = TRUE)
names(stdevOTC)[2] <- "stdevOTC"
nOTC <- aggregate(Sep ~ Site, data = df[df$tt == "OTC", ], FUN=length)
names(nOTC)[2] <- "nOTC"

dfC <- merge(meanC, stdevC, by = "Site")
dfC <- merge(dfC, nC, by = "Site")
dfC <- merge(df,dfC , by = "Site")
dfC <- dfC[, c("Site", "meanC", "stdevC", "nC")]
names(dfC)

dfOTC <- merge(meanOTC, stdevOTC, by = "Site")
dfOTC <- merge(dfOTC, nOTC, by = "Site")
dfOTC <- merge(df,dfOTC , by = "Site")
dfOTC <- dfOTC[, !(names(dfOTC) %in% c("Sep", "tt"))]

dat1 <- merge(dfOTC, dfC, by = "Site")
dat1<-unique(dat1)

OK<-escalc(measure= "SMD", m2i=meanC, sd2i=stdevC, n2i=nC, 
           m1i=meanOTC, sd1i=stdevOTC, n1i=nOTC, data=dat1)

overallresult <- rma(yi, vi, data=OK)
summary(overallresult)
forest.rma(overallresult, slab = OK$Site, addpred=T)
forest.rma(overallresult, slab = OK$Site, addpred=T, order=overallresult$yi, annotate=F)

mod.gradeq <-rma(yi,vi,mods = ~ factor(Vegetationtype), data=OK)
mod.gradeq
               
##############################
#unsept prepare
##############################
plots<-read.table("Metadata_general.csv", header=T, sep=";", dec='.')
hyf<-read.table("Hyphal_length_analyses_ITEX.csv", header=T, sep=";", dec=',')
names(hyf)[names(hyf) == "Hyphal_length_uns"] <- "Unsep"
names(hyf)[names(hyf) == "Hyphal_length_sept"] <- "Sep"
names(hyf)[names(hyf) == "Total_hyphal_length"] <- "Total"

unsephyf <- hyf[, c("Label", "Unsep", "Treatment")]
plots <- plots[, c("Label", "Site", "Vegetationtype")]


df<-merge(plots,unsephyf, by="Label", drop=F)
df <- df[, !(names(df) == "Label")]
names(df)[names(df) == "Treatment"] <- "tt"


df[] <- lapply(df, function(x) if (is.character(x) | is.integer(x)) as.factor(x) else x)


meanC <- aggregate(Unsep ~ Site, data = df[df$tt == "Control", ], FUN = mean , na.rm = TRUE)
names(meanC)[2] <- "meanC"
stdevC <- aggregate(Unsep ~ Site, data = df[df$tt == "Control", ], FUN= sd , na.rm = TRUE)
names(stdevC)[2] <- "stdevC"
nC <- aggregate(Unsep ~ Site, data = df[df$tt == "Control", ], FUN=length)
names(nC)[2] <- "nC"

meanOTC <- aggregate(Unsep ~ Site, data = df[df$tt == "OTC", ], FUN = mean , na.rm = TRUE)
names(meanOTC)[2] <- "meanOTC"
stdevOTC <- aggregate(Unsep ~ Site, data = df[df$tt == "OTC", ], FUN= sd , na.rm = TRUE)
names(stdevOTC)[2] <- "stdevOTC"
nOTC <- aggregate(Unsep ~ Site, data = df[df$tt == "OTC", ], FUN=length)
names(nOTC)[2] <- "nOTC"

dfC <- merge(meanC, stdevC, by = "Site")
dfC <- merge(dfC, nC, by = "Site")
dfC <- merge(df,dfC , by = "Site")
dfC <- dfC[, c("Site", "meanC", "stdevC", "nC")]
names(dfC)

dfOTC <- merge(meanOTC, stdevOTC, by = "Site")
dfOTC <- merge(dfOTC, nOTC, by = "Site")
dfOTC <- merge(df,dfOTC , by = "Site")
dfOTC <- dfOTC[, !(names(dfOTC) %in% c("Unsep", "tt"))]

dat1 <- merge(dfOTC, dfC, by = "Site")
dat1<-unique(dat1)

OK<-escalc(measure= "SMD", m2i=meanC, sd2i=stdevC, n2i=nC, 
           m1i=meanOTC, sd1i=stdevOTC, n1i=nOTC, data=dat1)

overallresult <- rma(yi, vi, data=OK)
summary(overallresult)
forest.rma(overallresult, slab = OK$Site, addpred=T)

mod.gradeq <-rma(yi,vi,mods = ~ factor(Vegetationtype), data=OK)
fitted.rma(mod.gradeq)
predict.rma(mod.gradeq)

mod.gradeq
influence(mod.gradeq)

################################
#prepare ratio
###############################

plots<-read.table("Metadata_general.csv", header=T, sep=";", dec='.')
hyf<-read.table("Hyphal_length_analyses_ITEX.csv", header=T, sep=";", dec=',')
names(hyf)[names(hyf) == "Hyphal_length_uns"] <- "Unsep"
names(hyf)[names(hyf) == "Hyphal_length_sept"] <- "Sep"
names(hyf)[names(hyf) == "Total_hyphal_length"] <- "Total"

hyf$Ratio<-(hyf$Unsep/(hyf$Sep+hyf$Unsep))*100
hyf <- na.omit(hyf)

Ratiohyf <- hyf[, c("Label", "Ratio", "Treatment")]
plots <- plots[, c("Label", "Site", "Vegetationtype")]


df<-merge(plots,Ratiohyf, by="Label")
df <- df[, !(names(df) == "Label")]
names(df)[names(df) == "Treatment"] <- "tt"


meanC <- aggregate(Ratio ~ Site, data = df[df$tt == "Control", ], FUN = mean , na.rm = TRUE)
names(meanC)[2] <- "meanC"
stdevC <- aggregate(Ratio ~ Site, data = df[df$tt == "Control", ], FUN= sd , na.rm = TRUE)
names(stdevC)[2] <- "stdevC"
nC <- aggregate(Ratio ~ Site, data = df[df$tt == "Control", ], FUN=length)
names(nC)[2] <- "nC"

meanOTC <- aggregate(Ratio ~ Site, data = df[df$tt == "OTC", ], FUN = mean , na.rm = TRUE)
names(meanOTC)[2] <- "meanOTC"
stdevOTC <- aggregate(Ratio ~ Site, data = df[df$tt == "OTC", ], FUN= sd , na.rm = TRUE)
names(stdevOTC)[2] <- "stdevOTC"
nOTC <- aggregate(Ratio ~ Site, data = df[df$tt == "OTC", ], FUN=length)
names(nOTC)[2] <- "nOTC"

dfC <- merge(meanC, stdevC, by = "Site")
dfC <- merge(dfC, nC, by = "Site")
dfC <- merge(df,dfC , by = "Site")
dfC <- dfC[, c("Site", "meanC", "stdevC", "nC")]
names(dfC)

dfOTC <- merge(meanOTC, stdevOTC, by = "Site")
dfOTC <- merge(dfOTC, nOTC, by = "Site")
dfOTC <- merge(df,dfOTC , by = "Site")
dfOTC <- dfOTC[, !(names(dfOTC) %in% c("Sep", "tt"))]

dat1 <- merge(dfOTC, dfC, by = "Site")
dat1$Ratio<-NULL
dat1<-unique(dat1)

OK<-escalc(measure= "SMD", m2i=meanC, sd2i=stdevC, n2i=nC, 
           m1i=meanOTC, sd1i=stdevOTC, n1i=nOTC, data=dat1)

overallresult <- rma(yi, vi, data=OK)
summary(overallresult)
forest.rma(overallresult, slab = OK$Site, addpred=F)
forest.rma(overallresult, slab = OK$Site, addpred=T, order=overallresult$yi, annotate=F, cex =1, psize = 1)   



OK <- OK[match(desired_order, OK$Site), ]

overallresult <- rma(yi, vi, data = OK, slab = OK$Site)

forest.rma(overallresult, addpred = TRUE, annotate=F,
           psize = 1, header = F,
           main = "Effectsize of warming on sept hyphal abundance ",
           cex = 1.5)


mod.gradeq <-rma(yi,vi,mods = ~ factor(Vegetationtype), data=OK)
mod.gradeq
