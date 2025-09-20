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

df$Site[df$Site == "Pga1"] <- "cc.it"
df$Site[df$Site == "Pga2"] <- "ps.it"
df$Site[df$Site == "Pga3"] <- "sa.it"
df$Site[df$Site == "fn.pa.no"] <- "pa.no"
df$Site[df$Site == "fn.th.no"] <- "th.no"
df$Site<-as.factor(df$Site)
df[] <- lapply(df, function(x) if (is.character(x) | is.integer(x)) as.factor(x) else x)


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



desired_order <- c("15.ch","20.ch","25.ch","nr.us","cc.it","vb.ch", "fi.no", 
                   "ko.gl", "wh.ch", "wl.ch", "yk.up", "th.no", "ki.fi", 
                   "ps.it", "sk.no", "sv.no", "ad.ad", "pa.no", "pc.fr", "sa.it",
                   "ra.sp", "ul.sp", "aw.us", "ul.no", "yk.low")

OK <- OK[match(desired_order, OK$Site), ]

overallresult <- rma(yi, vi, data = OK, slab = OK$Site)

forest.rma(overallresult, addpred = TRUE, annotate=F,
           psize = 1, header = F,
           main = "Effectsize of warming on sept hyphal abundance ",
           cex = 1.5)


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
#3 cont sites salix
plots <- plots[, c("Label", "Site", "Vegetationtype")]
#5cont???

df<-merge(plots,unsephyf, by="Label", drop=F)
df <- df[, !(names(df) == "Label")]
names(df)[names(df) == "Treatment"] <- "tt"



df$Site[df$Site == "Pga1"] <- "cc.it"
df$Site[df$Site == "Pga2"] <- "ps.it"
df$Site[df$Site == "Pga3"] <- "sa.it"
df$Site[df$Site == "fn.pa.no"] <- "pa.no"
df$Site[df$Site == "fn.th.no"] <- "th.no"
df$Site<-as.factor(df$Site)
#df$Unsep<-log1p(df$Unsep)

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
#################################
#analyses
################################
OK<-escalc(measure= "SMD", m2i=meanC, sd2i=stdevC, n2i=nC, 
           m1i=meanOTC, sd1i=stdevOTC, n1i=nOTC, data=dat1)
#*k = 25 studies
#*
#*tau^2 = estimate of var of true effect = var between study heterogeneity
#*est: 0.245 = effect size
#*pval = sign!
#*
#*test for heterogeneity, is there sign heterogeneity? p = 0.1 --> NO
#*there is no sign heterogeneity
#*I^2 = variation due to true variability = 23% good
#*

overallresult <- rma(yi, vi, data=OK)
summary(overallresult)
forest.rma(overallresult, slab = OK$Site, addpred=T)
forest.rma(overallresult, slab = OK$Site, addpred=T, order=overallresult$yi, annotate=F)

desired_order <- c("15.ch","20.ch","25.ch","nr.us","cc.it","vb.ch", "fi.no", 
                   "ko.gl", "wh.ch", "wl.ch", "yk.up", "th.no", "ki.fi", 
                   "ps.it", "sk.no", "sv.no", "ad.ad", "pa.no", "pc.fr", "sa.it",
                   "ra.sp", "ul.sp", "aw.us", "ul.no", "yk.low")

# reorder your dataset
OK <- OK[match(desired_order, OK$Site), ]

# refit model on reordered data
overallresult <- rma(yi, vi, data = OK)

overallresult <- rma(yi, vi, data = OK, slab = OK$Site)

forest.rma(overallresult, addpred = TRUE, annotate=F,
           psize = 1, header = F,
           main = "Effectsize of warming on AMF hyphal abundance ",
           cex = 1.5)


png("forest_plot.png", width = 4000, height = 5200, res = 600)
forest.rma(overallresult, 
           addpred = TRUE, 
           annotate = FALSE,
           psize = 1, 
           header = FALSE,
           main = "Effect size of warming on AMF hyphal abundance",
           cex = 1.5)

dev.off()

#no inf * probably no outlier or sign influence
#################################
#categorical moderator analyses
################################
#*calculate Qbetween
#*are there sign siffernces between levels of our moderator
#*is there a sign differance based on a category?
mod.gradeq <-rma(yi,vi,mods = ~ factor(Vegetationtype), data=OK)
fitted.rma(mod.gradeq)
predict.rma(mod.gradeq)
#categorisch op de vegetatie
mod.gradeq
influence(mod.gradeq)
#*test for moderators
#*no sign difference between levels.
#*we have one factor that is significant it is not real? (dont understand)
#*because there are no sign differences between levels
#*
#*test for residual heterogeneity
#*p is not significant, there is no unexplained variance
#####################################################


####################################################
# Save QM Test and write it into a text file
mod.gradeq <-rma(yi,vi,mods = ~ factor(Site), data=OK)
#categorisch op de vegetatie
mod.gradeq

Qgrade_collapsed_string <- paste(mod.gradeq[["QMdf"]])
Qgrade_type1 <- data.frame(CollapsedQMdf = Qgrade_collapsed_string)
Qgrade_type2 <- round(mod.gradeq$QM,2)
Qgrade_type3 <- round(mod.gradeq$QMp,3)
QgradeQ <- paste(
  "Qb(",Qgrade_collapsed_string,") =", 
  Qgrade_type2,
  ", p =", 
  Qgrade_type3,
  collapse = " "
)
cat(QgradeQ, "\n")
gradeQtest <- data.frame(Text = QgradeQ)
write.table(gradeQtest, file = "QgradeQ.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

#writing results for a paper