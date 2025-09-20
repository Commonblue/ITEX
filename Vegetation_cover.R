library(nlme);library(glmmTMB);library(ggplot2)

dfA<-read.table('Metadata_veg_prop.csv',header=T, sep=";", dec=".")
M<-read.table('Metadata_general.csv',header=T, sep=";", dec=",")
dfA[is.na(dfA)] <- 0

df<-merge(dfA,M,"Label")

str(df)
df$elevation<-as.numeric(df$elevation)
df$Time <- as.numeric(as.factor(df$Time))
df$Treatment <- as.factor(df$Treatment)


#shrubs--------------

ggplot(df, aes(x = Vegetationtype,
                   y = Shrub,
                   fill = Treatment)) +
  geom_boxplot(position = position_dodge()) +
  theme_minimal()

df<- subset(df,
                 !(Vegetationtype == "Alpine grassland"      & Shrub > 0.25)&
                   !(Vegetationtype == "Shrub snowbed"   & Treatment == "OTC" & Shrub > 0.72))




m0<-gls(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, method= "REML", data=df)


par(mfrow=c(3,4))
qqnorm(scale(resid(m0)));abline(0,1,col=2)
hist(resid(m0))
plot(resid(m0)~fitted(m0));abline(h=0,col=2)
plot(resid(m0)~df$Treatment);abline(h=0,col=2)
plot(resid(m0)~df$MAT);abline(h=0,col=2)
plot(resid(m0)~df$MAT);abline(h=0,col=2)
plot(resid(m0)~df$elevation);abline(h=0,col=2)
plot(resid(m0)~df$Lattitude);abline(h=0,col=2)
plot(resid(m0)~df$Time);abline(h=0,col=2)
shapiro.test(resid(m0))


m1<-lme(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2,random = ~1|Site, method= "REML", data=df)
anova(m0,m1)

m2<-lme(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype)))
m3<-lme(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment)))
#m4<-lme(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Site)))
m5<-lme(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment),varIdent(form=~1|Vegetationtype)))
#m6<-lme(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varPower(form=~fitted(.))))
m7<-lme(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~fitted(.))))
m8<-lme(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype), varExp(form=~fitted(.))))
#m9<-lme(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment), varPower(form=~fitted(.))))
#m10<-lme(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment),varIdent(form=~1|Vegetationtype),varPower(form=~fitted(.))))
m11<-lme(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~MAT)))
m12<-lme(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varPower(form=~MAP)))
m13<-lme(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~Lattitude)))
#m14<-lme(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~Elevation)))
m15<-lme(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~Time)))
anova(m1,m2,m3,m5,m7,m8,m11,m12,m13,m15)


par(mfrow=c(3,4))
qqnorm(scale(resid(m2,type="n")));abline(0,1,col=2)
hist(resid(m2,type="n"))
plot(resid(m2,type="n")~fitted(m2));abline(h=0,col=2)
plot(resid(m2,type="n")~df$Treatment);abline(h=0,col=2)
plot(resid(m2,type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m2,type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m2,type="n")~df$elevation);abline(h=0,col=2)
plot(resid(m2,type="n")~df$Lattitude);abline(h=0,col=2)
plot(resid(m2,type="n")~df$Time);abline(h=0,col=2)
shapiro.test(resid(m2,type="n"))



#m2
m2.ml<-lme(Shrub~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="ML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype)))
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:Lattitude)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:Time)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAP:Lattitude)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-elevation:Lattitude)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:MAP)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:elevation)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:MAT)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAT:Lattitude)
anova(m2.ml)



m2.reml<-update(m2.ml,method="REML")

par(mfrow=c(3,4))
qqnorm(scale(resid(m2.reml, type="n")));abline(0,1,col=2)
hist(resid(m2.reml, type="n"))
plot(resid(m2.reml, type="n")~fitted(m2.reml, type="n"));abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Treatment);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Vegetationtype);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$MAP);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Elevation);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Lattitude);abline(h=0,col=2)
shapiro.test(resid(m2.reml, type="n"))
summary(m2.reml)


mA<-gls(Shrub~Treatment*Vegetationtype, method= "REML", data=df)
summary(mA)


n <- nrow(df)
df$Shrub <- (df$Shrub * (n - 1) + 0.05) / n

mFr<- glmmTMB(Shrub~Treatment*Vegetationtype,beta_family(), data = df)
sim_res <- simulateResiduals(mFr)
par(mfrow=c(2,2))
plot(sim_res)        # overall diagnostic plots
plotResiduals(sim_res, form = df$Treatment)
plotResiduals(sim_res, form = df$Vegetationtype)
testDispersion(sim_res)   # overdispersion check
testZeroInflation(sim_res)  # zero-inflation check
summary(mFr)


#forbs--------------
ggplot(df, aes(x = Vegetationtype,
               y = Forb,
               fill = Treatment)) +
  geom_boxplot(position = position_dodge()) +
  theme_minimal()


df<- subset(df,
            !(Vegetationtype == "Heath"      & Treatment == "OTC" & Forb > 0.5)&
              !(Vegetationtype == "Heath"      & Treatment == "C" & Forb > 0.2)&
              !(Vegetationtype == "Shrub snowbed"      & Forb > 0.49)&
              !(Vegetationtype == "Wet meadow tundra"   & Treatment == "C" & Forb > 0.5))


m0<-gls(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, method= "REML", data=df)

par(mfrow=c(3,4))
qqnorm(scale(resid(m0)));abline(0,1,col=2)
hist(resid(m0))
plot(resid(m0)~fitted(m0));abline(h=0,col=2)
plot(resid(m0)~df$Treatment);abline(h=0,col=2)
plot(resid(m0)~df$MAT);abline(h=0,col=2)
plot(resid(m0)~df$MAT);abline(h=0,col=2)
plot(resid(m0)~df$elevation);abline(h=0,col=2)
plot(resid(m0)~df$Lattitude);abline(h=0,col=2)
plot(resid(m0)~df$Time);abline(h=0,col=2)
shapiro.test(resid(m0))


m1<-lme(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2,random = ~1|Site, method= "REML", data=df)
anova(m0,m1)

par(mfrow=c(3,4))
qqnorm(scale(resid(m1)));abline(0,1,col=2)
hist(resid(m1))
plot(resid(m1)~fitted(m1));abline(h=0,col=2)
plot(resid(m1)~df$Treatment);abline(h=0,col=2)
plot(resid(m1)~df$MAT);abline(h=0,col=2)
plot(resid(m1)~df$MAT);abline(h=0,col=2)
plot(resid(m1)~df$elevation);abline(h=0,col=2)
plot(resid(m1)~df$Lattitude);abline(h=0,col=2)
plot(resid(m1)~df$Time);abline(h=0,col=2)
shapiro.test(resid(m1))


m2<-lme(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype)))
m3<-lme(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment)))
#m4<-lme(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Site)))
m5<-lme(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment),varIdent(form=~1|Vegetationtype)))
#m6<-lme(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varPower(form=~fitted(.))))
m7<-lme(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~fitted(.))))
m8<-lme(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype), varExp(form=~fitted(.))))
#m9<-lme(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment), varPower(form=~fitted(.))))
#m10<-lme(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment),varIdent(form=~1|Vegetationtype),varPower(form=~fitted(.))))
m11<-lme(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~MAT)))
m12<-lme(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varPower(form=~MAP)))
m13<-lme(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~Lattitude)))
#m14<-lme(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~Elevation)))
m15<-lme(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~Time)))
anova(m1,m2,m3,m5,m7,m8,m11,m12,m13,m15)


par(mfrow=c(3,4))
qqnorm(scale(resid(m2,type="n")));abline(0,1,col=2)
hist(resid(m2,type="n"))
plot(resid(m2,type="n")~fitted(m2));abline(h=0,col=2)
plot(resid(m2,type="n")~df$Treatment);abline(h=0,col=2)
plot(resid(m2,type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m2,type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m2,type="n")~df$elevation);abline(h=0,col=2)
plot(resid(m2,type="n")~df$Lattitude);abline(h=0,col=2)
plot(resid(m2,type="n")~df$Time);abline(h=0,col=2)
shapiro.test(resid(m2,type="n"))


par(mfrow=c(3,4))
qqnorm(scale(resid(m7,type="n")));abline(0,1,col=2)
hist(resid(m7,type="n"))
plot(resid(m7,type="n")~fitted(m7,type="n"));abline(h=0,col=2)
plot(resid(m7,type="n")~df$Treatment);abline(h=0,col=2)
plot(resid(m7,type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m7,type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m7,type="n")~df$elevation);abline(h=0,col=2)
plot(resid(m7,type="n")~df$Lattitude);abline(h=0,col=2)
plot(resid(m7,type="n")~df$Time);abline(h=0,col=2)
shapiro.test(resid(m7,type="n"))

#m2, m7 not stable
m2.ml<-lme(Forb~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="ML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype)))
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Lattitude:elevation)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAT:Lattitude)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:elevation)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAP:elevation)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:MAP)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAP:Lattitude)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:MAT)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:Time)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:Lattitude)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Lattitude)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAT:elevation)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Time)
anova(m2.ml)




m2.reml<-update(m2.ml,method="REML")

par(mfrow=c(3,4))
qqnorm(scale(resid(m2.reml, type="n")));abline(0,1,col=2)
hist(resid(m2.reml, type="n"))
plot(resid(m2.reml, type="n")~fitted(m2.reml, type="n"));abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Treatment);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Vegetationtype);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$MAP);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Elevation);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Lattitude);abline(h=0,col=2)
shapiro.test(resid(m2.reml, type="n"))
library(effects)
plot(allEffects(m2.reml), multiline=T)
summary(m2.reml)
par(mfrow=c(1,1))

#Graminoid----------------------------------------
ggplot(df, aes(x = Vegetationtype,
               y = Graminoid,
               fill = Treatment)) +
  geom_boxplot(position = position_dodge()) +
  theme_minimal()


df<- subset(df,
              !(Vegetationtype == "Heath"      & Treatment == "OTC" & Graminoid > 0.25)&
              !(Vegetationtype == "Moss snowbed" & Treatment == "C"      & Graminoid > 0.15))


m0<-gls(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, method= "REML", data=df)

par(mfrow=c(3,4))
qqnorm(scale(resid(m0)));abline(0,1,col=2)
hist(resid(m0))
plot(resid(m0)~fitted(m0));abline(h=0,col=2)
plot(resid(m0)~df$Treatment);abline(h=0,col=2)
plot(resid(m0)~df$MAT);abline(h=0,col=2)
plot(resid(m0)~df$MAT);abline(h=0,col=2)
plot(resid(m0)~df$elevation);abline(h=0,col=2)
plot(resid(m0)~df$Lattitude);abline(h=0,col=2)
plot(resid(m0)~df$Time);abline(h=0,col=2)
shapiro.test(resid(m0))


m1<-lme(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2,random = ~1|Site, method= "REML", data=df)
anova(m0,m1)

par(mfrow=c(3,4))
qqnorm(scale(resid(m1)));abline(0,1,col=2)
hist(resid(m1))
plot(resid(m1)~fitted(m1));abline(h=0,col=2)
plot(resid(m1)~df$Treatment);abline(h=0,col=2)
plot(resid(m1)~df$MAT);abline(h=0,col=2)
plot(resid(m1)~df$MAT);abline(h=0,col=2)
plot(resid(m1)~df$elevation);abline(h=0,col=2)
plot(resid(m1)~df$Lattitude);abline(h=0,col=2)
plot(resid(m1)~df$Time);abline(h=0,col=2)
shapiro.test(resid(m1))


m2<-lme(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype)))
m3<-lme(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment)))
#m4<-lme(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Site)))
m5<-lme(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment),varIdent(form=~1|Vegetationtype)))
#m6<-lme(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varPower(form=~fitted(.))))
m7<-lme(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~fitted(.))))
m8<-lme(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype), varExp(form=~fitted(.))))
#m9<-lme(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment), varPower(form=~fitted(.))))
#m10<-lme(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment),varIdent(form=~1|Vegetationtype),varPower(form=~fitted(.))))
m11<-lme(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~MAT)))
m12<-lme(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varPower(form=~MAP)))
m13<-lme(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~Lattitude)))
#m14<-lme(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~Elevation)))
m15<-lme(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~Time)))
anova(m1,m2,m3,m5,m7,m8,m11,m12,m13,m15)


par(mfrow=c(3,4))
qqnorm(scale(resid(m2,type="n")));abline(0,1,col=2)
hist(resid(m2,type="n"))
plot(resid(m2,type="n")~fitted(m2,type="n"));abline(h=0,col=2)
plot(resid(m2,type="n")~df$Treatment);abline(h=0,col=2)
plot(resid(m2,type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m2,type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m2,type="n")~df$elevation);abline(h=0,col=2)
plot(resid(m2,type="n")~df$Lattitude);abline(h=0,col=2)
plot(resid(m2,type="n")~df$Time);abline(h=0,col=2)
shapiro.test(resid(m2,type="n"))


m2.ml<-lme(Graminoid~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="ML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype)))
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Lattitude:elevation)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:MAT)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:Time)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAT:Lattitude)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:elevation)
anova(m2.ml)

m2.reml<-update(m2.ml,method="REML")

par(mfrow=c(3,4))
qqnorm(scale(resid(m2.reml, type="n")));abline(0,1,col=2)
hist(resid(m2.reml, type="n"))
plot(resid(m2.reml, type="n")~fitted(m2.reml, type="n"));abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df_sep$Treatment);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df_sep$Vegetationtype);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df_sep$MAP);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df_sep$MAT);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df_sep$Elevation);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df_sep$Lattitude);abline(h=0,col=2)
shapiro.test(resid(m2.reml, type="n"))
library(effects)
plot(allEffects(m2.reml), multiline=T)
summary(m2.reml)
par(mfrow=c(1,1))


#Bryophyte----------------
ggplot(df, aes(x = Vegetationtype,
               y = Bryophyte,
               fill = Treatment)) +
  geom_boxplot(position = position_dodge()) +
  theme_minimal()

#from here
df<- subset(df,
            !(Vegetationtype == "Alpine grassland"       & Bryophyte > 0.1)&
              !(Vegetationtype == "Heath" & Treatment == "OTC"      & Bryophyte > 0.4)&
              !(Vegetationtype == "Moss snowbed" & Treatment == "C"      & Bryophyte > 0.8)&
              !(Vegetationtype == "Shrub snowbed" & Treatment == "OTC"      & Bryophyte > 0.3))


m0<-gls(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, method= "REML", data=df)

par(mfrow=c(3,4))
qqnorm(scale(resid(m0)));abline(0,1,col=2)
hist(resid(m0))
plot(resid(m0)~fitted(m0));abline(h=0,col=2)
plot(resid(m0)~df$Treatment);abline(h=0,col=2)
plot(resid(m0)~df$MAT);abline(h=0,col=2)
plot(resid(m0)~df$MAT);abline(h=0,col=2)
plot(resid(m0)~df$elevation);abline(h=0,col=2)
plot(resid(m0)~df$Lattitude);abline(h=0,col=2)
plot(resid(m0)~df$Time);abline(h=0,col=2)
shapiro.test(resid(m0))


m1<-lme(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2,random = ~1|Site, method= "REML", data=df)
anova(m0,m1)

par(mfrow=c(3,4))
qqnorm(scale(resid(m1)));abline(0,1,col=2)
hist(resid(m1))
plot(resid(m1)~fitted(m1));abline(h=0,col=2)
plot(resid(m1)~df$Treatment);abline(h=0,col=2)
plot(resid(m1)~df$MAT);abline(h=0,col=2)
plot(resid(m1)~df$MAT);abline(h=0,col=2)
plot(resid(m1)~df$elevation);abline(h=0,col=2)
plot(resid(m1)~df$Lattitude);abline(h=0,col=2)
plot(resid(m1)~df$Time);abline(h=0,col=2)
shapiro.test(resid(m1))


m2<-lme(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype)))
m3<-lme(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment)))
#m4<-lme(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Site)))
m5<-lme(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment),varIdent(form=~1|Vegetationtype)))
#m6<-lme(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varPower(form=~fitted(.))))
m7<-lme(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~fitted(.))))
m8<-lme(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype), varExp(form=~fitted(.))))
#m9<-lme(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment), varPower(form=~fitted(.))))
#m10<-lme(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment),varIdent(form=~1|Vegetationtype),varPower(form=~fitted(.))))
m11<-lme(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~MAT)))
m12<-lme(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varPower(form=~MAP)))
m13<-lme(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~Lattitude)))
#m14<-lme(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~Elevation)))
m15<-lme(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~Time)))
anova(m1,m2,m3,m5,m7,m8,m11,m12,m13,m15)


par(mfrow=c(3,4))
qqnorm(scale(resid(m2,type="n")));abline(0,1,col=2)
hist(resid(m2,type="n"))
plot(resid(m2,type="n")~fitted(m2,type="n"));abline(h=0,col=2)
plot(resid(m2,type="n")~df$Treatment);abline(h=0,col=2)
plot(resid(m2,type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m2,type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m2,type="n")~df$elevation);abline(h=0,col=2)
plot(resid(m2,type="n")~df$Lattitude);abline(h=0,col=2)
plot(resid(m2,type="n")~df$Time);abline(h=0,col=2)
shapiro.test(resid(m2,type="n"))


m2.ml<-lme(Bryophyte~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="ML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype)))
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:Time)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAP:MAT)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAT:Lattitude)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAP:Lattitude)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAP:elevation)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Lattitude:elevation)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:MAT)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAT:elevation)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:MAP)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:Lattitude)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:elevation)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAP)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Time)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-elevation)
anova(m2.ml)

m2.reml<-update(m2.ml,method="REML")
anova(m2.reml)


par(mfrow=c(3,4))
qqnorm(scale(resid(m2.reml, type="n")));abline(0,1,col=2)
hist(resid(m2.reml, type="n"))
plot(resid(m2.reml, type="n")~fitted(m2.reml, type="n"));abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Treatment);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Vegetationtype);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$MAP);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Elevation);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Lattitude);abline(h=0,col=2)
shapiro.test(resid(m2.reml, type="n"))
library(effects)
plot(allEffects(m2.reml), multiline=T)
summary(m2.reml)
par(mfrow=c(1,1))
#Lichen --------------------------------------
ggplot(df, aes(x = Vegetationtype,
               y = Lichen,
               fill = Treatment)) +
  geom_boxplot(position = position_dodge()) +
  theme_minimal()

#from here
df<- subset(df,
            !(Vegetationtype == "Alpine grassland"       & Lichen > 0.10)&
              !(Vegetationtype == "Heath"  & Treatment == "C"  & Lichen > 0.3)&
            !(Vegetationtype == "Shrub snowbed"       & Lichen > 0.15)&
              !(Site %in% c("15.ch","20.ch", "wl.ch", "ra.sp", "nr.us", "ul.sp",  "aw.us", "yk.low", "pc.fr")))


m0<-gls(Lichen~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, method= "REML", data=df)

par(mfrow=c(3,4))
qqnorm(scale(resid(m0)));abline(0,1,col=2)
hist(resid(m0))
plot(resid(m0)~fitted(m0));abline(h=0,col=2)
plot(resid(m0)~df$Treatment);abline(h=0,col=2)
plot(resid(m0)~df$MAT);abline(h=0,col=2)
plot(resid(m0)~df$MAT);abline(h=0,col=2)
plot(resid(m0)~df$elevation);abline(h=0,col=2)
plot(resid(m0)~df$Lattitude);abline(h=0,col=2)
plot(resid(m0)~df$Time);abline(h=0,col=2)
shapiro.test(resid(m0))
#not fixable


n <- nrow(df)
df$Lichen <- (df$Lichen * (n - 1) + 0.5) / n


m0<- glmmTMB(Lichen~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2,family = "beta",ziformula = ~ 1, dispformula = ~Vegetationtype, data = df)
m0<- glmmTMB(Lichen~Treatment+MAT+MAP+Lattitude+elevation+Time+Vegetationtype ,ziformula = ~ 1,family = "beta", data = df)

#deleting till something functions
#only one without issues
m0<- glmmTMB(Lichen~Treatment+Vegetationtype ,ziformula = ~ 1,family = "beta", data = df)
par(mfrow=c(3,4))


library(DHARMa)
sim_res <- simulateResiduals(m0)
par(mfrow=c(2,2))
plot(sim_res)        # 
plotResiduals(sim_res, form = df$Treatment)
plotResiduals(sim_res, form = df$Vegetationtype)
testDispersion(sim_res) 
testZeroInflation(sim_res)  
#ok
summary(m0)




#Bare soil-------------------------
ggplot(df, aes(x = Vegetationtype,
               y = Bare_Rock_Litter,
               fill = Treatment)) +
  geom_boxplot(position = position_dodge()) +
  theme_minimal()

df<- subset(df,
            !(Vegetationtype == "Alpine grassland"      & Bare_Rock_Litter > 0.5)&
              !(Vegetationtype == "Moss snowbed"      & Bare_Rock_Litter > 0.2)&
              !(Vegetationtype == "Moss snowbed"      & Treatment == "C" & Bare_Rock_Litter > 0.1)&
              !(Vegetationtype == "Alpine grassland"  & Treatment == "C"    & Bare_Rock_Litter > 0.3)&
              !(Vegetationtype == "Heath"   &  Bare_Rock_Litter > 0.1)&
              !(Site %in% c("15.ch", "aw.us", "fn.th.no", "ki.fi", "ko.gl", "yk.up", "yk.low")))




m0<-gls(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, method= "REML", data=df)


par(mfrow=c(3,4))
qqnorm(scale(resid(m0)));abline(0,1,col=2)
hist(resid(m0))
plot(resid(m0)~fitted(m0));abline(h=0,col=2)
plot(resid(m0)~df$Treatment);abline(h=0,col=2)
plot(resid(m0)~df$MAT);abline(h=0,col=2)
plot(resid(m0)~df$MAT);abline(h=0,col=2)
plot(resid(m0)~df$elevation);abline(h=0,col=2)
plot(resid(m0)~df$Lattitude);abline(h=0,col=2)
plot(resid(m0)~df$Time);abline(h=0,col=2)
shapiro.test(resid(m0))


m1<-lme(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2,random = ~1|Site, method= "REML", data=df)
anova(m0,m1)
par(mfrow=c(3,4))
qqnorm(scale(resid(m1)));abline(0,1,col=2)
hist(resid(m1))
plot(resid(m1)~fitted(m1));abline(h=0,col=2)
plot(resid(m1)~df$Treatment);abline(h=0,col=2)
plot(resid(m1)~df$MAT);abline(h=0,col=2)
plot(resid(m1)~df$MAT);abline(h=0,col=2)
plot(resid(m1)~df$elevation);abline(h=0,col=2)
plot(resid(m1)~df$Lattitude);abline(h=0,col=2)
plot(resid(m1)~df$Time);abline(h=0,col=2)
shapiro.test(resid(m1))


m2<-lme(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype)))
m3<-lme(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment)))
#m4<-lme(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Site)))
m5<-lme(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment),varIdent(form=~1|Vegetationtype)))
#m6<-lme(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varPower(form=~fitted(.))))
#m7<-lme(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~fitted(.))))
#m8<-lme(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype), varExp(form=~fitted(.))))
#m9<-lme(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment), varPower(form=~fitted(.))))
#m10<-lme(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment),varIdent(form=~1|Vegetationtype),varPower(form=~fitted(.))))
m11<-lme(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~MAT)))
m12<-lme(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varPower(form=~MAP)))
m13<-lme(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~Lattitude)))
#m14<-lme(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~Elevation)))
m15<-lme(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="REML", data=df,weights = varComb(varExp(form=~Time)))
anova(m1,m2,m3,m5,m11,m12,m13,m15)


par(mfrow=c(3,4))
qqnorm(scale(resid(m2, type="n")));abline(0,1,col=2)
hist(resid(m2, type="n"))
plot(resid(m2, type="n")~fitted(m2, type="n"));abline(h=0,col=2)
plot(resid(m2, type="n")~df$Treatment);abline(h=0,col=2)
plot(resid(m2, type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m2, type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m2, type="n")~df$elevation);abline(h=0,col=2)
plot(resid(m2, type="n")~df$Lattitude);abline(h=0,col=2)
plot(resid(m2, type="n")~df$Time);abline(h=0,col=2)
shapiro.test(resid(m2, type="n"))



#m2
m2.ml<-lme(Bare_Rock_Litter~Treatment*(MAT+MAP+Lattitude+elevation+Time+Vegetationtype) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site, method="ML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype)))
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:Lattitude)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:Time)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAT:elevation)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:elevation)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAP:Lattitude)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Lattitude:elevation)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:Vegetationtype)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:MAP)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment:MAT)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAP:elevation)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAT:Lattitude)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-elevation)
anova(m2.ml)


m2.reml<-update(m2.ml,method="REML")


par(mfrow=c(3,4))
qqnorm(scale(resid(m2.reml, type="n")));abline(0,1,col=2)
hist(resid(m2.reml, type="n"))
plot(resid(m2.reml, type="n")~fitted(m2.reml, type="n"));abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Treatment);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Vegetationtype);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$MAP);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Elevation);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Lattitude);abline(h=0,col=2)
shapiro.test(resid(m2.reml, type="n"))
library(effects)
plot(allEffects(m2.reml), multiline=T)
summary(m2.reml)
par(mfrow=c(1,1))

#visualisation---------------------------------
library(tidyr)
library(dplyr)

df_long <- df %>%
  pivot_longer(cols = c(Shrub, Forb, Graminoid, Lichen, Bryophyte, Bare_Rock_Litter),
               names_to = "Vegetation",
               values_to = "Cover")


ggplot(df_long, aes(x = Treatment, y = Cover, fill = Vegetation)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Vegetationtype) +  # if you want to separate vegetation groups
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(x = "Treatment",
       y = "Percent cover",
       fill = "Vegetation type")
