library(ggplot2);library(MASS);library(nlme)

df<-read.table('Hyphal_length_analyses_ITEX.csv', dec=",", sep=";", header=T)
A<-read.table('Metadata_general.csv', dec=",", sep=";", header=T)
df<-merge(df,A,"Label")
df[] <- lapply(df, function(x) if (is.character(x)) as.factor(x) else x)


sum(df$Hyphal_length_sept)
sum(df$Hyphal_length_uns)

df$Ratio<-(df$Hyphal_length_uns/(df$Hyphal_length_sep+df$Hyphal_length_uns))*100
df <- na.omit(df)


ggplot(df, aes(x = Vegetationtype,
               y = Ratio,
               fill = Treatment.x)) +
  geom_boxplot(position = position_dodge()) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()


tapply(log1p(df$Ratio),df$Treatment.x,mean);tapply(log1p(df$Ratio),df$Treatment.x,function(x) sd(x) / mean(x))
tapply(log1p(df$Ratio),df$Vegetationtype,mean);tapply(log1p(df$Ratio),df$Vegetationtype,function(x) sd(x) / mean(x))


m0<-gls(Ratio~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, method="REML", data=df)
par(mfrow=c(3,4))
qqnorm(scale(resid(m0)));abline(0,1,col=2)
hist(resid(m0))
plot(resid(m0)~fitted(m0));abline(h=0,col=2)
plot(resid(m0)~df$Treatment.x);abline(h=0,col=2)
plot(resid(m0)~df$Vegetationtype);abline(h=0,col=2)
plot(resid(m0)~df$MAT);abline(h=0,col=2)
plot(resid(m0)~df$MAP);abline(h=0,col=2)
plot(resid(m0)~df$elevation);abline(h=0,col=2)
plot(resid(m0)~df$Lattitude);abline(h=0,col=2)
plot(resid(m0)~df$Time);abline(h=0,col=2)
shapiro.test(resid(m0))
anova(m0)



m1<-lme(Ratio~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df)
shapiro.test(resid(m1))
anova(m0,m1)



m2<-lme(Ratio~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype)))
m3<-lme(Ratio~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment.x)))
#m4<-lme(Ratio~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varIdent(form=~1|Site.x)))
m5<-lme(Ratio~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment.x),varIdent(form=~1|Vegetationtype)))
m6<-lme(Ratio~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varPower(form=~fitted(.))))
m7<-lme(Ratio~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varExp(form=~fitted(.))))
m8<-lme(Ratio~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varExp(form=~MAT)))
m9<-lme(Ratio~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varPower(form=~MAP)))
m10<-lme(Ratio~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varExp(form=~Lattitude)))
m11<-lme(Ratio~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varPower(form=~elevation)))
m12<-lme(Ratio~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varPower(form=~Time)))
anova(m0,m1,m2,m3,m5,m6,m7,m8,m9, m10, m11,m12)
#m5,m15

par(mfrow=c(3,4))
qqnorm(scale(resid(m12, type="n")));abline(0,1,col=2)
hist(resid(m12, type="n"))
plot(resid(m12, type="n")~fitted(m12, type="n"));abline(h=0,col=2)
plot(resid(m12, type="n")~df$Treatment.x);abline(h=0,col=2)
plot(resid(m12, type="n")~df$Vegetationtype);abline(h=0,col=2)
plot(resid(m12, type="n")~df$MAP);abline(h=0,col=2)
plot(resid(m12, type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m12, type="n")~df$elevation);abline(h=0,col=2)
plot(resid(m12, type="n")~df$Lattitude);abline(h=0,col=2)
plot(resid(m12, type="n")~df$Time);abline(h=0,col=2)
shapiro.test(resid(m12, type="n"))
anova(m12)



m12.ml<-lme(Ratio~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="ML", data=df,weights = varComb(varPower(form=~Time)))
anova(m12.ml)
m12.ml<-update(m12.ml,~.-MAP:MAT)
anova(m12.ml)
m12.ml<-update(m12.ml,~.-Treatment.x:MAT)
anova(m12.ml)
m12.ml<-update(m12.ml,~.-Treatment.x:MAP)
anova(m12.ml)
m12.ml<-update(m12.ml,~.-Treatment.x:Lattitude)
anova(m12.ml)
m12.ml<-update(m12.ml,~.-Treatment.x:elevation)
anova(m12.ml)
m12.ml<-update(m12.ml,~.-MAP:elevation)
anova(m12.ml)
m12.ml<-update(m12.ml,~.-MAT:Lattitude)
anova(m12.ml)
m12.ml<-update(m12.ml,~.-MAT:elevation)
anova(m12.ml)
m12.ml<-update(m12.ml,~.-Treatment.x:Time)
anova(m12.ml)
m12.ml<-update(m12.ml,~.-Lattitude:elevation)
anova(m12.ml)
m12.ml<-update(m12.ml,~.-Time)
anova(m12.ml)


m12.reml<-update(m12.ml,method="REML")

par(mfrow=c(3,4))
qqnorm(scale(resid(m12.reml, type="n")));abline(0,1,col=2)
hist(resid(m12.reml, type="n"))
plot(resid(m12.reml, type="n")~fitted(m12.reml, type="n"));abline(h=0,col=2)
plot(resid(m12.reml, type="n")~df$Treatment.x);abline(h=0,col=2)
plot(resid(m12.reml, type="n")~df$Vegetationtype);abline(h=0,col=2)
plot(resid(m12.reml, type="n")~df$MAP);abline(h=0,col=2)
plot(resid(m12.reml, type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m12.reml, type="n")~df$elevation);abline(h=0,col=2)
plot(resid(m12.reml, type="n")~df$Lattitude);abline(h=0,col=2)
shapiro.test(resid(m12.reml, type="n"))
anova(m12.reml)

summary(m12.reml)
library(effects)
plot(allEffects(m12.reml), multiline=T)


