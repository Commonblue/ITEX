library(ggplot2);library(MASS);library(nlme)

df<-read.table('Hyphal_length_analyses_ITEX.csv', dec=",", sep=";", header=T)
A<-read.table('Metadata_general.csv', dec=",", sep=";", header=T)
df<-merge(df,A,"Label")
df[] <- lapply(df, function(x) if (is.character(x)) as.factor(x) else x)
str(df)


par(mfrow=c(1,1))
boxplot(Hyphal_length_uns~Treatment.x*Vegetationtype, data=df)

ggplot(df, aes(x = Vegetationtype,
               y = log1p(Hyphal_length_uns),
               fill = Treatment.x)) +
  geom_boxplot(position = position_dodge()) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()



tapply(log1p(df$Hyphal_length_uns),df$Treatment.x,mean);tapply(log1p(df$Hyphal_length_uns),df$Treatment.x,function(x) sd(x) / mean(x))
tapply(log1p(df$Hyphal_length_uns),df$Vegetationtype,mean);tapply(log1p(df$Hyphal_length_uns),df$Vegetationtype,function(x) sd(x) / mean(x))


m0<-gls(log1p(Hyphal_length_uns)~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, method="REML", data=df)
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



m1<-lme(log1p(Hyphal_length_uns)~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df)
par(mfrow=c(3,4))
qqnorm(scale(resid(m1)));abline(0,1,col=2)
hist(resid(m1))
plot(resid(m1)~fitted(m1));abline(h=0,col=2)
plot(resid(m1)~df$Treatment.x);abline(h=0,col=2)
plot(resid(m1)~df$Vegetationtype);abline(h=0,col=2)
plot(resid(m1)~df$MAT);abline(h=0,col=2)
plot(resid(m1)~df$MAP);abline(h=0,col=2)
plot(resid(m1)~df$elevation);abline(h=0,col=2)
plot(resid(m1)~df$Lattitude);abline(h=0,col=2)
plot(resid(m1)~df$Time);abline(h=0,col=2)
shapiro.test(resid(m1))
anova(m0,m1)



m2<-lme(log1p(Hyphal_length_uns)~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype)))
m3<-lme(log1p(Hyphal_length_uns)~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment.x)))
#m4<-lme(log1p(Hyphal_length_uns)~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varIdent(form=~1|Site.x)))
m5<-lme(log1p(Hyphal_length_uns)~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment.x),varIdent(form=~1|Vegetationtype)))
m6<-lme(log1p(Hyphal_length_uns)~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varPower(form=~fitted(.))))
m7<-lme(log1p(Hyphal_length_uns)~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varExp(form=~fitted(.))))
m8<-lme(log1p(Hyphal_length_uns)~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varExp(form=~MAT)))
m9<-lme(log1p(Hyphal_length_uns)~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varPower(form=~MAP)))
m10<-lme(log1p(Hyphal_length_uns)~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varExp(form=~Lattitude)))
m11<-lme(log1p(Hyphal_length_uns)~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varPower(form=~elevation)))
m12<-lme(log1p(Hyphal_length_uns)~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="REML", data=df,weights = varComb(varPower(form=~Time)))
anova(m0,m1,m2,m3,m5,m6,m7,m8,m9,m10,m11,m12)
#m5,m15

par(mfrow=c(3,4))
qqnorm(scale(resid(m5, type="n")));abline(0,1,col=2)
hist(resid(m5, type="n"))
plot(resid(m5, type="n")~fitted(m5, type="n"));abline(h=0,col=2)
plot(resid(m5, type="n")~df$Treatment.x);abline(h=0,col=2)
plot(resid(m5, type="n")~df$Vegetationtype);abline(h=0,col=2)
plot(resid(m5, type="n")~df$MAP);abline(h=0,col=2)
plot(resid(m5, type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m5, type="n")~df$elevation);abline(h=0,col=2)
plot(resid(m5, type="n")~df$Lattitude);abline(h=0,col=2)
plot(resid(m5, type="n")~df$Time);abline(h=0,col=2)
shapiro.test(resid(m5, type="n"))
anova(m5)


#m5
m5.ml<-lme(log1p(Hyphal_length_uns)~Treatment.x*(Vegetationtype+MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2, random = ~1|Site.x, method="ML", data=df,weights = varComb(varIdent(form=~1|Treatment.x),varIdent(form=~1|Vegetationtype)))
anova(m5.ml)
m5.ml<-update(m5.ml,~.-Treatment.x:elevation)
anova(m5.ml)
m5.ml<-update(m5.ml,~.-Lattitude:elevation)
anova(m5.ml)
m5.ml<-update(m5.ml,~.-MAP:MAT)
anova(m5.ml)
m5.ml<-update(m5.ml,~.-Treatment.x:MAP)
anova(m5.ml)
m5.ml<-update(m5.ml,~.-MAT:Lattitude)
anova(m5.ml)
m5.ml<-update(m5.ml,~.-Treatment.x:MAT)
anova(m5.ml)
m5.ml<-update(m5.ml,~.-Treatment.x:Lattitude)
anova(m5.ml)
m5.ml<-update(m5.ml,~.-MAT:elevation)
anova(m5.ml)
m5.ml<-update(m5.ml,~.-Treatment.x:Time)
anova(m5.ml)
m5.ml<-update(m5.ml,~.-MAP:elevation)
anova(m5.ml)
m5.ml<-update(m5.ml,~.-MAP:Lattitude)
anova(m5.ml)
m5.ml<-update(m5.ml,~.-Lattitude)
anova(m5.ml)


m5.reml<-update(m5.ml,method="REML")

par(mfrow=c(3,4))
qqnorm(scale(resid(m5.reml, type="n")));abline(0,1,col=2)
hist(resid(m5.reml, type="n"))
plot(resid(m5.reml, type="n")~fitted(m5.reml, type="n"));abline(h=0,col=2)
plot(resid(m5.reml, type="n")~df$Treatment.x);abline(h=0,col=2)
plot(resid(m5.reml, type="n")~df$Vegetationtype);abline(h=0,col=2)
plot(resid(m5.reml, type="n")~df$MAP);abline(h=0,col=2)
plot(resid(m5.reml, type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m5.reml, type="n")~df$elevation);abline(h=0,col=2)
plot(resid(m5.reml, type="n")~df$Lattitude);abline(h=0,col=2)
shapiro.test(resid(m5.reml, type="n"))
anova(m5.reml)

summary(m5.reml)
library(effects)
plot(allEffects(m5.reml), multiline=T)


