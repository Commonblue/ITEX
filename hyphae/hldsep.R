library(ggplot2);library(MASS);library(nlme)

df<-read.table('Hyphal_length_analyses_ITEX.csv', dec=",", sep=";", header=T)
A<-read.table('Metadata_general.csv', dec=",", sep=";", header=T)
df<-merge(df,A,"Label")
df[] <- lapply(df, function(x) if (is.character(x)) as.factor(x) else x)

#Unsep------------------

boxplot(Hyphal_length_sept~Treatment.x*Vegetationtype, data=df)


ggplot(df, aes(x = Vegetationtype,
               y = log1p(Hyphal_length_sept),
               fill = Treatment.x)) +
  geom_boxplot(position = position_dodge()) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()

tapply(log1p(df$Hyphal_length_sept),df$Treatment.x,mean);tapply(log1p(df$Hyphal_length_sept),df$Treatment.x,function(x) sd(x) / mean(x))
tapply(log1p(df$Hyphal_length_sept),df$Vegetationtype,mean);tapply(log1p(df$Hyphal_length_sept),df$Vegetationtype,function(x) sd(x) / mean(x))

par(mfrow=c(2,2))
qqnorm(log1p(df$Hyphal_length_sept))
qqline(log1p(df$Hyphal_length_sept))
hist(log1p(df$Hyphal_length_sept))



m0<-gls(log1p(Hyphal_length_sept)~Treatment.x*(MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2 + (Treatment.x*Vegetationtype*Time), method="REML", data=df)
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



m1<-lme(log1p(Hyphal_length_sept)~Treatment.x*(MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2 + (Treatment.x*Vegetationtype*Time), random = ~1|Site.x, method="REML", data=df)
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



m2<-lme(log1p(Hyphal_length_sept)~Treatment.x*(MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2 + (Treatment.x*Vegetationtype*Time), random = ~1|Site.x, method="REML", data=df,weights = varComb(varIdent(form=~1|Vegetationtype)))
m3<-lme(log1p(Hyphal_length_sept)~Treatment.x*(MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2 + (Treatment.x*Vegetationtype*Time), random = ~1|Site.x, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment.x)))
#m4<-lme(log1p(Hyphal_length_sept)~Treatment.x*(MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2 + (Treatment.x*Vegetationtype*Time), random = ~1|Site.x, method="REML", data=df,weights = varComb(varIdent(form=~1|Site.x)))
m5<-lme(log1p(Hyphal_length_sept)~Treatment.x*(MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2 + (Treatment.x*Vegetationtype*Time), random = ~1|Site.x, method="REML", data=df,weights = varComb(varIdent(form=~1|Treatment.x),varIdent(form=~1|Vegetationtype)))
m6<-lme(log1p(Hyphal_length_sept)~Treatment.x*(MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2 + (Treatment.x*Vegetationtype*Time), random = ~1|Site.x, method="REML", data=df,weights = varComb(varPower(form=~fitted(.))))
m7<-lme(log1p(Hyphal_length_sept)~Treatment.x*(MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2 + (Treatment.x*Vegetationtype*Time), random = ~1|Site.x, method="REML", data=df,weights = varComb(varExp(form=~fitted(.))))
m8<-lme(log1p(Hyphal_length_sept)~Treatment.x*(MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2 + (Treatment.x*Vegetationtype*Time), random = ~1|Site.x, method="REML", data=df,weights = varComb(varExp(form=~MAT)))
m9<-lme(log1p(Hyphal_length_sept)~Treatment.x*(MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2 + (Treatment.x*Vegetationtype*Time), random = ~1|Site.x, method="REML", data=df,weights = varComb(varPower(form=~MAP)))
m10<-lme(log1p(Hyphal_length_sept)~Treatment.x*(MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2 + (Treatment.x*Vegetationtype*Time), random = ~1|Site.x, method="REML", data=df,weights = varComb(varExp(form=~Lattitude)))
m11<-lme(log1p(Hyphal_length_sept)~Treatment.x*(MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2 + (Treatment.x*Vegetationtype*Time), random = ~1|Site.x, method="REML", data=df,weights = varComb(varPower(form=~elevation)))
m12<-lme(log1p(Hyphal_length_sept)~Treatment.x*(MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2 + (Treatment.x*Vegetationtype*Time), random = ~1|Site.x, method="REML", data=df,weights = varComb(varPower(form=~Time)))
anova(m0,m1,m2,m3,m5,m6,m7,m8,m9, m10, m11,m2)
#m2

par(mfrow=c(3,4))
qqnorm(scale(resid(m2, type="n")));abline(0,1,col=2)
hist(resid(m2, type="n"))
plot(resid(m2, type="n")~fitted(m2, type="n"));abline(h=0,col=2)
plot(resid(m2, type="n")~df$Treatment.x);abline(h=0,col=2)
plot(resid(m2, type="n")~df$Vegetationtype);abline(h=0,col=2)
plot(resid(m2, type="n")~df$MAP);abline(h=0,col=2)
plot(resid(m2, type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m2, type="n")~df$elevation);abline(h=0,col=2)
plot(resid(m2, type="n")~df$Lattitude);abline(h=0,col=2)
plot(resid(m2, type="n")~df$Time);abline(h=0,col=2)
shapiro.test(resid(m2, type="n"))

par(mfrow=c(3,4))
qqnorm(scale(resid(m10, type="n")));abline(0,1,col=2)
hist(resid(m10, type="n"))
plot(resid(m10, type="n")~fitted(m10, type="n"));abline(h=0,col=2)
plot(resid(m10, type="n")~df$Treatment.x);abline(h=0,col=2)
plot(resid(m10, type="n")~df$Vegetationtype);abline(h=0,col=2)
plot(resid(m10, type="n")~df$MAP);abline(h=0,col=2)
plot(resid(m10, type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m10, type="n")~df$elevation);abline(h=0,col=2)
plot(resid(m10, type="n")~df$Lattitude);abline(h=0,col=2)
plot(resid(m10, type="n")~df$Time);abline(h=0,col=2)
shapiro.test(resid(m10, type="n"))



m2.ml<-lme(log1p(Hyphal_length_sept)~Treatment.x*(MAT+MAP+Lattitude+elevation+Time) +  (MAT+MAP+Lattitude+elevation)^2 + (Treatment.x*Vegetationtype*Time), random = ~1|Site.x, method="ML", data=df,weights =  varComb(varExp(form=~Lattitude),varIdent(form=~1|Vegetationtype)))
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment.x:Time:Vegetationtype)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAP:elevation)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment.x:elevation)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment.x:MAP)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment.x:Lattitude)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Treatment.x:Vegetationtype)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAP:MAT)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Lattitude:MAP)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-Time:Vegetationtype)
anova(m2.ml)
m2.ml<-update(m2.ml,~.-MAP)
anova(m2.ml)

m2.reml<-update(m2.ml,method="REML")


par(mfrow=c(3,4))
qqnorm(scale(resid(m2.reml, type="n")));abline(0,1,col=2)
hist(resid(m2.reml, type="n"))
plot(resid(m2.reml, type="n")~fitted(m2.reml, type="n"));abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Treatment.x);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Vegetationtype);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$MAP);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$MAT);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$elevation);abline(h=0,col=2)
plot(resid(m2.reml, type="n")~df$Lattitude);abline(h=0,col=2)
shapiro.test(resid(m2.reml, type="n"))
anova(m2.reml)

summary(m2.reml)
library(effects)
plot(allEffects(m2.reml), multiline=T)



