#### Sex and genotype experiments
setwd("C:/Users/Michael/Desktop/project/analysis/data")
SP <- read.csv("sexpheno compilied.csv")

# packages
library(ggplot2)
library(lme4)
library(multcomp)
library(aod)
library(gtools)
library(lsmeans)
library(Rmisc)
library(chemCal)
library(plyr)
library(dplyr)

## check and fix to factors
str(SP)
SP$pop <- factor(SP$pop)
SP <- na.omit(SP)
str(SP)
summary(SP)



########################################################################################
#########################################################################################
###### calculate a proportion of ebony females out of total ebony ############

SP$propEBF <- SP$comptotaleeF/SP$comptotalEB 

## GLMM for the affect of treatment on the proportion of Ebony females

GLMM1 <- glmer(propEBF ~ Treatment + (1 | pop) + (1 | Block), data=SP,
               weights=comptotalEB, family=binomial(link=logit))

## but 1 is not different by mucg and makes more sense 
res <- resid(GLMM1)
ran <- ranef(GLMM1)
fit <- fitted(GLMM1)
par(mfrow=c(1,4))
hist(ran$pop[,1], main='')
hist(ran$Block[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(GLMM1)
anova(GLMM1)

wald.test(b=fixef(GLMM1), Sigma=vcov(GLMM1), Terms = 2, df=8)

## no affect of treatment here

# ls means
lsEBF <- lsmeans(GLMM1, "Treatment")
lsEBF

lsEBFtreat=data.frame(treatment=c("control","SA"), mean=c(0.04072487, 0.07968419), lower=c(-0.00151913, 0.03745035), upper=c(0.07993661, 0.121918))
lsEBFtreat$invmean <- inv.logit(lsEBFtreat$mean)
lsEBFtreat$invlow <- inv.logit(lsEBFtreat$lower)
lsEBFtreat$invhigh <- inv.logit(lsEBFtreat$upper)
lsEBFtreat

###### ls mean plot propEBF ~ treatment #####################
ggplot() + 
  geom_pointrange(data=lsEBFtreat, mapping=aes(x=treatment, y=invmean, ymin=invlow, ymax=invhigh), width=0.2, size=1, color="blue", fill="white", shape=22) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black")) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_y_continuous(name="Proportion of ebony females") + 
  scale_x_discrete(name="Treatment", limits=c("SA", "control"), labels=c("Selected", "Control")) + 
  coord_cartesian(ylim=c(0.48, 0.54))


###### means and se ############
Tmeans <- tapply(SP$propEBF, SP$Treatment, mean)
Tmeans
TSD <- tapply(SP$propEBF, SP$Treatment, sd)
TSD
SASE <- 0.1229896/sqrt(118)
SASE
CSASE <- 0.1240056/sqrt(110)
CSASE
# for population
Pmeans <- tapply(SP$propEBF, SP$pop, mean)
Pmeans
PSD <- tapply(SP$propEBF, SP$pop, sd)
PSD
SE1 <- 0.1297759/sqrt(24)
SE1
SE2 <- 0.1214975/sqrt(22)
SE2
SE3 <- 0.1145658/sqrt(24)
SE3
SE4 <- 0.1310569/sqrt(24)
SE4
SE5 <- 0.1162799/sqrt(24)
SE5
SE6 <- 0.1151343/sqrt(20)
SE6
SE7 <- 0.1435755/sqrt(18)
SE7
SE8 <- 0.1280776/sqrt(24)
SE8
SE9 <- 0.1165792/sqrt(24)
SE9
SE10 <- 0.1153703/sqrt(24)
SE10
## for block
Bmeans <- tapply(SP$propEBF, SP$Block, mean)
Bmeans
BSD <- tapply(SP$propEBF, SP$Block, sd)
BSD
B1SE <- 0.1109837/sqrt(112)
B1SE
B2SE <- 0.1195541/sqrt(116)
B2SE

### plots
# 1: prop EBF  ~ treatment
WEB.treat <- ggplot(aes(y=propEBF, x=Treatment), data=SP) + 
  geom_boxplot(fill=c("orchid2", "palegreen2"), notch = TRUE) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  stat_summary(fun.y=mean, geom="point", shape=1, size=4) + 
  theme(axis.title.y = element_text(size=18), 
        axis.title.x = element_text(size=18), 
        axis.text.x=element_text(size=15, colour="black"), 
        axis.text.y=element_text(size=15)) + 
  scale_y_continuous(name="proportion of ebony females/total ebony") +
  scale_x_discrete(name="Selection Treatment", limits=c("SA", "Control"), labels=c("Selected", "Control")) + 
  coord_cartesian(ylim=c(0, 1))
WEB.treat

# with pop and fill=treatment
WEB.pop <- ggplot(aes(y=propEBF, x=pop, fill=Treatment), data=SP) + 
  geom_boxplot() + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black"), 
        legend.position="none") +
  scale_y_continuous(name="Proportion of ebony females") +
  scale_x_discrete(name="Population") + 
  scale_fill_manual(values=c("palegreen2", "mediumorchid2")) + 
  coord_cartesian(ylim=c(0, 1))
WEB.pop





###################################################################
##### now to look at the prop of WTF ########################### 



SP$propWTF <- SP$comptotalWTF/SP$comptotalWT 

## GLMM for the affect of treatment on the proportion of Ebony females

GLMM2 <- glmer(propWTF ~ Treatment + (1 | pop) + (1 | Block), data=SP,
               weights=comptotalWT, family=binomial(link=logit))

## but 1 is not different by mucg and makes more sense 
res <- resid(GLMM2)
ran <- ranef(GLMM2)
fit <- fitted(GLMM2)
par(mfrow=c(1,4))
hist(ran$pop[,1], main='')
hist(ran$Block[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(GLMM2)
anova(GLMM2)

wald.test(b=fixef(GLMM2), Sigma=vcov(GLMM2), Terms = 2, df=8)

lsWTF <- lsmeans(GLMM2, "Treatment")
lsWTF

lsWTFtreat=data.frame(treatment=c("control","SA"), mean=c(0.09647805, -0.02824233), lower=c(0.05182007, -0.07130551), upper=c(0.141136, 0.01482085))
lsWTFtreat$invmean <- inv.logit(lsWTFtreat$mean)
lsWTFtreat$invlow <- inv.logit(lsWTFtreat$lower)
lsWTFtreat$invhigh <- inv.logit(lsWTFtreat$upper)
lsWTFtreat
## nice lsmeans plot
ggplot() + 
  geom_pointrange(data=lsWTFtreat, mapping=aes(x=treatment, y=invmean, ymin=invlow, ymax=invhigh), width=0.2, size=1, color="blue", fill="white", shape=22) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_y_continuous(name="Proportion of Wildtype females") + 
  scale_x_discrete(name="Treatment", limits=c("SA", "control"), labels=c("Selected", "Control")) + 
  coord_cartesian(ylim=c(0.48, 0.54))


### means and se
# for treatment
Tmeans <- tapply(SP$propWTF, SP$Treatment, mean)
Tmeans
TSD <- tapply(SP$propWTF, SP$Treatment, sd)
TSD
SASE <- 0.09940227/sqrt(118)
SASE
CSASE <- 0.12677454/sqrt(110)
CSASE
# for population
Pmeans <- tapply(SP$propWTF, SP$pop, mean)
Pmeans
PSD <- tapply(SP$propWTF, SP$pop, sd)
PSD
SE1 <- 0.14745615/sqrt(24)
SE1
SE2 <- 0.10488108/sqrt(22)
SE2
SE3 <- 0.09904769/sqrt(24)
SE3
SE4 <- 0.09901405/sqrt(24)
SE4
SE5 <- 0.16172889/sqrt(24)
SE5
SE6 <- 0.09007699/sqrt(20)
SE6
SE7 <- 0.07619439/sqrt(18)
SE7
SE8 <- 0.12426938/sqrt(24)
SE8
SE9 <- 0.09820765/sqrt(24)
SE9
SE10 <- 0.09426901/sqrt(24)
SE10
## for block
Bmeans <- tapply(SP$propWTF, SP$Block, mean)
Bmeans
BSD <- tapply(SP$propWTF, SP$Block, sd)
BSD
B1SE <- 0.1298255/sqrt(112)
B1SE
B2SE <- 0.1166414/sqrt(116)
B2SE

## plots

# 1: propwtf  ~ treatment
WTF.treat <- ggplot(aes(y=propWTF, x=Treatment), data=SP) + 
  geom_boxplot(fill=c("orchid2", "palegreen2"), notch = TRUE) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  stat_summary(fun.y=mean, geom="point", shape=1, size=4) + 
  theme(axis.title.y = element_text(size=18), 
        axis.title.x = element_text(size=18), 
        axis.text.x=element_text(size=15, colour="black"), 
        axis.text.y=element_text(size=15)) + 
  scale_y_continuous(name="proportion of WT females (over total WT)") +
  scale_x_discrete(name="Selection Treatment", limits=c("SA", "Control"), labels=c("Selected", "Control")) + 
  coord_cartesian(ylim=c(0, 1))
WTF.treat


# with pop and fill=treatment
WTF.pop <- ggplot(aes(y=propWTF, x=pop, fill=Treatment), data=SP) + 
  geom_boxplot() + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black"), 
        legend.position="none") + 
  scale_y_continuous(name="Proportion of wildtype females") +
  scale_x_discrete(name="Population") + 
  scale_fill_manual(values=c("palegreen2", "mediumorchid2")) + 
  coord_cartesian(ylim=c(0, 1))
WTF.pop
