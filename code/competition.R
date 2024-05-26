## read in data
setwd("C:/Users/Michael/Desktop/project/analysis/data")
comp <-read.csv("fights.csv")

## packages
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

##check data and convert the relavent intergers into factors
str(comp)
comp$population <- factor(comp$population)
comp$day <- factor(comp$day)
str(comp)
summary(comp)
comp$lgtime <- log(comp$time)
comp <- na.omit(comp)



## some preliminary plots to eyeball data
par(mfrow=c(1,1))
hist(comp$time)
#log time to make in more normal
hist(comp$lgtime)

### calculating means for time
# for treatment
Tmeans <- tapply(comp$time, comp$treatment, mean)
Tmeans
TSD <- tapply(comp$time, comp$treatment, sd)
TSD
SASE <- 434.7038/sqrt(125)
SASE
CSASE <- 516.7232/sqrt(121)
CSASE
# for phenotypes
Pmeans <- tapply(comp$time, comp$success.male, mean)
Pmeans
PSD <- tapply(comp$time, comp$success.male, sd)
PSD
EBSE <- 464.1863/sqrt(106)
EBSE
WTSE <- 487.4724/sqrt(140)
WTSE
# for populations
Popmeans <- tapply(comp$time, comp$population, mean)
Popmeans
PopSD <- tapply(comp$time, comp$population, sd)
PopSD
SE1 <- 488.1812/sqrt(25)
SE1
SE2 <- 211.4838/sqrt(25)
SE2
SE3 <- 505.4976/sqrt(25)
SE3
SE4 <- 509.8092/sqrt(25)
SE4
SE5 <- 286.6308/sqrt(25)
SE5
SE6 <- 313.1094/sqrt(24)
SE6
SE7 <- 701.9115/sqrt(24)
SE7
SE8 <- 523.776/sqrt(25)
SE8
SE9 <- 424.1952/sqrt(25)
SE9
SE10 <- 531.1302/sqrt(23)
SE10
# treatment x phenotype
comp$TP <- with(comp, treatment:success.male)
TPmeans <- tapply(comp$time, comp$TP, mean)
TPmeans
TPSD <- tapply(comp$time, comp$TP, sd)
TPSD
CSAEB <- 516.8160/sqrt(57)
CSAEB
CSAWT <- 517.2019/sqrt(64)
CSAWT
SAEB <- 396.8203/sqrt(49)
SAEB
SAWT <- 459.0588/sqrt(76)
SAWT
# for day
Dmeans <- tapply(comp$time, comp$day, mean)
Dmeans
DSD <- tapply(comp$time, comp$day, sd)
DSD
DSE1 <- 304.8671/sqrt(30)
DSE1
DSE2 <- 612.218/sqrt(56)
DSE2
DSE3 <- 408.6051/sqrt(40)
DSE3
DSE4 <- 410.2576/sqrt(40)
DSE4
DSE5 <- 442.9746/sqrt(40)
DSE5
DSE6 <- 498.801/sqrt(40)
DSE6

### but is log time better?
### calculating means for time
# for treatment
Tmeans <- tapply(comp$lgtime, comp$treatment, mean)
Tmeans
TSD <- tapply(comp$lgtime, comp$treatment, sd)
TSD
SASE <- 0.5403448/sqrt(125)
SASE
CSASE <- 0.5412490/sqrt(121)
CSASE

# for phenotypes
Pmeans <- tapply(comp$lgtime, comp$success.male, mean)
Pmeans
PSD <- tapply(comp$lgtime, comp$success.male, sd)
PSD
EBSE <- 0.5370911/sqrt(106)
EBSE
WTSE <- 0.5444524/sqrt(140)
WTSE

# for populations
Popmeans <- tapply(comp$lgtime, comp$population, mean)
Popmeans
PopSD <- tapply(comp$lgtime, comp$population, sd)
PopSD
SE1 <- 0.5044386/sqrt(25)
SE1
SE2 <- 0.5168098/sqrt(25)
SE2
SE3 <- 0.4734863/sqrt(25)
SE3
SE4 <- 0.6452218/sqrt(25)
SE4
SE5 <- 0.4392317/sqrt(25)
SE5
SE6 <- 0.453466/sqrt(24)
SE6
SE7 <- 0.5884982/sqrt(24)
SE7
SE8 <- 0.5493233/sqrt(25)
SE8
SE9 <- 0.4564127/sqrt(25)
SE9
SE10 <- 0.6293858/sqrt(23)
SE10

# treatment x phenotype
comp$TP <- with(comp, treatment:success.male)
TPmeans <- tapply(comp$lgtime, comp$TP, mean)
TPmeans
TPSD <- tapply(comp$lgtime, comp$TP, sd)
TPSD
CSAEB <- 0.5600772/sqrt(57)
CSAEB
CSAWT <- 0.5209513/sqrt(64)
CSAWT
SAEB <- 0.5137642/sqrt(49)
SAEB
SAWT <- 0.5596357/sqrt(76)
SAWT

# for day
Dmeans <- tapply(comp$lgtime, comp$day, mean)
Dmeans
DSD <- tapply(comp$lgtime, comp$day, sd)
DSD
DSE1 <- 0.4791716/sqrt(30)
DSE1
DSE2 <- 0.6224479/sqrt(56)
DSE2
DSE3 <- 0.5557697/sqrt(40)
DSE3
DSE4 <- 0.5539011/sqrt(40)
DSE4
DSE5 <- 0.4126920/sqrt(40)
DSE5
DSE6 <- 0.4987733/sqrt(40)
DSE6



#########################################################
####### linear mixed model by time and male phenotype ############
########################################################


LMM1 <- lmer(lgtime ~ treatment*success.male + (1 | population) + 
               (1 | day), data=comp)

res <- resid(LMM1)
ran <- ranef(LMM1)
fit <- fitted(LMM1)
par(mfrow=c(1,4))
hist(ran$population[,1], main='')
hist(ran$day[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

## seems fine...

summary(LMM1)
anova(LMM1)
wald.test(b=fixef(LMM1), Sigma=vcov(LMM1), Terms = 2, df=8)
wald.test(b=fixef(LMM1), Sigma=vcov(LMM1), Terms = 3, df=244)
wald.test(b=fixef(LMM1), Sigma=vcov(LMM1), Terms = 4, df=8)

### so there is no effect of the treatment, success.male or an interaction
### of the two on the time to compulation


## lsmeans
# for treatment
lstreat <- lsmeans(LMM1, "treatment")
lstreat

# for phenotype
lspheno <- lsmeans(LMM1, "success.male")
lspheno

# for interaction
lsinter <- lsmeans(LMM1, list(pairwise ~ treatment|success.male))
lsinter

######## nice lsmeans plots despite nothing is happening ###########

## treatment
LSMTT=data.frame(treatment=c("control","SA"), mean=c(6.592805, 6.506757), lower=c(6.506324, 6.420513), upper=c(6.679286, 6.593001))
ggplot() + 
  geom_pointrange(data=LSMTT, mapping=aes(x=treatment, y=mean, ymin=upper, ymax=lower), width=0.2, size=1, color="blue", fill="white", shape=22) + 
  theme(axis.title.y = element_text(size=18), 
        axis.title.x = element_text(size=18), 
        axis.text.x=element_text(size=15, colour="black"), 
        axis.text.y=element_text(size=15)) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_y_continuous(name="log time till mating occured (s)") + 
  scale_x_discrete(name="Treatment", limits=c("SA", "control"), labels=c("Selected", "Control")) + 
  coord_cartesian(ylim=c(6.4, 6.8))

## phenotype
LSMTP=data.frame(phenotype=c("Ebony","Wildtype"), mean=c(6.508641, 6.59092), lower=c(6.42962, 6.51637), upper=c(6.587662, 6.665443))
ggplot() + 
  geom_pointrange(data=LSMTP, mapping=aes(x=phenotype, y=mean, ymin=upper, ymax=lower), width=0.2, size=1, color="blue", fill="white", shape=22) + 
  theme(axis.title.y = element_text(size=18), 
        axis.title.x = element_text(size=18), 
        axis.text.x=element_text(size=15, colour="black"), 
        axis.text.y=element_text(size=15)) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_y_continuous(name="log time till mating occured (s)") + 
  scale_x_discrete(name="Successful Male Phenotype") + 
  coord_cartesian(ylim=c(6.4, 6.8))

## interaction
LSMINT=data.frame(combo=c("Selected-Ebony","Selected-Wildtype", "Control-Ebony", "Control-Wildtype"), mean=c(6.499671, 6.513843, 6.517612, 6.667997), lower=c(6.395538, 6.419993, 6.417184, 6.570666), upper=c(6.603804, 6.607693, 6.61804, 6.765328))
ggplot() + 
  geom_pointrange(data=LSMINT, mapping=aes(x=combo, y=mean, ymin=upper, ymax=lower), width=0.2, size=1, color="blue", fill="white", shape=22) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_y_continuous(name="log time till mating occured (s)") + 
  scale_x_discrete(name="Phenotype/treatment interaction", limits=c("Selected-Ebony", "Selected-Wildtype", "Control-Ebony", "Control-Wildtype")) + 
  coord_cartesian(ylim=c(6.39, 6.8))


########################################################
####### plots of log time to mating####################

dev.off()
## first simple plots of factors versus the lgtime
# 1: lgtime ~ treatment
bp2.1 <- ggplot(aes(y=lgtime, x=treatment), data=comp) + 
  geom_boxplot(fill=c("orchid2", "palegreen2"), notch=TRUE) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  stat_summary(fun.y=mean, geom="point", shape=1, size=4) + 
  theme(axis.title.y = element_text(size=18), 
        axis.title.x = element_text(size=18), 
        axis.text.x=element_text(size=15, colour="black"), 
        axis.text.y=element_text(size=15)) + 
  scale_y_continuous(name="log time till mating occured (s)") +
  scale_x_discrete(name="Selection Treatment", limits=c("SA", "control"), labels=c("Selected", "Control"))
bp2.1

# 2: lgtime ~ success.male
bp2.2 <- ggplot(aes(y=lgtime, x=success.male), data=comp) + 
  geom_boxplot(fill=c("grey70", "goldenrod3"), notch=TRUE) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  stat_summary(fun.y=mean, geom="point", shape=1, size=4) + 
  theme(axis.title.y = element_text(size=18), 
        axis.title.x = element_text(size=18), 
        axis.text.x=element_text(size=15, colour="black"), 
        axis.text.y=element_text(size=15)) + 
  scale_y_continuous(name="log time till mating occured (s)") +
  scale_x_discrete(name="Phenotype of the Successful Male", labels=c("Ebony", "Wildtype"))
bp2.2

## not much variatipon from pop or day so leave out day (3.95% var) and just do pop?

## all pops grouped by treatment
bp2.5 <- ggplot(aes(y=lgtime, x=population, fill=treatment), data=comp) + 
  geom_boxplot() + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  stat_summary(fun.y=mean, geom="point", shape=1, size=4) + theme(axis.title.y = element_text(size=18), 
                                                                  axis.title.x = element_text(size=18), 
                                                                  axis.text.x=element_text(size=15, colour="black"), 
                                                                  axis.text.y=element_text(size=15)) + 
  scale_y_continuous(name="log time till mating occured (s)") +
  scale_x_discrete(name="Population") + 
  scale_fill_manual(values=c("palegreen2", "orchid2"))
bp2.5

## interaction of treatment and success.male

## intercation with fills for phenotype
bp2.7 <- ggplot(aes(y=lgtime, x=treatment, fill=success.male), data=comp) + 
  geom_boxplot(notch=TRUE)  + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black"))
bp2.7 + theme(axis.title.y = element_text(size=18), 
              axis.title.x = element_text(size=18), 
              axis.text.x=element_text(size=15, colour="black"), 
              axis.text.y=element_text(size=15)) + 
  scale_y_continuous(name="log time till mating occurred (s)") +
  scale_x_discrete(name="Selection treatment", limits=c("SA", "control"), labels=c("Selected", "Control")) + 
  scale_fill_manual(values=c("grey70", "goldenrod3"))

## everything together fill=interaction
comp$succes.maletreat <- interaction(comp$success.male, comp$treatment)
bp2.8 <- ggplot(aes(y=lgtime, x=population, fill=succes.maletreat), data=comp) + 
  geom_boxplot()  + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black"))
bp2.8 + theme(axis.title.y = element_text(size=18), 
              axis.title.x = element_text(size=18), 
              axis.text.x=element_text(size=15, colour="black"), 
              axis.text.y=element_text(size=15)) + 
  scale_y_continuous(name="log time till mating occurred (s)") +
  scale_x_discrete(name="Population")


# without the colours -> we know what pops are SA and CSA
bp2.9 <- ggplot(aes(y=lgtime, x=population, fill=success.male), data=comp) + 
  geom_boxplot()  + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black"))
bp2.9 + theme(axis.title.y = element_text(size=20), 
              axis.title.x = element_text(size=20), 
              axis.text.x=element_text(size=16, colour="black"), 
              axis.text.y=element_text(size=16, colour="black"), 
              legend.position="none") +
  scale_y_continuous(name="log time till mating occurred (s)") +
  scale_x_discrete(name="Population") + 
  scale_fill_manual(values=c("grey75", "darkgoldenrod3"))

## are these plots ok? which ones are best to use?
## any further things to add?


#############################################
#######  extra calcs and plots for differences ###########
# population x phenotype
comp$PP <- with(comp, population:success.male)

PPmeans <- tapply(comp$lgtime, comp$PP, mean)
PPmeans
PopSD <- tapply(comp$lgtime, comp$PP, sd)
PopSD

## now read in the proportional differences dara frame
timediff <- read.csv("time difference.csv")
str(timediff)
timediff$mean.time.EB <- exp(timediff$mean.time.EB)
timediff$mean.time.WT <- exp(timediff$mean.time.WT)
timediff$prop.difference.e.wt <- timediff$mean.time.EB/timediff$mean.time.WT

par(mfrow=c(1,1))
hist(timediff$prop.difference.e.wt)

test <- t.test(prop.difference.e.wt ~ treatment, data=timediff, var.equal=TRUE)
test

summarySE(timediff, measurevar=c("prop.difference.e.wt"), groupvars=c("treatment"))

####### relative difference between conrol and sa pops for time to mating ###############
d=data.frame(treatment=c("control","SA"), mean=c(0.891749, 0.9797992), lower=c(0.8137173, 0.9058728), upper=c(0.9697807, 1.053726))
ggplot() + 
  geom_pointrange(data=d, mapping=aes(x=treatment, y=mean, ymin=upper, ymax=lower), width=0.2, size=1, color="blue", fill="white", shape=22) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_y_continuous(name="relative time to mating (ebony/wildtype)") + 
  scale_x_discrete(name="Treatment", limits=c("SA", "control"), labels=c("Selected", "Control"))








####################################################################
#### now to the more interesting part of analysing the effect of treatment
#### on the likely success of each male phenotype
### first the success.male catogory need to be converted into a binary vector
### so that it may be used as a respnse variable
#####################################################################

############### male success phentype analysis ####################
comp$phenobinary <- ifelse(comp$success.male %in% c("E"), 0, 1)

## E=0 and WT=1

### do some plots

plot(phenobinary ~ treatment, data=comp)
plot(phenobinary ~ day, data=comp)
plot(phenobinary ~ time, data=comp)

## these are not a great deal of help -> come back to later

### get some means for reference
summarySE(comp, measurevar=c("phenobinary"), groupvars=c("treatment"))
summarySE(comp, measurevar=c("phenobinary"), groupvars=c("day"))
sumbars <- summarySE(comp, measurevar=c("phenobinary"), groupvars=c("population", "treatment"))
sumbars ## for later to make the charts

### some general chisq stuff to descrive the data
### (some has been worked out by hand before)


### bias in the matings at all?
## a chi squre on the ratios of eb:WT between treatments
success.ratio <- matrix(c(49, 76, 57, 63), nrow=2)
success.ratio
chisq.test(success.ratio)


####### random factors and binary data = GLMM #############

BGLMM1 <- glmer(phenobinary ~ treatment + lgtime + (1 | population) + (1 | day),
                data=comp, family=binomial(link=logit))
BGLMM2 <- glmer(phenobinary ~ treatment + (1 | population) + (1 | day),
                data=comp, family=binomial(link=logit))

anova(BGLMM1, BGLMM2)

## all other combinations of random factors don't appear to work
## supports removal of lgtime 

# plot
res <- resid(BGLMM2)
ran <- ranef(BGLMM2)
fit <- fitted(BGLMM2)
par(mfrow=c(1,4))
hist(ran$population[,1], main='')
hist(ran$day[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
# these don't show muuch but they show that the binary thing is working

summary(BGLMM2)
anova (BGLMM2)

wald.test(b=fixef(BGLMM2), Sigma=vcov(BGLMM2), Terms = 2, df=8)

## lsmeans possible?
lsbintreat <- lsmeans(BGLMM2, "treatment")
lsbintreat

## beter ls plots
bintreat=data.frame(treatment=c("control","SA"), mean=c(0.1204685, 0.4453519), lower=c(-0.07134404, 0.2515994), upper=c(0.3122774, 0.6391044))
bintreat$invmean <- inv.logit(bintreat$mean)
bintreat$invlow <- inv.logit(bintreat$lower)
bintreat$invhigh <- inv.logit(bintreat$upper)
bintreat

################# point range for lsmean prop male success #######################

ggplot() + 
  geom_pointrange(data=bintreat, mapping=aes(x=treatment, y=invmean, ymin=invlow, ymax=invhigh), width=0.2, size=1, color="blue", fill="white", shape=22) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_y_continuous(name="Proportion of matings by wildtype males") + 
  scale_x_discrete(name="Treatment", limits=c("SA", "control"), labels=c("Selected", "Control")) + 
  coord_cartesian(ylim=c(0.46, 0.68))



####### better bar charts from brian surggestions #####################
write.csv(sumbars, "sumbars.csv")

bars <- read.csv("sumbars-use.csv")
summary(bars)
bars$population <- factor(bars$population)
summary(bars)
str(bars)


plot2 <- ggplot(data=bars, aes(y=phenobinary, x=population, fill=treatment)) + 
  geom_bar(stat="identity", colour="black", position="dodge") + 
  geom_errorbar(aes(ymin=phenobinary-se, ymax=phenobinary+se),
                width=.5,                    
                position=position_dodge(.2)) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black"), 
        legend.position="none") +
  scale_y_continuous(name="Proportion of matings by Wildtype males") + 
  scale_x_discrete(name="Population", labels=c("1", "2", "3", "4", "5", " ", " ", " ", "6", "7", "8", "9", "10")) + 
  scale_fill_manual(values=c("palegreen2", "mediumorchid2")) + 
  coord_cartesian(ylim=c(0, 1))
plot2



## no difference between SA and CSA populations





############### EXTRA STUFF THAT MIGHT NOT BE NEEDED (probs not) #####################

#### bar plots for no. of male success
# no. of males mating with treatment 
ggplot(aes(x=treatment, fill=success.male), data=comp) + 
  geom_bar(stat="bin", position=position_dodge(), colour="black") + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  theme(axis.title.y = element_text(size=18), 
        axis.title.x = element_text(size=18), 
        axis.text.x=element_text(size=15, colour="black"), 
        axis.text.y=element_text(size=15)) + 
  scale_y_continuous(name="Number of mating successes") + 
  scale_x_discrete(name="Treatment", limits=c("SA", "control"), labels=c("Selected", "Control")) + 
  scale_fill_manual(values=c("grey70", "goldenrod3"))


#### no of males per block
ggplot(aes(x=day, fill=success.male), data=comp) + 
  geom_bar(stat="bin", position=position_dodge(), colour="black") + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  theme(axis.title.y = element_text(size=18), 
        axis.title.x = element_text(size=18), 
        axis.text.x=element_text(size=15, colour="black"), 
        axis.text.y=element_text(size=15)) + 
  scale_y_continuous(name="Number of mating successes") + 
  scale_fill_manual(values=c("grey70", "goldenrod3"))



