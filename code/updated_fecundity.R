######################################################################

getwd()

##packages
library(ggplot2)
library(nlme)
library(lme4)
library(multcomp)
library(lsmeans)

### read in the data
# calling this data set 'egg'
egg <- read.csv("../data/eggs.csv")
##check
summary(egg)
str(egg)
## 180 obs of 8 variables

## change population to a factor
egg$population <- factor(egg$population)
str(egg)
summary(egg)

## remove zeros -> no laying, duff flies
str(egg)
## log total? - no harm in trying
egg$lgtotal <- log(egg$total)


#### average number of eggs laid
# for treatment
Tmeans <- tapply(egg$total, egg$treatment, mean)
Tmeans
TSD <- tapply(egg$total, egg$treatment, sd)
TSD
SASE <- 18.12243/sqrt(90)
SASE
CSASE <- 14.87208/sqrt(90)
CSASE
# for phenotypes
Pmeans <- tapply(egg$total, egg$phenotype, mean)
Pmeans
PSD <- tapply(egg$total, egg$phenotype, sd)
PSD
EBSE <- 14.16403/sqrt(90)
EBSE
WTSE <- 18.60524/sqrt(90)
WTSE
# for populations
Popmeans <- tapply(egg$total, egg$population, mean)
Popmeans
PopSD <- tapply(egg$total, egg$population, sd)
PopSD
SE1 <- 16.77845/sqrt(18)
SE1
SE2 <- 19.17455/sqrt(18)
SE2
SE3 <- 17.3868/sqrt(18)
SE3
SE4 <- 11.03767/sqrt(18)
SE4
SE5 <- 16.44649/sqrt(18)
SE5
SE6 <- 11.86319/sqrt(18)
SE6
SE7 <- 16.37879/sqrt(18)
SE7
SE8 <- 13.09006/sqrt(18)
SE8
SE9 <- 14.64080/sqrt(18)
SE9
SE10 <- 18.13205/sqrt(18)
SE10
# treatment x phenotype
egg$TP <- with(egg, treatment:phenotype)
TPmeans <- tapply(egg$total, egg$TP, mean)
TPmeans
TPSD <- tapply(egg$total, egg$TP, sd)
TPSD
CSAEB <- 13.70405/sqrt(45)
CSAEB
CSAWT <- 15.93189/sqrt(45)
CSAWT
SAEB <- 14.61499/sqrt(45)
SAEB
SAWT <- 21.05328/sqrt(45)
SAWT

### models: poisson probs cause they are counts, use total or average?
## generalised linear model likely
## 2 ways thus can be done
## 1: averaging out the replication between the two days
## 2: including day as a random factor in the model and nesting the things

## trying model with total rather than average
#GLMM1 <- glmer(total ~ treatment*phenotype + (1 | population) +(1 | vial), 
#              data=egg,  
#              family = poisson(link=log))
## this midel includes a random effect for vial (good in other analysis but not here)

TGLMM <- glmer(total ~ treatment*phenotype + (1 | population), 
               data=egg,  
               family = poisson(link=log))

res <- resid(TGLMM)
ran <- ranef(TGLMM)
fit <- fitted(TGLMM)
par(mfrow=c(1,4))
hist(ran$population[,1], main='')
hist(ran$population[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

anova(TGLMM)
summary(TGLMM)

wald.test(b=fixef(TGLMM), Sigma=vcov(TGLMM), Terms = 2, df=8)
wald.test(b=fixef(TGLMM), Sigma=vcov(TGLMM), Terms = 3, df=178)
wald.test(b=fixef(TGLMM), Sigma=vcov(TGLMM), Terms = 4, df=8)



############# lsmeans ##########
# for treatment
lstreat <- lsmeans(TGLMM, "treatment")
lstreat
## undo logs for SA
exp(3.189359)
exp(0.1160505)
##undo logs for CSA
exp(3.127488)
exp(0.166064)

# for phenotype
lspheno <- lsmeans(TGLMM, "phenotype")
lspheno
## undo logs for EB
exp(3.2311)
exp(0.08325038)
##undo logs for WT
exp(3.084737)
exp(0.0836548)

#### ls plot for pheno
LSMFP=data.frame(phenotype=c("Ebony","Wildtype"), mean=c(25.30748, 21.86172), lower=c(24.22067, 20.77447), upper=c(26.39429, 22.94897))
ggplot() + 
  geom_pointrange(data=LSMFP, mapping=aes(x=phenotype, y=mean, ymin=upper, ymax=lower), width=0.2, size=1, color="blue", fill="white", shape=22) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_y_continuous(name="Total number of eggs laid") + 
  scale_x_discrete(name="Phenotype") + 
  coord_cartesian(ylim=c(20, 28))

# for interaction
lsinter <- lsmeans(TGLMM, list(pairwise ~ treatment|phenotype))
lsinter
## undo logs for SAEB
exp(3.263762)
exp(0.1176309)
##undo logs for CSAEB
exp(3.200459)
exp(0.1178333)
## undo logs for SAWT
exp(3.114956)
exp(0.1181773)
##undo logs for CSAWT
exp(3.054517)
exp(0.1184308)

#### ls plot for inter
## interaction
LSMFINT=data.frame(combo=c("Selected-Ebony","Selected-Wildtype", "Control-Ebony", "Control-Wildtype"), mean=c(26.14772, 22.53244, 24.54379, 21.21094), lower=c(25.02289, 21.407, 23.41873, 20.08521), upper=c(27.27255, 23.65788, 25.66885, 22.33667))
ggplot() + 
  geom_pointrange(data=LSMFINT, mapping=aes(x=combo, y=mean, ymin=upper, ymax=lower), width=0.2, size=1, color="blue", fill="white", shape=22) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_y_continuous(name="Total number of eggs laid") + 
  scale_x_discrete(name="Phenotype/treatment interaction", limits=c("Selected-Ebony", "Selected-Wildtype", "Control-Ebony", "Control-Wildtype")) + 
  coord_cartesian(ylim=c(20, 28))





############################################
######### PLOTS ######################### 
#############################

## box plots of all populations by treatments 
dev.off()

# simple box plots
# 1: total ~ treatment
ggplot(aes(y=total, x=treatment), data=egg) + 
  stat_boxplot(geom="errorbar") + 
  geom_boxplot(fill=c("orchid2", "palegreen2"), notch=TRUE) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  stat_summary(fun.y=mean, geom="point", shape=1, size=4) + 
  theme(axis.title.y = element_text(size=18), 
        axis.title.x = element_text(size=18), 
        axis.text.x=element_text(size=15, colour="black"), 
        axis.text.y=element_text(size=15)) + 
  scale_y_continuous(name="Total number of eggs laid") +
  scale_x_discrete(name="Selection Treatment", limits=c("selected", "control"), labels=c("Selected", "Control"))

# 2: total ~ phenotype
bp2 <-ggplot(aes(y=total, x=phenotype), data=egg) + 
  geom_boxplot(fill=c("grey70", "goldenrod3"), notch=TRUE) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  stat_summary(fun.y=mean, geom="point", shape=1, size=4)
bp2 + theme(axis.title.y = element_text(size=18), 
            axis.title.x = element_text(size=18), 
            axis.text.x=element_text(size=15, colour="black"), 
            axis.text.y=element_text(size=15)) + 
  scale_y_continuous(name="Total number of eggs laid") +
  scale_x_discrete(name="Female Phenotype", labels=c("Ebony", "Wildtype"))


# all pops showing grouped by treatment
bp3 <- ggplot(aes(y=total, x=population, fill=treatment), data=egg) + 
  geom_boxplot() + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  stat_summary(fun.y=mean, geom="point", shape=1, size=4)

bp3 + theme(axis.title.y = element_text(size=18), 
            axis.title.x = element_text(size=18), 
            axis.text.x=element_text(size=15, colour="black"), 
            axis.text.y=element_text(size=15)) + 
  scale_y_continuous(name="Total number of eggs laid") +
  scale_x_discrete(name="Population") + 
  scale_fill_manual(values=c("palegreen2", "orchid2"))

## box plots for interaction of phenotype and treatment

# using fills for phenotype
bp5 <- ggplot(aes(y=total, x=treatment, fill=phenotype), data=egg) + 
  geom_boxplot(notch=TRUE)  + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black"))
bp5 + theme(axis.title.y = element_text(size=18), 
            axis.title.x = element_text(size=18), 
            axis.text.x=element_text(size=15, colour="black"), 
            axis.text.y=element_text(size=15)) + 
  scale_y_continuous(name="Total number of eggs laid") +
  scale_x_discrete(name="Selection treatment", limits=c("selected", "control"), labels=c("Selected", "Control")) + 
  scale_fill_manual(values=c("grey70", "goldenrod3"))

# trying to combine both
egg$phenotreat <- interaction(egg$phenotype, egg$treatment)
bp6 <- ggplot(aes(y=total, x=population, fill=phenotreat), data=egg) + 
  geom_boxplot() + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black"))
bp6 + theme(axis.title.y = element_text(size=18), 
            axis.title.x = element_text(size=18), 
            axis.text.x=element_text(size=15, colour="black"), 
            axis.text.y=element_text(size=15)) + 
  scale_y_continuous(name="Total number of eggs laid") +
  scale_x_discrete(name="Population")

# this is a bit messy

# try with just population and phenotype
bp7 <- ggplot(aes(y=total, x=population, fill=phenotype), data=egg) + 
  geom_boxplot() + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black"))
bp7 + theme(axis.title.y = element_text(size=20), 
            axis.title.x = element_text(size=20), 
            axis.text.x=element_text(size=16, colour="black"), 
            axis.text.y=element_text(size=16, colour="black"), 
            legend.position="none") +
  scale_y_continuous(name="Total number of eggs laid") +
  scale_x_discrete(name="Population") + 
  scale_fill_manual(values=c("grey75", "darkgoldenrod3"))

## are these plots ok? which ones are best to use?
## any further things to add?


###  extra calcs and plots for differences
# population x phenotype
egg$PP <- with(egg, population:phenotype)

PPmeans <- tapply(egg$total, egg$PP, mean)
PPmeans
PopSD <- tapply(egg$total, egg$PP, sd)
PopSD

## now read in the proportional differences dara frame
totaldiff <- read.csv("fecundity difference.csv")
str(totaldiff)

par(mfrow=c(1,1))
hist(totaldiff$prop.difference.e.wt)

test <- t.test(prop.difference.e.wt ~ treatment, data=totaldiff, var.equal=TRUE)
test

TDmean <- tapply(totaldiff$prop.difference.e.wt, totaldiff$treatment, mean)
TDmean
TDSD <- tapply(totaldiff$prop.difference.e.wt, totaldiff$treatment, sd)
TDSD
SESA <- 1.2941518/sqrt(5)
SESA
SECSA <-0.4095912/sqrt(5)
SECSA

### now tomplot this - point range?


d=data.frame(Treatment=c("control","SA"), mean=c(1.23712, 1.7511), lower=c(1.053945, 1.172338), upper=c(1.420295, 2.329862))
ggplot() + 
  geom_pointrange(data=d, mapping=aes(x=Treatment, y=mean, ymin=upper, ymax=lower), width=0.2, size=1, color="blue", fill="white", shape=22) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_y_continuous(name="relative total eggs laid (ebony/wildtype)") + 
  scale_x_discrete(limits=c("SA", "control"), labels=c("Selected", "Control"))

