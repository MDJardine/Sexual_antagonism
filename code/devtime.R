#### Sex and genotype experiments
setwd("C:/Users/Michael/Documents/masters/project/analysis/data")
sexgeno1 <- read.csv("sexgenotype.csv")

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
sexgeno$pop <- factor(sexgeno$pop)
sexgeno <- na.omit(sexgeno)
str(sexgeno)
summary(sexgeno)

## some histograms of the spread of total data per vial
par(mfrow=c(1,1))
hist(sexgeno$total)
par(mfrow=c(2,2))
hist(sexgeno$totF)
hist(sexgeno$totM)
hist(sexgeno$totee)
hist(sexgeno$totWT)

hist(sexgeno$totF)
hist(sexgeno$totM)
hist(sexgeno$totee)
hist(sexgeno$totWT)


plot(total~day, data=sexgeno)

plot(totF~day, data=sexgeno1, col=Treatment, ylab="Total female flies eclosed", xlab="Days since egg laying", pch=16)
plot(totM~day, data=sexgeno1, col=Treatment,  ylab="Total male flies eclosed", xlab="Days since egg laying", pch=16)
plot(totee~day, data=sexgeno1,  col=Treatment, ylab="Total 'ebony' flies eclosed", xlab="Days since egg laying", pch=16)
plot(totWT~day, data=sexgeno1, col=Treatment, ylab="Total Wildtype flies eclosed", xlab="Days since egg laying", pch=16)



totaperday <- ggplot(sexgeno, aes(x=day, y=total)) + 
  geom_bar(position=position_dodge(), stat="identity") + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black"), 
        legend.position="none") + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_y_continuous(name="Total flies ecolosed per day") + 
  scale_x_continuous(name="Day since laying") + 


## this leads to lots of zeros as most days there were very few flies
## could we delete some days of data and therefore have a cleaner data set?
## which days to delete -> look at plots of total vs. day to see when the drop off is


## seems that anything over day 17 has too few flies, remove these days

sexgeno = sexgeno[which(sexgeno$day < 17),]
sexgeno = sexgeno[which(sexgeno$day > 9),]
str(sexgeno)
summary(sexgeno)

## seems to have worked -> plot again

# histograms
par(mfrow=c(1,1))
hist(sexgeno$total)
par(mfrow=c(2,2))
hist(sexgeno$totF)
hist(sexgeno$totM)
hist(sexgeno$totee)
hist(sexgeno$totWT)
## still lots of zeros and not normal but an awful lot better than before

### from the current data set to get totals of sex and genotype to use in chisq tests

## sum of total sex and phenotype
?colSums
sum(sexgeno$totF, na.rm=FALSE)
sum(sexgeno$totM, na.rm=FALSE)
sum(sexgeno$totee, na.rm=FALSE)
sum(sexgeno$totWT, na.rm=FALSE)

# sum of each sex and phenotype
sum(sexgeno$eeF, na.rm=FALSE)
sum(sexgeno$eeM, na.rm=FALSE)
sum(sexgeno$WTF, na.rm=FALSE)
sum(sexgeno$WTM, na.rm=FALSE)

## now the SA for quick chisq tests
sexgenoSA <- sexgeno[which(sexgeno$Treatment=='SA'),]
sum(sexgenoSA$totF, na.rm=FALSE)
sum(sexgenoSA$totM, na.rm=FALSE)
sum(sexgenoSA$totee, na.rm=FALSE)
sum(sexgenoSA$totWT, na.rm=FALSE)

# sex and genotype
sum(sexgenoSA$eeF, na.rm=FALSE)
sum(sexgenoSA$eeM, na.rm=FALSE)
sum(sexgenoSA$WTF, na.rm=FALSE)
sum(sexgenoSA$WTM, na.rm=FALSE)

## and for the c-sa
sexgenoCSA <- sexgeno[which(sexgeno$Treatment=='Control'),]
sum(sexgenoCSA$totF, na.rm=FALSE)
sum(sexgenoCSA$totM, na.rm=FALSE)
sum(sexgenoCSA$totee, na.rm=FALSE)
sum(sexgenoCSA$totWT, na.rm=FALSE)

# SEX AND GENOTYPE
sum(sexgenoCSA$eeF, na.rm=FALSE)
sum(sexgenoCSA$eeM, na.rm=FALSE)
sum(sexgenoCSA$WTF, na.rm=FALSE)
sum(sexgenoCSA$WTM, na.rm=FALSE)

### comparing the two treatments

# sex ratio
sex <- matrix(c(2650, 2754, 2603, 2573), nrow=2)
sex
chisq.test(sex)

# phenotype ratio
pheno <- matrix(c(2235, 2603, 3018, 2724), nrow=2)
pheno
chisq.test(pheno)

# sex amd pheno
SAP <- matrix(c(1162, 1328, 1073, 1275, 1488, 1426, 1530, 1298), nrow=2)
SAP
chisq.test(SAP)

### there appears to be no difference in the sex ratio between the two treatments
### difference in pheno -> more WT in SA?
### difference in sex/pheno combos between SA and cSA




######################################################################
####### DEVELOPMENT TIME AND CUMULATION PLOTS ##################


### calculate cumulative eclosion per day per vial
propSP <- ddply(sexgeno, "replicate", plyr::mutate, 
                eeF_prop = cumsum(eeF)/sum(eeF),
                eeM_prop = cumsum(eeM)/sum(eeM),
                WTF_prop = cumsum(WTF)/sum(WTF),
                WTM_prop = cumsum(WTM)/sum(WTM),
                totF_prop = cumsum(totF)/sum(totF),
                totM_prop = cumsum(totM)/sum(totM),
                totee_prop = cumsum(totee)/sum(totee),
                totWT_prop = cumsum(totWT)/sum(totWT),
                total_prop = cumsum(total)/sum(total))

propSP <- na.omit(propSP)

### claculate the day at 50% eclosed in new data frame
models <- propSP %>% group_by(replicate) %>% do(Block = .$Block[1], Treatment = .$Treatment[1], pop = .$pop[1], 
                                                eeF_mod = inverse.predict(lm(eeF_prop~day, data = .), 0.5)$Prediction,
                                                eeM_mod = inverse.predict(lm(eeM_prop~day, data = .), 0.5)$Prediction,
                                                WTF_mod = inverse.predict(lm(WTF_prop~day, data = .), 0.5)$Prediction,
                                                WTM_mod = inverse.predict(lm(WTM_prop~day, data = .), 0.5)$Prediction,
                                                totF_mod = inverse.predict(lm(totF_prop~day, data = .), 0.5)$Prediction,
                                                totM_mod = inverse.predict(lm(totM_prop~day, data = .), 0.5)$Prediction,
                                                totee_mod = inverse.predict(lm(totee_prop~day, data = .), 0.5)$Prediction,
                                                totWT_mod = inverse.predict(lm(totWT_prop~day, data = .), 0.5)$Prediction,
                                                total = inverse.predict(lm(eeM_prop~day, data = .), 0.5)$Prediction)

glimpse(models)

### fixing this data frame
mod <- data.frame(matrix(unlist(models), nrow =226))
str(mod)
names(mod)[names(mod)=='X4'] <- 'pop'
names(mod)[names(mod)=='X3'] <- 'treat'
names(mod)[names(mod)=='X2'] <- 'block'
names(mod)[names(mod)=='X5'] <- 'eeF_mod'
names(mod)[names(mod)=='X6'] <- 'eeM_mod'
names(mod)[names(mod)=='X7'] <- 'WTF_mod'
names(mod)[names(mod)=='X8'] <- 'WTM_mod'
names(mod)[names(mod)=='X9'] <- 'totF_mod'
names(mod)[names(mod)=='X10'] <- 'totM_mod'
names(mod)[names(mod)=='X11'] <- 'totee_mod'
names(mod)[names(mod)=='X12'] <- 'totWT_mod'
names(mod)[names(mod)=='X13'] <- 'total_mod'

##  export to fix finally in excel
write.csv(mod, "mod.csv")
## back from excel 

######## START HERE ############################ 
########################################################
############# DEVELOPMENT TIME ANALYSIS AND AVERGAE DAY PLOTS #############
dev <- read.csv("fixmod.csv")
dev$treat <- factor(dev$treat)
dev$block <- factor(dev$block)
dev$pop <- factor(dev$pop)
str(dev)


#### CALCULATE AVERGAES FOR BAR PLOT OF ECLOSION DAY

######## EBONY###############
## means and se for ebony 
eeFbar <- summarySE(dev, measurevar=c("eeF_mod"), groupvars=c("pop", "treat"))
eeFbar
eeMbar <- summarySE(dev, measurevar=c("eeM_mod"), groupvars=c("pop", "treat"))
eeMbar
str(eeFbar)
eeFbar$eeM_mod <- eeMbar$eeM_mod
eeFbar$mSE <- eeMbar$se
eeFbar$mSD <- eeMbar$sd
eeFbar$mCI <- eeMbar$ci
str(eeFbar)

#### created new df now export to fix
write.csv(eeFbar, "eebars.csv")

## get it back and alter the 
eedaybars <- read.csv("eebars-use.csv")
str(eedaybars)
eedaybars$pop <- factor(eedaybars$pop)
str(eedaybars)

### now plot mean day eclosed for pop and sex
eedayplot <- ggplot(eedaybars, aes(x=pop, y=day, fill=sex)) + 
  geom_bar(position=position_dodge(), stat="identity", colour="black") +
  geom_errorbar(aes(ymin=day-se, ymax=day+se),
                width=.5,                    
                position=position_dodge(.9)) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black"), 
        legend.position="none") + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_y_continuous(name="Mean predicted ebony development time (days)") + 
  coord_cartesian(ylim=c(9,13)) + 
  scale_x_discrete(name="Population", labels=c("1", "2", "3", "4", "5", " ", " ", " ", "6", "7", "8", "9", "10")) + 
  scale_fill_manual(values=c("red1", "cornflowerblue"))
eedayplot

##### NOW THE EBONY DIFFERENCE plot #########
eeFbar$reltime <- eeFbar$eeF_mod/eeFbar$eeM_mod
str(eeFbar)
eereldif <- summarySE(eeFbar, measurevar=c("reltime"), groupvars=c("treat"))
eereldif

eetimedif=data.frame(treatment=c("control","SA"), mean=c(0.9760862, 0.9819886), lower=c(0.96303375, 0.9733188), upper=c(0.9891349, 0.9906584))
ggplot() + 
  geom_pointrange(data=eetimedif, mapping=aes(x=treatment, y=mean, ymin=upper, ymax=lower), width=0.2, size=1, color="blue", fill="white", shape=22) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black")) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_y_continuous(name="Mean predicted ebi development time (days)") + 
  scale_x_discrete(name="Treatment", limits=c("SA", "control"), labels=c("Selected", "Control")) + 
  coord_cartesian(ylim=c(0.94, 0.995))


######## WILDTYPE ###############
## means and se for wildtype 
WTFbar <- summarySE(dev, measurevar=c("WTF_mod"), groupvars=c("pop", "treat"))
WTFbar
WTMbar <- summarySE(dev, measurevar=c("WTM_mod"), groupvars=c("pop", "treat"))
WTMbar
str(WTFbar)
WTFbar$WTM_mod <- WTMbar$WTM_mod
WTFbar$mSE <- WTMbar$se
WTFbar$mSD <- WTMbar$sd
WTFbar$mCI <- WTMbar$ci
str(eeFbar)

#### created new df now export to fix
write.csv(WTFbar, "WTbars.csv")

## get it back and alter the 
WTdaybars <- read.csv("WTbars-use.csv")
str(WTdaybars)
WTdaybars$pop <- factor(WTdaybars$pop)
str(WTdaybars)

### now plot mean day eclosed for pop and sex
WTdayplot <- ggplot(WTdaybars, aes(x=pop, y=day, fill=sex)) + 
  geom_bar(position=position_dodge(), stat="identity", colour="black") +
  geom_errorbar(aes(ymin=day-se, ymax=day+se),
                width=.5,                    
                position=position_dodge(.9)) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black"), 
        legend.position="none") +  
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_y_continuous(name="Mean predicted Wildtype development time (days)") + 
  coord_cartesian(ylim=c(9,13)) + 
  scale_x_discrete(name="Population", labels=c("1", "2", "3", "4", "5", " ", " ", " ", "6", "7", "8", "9", "10")) + 
  scale_fill_manual(values=c("red1", "cornflowerblue"))
WTdayplot



##### NOW THE WT DIFFERENCE plot #########

WTFbar$reltime <- WTFbar$WTF_mod/WTFbar$WTM_mod
str(WTFbar)
WTreldif <- summarySE(WTFbar, measurevar=c("reltime"), groupvars=c("treat"))
WTreldif


WTtimedif=data.frame(treatment=c("control","SA"), mean=c(0.9538802, 0.9611793), lower=c(0.9427964, 0.9495046), upper=c(0.964964, 0.972854))
ggplot() + 
  geom_pointrange(data=WTtimedif, mapping=aes(x=treatment, y=mean, ymin=upper, ymax=lower), width=0.2, size=1, color="blue", fill="white", shape=22) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black")) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        panel.background=element_blank(), axis.line=element_line(colour="black")) + 
  scale_y_continuous(name="Wild-type female/male relative development time") + 
  scale_x_discrete(name="Treatment", limits=c("SA", "control"), labels=c("Selected", "Control")) + 
  coord_cartesian(ylim=c(0.94, 0.995))



################# MODELS FOR DEV TIME ####################################################
timemod <- read.csv("daymodeldata.csv")
str(timemod)
timemod$block <- factor(timemod$block)
timemod$pop <- factor(timemod$pop)
str(timemod)

#### checks of data
hist(timemod$ebony.day50)
hist(timemod$WT.day50)
#### one data point at 0.6 so remove
timemod = timemod[which(timemod$ebony.day50 >6),]
timemod = timemod[which(timemod$WT.day50 >6),]
#### hists look normal enough now

#### NEED AVERGAE AND STANDARD ERROR DATA FROM HERE?#####
#sex
summarySE(timemod, measurevar=c("ebony.day50"), groupvars=c("sex"))
summarySE(timemod, measurevar=c("WT.day50"), groupvars=c("sex"))

# treat
summarySE(timemod, measurevar=c("ebony.day50"), groupvars=c("treat"))
summarySE(timemod, measurevar=c("WT.day50"), groupvars=c("treat"))

# pop
summarySE(timemod, measurevar=c("ebony.day50"), groupvars=c("pop"))
summarySE(timemod, measurevar=c("WT.day50"), groupvars=c("pop"))

# block
summarySE(timemod, measurevar=c("ebony.day50"), groupvars=c("block"))
summarySE(timemod, measurevar=c("WT.day50"), groupvars=c("block"))

# sex and treat
summarySE(timemod, measurevar=c("ebony.day50"), groupvars=c("sex", "treat"))
summarySE(timemod, measurevar=c("WT.day50"), groupvars=c("sex", "treat"))


##############################################################
##### model1 ebony ###############################

eeday1 <- lmer(ebony.day50 ~ treat*sex + (1 | pop) + (1 | block) + (1 | vial), data=timemod) 

res <- resid(eeday1)
ran <- ranef(eeday1)
fit <- fitted(eeday1)
par(mfrow=c(1,4))
hist(ran$pop[,1], main='')
hist(ran$block[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(eeday1)
anova(eeday1)

wald.test(b=fixef(eeday1), Sigma=vcov(eeday1), Terms = 2, df=8)
wald.test(b=fixef(eeday1), Sigma=vcov(eeday1), Terms = 3, df=8)
wald.test(b=fixef(eeday1), Sigma=vcov(eeday1), Terms = 3, df=223)
wald.test(b=fixef(eeday1), Sigma=vcov(eeday1), Terms = 4, df=8)

#### an affect of sex (males longer) and treatment (Sa longer) but no interaction

#### lsmeans for ebony 50% day
lsmeans(eeday1, "treat")
lsmeans(eeday1, "sex")
lsmeans(eeday1, list(pairwise ~ treat|sex))

##### model2 wildtype ################### 

WTday1 <- lmer(WT.day50 ~ treat*sex + (1 | pop) + (1 | block) + (1 | vial), data=timemod) 

res <- resid(WTday1)
ran <- ranef(WTday1)
fit <- fitted(WTday1)
par(mfrow=c(1,4))
hist(ran$pop[,1], main='')
hist(ran$block[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(WTday1)
anova(WTday1)

wald.test(b=fixef(WTday1), Sigma=vcov(WTday1), Terms = 2, df=8)
wald.test(b=fixef(WTday1), Sigma=vcov(WTday1), Terms = 3, df=223)
wald.test(b=fixef(WTday1), Sigma=vcov(WTday1), Terms = 4, df=8)

#### same as for ebony but the treatment if not significant but almost so

#### ls means for wt 50% day #######
lsmeans(WTday1, "treat")
lsmeans(WTday1, "sex")
lsmeans(WTday1, list(pairwise ~ treat|sex))
