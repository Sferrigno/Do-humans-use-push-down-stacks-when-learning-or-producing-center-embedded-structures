# REDACTED - Creating recursive, center-embedded sequences with an iterative queue.

############################################################
#### Load Libs, Import & Prep Data Frame 
############################################################
###Load Libraries 
library(ggplot2)
library(rstudioapi)
library(plyr)
library(Matrix)
library(lme4)
library(emmeans)
library(sjPlot)
library(lmerTest)

### Set theme for plots 
SAG_theme <- theme(legend.title=element_text( size = 14, face="plain"), 
                   legend.text=element_text(size = 12),
                   axis.title.x = element_text(size=14),
                   axis.text.x=element_text(colour="black", size = 14), 
                   axis.title.y = element_text(size = 14, angle = 90),
                   axis.text.y  = element_text(size = 12),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
                   axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"),
                   plot.title = element_text(hjust =.5,size = 18 ),
                   strip.background = element_blank(), strip.placement = "outside",
                   strip.text.x = element_text(size = 14)
)
theme_set(SAG_theme)

#Get Data
current_path <- getActiveDocumentContext()$path   
setwd(dirname(current_path))                      
d <- read.csv('Data_ForReview.csv', na = c("", "NA"))


############################################################
#### Phase 1: Training Results
############################################################
# Make dataframe with just Training Trials 
dTrainingTrials <- d[ d$round == "TrainingList1" | d$round == "TrainingList2", ]

#Summary by condition*array order (for each participant)
TrainingSummaryCondArraySubs<- ddply(dTrainingTrials, c( "Condition", "participant", "round"), summarise,
                            TrialsToCrit = max(trialRound, na.rm = T)
)

#Summary by condition (for each participant)
TrainingSummarySubs<- ddply(TrainingSummaryCondArraySubs, c( "Condition", "participant"), summarise,
                            TrialsToCrit = mean(TrialsToCrit, na.rm = T)
)

#Plot Trials to Criterion  Data 
boxplot(TrialsToCrit~Condition,data=TrainingSummarySubs, main="Trials to Criterion",
        xlab="Condition", ylab="Average Trials to Criterion") 

#Find and save outliers (12 total, 5 center-embedded, 7 cross-serial)
out <- boxplot.stats(TrainingSummarySubs$TrialsToCrit)$out
out_ind <- which(TrainingSummarySubs$TrialsToCrit %in% c(out))
TrainingOutliers <- TrainingSummarySubs[out_ind, ]

#save data with outliers removed
TrainingSummarySubsOutliersRemoved  <- TrainingSummarySubs[-which(TrainingSummarySubs$participant %in% TrainingOutliers$participant),]
TrainingSummarySubsEachArrayOutliersRemoved <- TrainingSummaryCondArraySubs[-which(TrainingSummaryCondArraySubs$participant %in% TrainingOutliers$participant),]

#Get Mean and SD of trials to criterion (outliers removed)
TrainingSummaryByRound<- ddply(TrainingSummarySubsEachArrayOutliersRemoved, c( "round"), summarise,
                            MeanTrialsToCrit = mean(TrialsToCrit, na.rm = T),
                            StdevTrialsToCrit = sd(TrialsToCrit, na.rm = T)
)

# TTEST No effect of condition (as specified in the Preregistration: Outliers excluded: t = 1.1834, df = 85.646, p-value = 0.2399; All Subs included: t = 0.28567, df = 97.563, p-value = 0.7757) 
t.test(TrialsToCrit ~ Condition, data = TrainingSummarySubsOutliersRemoved)
t.test(TrialsToCrit ~ Condition, data = TrainingSummarySubs)
mean(TrainingSummarySubsOutliersRemoved$TrialsToCrit)
sd(TrainingSummarySubsOutliersRemoved$TrialsToCrit)

# LMM (array*condition): Analysis in paper (qualitatively similar results as both ttests above, this is the maximal structure that is identifiable)
#Sig effect of Round, No effect of condition or interaction
#Maximal data structures including round and condition were unidentifiable.
TrialsToCritNull<- lmer(TrialsToCrit ~ 1 + (1|participant), data = TrainingSummarySubsEachArrayOutliersRemoved)
TrialsToCritA<- lmer(TrialsToCrit ~ round + (1 |participant), data = TrainingSummarySubsEachArrayOutliersRemoved)
TrialsToCritB<- lmer(TrialsToCrit ~ round + Condition + (1 | participant), data = TrainingSummarySubsEachArrayOutliersRemoved)

anova(TrialsToCritNull, TrialsToCritA)
anova(TrialsToCritA, TrialsToCritB)
# Best fitting model is with round only

summary(TrialsToCritA)
confint(TrialsToCritA,method="Wald")
tab_model(TrialsToCritA, title = "Table S1 - Phase 1 - Trials to CriterionResults", dv.labels = c("Trials to Criterion"), file = "S1.doc")


############################################################
#### Phase 2: Novel Combination of 4-item arrays
############################################################

# Get Phase 2 test trials
dTest <- d[d$round == "test", ]

#Get Only Nondifferentially rewarded  generalization (test) trials (and remove extra rows)
dPhase2 <- dTest[dTest$Test_TrialType == "test", ]
dPhase2 <- subset(dPhase2, !is.na(dPhase2$FullyCorrectOrderAndColor))

### Remove exclusions ###
#Remove Subjects who touched the same image on 50% or more of the transfer trial (there are none of these)
dPhase2 <- subset(dPhase2, ExcludeTestGen != 1)

#Remove trials with >2 same item touches (There are 1 of these trials)
dPhase2 <- dPhase2[is.na(dPhase2$MaxSameItem) | dPhase2$MaxSameItem == 2 , ]


####Phase 2 Stats
#Average accuracy for each subject
dTransferTrialsSum <- ddply(dPhase2 , c( "Condition", "participant"), summarise,
                            Accuracy = mean(FullyCorrectOrderAndColor, na.rm = T))

#TTESTs Condition compared to all red then blue
t.test(subset(dTransferTrialsSum$Accuracy, dTransferTrialsSum$Condition == "Center"), mu = 0.5)
t.test(subset(dTransferTrialsSum$Accuracy, dTransferTrialsSum$Condition == "Crossed"), mu = 0.5)

#TTEST Between Conditions
t.test(subset(dTransferTrialsSum$Accuracy, dTransferTrialsSum$Condition == "Crossed"),subset(dTransferTrialsSum$Accuracy, dTransferTrialsSum$Condition == "Center"), var.equal = TRUE)

# Qualitatively similar results using mixed effects logistic regression (See SI)
Test.TrialsMAX.glm <- glmer(FullyCorrectOrderAndColor ~ Condition  + (1 + Condition|participant), data = dPhase2, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
# Above Maximal model with randoms slopes is singular
Test.Trials.glm <- glmer(FullyCorrectOrderAndColor ~ Condition  + (1 |participant), data = dPhase2, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
Test.Trials.glm.Null <- glmer(FullyCorrectOrderAndColor ~ 1 + (1 |participant), data = dPhase2, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
anova(Test.Trials.glm.Null, Test.Trials.glm)

summary(Test.Trials.glm)
confint(Test.Trials.glm,method="Wald")
tab_model(Test.Trials.glm, title = "Table S2 - Phase 2 - Novel combination of base pairs ", dv.labels = c("Correct Generalizations"), file = "S2.doc")



############################################################
#### Phase 3: Generalization trials. Completely novel 4-item, 6-item, & 8-item arrays
############################################################

### Get Generalization  Trials
d4Gen <- d[d$round == "genList4", ]
d6Gen <- d[d$round == "genList6", ]
d8Gen <- d[d$round == "genList8", ]

#Remove extra rows that dont have the full order
d4Gen <- subset(d4Gen, !is.na(d4Gen$Order1))
d6Gen <- subset(d6Gen, !is.na(d6Gen$Order1))
d8Gen <- subset(d8Gen, !is.na(d8Gen$Order1))

### Remove exclusions ###
#Remove Subjects who touched the same image on 50% or more of the transfer trial
d4Gen <- subset(d4Gen, ExcludeGen4 == 0) #7 total (1 Center-embedded) 
d6Gen <- subset(d6Gen, ExcludeGen6 == 0) #6 total (1 center-embedded; 2 cross-serial)
d8Gen <- subset(d8Gen, ExcludeGen8 == 0) #9 total (5 center-embedded; 1 cross-serial)

#Remove trials with >2 same item touches
d4Gen <-d4Gen[is.na(d4Gen$MaxSameItem) | d4Gen$MaxSameItem == 2 , ] # None
d6Gen <-d6Gen[is.na(d6Gen$MaxSameItem) | d6Gen$MaxSameItem == 2 , ] # 2 total trials
d8Gen <-d8Gen[is.na(d8Gen$MaxSameItem) | d8Gen$MaxSameItem == 2 , ] # 2 total trials

#reCombine Generalization lists
dGen <- rbind(d4Gen, d6Gen, d8Gen)
#reCombine Generalization lists & Phase 2 list
dGenPlusPhase2 <- rbind(dPhase2,d4Gen, d6Gen, d8Gen)

#Make summary tables
dGenSum <- ddply(dGen , c( "Condition", "round", "participant"), summarise,
                 Accuracy = mean(FullyCorrectOrderAndColor, na.rm = T))

dGenSumOverall <- ddply(dGenSum , c( "Condition", "round"), summarise,
                        OverallAccuracy = mean(Accuracy, na.rm = T),
                        SD = sd(Accuracy, na.rm = T),
                        SE = sd(Accuracy, na.rm = T)/ sqrt(length(participant)),
                        NumParticipants = length(participant)
)
dGenSumOverall

#Make summary tables w/ phase 2
dGenPlusPhase2Sum <- ddply(dGenPlusPhase2 , c( "Condition", "round", "participant"), summarise,
                 Accuracy = mean(FullyCorrectOrderAndColor, na.rm = T))

dGenPlusPhase2SumOverall <- ddply(dGenPlusPhase2Sum , c( "Condition", "round"), summarise,
                        OverallAccuracy = mean(Accuracy, na.rm = T),
                        SD = sd(Accuracy, na.rm = T),
                        SE = sd(Accuracy, na.rm = T)/ sqrt(length(participant)),
                        NumParticipants = length(participant)
)
dGenPlusPhase2SumOverall


### 4-item Generalization  Trials
d4GenSum <- ddply(d4Gen , c( "Condition",  "participant"), summarise,
                 Accuracy = mean(FullyCorrectOrderAndColor, na.rm = T))

#TTESTs Condition compared to all red then blue
t.test(subset(d4GenSum$Accuracy, d4GenSum$Condition == "Center"), mu = 0.5)
t.test(subset(d4GenSum$Accuracy, d4GenSum$Condition == "Crossed"), mu = 0.5)
#between groups
t.test(subset(d4GenSum$Accuracy, d4GenSum$Condition == "Crossed"),subset(d4GenSum$Accuracy, d4GenSum$Condition == "Center"))


### 6-item Generalization  Trials
d6GenSum <- ddply(d6Gen , c( "Condition",  "participant"), summarise,
                  Accuracy = mean(FullyCorrectOrderAndColor, na.rm = T))

#TTESTs Condition compared to all red then blue
t.test(subset(d6GenSum$Accuracy, d6GenSum$Condition == "Center"), mu = (1/6))
t.test(subset(d6GenSum$Accuracy, d6GenSum$Condition == "Crossed"), mu = (1/6))


### 8-item Generalization Trials
d8GenSum <- ddply(d8Gen , c( "Condition",  "participant"), summarise,
                  Accuracy = mean(FullyCorrectOrderAndColor, na.rm = T))

#TTESTs Condition compared to all red then blue
t.test(subset(d8GenSum$Accuracy, d8GenSum$Condition == "Center"), mu = (1/24))
t.test(subset(d8GenSum$Accuracy, d8GenSum$Condition == "Crossed"), mu = (1/24))


###############
# Logistic Regression Accuracy by condition * array length, Random effects term of Participant (PREREGISTERED, see below for maximal model - qualitatively similar results)
###############
LogReg<- glmer(FullyCorrectOrderAndColor ~ Condition * round +  (1 | participant), data = dGen, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
summary(LogReg)
confint(LogReg ,method="Wald")
tab_model(LogReg, title = "Table S1: Phase 3 - Logistic Regression Results", dv.labels = c("Correct Generalizations"))

#Pairwise tests - Overall effect of Condition
condition.emm.s <- emmeans(LogReg, ~ Condition)
pairs(condition.emm.s)

#Pairwise tests - Overall effect of Array Length
round.emm.s <- emmeans(LogReg, ~ round)
pairs(round.emm.s)

#Pairwise tests -  Array Length * Condition
conditionround.emm.s <- emmeans(LogReg, ~ Condition * round)
pairs(conditionround.emm.s)
#Make table for SI
xtable::xtable(conditionround.emm.s, type = "response")

###############
# Maximal random effects models - Logistic Regression Accuracy by condition * array length
###############

#Mixed effects logistic regressions
MaxLogRegNull<- glmer(FullyCorrectOrderAndColor ~ 1 + (1  | participant), data = dGen, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
MaxLogRegA<- glmer(FullyCorrectOrderAndColor ~   round + (1 + round | participant), data = dGen, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
MaxLogRegB<- glmer(FullyCorrectOrderAndColor ~ Condition + round + (1 + Condition + round | participant), data = dGen, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
MaxLogRegC<- glmer(FullyCorrectOrderAndColor ~  Condition * round + (1 + Condition + round | participant), data = dGen, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
MaxLogRegD<- glmer(FullyCorrectOrderAndColor ~  Condition * round + (1 + Condition * round | participant), data = dGen, family=binomial,  control = glmerControl(optimizer = "bobyqa") )

#Model comparisons
anova(MaxLogRegNull, MaxLogRegA, MaxLogRegB,  MaxLogRegC, MaxLogRegD)
### Model C is the best fit. Interaction term as fixed effects, but additional variance from random effects
summary(MaxLogRegC)
tab_model(MaxLogRegC, title = "Table S3: Phase 3 - Completely novel 4-item, 6-item, & 8-item arrays", dv.labels = c("Correct Generalizations"), file = "S3.doc")

#Pairwise tests - Overall effect of Condition
Maxcondition.emm.s  <- emmeans(MaxLogRegC, ~ Condition)
pairs(Maxcondition.emm.s )

#Pairwise tests - Overall effect of Array Length
Maxround.emm.s <- emmeans(MaxLogRegC, ~ round)
pairs(Maxround.emm.s)

#Pairwise tests -  Array Length * Condition
Maxconditionround.emm.s <- emmeans(MaxLogRegC, ~ Condition * round)
pairs(Maxconditionround.emm.s)

##### Figure 4
#Reorder conditions
dGenPlusPhase2SumOverall$round <- as.character(dGenPlusPhase2SumOverall$round)
dGenPlusPhase2SumOverall$round <- factor(dGenPlusPhase2SumOverall$round, levels=c("test","genList4", "genList6", "genList8"))

g = ggplot (data = dGenPlusPhase2SumOverall , aes(x=Condition, y = OverallAccuracy, fill = Condition)) +
  facet_wrap(~ round, strip.position = "top", nrow = 1) +
  scale_colour_manual(values=c("#009E73", "red")) +
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=c("#009E73", "red")) +
  geom_errorbar(aes(ymin=OverallAccuracy-SE, ymax=OverallAccuracy+SE), width=.01) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(y = "Proportion Correct") +
  labs(x = "Condition") +
  ggtitle ("Generalization trials"
)

#Make chance lines dataframe and add to graph
chance_hlines <- data.frame(round = unique(dGenPlusPhase2SumOverall$round),  # Create data for lines
                            hline = c(1/2, 1/6, 1/24, 1/2),                        
                            
                            hline2 = c(1/12, 1/120, 1/1680, 1/12)
                            )
g + geom_hline(data = chance_hlines, aes(yintercept = hline),linetype="dotted") +
    geom_hline(data = chance_hlines, aes(yintercept = hline2), linetype="dashed")


############################################################
#### Response Time analyses by Item # (correct sequences in the further training trials: 6- & and 8-item arrays)
############################################################

#Get correct trials only
dcorrect <-subset(d, d$TrialAccuracyRepeated == 1)

#Split by length
SixItemsRT <- dcorrect [dcorrect$round == "6Test", ]
EightItemsRT <- dcorrect [dcorrect$round == "8Test",]

###Remove Exclusions
#Remove excluded subjects (low accuracy)
SixItemsRT <- subset(SixItemsRT, SixItemsRT$Exclude6item == 0) 
EightItemsRT <- subset(EightItemsRT, EightItemsRT$Exclude8item == 0) 

#Test with different Exclusion criteria (from review, toggle below ON and above OFF to run)
#SixItemsRT <- subset(SixItemsRT, SixItemsRT$Exclude6itemTest_80Acc == 0) 
#EightItemsRT <- subset(EightItemsRT, EightItemsRT$Exclude8itemTest_80Acc == 0) 

#Remove long trials >2 SD from their own average
SixItemsRT <- subset(SixItemsRT, SixItemsRT$CorrectTrialTotalRTRepeated < SixItemsRT$OverallRTLimit6) 
EightItemsRT <- subset(EightItemsRT, EightItemsRT$CorrectTrialTotalRTRepeated < EightItemsRT$OverallRTLimit8) 

#Find & Remove individual touches that are >2 sd away from the subjects own mean
SixItemsRT <- subset(SixItemsRT, SixItemsRT$rt < SixItemsRT$X6RTLimit) 
EightItemsRT <- subset(EightItemsRT, EightItemsRT$rt < EightItemsRT$X8RTLimit) 


##### 6 item arrays
#Get just the end items (with and without the first item in the 2nd half)
SixItemsRTEndItems<- SixItemsRT [ SixItemsRT$Touch == "5" | SixItemsRT$Touch == "6" ,]
SixItemsRTEndItemsW4<- SixItemsRT [ SixItemsRT$Touch == "4" |SixItemsRT$Touch == "5" | SixItemsRT$Touch == "6" ,]

### Mixed effects model (with participant as random effects term)
rt6 <- lmer(rt ~ Condition * Touch + (1 | participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
summary(rt6)

### Mixed effects models
MaxRegRT6Null <- lmer(rt ~ 1 + (1 | participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
MaxRegRT6A <- lmer(rt ~  Condition + (1 + Condition| participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
#Above doesn't converge, removing random effects of condition
MaxRegRT6A <- lmer(rt ~  Condition + (1 | participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
MaxRegRT6B <- lmer(rt ~ Condition + Touch + (1 + Condition + Touch| participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
#Above doesn't converge, removing random effects of condition
MaxRegRT6B <- lmer(rt ~ Condition + Touch + (1 +  Touch| participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
MaxRegRT6C <- lmer(rt ~ Condition * Touch + (1 + Condition * Touch| participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
#Above doesn't converge,removing random effects starting with interaction
MaxRegRT6C <- lmer(rt ~ Condition * Touch + (1 + Condition + Touch| participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
#Above doesn't converge,removing random effects starting with interaction
MaxRegRT6C <- lmer(rt ~ Condition * Touch + (1 + Touch| participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))

###Model Comparisons
anova(MaxRegRT6Null, MaxRegRT6A)
anova(MaxRegRT6A, MaxRegRT6B)
anova(MaxRegRT6B, MaxRegRT6C)
#Best fitting model is C 
summary(MaxRegRT6C)


#Get Data Ready to Graph
SixItemSumRT <- ddply(SixItemsRTEndItems , c( "Condition", "Touch", "participant"), summarise,
                         RT = mean(rt, na.rm = T))

SixItemSumRTOverall <- ddply(SixItemSumRT , c( "Condition", "Touch"), summarise,
                                OverallRT = mean(RT, na.rm = T),
                                SD = sd(RT, na.rm = T),
                                SE = sd(RT, na.rm = T)/ sqrt(length(participant)),
                                NumParticipants = length(participant)
)

#Get Data Ready to Graph
SixItemSumRTW4 <- ddply(SixItemsRTEndItemsW4 , c( "Condition", "Touch", "participant"), summarise,
                      RT = mean(rt, na.rm = T))

SixItemSumRTOverallW4 <- ddply(SixItemSumRTW4 , c( "Condition", "Touch"), summarise,
                             OverallRT = mean(RT, na.rm = T),
                             SD = sd(RT, na.rm = T),
                             SE = sd(RT, na.rm = T)/ sqrt(length(participant)),
                             NumParticipants = length(participant)
)


#Plot 6 item stuff    The Linear fit in the graph just uses subject's averages (and ignores trials per subject - see lmer above & SI for full model
g = ggplot (data = SixItemSumRTOverallW4, aes(x=Touch, y = OverallRT, colour = Condition)) +
  scale_colour_manual(values=c("#009E73", "red")) +
  geom_smooth(data = SixItemSumRT, aes(x=Touch, y = RT, colour = Condition), method = "lm", se=TRUE) +
  geom_point( stat="identity", size = 3) +
  scale_fill_manual(values=c("#009E73", "red")) +
  geom_errorbar(aes(ymin=OverallRT-SE, ymax=OverallRT+SE), width=.01) +
  #scale_x_continuous(limits = c(4,5,6)) +
  coord_cartesian(ylim=c(0, 1500), xlim=c(4,6)) +
  scale_x_continuous(breaks=seq(4, 6, 1))+
  labs(y = "RT (ms)") +
  labs(x = "Touch Number")+
  ggtitle ("6-item lists")
g


##### 8 item Stuff
#Get just the end items (with and without the first item in the 2nd half)
EightItemsRTEndItems<- EightItemsRT [EightItemsRT$Touch == "6" | EightItemsRT$Touch == "7" | EightItemsRT$Touch == "8",]
EightItemsRTEndItemsW5<- EightItemsRT [EightItemsRT$Touch == "5" |EightItemsRT$Touch == "6" | EightItemsRT$Touch == "7" | EightItemsRT$Touch == "8",]


### Mixed effects model (with participant as random effects term) 
rt8 <- lmer(rt ~ Condition * Touch + (1 + Touch| participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
summary(rt8)

### Mixed effects models 
MaxRegRT8Null <- lmer(rt ~ 1 + (1 | participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
MaxRegRT8A <- lmer(rt ~  Condition + (1 + Condition| participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
#Above doesn't converge, removing random effects of condition
MaxRegRT8A <- lmer(rt ~  Condition + (1 | participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
MaxRegRT8B <- lmer(rt ~ Condition + Touch + (1 + Condition + Touch| participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
#Above doesn't converge, removing random effects of condition
MaxRegRT8B <- lmer(rt ~ Condition + Touch + (1 +  Touch| participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
MaxRegRT8C <- lmer(rt ~ Condition * Touch + (1 + Condition * Touch| participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
#Above doesn't converge,removing random effects starting with interaction
MaxRegRT8C <- lmer(rt ~ Condition * Touch + (1 + Condition + Touch| participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
#Above doesn't converge,removing random effects starting with interaction
MaxRegRT8C <- lmer(rt ~ Condition * Touch + (1 + Touch| participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))

###Model Comparisons
anova(MaxRegRT8Null, MaxRegRT8A)
anova(MaxRegRT8A, MaxRegRT8B)
anova(MaxRegRT8B, MaxRegRT8C)
#Best fitting model is C (w/ interaction)
summary(MaxRegRT8C)

#Get Data Ready to Graph
EightItemSumRT <- ddply(EightItemsRTEndItems , c( "Condition", "Touch", "participant"), summarise,
                      RT = mean(rt, na.rm = T))

EightItemSumRTOverall <- ddply(EightItemSumRT , c( "Condition", "Touch"), summarise,
                             OverallRT = mean(RT, na.rm = T),
                             SD = sd(RT, na.rm = T),
                             SE = sd(RT, na.rm = T)/ sqrt(length(participant)),
                             NumParticipants = length(participant)
)
EightItemSumRTW5 <- ddply(EightItemsRTEndItemsW5 , c( "Condition", "Touch", "participant"), summarise,
                        RT = mean(rt, na.rm = T))

EightItemSumRTOverallW5 <- ddply(EightItemSumRTW5 , c( "Condition", "Touch"), summarise,
                               OverallRT = mean(RT, na.rm = T),
                               SD = sd(RT, na.rm = T),
                               SE = sd(RT, na.rm = T)/ sqrt(length(participant)),
                               NumParticipants = length(participant)
)


#Plot 8 item stuff    The Linear fit in the graph just uses subject's averages (and ignores trials per subject - see lmer above & SI for full model
g = ggplot (data = EightItemSumRTOverallW5, aes(x=Touch, y = OverallRT, colour = Condition)) +
  scale_colour_manual(values=c("#009E73", "red")) +
  geom_smooth(data = EightItemSumRT, aes(x=Touch, y = RT, colour = Condition), method = "lm", se=TRUE) +
  geom_point( stat="identity", size = 3) +
  scale_fill_manual(values=c("#009E73", "red")) +
  geom_errorbar(aes(ymin=OverallRT-SE, ymax=OverallRT+SE), width=.01) +
  #scale_x_continuous(limits = c(4,5,6)) +
  coord_cartesian(ylim=c(0, 1500), xlim=c(5,8)) +
  scale_x_continuous(breaks=seq(5, 8, 1))+
  labs(y = "RT (ms)") +
  labs(x = "Touch Number")+
  ggtitle ("8-item lists")
g

tab_model(MaxRegRT6C, MaxRegRT8C,  title = "Table S7: Response time by item number (80% Accuracy Criterion)", dv.labels = c("6-item Array", "8-item Array"), file = "80Acc RT Graph.doc")


############################################################
####Further Training Accuracy by condition (for each length) Not included in main text  (see SI)
############################################################
###Get Further Training Trials
FourItemsAcc <- d [d$round == "4Test", ]
SixItemsAcc <- d [d$round == "6Test", ]
EightItemsAcc <- d [d$round == "8Test",]

#Remove extra rows (Just keep rows with OverallTrial Accuracy)
FourItemsAcc <- subset(FourItemsAcc, !is.na(FourItemsAcc$TrialAccuracy))
SixItemsAcc <- subset(SixItemsAcc, !is.na(SixItemsAcc$TrialAccuracy))
EightItemsAcc <- subset(EightItemsAcc, !is.na(EightItemsAcc$TrialAccuracy))

#Remove excluded subjects (low accuracy)
FourItemsAcc <- subset(FourItemsAcc, FourItemsAcc$Exclude4item == 0) 
SixItemsAcc <- subset(SixItemsAcc, SixItemsAcc$Exclude6item == 0) 
EightItemsAcc <- subset(EightItemsAcc, EightItemsAcc$Exclude8item == 0) 

#reCombine Generalization lists
dFurtherTraining <- rbind(FourItemsAcc, SixItemsAcc, EightItemsAcc)

#Make summary tables
dFurtherTrainingSum <- ddply(dFurtherTraining , c( "Condition", "round", "participant"), summarise,
                             Accuracy = mean(TrialAccuracy, na.rm = T))

dFurtherTrainingSumOverall <- ddply(dFurtherTrainingSum , c( "Condition", "round"), summarise,
                                    OverallAccuracy = mean(Accuracy, na.rm = T),
                                    SD = sd(Accuracy, na.rm = T),
                                    SE = sd(Accuracy, na.rm = T)/ sqrt(length(participant)),
                                    NumParticipants = length(participant)
)
dFurtherTrainingSumOverall


### 4-item Further Training  Trials
FourItemsAccSum <- ddply(FourItemsAcc , c( "Condition",  "participant"), summarise,
                         Accuracy = mean(TrialAccuracy, na.rm = T))

#TTESTs Condition compared to all red then blue
t.test(subset(FourItemsAccSum$Accuracy, FourItemsAccSum$Condition == "Center"), mu = 0.5)
t.test(subset(FourItemsAccSum$Accuracy, FourItemsAccSum$Condition == "Crossed"), mu = 0.5)

### 6-item Generalization  Trials
SixItemsAccSum <- ddply(SixItemsAcc , c( "Condition",  "participant"), summarise,
                        Accuracy = mean(TrialAccuracy, na.rm = T))

#TTESTs Condition compared to all red then blue
t.test(subset(SixItemsAccSum$Accuracy, SixItemsAccSum$Condition == "Center"), mu = (1/6))
t.test(subset(SixItemsAccSum$Accuracy, SixItemsAccSum$Condition == "Crossed"), mu = (1/6))

### 8-item Generalization Trials
EightItemsAccSum <- ddply(EightItemsAcc , c( "Condition",  "participant"), summarise,
                          Accuracy = mean(TrialAccuracy, na.rm = T))

#TTESTs Condition compared to all red then blue
t.test(subset(EightItemsAccSum$Accuracy, EightItemsAccSum$Condition == "Center"), mu = (1/24))
t.test(subset(EightItemsAccSum$Accuracy, EightItemsAccSum$Condition == "Crossed"), mu = (1/24))

### Logistic Regression - Accuracy by condition * array length, Random effects term of Participant (Preregistered)
LogRegFurtherTraining<- glmer(TrialAccuracy ~ Condition * round +  (1 | participant), data = dFurtherTraining, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
summary(LogRegFurtherTraining)
tab_model(LogRegFurtherTraining)

#With maximal random effects 
MaxLogRegFurtherTrainingNull<- glmer(TrialAccuracy ~ 1 + (1  | participant), data = dFurtherTraining, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
MaxLogRegFurtherTrainingA<- glmer(TrialAccuracy ~   round + (1 +  round | participant), data = dFurtherTraining, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
MaxLogRegFurtherTrainingB<- glmer(TrialAccuracy ~ Condition + round + (1 + Condition + round | participant), data = dFurtherTraining, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
MaxLogRegFurtherTrainingC<- glmer(TrialAccuracy ~  Condition * round + (1 + Condition * round | participant), data = dFurtherTraining, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
#Model C has singular fit, removing random effects starting with interaction
MaxLogRegFurtherTrainingC<- glmer(TrialAccuracy ~  Condition * round + (1 + Condition + round | participant), data = dFurtherTraining, family=binomial,  control = glmerControl(optimizer = "bobyqa") )

###Model Comparisons
anova(MaxLogRegFurtherTrainingNull, MaxLogRegFurtherTrainingA)
anova(MaxLogRegFurtherTrainingA, MaxLogRegFurtherTrainingB)
anova(MaxLogRegFurtherTrainingB, MaxLogRegFurtherTrainingC)
#Best fitting model is B (no interaction)

tab_model(MaxLogRegFurtherTrainingC, title = "Table S8: Further training trials", dv.labels = c("Trial Accuracy"), file = "S8.doc")

# Figure S1
g = ggplot (data = dFurtherTrainingSumOverall , aes(x=Condition, y = OverallAccuracy, fill = Condition)) +
  facet_wrap(~ round, strip.position = "top") +
  scale_colour_manual(values=c("#009E73", "red")) +
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=c("#009E73", "red")) +
  geom_errorbar(aes(ymin=OverallAccuracy-SE, ymax=OverallAccuracy+SE), width=.01) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(y = "Proportion Correct") +
  labs(x = "Condition") +
  ggtitle ("Further Training Trials")

#Make chance lines dataframe and add to graph
chance_hlines <- data.frame(round = unique(dFurtherTrainingSumOverall$round),  # Create data for lines
                            hline = c(1/2, 1/6, 1/24),
                            hline2 = c(1/12, 1/120, 1/1680)
)
g + geom_hline(data = chance_hlines, aes(yintercept = hline),linetype="dotted") +
  geom_hline(data = chance_hlines, aes(yintercept = hline2), linetype="dashed")


############################################################
####Accuracy by condition * Touch (for each length) Not included in main text  (see SI)
############################################################
#Split by length
SixItemsAcc <- d [d$round == "6Test", ]
EightItemsAcc <- d [d$round == "8Test",]

###Need to remove exclusions here (as specified in the prereg.) 
#Remove excluded subjects (low accuracy)
SixItemsAcc <- subset(SixItemsAcc, SixItemsAcc$Exclude6item == 0) 
EightItemsAcc <- subset(EightItemsAcc, EightItemsAcc$Exclude8item == 0) 

# Get 2nd half of 6- and 8-item trials
SixItemsAccEndItems<- SixItemsAcc [ SixItemsAcc$Touch == "4" |SixItemsAcc$Touch == "5" | SixItemsAcc$Touch == "6" ,]
EightItemsAccEndItems<- EightItemsAcc  [EightItemsAcc$Touch == "5" |EightItemsAcc$Touch == "6" | EightItemsAcc$Touch == "7" | EightItemsAcc$Touch == "8",]

### Logistic Regression predicting Accuracy by condition and touch # - 6-item array
Acc6Null <- glmer(TouchAccuracy ~  1 + (1 | participant), data = SixItemsAccEndItems, family=binomial, control = glmerControl(optimizer="bobyqa"))
Acc6A <- glmer(TouchAccuracy ~  Touch + (1 + Touch| participant), data = SixItemsAccEndItems, family=binomial, control = glmerControl(optimizer="bobyqa"))
Acc6B <- glmer(TouchAccuracy ~ Condition + Touch + (1 + Condition + Touch| participant), data = SixItemsAccEndItems, family=binomial, control = glmerControl(optimizer="bobyqa"))
#above doesn't converge
Acc6B <- glmer(TouchAccuracy ~ Condition + Touch + (1 + Touch| participant), data = SixItemsAccEndItems, family=binomial, control = glmerControl(optimizer="bobyqa"))
Acc6C <- glmer(TouchAccuracy ~ Condition * Touch + (1 + Condition * Touch| participant), data = SixItemsAccEndItems, family=binomial, control = glmerControl(optimizer="bobyqa"))
#above is singular (random interaction term removed)
Acc6C <- glmer(TouchAccuracy ~ Condition * Touch + (1 + Condition + Touch| participant), data = SixItemsAccEndItems, family=binomial, control = glmerControl(optimizer="bobyqa"))
#above doesn't converge
Acc6C <- glmer(TouchAccuracy ~ Condition * Touch + (1 + Touch| participant), data = SixItemsAccEndItems, family=binomial, control = glmerControl(optimizer="bobyqa"))

anova(Acc6Null,Acc6A)
anova(Acc6A,Acc6B)
anova(Acc6B,Acc6C)
### Model B is the best fit. (no interaction term)

###  Logistic Regression predicting Accuracy by condition and touch # - 8-item array
Acc8Null <- glmer(TouchAccuracy ~ 1 + (1 | participant), data = EightItemsAccEndItems, family=binomial, control = glmerControl(optimizer="bobyqa"))
Acc8A <- glmer(TouchAccuracy ~  Touch + (1 + Touch | participant), data = EightItemsAccEndItems, family=binomial, control = glmerControl(optimizer="bobyqa"))
Acc8B <- glmer(TouchAccuracy ~ Condition + Touch + (1 + Condition + Touch | participant), data = EightItemsAccEndItems, family=binomial, control = glmerControl(optimizer="bobyqa"))
Acc8C <- glmer(TouchAccuracy ~ Condition * Touch + (1 + Condition * Touch | participant), data = EightItemsAccEndItems, family=binomial, control = glmerControl(optimizer="bobyqa"))
summary(Acc8)

anova(Acc8Null,Acc8A)
anova(Acc8A,Acc8B)
anova(Acc8B,Acc8C)
### Model B is the best fit. (no interaction term)

tab_model(Acc6B, Acc8B,  title = "Table S9: Accuracy by Item Number - Second half of list", dv.labels = c("6-item Array", "8-item Array"), file = "S9.doc")



############################################################
#### OVERALL Response Time analyses - not included in main text because it is part of the larger logistic model of RT (see SI for the below analyses)
############################################################

# Get Just Correct Trials
dcorrect <- subset(d, d$TrialAccuracyRepeated == 1)

#Separate into conditions
FourItems <- dcorrect[dcorrect$round == "4Test", ]
SixItems <- dcorrect[dcorrect$round == "6Test", ]
EightItems <- dcorrect[dcorrect$round == "8Test",]

#Remove excluded subs for each length  (Ones Accuracy > 2sd from the mean as specified in the prereg.)  
FourItemsE <- subset(FourItems, FourItems$Exclude4item == 0) 
SixItems <- subset(SixItems, SixItems$Exclude6item == 0) 
EightItems <- subset(EightItems, EightItems$Exclude8item == 0) 

##Remove any trial that has a total response time that is  >2 SD of the mean  
FourItems <- subset(FourItems, FourItems$CorrectTrialTotalRTRepeated < FourItems$OverallRTLimit4) 
SixItems <- subset(SixItems, SixItems$CorrectTrialTotalRTRepeated < SixItems$OverallRTLimit6) 
EightItems <- subset(EightItems, EightItems$CorrectTrialTotalRTRepeated < EightItems$OverallRTLimit8) 

#Recombine dataframes
dRT <- rbind(FourItems, SixItems, EightItems )


### Overall RTs for Second Half of the list
FourItemsEndItems<- FourItems [FourItems$Touch == "3" | FourItems$Touch == "4" ,]

rt4Null <- lmer(rt ~ 1 + (1 |participant), data = FourItemsEndItems, control = lmerControl(optimizer="bobyqa"))
rt4A <- lmer(rt ~ Condition + (1 + Condition|participant), data = FourItemsEndItems, control = lmerControl(optimizer="bobyqa"))
#Above does not converge, removing condition from random effects
rt4A <- lmer(rt ~ Condition + (1 |participant), data = FourItemsEndItems, control = lmerControl(optimizer="bobyqa"))
anova(rt4Null, rt4A)
#Null model is best fit

SixItemsEndItems<- SixItems [SixItems$Touch == "4" | SixItems$Touch == "5" | SixItems$Touch == "6" ,]
rt6Null <- lmer(rt ~ 1 + (1 |participant), data = SixItemsEndItems, control = lmerControl(optimizer="bobyqa"))
rt6A <- lmer(rt ~ Condition + (1 + Condition|participant), data = SixItemsEndItems, control = lmerControl(optimizer="bobyqa"))
anova(rt6Null, rt6A)
#Model A  has the best fit


EightItemsEndItems<- EightItems [EightItems$Touch == "5" | EightItems$Touch == "6" | EightItems$Touch == "7" | EightItems$Touch == "8",]
rt8Null <- lmer(rt ~ 1 +(1 |participant), data = EightItemsEndItems, control = lmerControl(optimizer="bobyqa"))
rt8A <- lmer(rt ~ Condition + (1 + Condition|participant), data = EightItemsEndItems, control = lmerControl(optimizer="bobyqa"))
#Above doesnt converge, removing condition from random effects
rt8A <- lmer(rt ~ Condition + (1 |participant), data = EightItemsEndItems, control = lmerControl(optimizer="bobyqa"))
anova(rt8Null, rt8A)
#Model A  has the best fit

# 6 and 8 item difference is significant
tab_model(rt4Null, rt6A, rt8A, title = "Table S10: Overall Response Time - Second half of lists", dv.labels = c("4-item List", "6-item List", "8-item List"), file = "S10.doc")


### Overall RTs for First Half of the list
FourItemsFirstItems<- FourItems [FourItems$Touch == "1" | FourItems$Touch == "2" ,]
rt4Null <- lmer(rt ~ 1 + (1 |participant), data = FourItemsFirstItems, control = lmerControl(optimizer="bobyqa"))
rt4A <- lmer(rt ~ Condition + (1 + Condition |participant), data = FourItemsFirstItems, control = lmerControl(optimizer="bobyqa"))
#Above did not converge
rt4A <- lmer(rt ~ Condition + (1 |participant), data = FourItemsFirstItems, control = lmerControl(optimizer="bobyqa"))

anova(rt4Null, rt4A)
#Null model is best fit

SixItemsFirstItems<- SixItems [SixItems$Touch == "1" | SixItems$Touch == "2" | SixItems$Touch == "3" ,]
rt6Null <- lmer(rt ~ 1 + (1 |participant), data = SixItemsFirstItems, control = lmerControl(optimizer="bobyqa"))
rt6A <- lmer(rt ~ Condition + (1 + Condition |participant), data = SixItemsFirstItems, control = lmerControl(optimizer="bobyqa"))

anova(rt6Null, rt6A)
#Null model is best fit


EightItemsFirstItems<- EightItems [EightItems$Touch == "1" | EightItems$Touch == "2" | EightItems$Touch == "3" | EightItems$Touch == "4",]
rt8Null <- lmer(rt ~ 1 + (1 |participant), data = EightItemsFirstItems, control = lmerControl(optimizer="bobyqa"))
rt8A <- lmer(rt ~ Condition + (1 + Condition |participant), data = EightItemsFirstItems, control = lmerControl(optimizer="bobyqa"))
#Above is singular, removing condition from random effects
rt8A <- lmer(rt ~ Condition + (1  |participant), data = EightItemsFirstItems, control = lmerControl(optimizer="bobyqa"))

anova(rt8Null, rt8A)
#Null model is best fit

# No Diffs in first half of lists
tab_model(rt4Null, rt6Null, rt8Null, title = "Table S11: Overall Response Time - First half of lists", dv.labels = c("4-item List", "6-item List", "8-item List"), file = "S11.doc")





############################################################
#### NO OUTLIERS REMOVED - Phase 3: Generalization trials. Completely novel 4-item, 6-item, & 8-item arrays
############################################################

### Get Generalization  Trials
d4Gen <- d[d$round == "genList4", ]
d6Gen <- d[d$round == "genList6", ]
d8Gen <- d[d$round == "genList8", ]

#Remove extra rows that dont have the full order
d4Gen <- subset(d4Gen, !is.na(d4Gen$Order1))
d6Gen <- subset(d6Gen, !is.na(d6Gen$Order1))
d8Gen <- subset(d8Gen, !is.na(d8Gen$Order1))

### Do Not Remove exclusions ###
#Remove Subjects who touched the same image on 50% or more of the transfer trial
#d4Gen <- subset(d4Gen, ExcludeGen4 == 0) #7 total (1 Center-embedded) 
#d6Gen <- subset(d6Gen, ExcludeGen6 == 0) #6 total (1 center-embedded; 2 cross-serial)
#d8Gen <- subset(d8Gen, ExcludeGen8 == 0) #9 total (5 center-embedded; 1 cross-serial)

#Remove trials with >2 same item touches
#d4Gen <-d4Gen[is.na(d4Gen$MaxSameItem) | d4Gen$MaxSameItem == 2 , ] # None
#d6Gen <-d6Gen[is.na(d6Gen$MaxSameItem) | d6Gen$MaxSameItem == 2 , ] # 2 total trials
#d8Gen <-d8Gen[is.na(d8Gen$MaxSameItem) | d8Gen$MaxSameItem == 2 , ] # 2 total trials

#reCombine Generalization lists
dGen <- rbind(d4Gen, d6Gen, d8Gen)
#reCombine Generalization lists & Phase 2 list
dGenPlusPhase2 <- rbind(d4Gen, d6Gen, d8Gen)

#Make summary tables
dGenSum <- ddply(dGen , c( "Condition", "round", "participant"), summarise,
                 Accuracy = mean(FullyCorrectOrderAndColor, na.rm = T))

dGenSumOverall <- ddply(dGenSum , c( "Condition", "round"), summarise,
                        OverallAccuracy = mean(Accuracy, na.rm = T),
                        SD = sd(Accuracy, na.rm = T),
                        SE = sd(Accuracy, na.rm = T)/ sqrt(length(participant)),
                        NumParticipants = length(participant)
)
dGenSumOverall


#Mixed effects logistic regressions
MaxLogRegNull<- glmer(FullyCorrectOrderAndColor ~ 1 + (1  | participant), data = dGen, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
MaxLogRegA<- glmer(FullyCorrectOrderAndColor ~   round + (1 + round | participant), data = dGen, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
MaxLogRegB<- glmer(FullyCorrectOrderAndColor ~ Condition + round + (1 + Condition + round | participant), data = dGen, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
MaxLogRegC<- glmer(FullyCorrectOrderAndColor ~  Condition * round + (1 + Condition * round | participant), data = dGen, family=binomial,  control = glmerControl(optimizer = "bobyqa") )
#Above model doesn't converge. Removing random effects interaction term 
MaxLogRegC<- glmer(FullyCorrectOrderAndColor ~  Condition * round + (1 + Condition + round | participant), data = dGen, family=binomial,  control = glmerControl(optimizer = "bobyqa") )



#Model comparisons
anova(MaxLogRegNull, MaxLogRegA, MaxLogRegB,  MaxLogRegC)
### Model C is the best fit. Include both main effects and interaction as fixed effects

summary(MaxLogRegC)
tab_model(MaxLogRegC, title = "Table S12: Phase 3 - Logistic Regression results. (No outliers excluded)", dv.labels = c("Correct Generalizations"), file = "S12.doc")

##### Figure S2
g = ggplot (data = dGenSumOverall , aes(x=Condition, y = OverallAccuracy, fill = Condition)) +
  facet_wrap(~ round, strip.position = "top", nrow = 1) +
  scale_colour_manual(values=c("#009E73", "red")) +
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=c("#009E73", "red")) +
  geom_errorbar(aes(ymin=OverallAccuracy-SE, ymax=OverallAccuracy+SE), width=.01) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(y = "Proportion Correct") +
  labs(x = "Condition") +
  ggtitle ("Fig. S2     Phase 3: Novel Arrays (no outliers removed)"
  )

#Make chance lines dataframe and add to graph
chance_hlines <- data.frame(round = unique(dGen$round),  # Create data for lines
                            hline = c(1/2, 1/6, 1/24),                        
                            
                            hline2 = c(1/12, 1/120, 1/1680)
)
g + geom_hline(data = chance_hlines, aes(yintercept = hline),linetype="dotted") +
  geom_hline(data = chance_hlines, aes(yintercept = hline2), linetype="dashed")






############################################################
#### NO OUTLIERS REMOVED - Response Time analyses by Item # (correct sequences in the further training trials: 6- & and 8-item arrays)
############################################################

#Get correct trials only
dcorrect <-subset(d, d$TrialAccuracyRepeated == 1)

#Split by length
SixItemsRT <- dcorrect [dcorrect$round == "6Test", ]
EightItemsRT <- dcorrect [dcorrect$round == "8Test",]

### DO NOT Remove Exclusions
#Remove excluded subjects (low accuracy)
#SixItemsRT <- subset(SixItemsRT, SixItemsRT$Exclude6item == 0) 
#EightItemsRT <- subset(EightItemsRT, EightItemsRT$Exclude8item == 0) 

#Remove long trials >2 SD from their own average
#SixItemsRT <- subset(SixItemsRT, SixItemsRT$CorrectTrialTotalRTRepeated < SixItemsRT$OverallRTLimit6) 
#EightItemsRT <- subset(EightItemsRT, EightItemsRT$CorrectTrialTotalRTRepeated < EightItemsRT$OverallRTLimit8) 

#Find & Remove individual touches that are >2 sd away from the subjects own mean
#SixItemsRT <- subset(SixItemsRT, SixItemsRT$rt < SixItemsRT$X6RTLimit) 
#EightItemsRT <- subset(EightItemsRT, EightItemsRT$rt < EightItemsRT$X8RTLimit) 


##### 6 item arrays
#Get just the end items (with and without the first item in the 2nd half)
SixItemsRTEndItems<- SixItemsRT [ SixItemsRT$Touch == "5" | SixItemsRT$Touch == "6" ,]
SixItemsRTEndItemsW4<- SixItemsRT [ SixItemsRT$Touch == "4" |SixItemsRT$Touch == "5" | SixItemsRT$Touch == "6" ,]

### Mixed effects model (with participant as random effects term)
rt6 <- lmer(rt ~ Condition * Touch + (1 | participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
summary(rt6)

### Mixed effects models
MaxRegRT6Null <- lmer(rt ~ 1 + (1 | participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
MaxRegRT6A <- lmer(rt ~  Touch + (1 + Touch| participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
MaxRegRT6B <- lmer(rt ~ Condition + Touch + (1 + Condition + Touch| participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
#Above doesn't converge, removing random effects of Condition
MaxRegRT6B <- lmer(rt ~ Condition + Touch + (1 +  Touch| participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
MaxRegRT6C <- lmer(rt ~ Condition * Touch + (1 + Condition * Touch| participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
#Above doesn't converge,removing random effects starting with interaction
MaxRegRT6C <- lmer(rt ~ Condition * Touch + (1 + Condition + Touch| participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
#Above doesn't converge,removing random effects starting with interaction
MaxRegRT6C <- lmer(rt ~ Condition * Touch + (1 + Touch| participant), data = SixItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
#Above doesn't converge,removing random effects starting with interaction


###Model Comparisons
anova(MaxRegRT6Null, MaxRegRT6A)
anova(MaxRegRT6A, MaxRegRT6B)
anova(MaxRegRT6A, MaxRegRT6C)
#Best fitting model is C 
summary(MaxRegRT6C)


#Get Data Ready to Graph
SixItemSumRT <- ddply(SixItemsRTEndItems , c( "Condition", "Touch", "participant"), summarise,
                      RT = mean(rt, na.rm = T))

SixItemSumRTOverall <- ddply(SixItemSumRT , c( "Condition", "Touch"), summarise,
                             OverallRT = mean(RT, na.rm = T),
                             SD = sd(RT, na.rm = T),
                             SE = sd(RT, na.rm = T)/ sqrt(length(participant)),
                             NumParticipants = length(participant)
)
#Get Data Ready to Graph
SixItemSumRTW4 <- ddply(SixItemsRTEndItemsW4 , c( "Condition", "Touch", "participant"), summarise,
                        RT = mean(rt, na.rm = T))

SixItemSumRTOverallW4 <- ddply(SixItemSumRTW4 , c( "Condition", "Touch"), summarise,
                               OverallRT = mean(RT, na.rm = T),
                               SD = sd(RT, na.rm = T),
                               SE = sd(RT, na.rm = T)/ sqrt(length(participant)),
                               NumParticipants = length(participant)
)


#Plot 6 item stuff    The Linear fit in the graph just uses subject's averages (and ignores trials per subject - see lmer above & SI for full model
g = ggplot (data = SixItemSumRTOverallW4, aes(x=Touch, y = OverallRT, colour = Condition)) +
  scale_colour_manual(values=c("#009E73", "red")) +
  geom_smooth(data = SixItemSumRT, aes(x=Touch, y = RT, colour = Condition), method = "lm", se=TRUE) +
  geom_point( stat="identity", size = 3) +
  scale_fill_manual(values=c("#009E73", "red")) +
  geom_errorbar(aes(ymin=OverallRT-SE, ymax=OverallRT+SE), width=.01) +
  #scale_x_continuous(limits = c(4,5,6)) +
  coord_cartesian(ylim=c(0, 1500), xlim=c(4,6)) +
  scale_x_continuous(breaks=seq(4, 6, 1))+
  labs(y = "RT (ms)") +
  labs(x = "Touch Number")+
  ggtitle ("6-item lists")
g


##### 8 item Stuff
#Get just the end items (with and without the first item in the 2nd half)
EightItemsRTEndItems<- EightItemsRT [EightItemsRT$Touch == "6" | EightItemsRT$Touch == "7" | EightItemsRT$Touch == "8",]
EightItemsRTEndItemsW5<- EightItemsRT [EightItemsRT$Touch == "5" |EightItemsRT$Touch == "6" | EightItemsRT$Touch == "7" | EightItemsRT$Touch == "8",]


### Mixed effects model (with participant as random effects term) 
rt8 <- lmer(rt ~ Condition * Touch + (1 + Touch| participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
summary(rt8)

### Mixed effects models 
MaxRegRT8Null <- lmer(rt ~ 1 + (1 | participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
MaxRegRT8A <- lmer(rt ~  Condition + (1 + Condition| participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
#Above doesn't converge, removing random effects of condition
MaxRegRT8A <- lmer(rt ~  Condition + (1 | participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
MaxRegRT8B <- lmer(rt ~ Condition + Touch + (1 + Condition + Touch| participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
MaxRegRT8C <- lmer(rt ~ Condition * Touch + (1 + Condition * Touch| participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
#Above has a singular fit, removing random effects starting with interaction
MaxRegRT8C <- lmer(rt ~ Condition * Touch + (1 + Condition + Touch| participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
#Above doesn't converge,removing random effects starting with interaction
MaxRegRT8C <- lmer(rt ~ Condition * Touch + (1 + Touch| participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))

###Model Comparisons
anova(MaxRegRT8Null, MaxRegRT8A)
anova(MaxRegRT8A, MaxRegRT8B)

#Removing condition from random effects for model comparison
MaxRegRT8B1 <- lmer(rt ~ Condition + Touch + (1 +  Condition + Touch| participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
MaxRegRT8B2 <- lmer(rt ~ Condition + Touch + (1 + Touch| participant), data = EightItemsRTEndItems, control = lmerControl(optimizer="bobyqa"))
anova(MaxRegRT8B1, MaxRegRT8B2) # no sign. difference between models with and without the random effects term of Condition

anova(MaxRegRT8B2, MaxRegRT8C)
#Best fitting model is C (w/ interaction)
summary(MaxRegRT8C)

#Get Data Ready to Graph
EightItemSumRT <- ddply(EightItemsRTEndItems , c( "Condition", "Touch", "participant"), summarise,
                        RT = mean(rt, na.rm = T))

EightItemSumRTOverall <- ddply(EightItemSumRT , c( "Condition", "Touch"), summarise,
                               OverallRT = mean(RT, na.rm = T),
                               SD = sd(RT, na.rm = T),
                               SE = sd(RT, na.rm = T)/ sqrt(length(participant)),
                               NumParticipants = length(participant)
)
EightItemSumRTW5 <- ddply(EightItemsRTEndItemsW5 , c( "Condition", "Touch", "participant"), summarise,
                          RT = mean(rt, na.rm = T))

EightItemSumRTOverallW5 <- ddply(EightItemSumRTW5 , c( "Condition", "Touch"), summarise,
                                 OverallRT = mean(RT, na.rm = T),
                                 SD = sd(RT, na.rm = T),
                                 SE = sd(RT, na.rm = T)/ sqrt(length(participant)),
                                 NumParticipants = length(participant)
)


#Plot 8 item stuff    The Linear fit in the graph just uses subject's averages (and ignores trials per subject - see lmer above & SI for full model
g = ggplot (data = EightItemSumRTOverallW5, aes(x=Touch, y = OverallRT, colour = Condition)) +
  scale_colour_manual(values=c("#009E73", "red")) +
  geom_smooth(data = EightItemSumRT, aes(x=Touch, y = RT, colour = Condition), method = "lm", se=TRUE) +
  geom_point( stat="identity", size = 3) +
  scale_fill_manual(values=c("#009E73", "red")) +
  geom_errorbar(aes(ymin=OverallRT-SE, ymax=OverallRT+SE), width=.01) +
  #scale_x_continuous(limits = c(4,5,6)) +
  coord_cartesian(ylim=c(0, 1500), xlim=c(5,8)) +
  scale_x_continuous(breaks=seq(5, 8, 1))+
  labs(y = "RT (ms)") +
  labs(x = "Touch Number")+
  ggtitle ("8-item lists")
g

tab_model(MaxRegRT6C, MaxRegRT8C,  title = "Table S13: Response time by item number", dv.labels = c("6-item Array", "8-item Array"), file = "S13.doc")
