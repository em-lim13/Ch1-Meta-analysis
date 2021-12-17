#Meta-analysis using metafor package for invasive vs non-invasive crab predatory impacts

#### SET UP #### 
setwd("~/Documents/My Documents/CAISN II-1.3/Projects/Invasion Meta-analysis")
dir.data <- "~/Documents/My Documents/CAISN II-1.3/Projects/Invasion Meta-analysis"
datafile <- "CrabMA.Final.csv"

setwd(dir.data)
alldata <- read.csv(datafile, header = TRUE)
spalding <- read.csv("Spalding2007.ProvincesList.csv", header=T)

### LOAD PACKAGES ###
library(metafor)
library(plyr)
library(car)

#### RUN ESCALC TO GENERATE EFFECT SIZES ####

#These can be used to sequentially remove large influential studies
#alldata <- subset(alldata, article!="106") #Fernandes et al 1999
#alldata <- subset(alldata, article!="61") #Griffen et al 2008 
#alldata <- subset(alldata, article!="14") #Tyrell et al. 2006
#alldata <- subset(alldata, article!="55") #Matassa & Trussell 2011

#alldata <- subset(alldata, crab.species!="Carcinus maenas") #remove all green crab experiments


effectsizes <- escalc(measure = 'SMD', n1i = n.wo, n2i = n.w, m1i = x.wo, m2i = x.w, sd1i = sd.wo, sd2i = sd.w, 
                      data=alldata, vtype = 'UB', append=TRUE)


adjusted_effect_size <- c() #loop to reverse sign on all survivorship studies (measure = 2) effect sizes 
                                  #(assuming crabs have a negative effect, then a -ve effect size always 
                                  #indicates a decline in responding species abundance in the presence of crabs) 
                                  
for (i in 1:length(effectsizes$yi)) {
  if (effectsizes[i,]$measure == 2) {
    adjusted_effect_size <- c(adjusted_effect_size, effectsizes[i,]$yi*(-1))
  } else {
    adjusted_effect_size <- c(adjusted_effect_size, effectsizes[i,]$yi*(1))
  }
}
rev_effects <- data.frame(effectsizes, adjusted_effect_size)

#SMD is standardized mean difference (Hedges). 
#vtype options are LS (large sample) and UB (unbiased), I'm not sure which is appropriate although it doesn't 
#seem to change the output? 

#Change intergers in dataset to categorical (factor)
rev_effects$crab.status <- as.factor(rev_effects$crab.status)
rev_effects$crab.sex <- as.factor(rev_effects$crab.sex)
rev_effects$crab.region <- as.factor(rev_effects$crab.region)
rev_effects$prey.status <- as.factor(rev_effects$prey.status)
rev_effects$prey.region <- as.factor(rev_effects$prey.region)
rev_effects$exp.region <- as.factor(rev_effects$exp.region)
rev_effects$lab.field <- as.factor(rev_effects$lab.field)
rev_effects$openess <- as.factor(rev_effects$openess)
rev_effects$measure <- as.factor(rev_effects$measure)
rev_effects$effect <- as.factor(rev_effects$effect)
rev_effects$independent <- as.factor(rev_effects$independent)
rev_effects$article <- as.factor(rev_effects$article)
rev_effects$experiment <- as.factor(rev_effects$experiment)

#Append a column for latitude where all values positive
rev_effects["latitude.pos"] <-  NA
rev_effects$latitude.pos <- abs(rev_effects$latitude)

#Append a column of individual weights 
rev_effects["inv.var"] <- NA
rev_effects$inv.var <- 1/rev_effects$vi

#Remove the influential experiments (n=6)
#rev_effects <- rev_effects[-c(43, 137, 141, 230, 295, 514, 515, 550),]
#rev_effects <- rev_effects[-c(44, 89, 286, 310, 311, 377),]
#rev_effects <- rev_effects[-c(43, 87, 281, 305, 306, 372),]

#Create seperate dataframes for direct and indirect effects 

direct.df <- data.frame(subset(rev_effects, effect=="1"))
indirect.df <- data.frame(subset(rev_effects, effect=="2"))

#### BASIC STATISTICS FOR FINAL DATASET ####

#Do not run every time! 

#Number of experiments from each study and average   #NOTE at some point show this code to Sean - surely there is a faster way to do this? 
table1 <- table(rev_effects$article)
article.data <- as.data.frame(table1)
mean(article.data$Freq)
sd(article.data$Freq)

#Number of non-independent experiments (shared control)
table(rev_effects$independent)

#Range of publication year: 
summary(alldata$year)
hist(alldata$year)

#Number of studies
count(rev_effects$article)

#Range of crab species and stacked plot of status by species
crabspecies <- count(alldata, vars='crab.species')
barplot(crabspecies$freq, names.arg = crabspecies$crab.species, cex.names = 0.4, las=3)

crabspecies.table <- table(rev_effects$crab.status, rev_effects$crab.species)
barplot(crabspecies.table, cex.names = 0.4, las=3)

#Range of crab families
crabfamily <- count(alldata, vars='crab.family')
barplot(crabfamily$freq, names.arg = crabfamily$crab.family, cex.names = 1, las=3)

#Number of invasive (1) vs native (2) crabs
count(alldata$crab.status)

#Origin of crab species 
count(alldata$crab.region)
alldata$crabprov <- as.character(spalding[match(alldata$crab.region, spalding$number), 'province']) #Create a new
  #column in alldata that matches the province number with the name (from Spalding) 
hist(alldata$latitude)

crabprov <- count(alldata, vars='crabprov')
barplot(crabprov$freq, names.arg = crabprov$crabprov, cex.names =1, las =3) 

#Number of prey (specific type = 1)
prey1 <- count(alldata, vars='preytype.1')
barplot(prey1$freq, names.arg = prey1$preytype.1, cex.names=0.8, las=3)

#Number of prey (general type = 3)
prey3 <- count(alldata, vars='preytype.3')
par(mar=c(9.1,4.1,2.1,2.1))
barplot(prey3$freq, names.arg = prey3$preytype.3, cex.names=1, las=2, ylim=c(0,275))

#Number of invasive (1), native prey (2) and unknown (3)
count(alldata$prey.status)

#Origin of prey species 
count(alldata$prey.region)
alldata$preyprov <- as.character(spalding[match(alldata$prey.region, spalding$number), 'province']) #Create a new
#column in alldata that matches the province number with the name (from Spalding) 

preyprov <- count(alldata, vars='preyprov')
barplot(preyprov$freq, names.arg = preyprov$preyprov, cex.names =1, las =3)

#Experimental region (location of study)
count(alldata$exp.region)
alldata$expprov <- as.character(spalding[match(alldata$exp.region, spalding$number), 'province']) #Create a new
#column in alldata that matches the province number with the name (from Spalding) 

expprov <- count(alldata, vars='expprov')
barplot(expprov$freq, names.arg = expprov$expprov, cex.names =1, las =3)

#Number of lab vs field experiments
count(alldata$lab.field)

#How 'open' are the experiments 
open <- count(alldata$openess)
par(mar=c(11.1,4.1,1.1,2.1))
barplot(open$freq, names.arg = c("field: natural experiment", "field: natural enclosure", 
                                 "field: stocked enclosure", "lab: mesocosm", 
                                 "lab: predation experiment"), cex.names =1, las = 2, ylim=c(0,275))

# Experiment duration (days) NOTE: Histogram looks horrible due to huge spread
summary(alldata$duration)

# Field density compared to experimental density NOTE: other/better ways to show this both graphically 
    #and statisitically needed. Will probably have to clean up data as well (too many blanks)
hist(alldata$dens.diff, breaks=25)
plot(alldata$dens.diff, alldata$point)

#If prey mortality was measured (=1), if prey survivorship was measured (=2)
count(alldata$measure)

#If effect is direct (=1), if indirect (=2)
count(alldata$effect)
count(alldata$effect)/(nrow(alldata))

#Histogram of latitudes
summary(rev_effects$latitude)
hist(rev_effects$latitude, breaks=25)

sum(alldata$latitude >23)
sum(alldata$latitude >-23 & alldata$latitude <23)
sum(alldata$latitude <(-23))

#### RANDOM-EFFECTS MODEL FOR OVERALL EFFECT SIZE #### 

no.mods <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = NULL, weighted = TRUE, subset = NULL, 
                   method = "REML") #no moderators, no subsets. 

#alternate method of weighting
no.mods2 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = NULL, weights = inv.var, weighted = FALSE, subset = NULL, 
                   method = "REML") #no moderators, no subsets. 

#NOTE: To make confint(no.mods) work (see the CIs on heterogeneity estimates) you may need to increase the
#upper bound by increasing tau2.max via control to 5000. Only did this when necessary.
#control = list(tau2.max = 5000

summary(no.mods)

#### MIXED/RANDOM EFFECTS MODELS - DIRECT vs. INDIRECT (mods) #####

#Test for difference between direct vs indirect 
effect.qm <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects,
                      mods = ~ factor(effect),
                      weighted = TRUE, level = 95, method = "REML") #without dummy variable (for Qm test)

effect.qt <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects,
                      mods = ~ factor(effect) - 1,
                      weighted = TRUE,  level = 95, method = "REML") #without dummy variable (for effect sizes and CIs)


#Random effects models for DIRECT & INDIRECT (subsets)

direct.rma <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.df, weighted = TRUE, method = "REML", level = 95)

indirect.rma <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.df, weighted = TRUE, method = "REML", level = 95) 

#### MIXED-EFFECTS MODELS - DIRECT INTERACTIONS ONLY (mods = crab status, openess, prey type, latitude) ####

#With intercept (for Qm stats)
mod6.1qm <-rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.df, 
                 mods = ~ factor(crab.status), 
                 weighted = TRUE, level = 95, method = "REML") #DIRECT, CRAB STATUS


mod6.2qm <-rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.df, 
                 mods = ~ factor(openess), 
                 weighted = TRUE, level = 95, method = "REML") #DIRECT, OPENESS


mod6.3qm <-rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.df, 
                 mods = ~ factor(preytype.3), 
                 weighted = TRUE, level = 95, method = "REML") #DIRECT, PREY TYPE


#Without intercept (for Figure 2, effect sizes, and CIs)
mod6.1 <-rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.df, 
                 mods = ~ factor(crab.status) - 1, 
                 weighted = TRUE, level = 95, method = "REML") #DIRECT, CRAB STATUS


mod6.2 <-rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.df, 
                 mods = ~ factor(openess) - 1, 
                 weighted = TRUE, level = 95, method = "REML") #DIRECT, OPENESS


mod6.3 <-rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.df, 
                 mods = ~ factor(preytype.3) - 1, 
                 weighted = TRUE, level = 95, method = "REML") #DIRECT, PREY TYPE

#### MIXED-EFFECT MODELS - INDIRECT INTERACTIONS ONLY (mods = crab status, openess, prey type, latitude) #####

#With intercept (for Qm stats)
mod7.1qm <-rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.df, 
                 mods = ~ factor(crab.status), 
                 weighted = TRUE, level = 95, method = "REML") #INDIRECT, CRAB STATUS


mod7.2qm <-rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.df, 
                 mods = ~ factor(openess), 
                 weighted = TRUE, level = 95, method = "REML") #INDIRECT, OPENESS 

mod7.3qm <-rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.df, 
                 mods = ~ factor(preytype.3), 
                 weighted = TRUE, level = 95, method = "REML") #INDIRECT, PREY TYPE


#Without intercept (for Figure 2, effect sizes, and CIs)
mod7.1 <-rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.df, 
                 mods = ~ factor(crab.status) - 1, 
                 weighted = TRUE, level = 95, method = "REML") #INDIRECT, CRAB STATUS

mod7.2 <-rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.df, 
                 mods = ~ factor(openess) - 1, 
                 weighted = TRUE, level = 95, method = "REML") #INDIRECT, OPENESS 

mod7.3 <-rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.df, 
                 mods = ~ factor(preytype.3) - 1, 
                 weighted = TRUE, level = 95, method = "REML") #INDIRECT, PREY TYPE

##### FIGURE 2 GRAPHING ##### 
#Individual plots done, combine in Pixelmator 

library(forestplot)

#Direct effects
DirectStatus <- coef(summary(mod6.1))
temprow <- matrix(c(rep.int(NA, length(DirectStatus))), nrow=1, ncol = length(DirectStatus))
newrow <- data.frame(temprow)
colnames(newrow) <- colnames(DirectStatus)
DirectStatus <- rbind(DirectStatus, newrow)
DirectOpen <- coef(summary(mod6.2))
temprow <- matrix(c(rep.int(NA, length(DirectOpen))), nrow=1, ncol = length(DirectOpen))
newrow <- data.frame(temprow)
colnames(newrow) <- colnames(DirectOpen)
DirectOpen <- rbind(DirectOpen, newrow)
DirectPrey <- coef(summary(mod6.3))  

DirectColumn <-rbind(DirectStatus, DirectOpen, DirectPrey)
DirectColumn

DirectCategories <- c("Non-native crab consumer", "Native crab consumer", NA, 
                      "Natural experiments","Unstocked enclosures", "Stocked enclosures","Mesocosms","Predation experiment", NA, 
                      "Mobile infauna", "Mobile epifauna", "Primary producers", "Sessile invertebrates", "Vertebrates")

DirectLabels <- cbind(
  c(DirectCategories),
  c("", "", NA, 
    "", "","", "", "", NA, 
    "", "","", "", ""),
  c("(166)", "(422)", NA, 
    "(15)", "(187)","(116)", "(92)", "(178)", NA, 
    "(222)", "(166)","(44)", "(150)", "(6)") 
)
DirectLabels

#Indirect effects
IndirectStatus <- coef(summary(mod7.1))
temprow <- matrix(c(rep.int(NA, length(IndirectStatus))), nrow=1, ncol = length(IndirectStatus))
newrow <- data.frame(temprow)
colnames(newrow) <- colnames(IndirectStatus)
IndirectStatus <- rbind(IndirectStatus, newrow)
IndirectOpen <- coef(summary(mod7.2))
temprow <- matrix(c(rep.int(NA, length(IndirectOpen))), nrow=2, ncol = length(IndirectOpen))
newrow <- data.frame(temprow)
colnames(newrow) <- colnames(IndirectOpen)
IndirectOpen <- rbind(IndirectOpen, newrow)
IndirectPrey <- coef(summary(mod7.3))

IndirectColumn <-rbind(IndirectStatus, IndirectOpen, IndirectPrey)
IndirectColumn

IndirectLabels <- cbind(
  c("", "", NA, 
    "", "", "", "", NA, NA, 
    "", "", "", "", ""),
  c("(75)", "(66)", NA, 
    "(7)", "(16)", "(33)", "(85)", NA, NA, 
    "(6)", "(2)", "(60)", "(67)", "(6)")
)
IndirectLabels


tiff("Figure2-ForestPlot-Direct.tiff", width = 105, height = 150, units="mm",
     pointsize=10, res=600, compression="lzw") #exports figure to working directory with proper formatting 

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 1)))

pushViewport(viewport(layout.pos.col = 1, width = 11, height = 12))
forestplot(DirectLabels, DirectColumn$estimate, DirectColumn$ci.lb, DirectColumn$ci.ub, 
           txt_gp = fpTxtGp(label = list(gpar(fontfamily = "Arial"), cex = 1.0),
                            ticks = gpar(fontfamily = "", cex = 1),
                            xlab = gpar(fontfamily = "Arial", cex = 8)),
           zero = 0, clip = c(-3, 3), xticks = c(-3, -2, -1, 0, 1), graph.pos = 3, lineheight = unit(0.85, "cm"), 
           graphwidth = unit(3.5, "cm"),
           ci.vertices = TRUE, boxsize = 0.15, lwd.zero = 1.5, lwd.xaxis = 1.5, 
           col = fpColors(box = "black", lines = "black", summary = "black", zero = "darkgrey", text = "black", axes = "black", hrz_lines = "black"))
popViewport()
popViewport()

dev.off()


tiff("Figure2-ForestPlot-Indirect.tiff", width = 55, height = 150, units="mm",
     pointsize=10, res=600, compression="lzw") #exports figure to working directory with proper formatting 

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 1)))

pushViewport(viewport(layout.pos.col = 1))
forestplot(IndirectLabels, IndirectColumn$estimate, IndirectColumn$ci.lb, IndirectColumn$ci.ub,
           txt_gp = fpTxtGp(label = list(gpar(fontfamily = "Arial"), cex = 1.0),
                            ticks = gpar(fontfamily = "", cex = 1),
                            xlab = gpar(fontfamily = "Arial", cex = 8)),
           zero = 0, clip = c(-2, 2), xticks = c(-2, -1, 0, 1, 2), graph.pos = 2, lineheight = unit(0.85, "cm"), 
           graphwidth = unit(3.5, "cm"),
           ci.vertices = TRUE, boxsize = 0.15, lwd.zero = 1.5, lwd.xaxis = 1.5,
           col = fpColors(box = "black", lines = "black", summary = "black", zero = "darkgrey", text = "black", axes = "black", hrz_lines = "black"))
popViewport()
popViewport()

dev.off()

#### INTERACTION MODELS #####

#Linear hypothesis and resulting chi-squared value (equivalent to Qm) tell you if comparison is significant ONLY 

#DIRECT
mod14.1 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.df, mods = ~ crab.status:openess -1, 
                   weighted = TRUE, level = 95, method = "REML")

linearHypothesis(mod14.1, c("crab.status1:openess1 - crab.status2:openess1 = 0"))
linearHypothesis(mod14.1, c("crab.status1:openess2 - crab.status2:openess2 = 0"))
linearHypothesis(mod14.1, c("crab.status1:openess3 - crab.status2:openess3 = 0"))
linearHypothesis(mod14.1, c("crab.status1:openess4 - crab.status2:openess4 = 0"))
linearHypothesis(mod14.1, c("crab.status1:openess5 - crab.status2:openess5 = 0"))

mod14.2 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.df, mods = ~ crab.status:preytype.3 -1, 
                   weighted = TRUE, level = 95, method = "REML")

linearHypothesis(mod14.2, c("crab.status1:preytype.3infauna - crab.status2:preytype.3infauna = 0"))
linearHypothesis(mod14.2, c("crab.status1:preytype.3mobile invertebrates - crab.status2:preytype.3mobile invertebrates = 0"))
linearHypothesis(mod14.2, c("crab.status1:preytype.3primary producer - crab.status2:preytype.3primary producer = 0"))
linearHypothesis(mod14.2, c("crab.status1:preytype.3sessile invertebrates - crab.status2:preytype.3sessile invertebrates = 0"))
linearHypothesis(mod14.2, c("crab.status1:preytype.3vertebrates - crab.status2:preytype.3vertebrates = 0"))

mod14.3 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.df, mods = ~ crab.status:latitude.pos -1, 
                   weighted = TRUE, level = 95, method = "REML") 

linearHypothesis(mod14.3, c("crab.status1:latitude.pos - crab.status2:latitude.pos = 0"))

#To look at the change per degree of latitude you need to use the qm model where you haven't forced the intercept. 
#Also how you determine if the effects are statisctically different from zero (Qm)
#The predict fcn can be used to get effect sizes for a range of latitudes.
mod14.3qm <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.df, mods = ~ crab.status:latitude.pos, 
                   weighted = TRUE, level = 95, method = "REML") 

predict.rma(mod14.3qm, newmods = c(0,0), addx = TRUE)  #effect size when both crab and latitude = 0, equivalent to intercept 
predict.rma(mod14.3qm, newmods = c(10,0), addx = TRUE) #effect size when invasive:latitude = 10 
predict.rma(mod14.3qm, newmods = c(20,0), addx = TRUE)
predict.rma(mod14.3qm, newmods = c(30,0), addx = TRUE)
predict.rma(mod14.3qm, newmods = c(40,0), addx = TRUE)
predict.rma(mod14.3qm, newmods = c(50,0), addx = TRUE)

predict.rma(mod14.3qm, newmods = c(0,0), addx = TRUE)  #effect size when both crab and latitude = 0, equivalent to intercept 
predict.rma(mod14.3qm, newmods = c(0,10), addx = TRUE) #effect size when native:latitude = 10 
predict.rma(mod14.3qm, newmods = c(0,20), addx = TRUE)
predict.rma(mod14.3qm, newmods = c(0,30), addx = TRUE)
predict.rma(mod14.3qm, newmods = c(0,40), addx = TRUE)
predict.rma(mod14.3qm, newmods = c(0,50), addx = TRUE)

table(direct.df$openess, direct.df$crab.status)
table(direct.df$preytype.3, direct.df$crab.status)
table(direct.df$crab.status)

#INDIRECT
mod15.1 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.df, mods = ~ crab.status:openess -1, 
                   weighted = TRUE, level = 95, method = "REML")

linearHypothesis(mod15.1, c("crab.status1:openess1 - crab.status2:openess1 = 0"))
linearHypothesis(mod15.1, c("crab.status1:openess2 - crab.status2:openess2 = 0"))
linearHypothesis(mod15.1, c("crab.status1:openess3 - crab.status2:openess3 = 0"))
linearHypothesis(mod15.1, c("crab.status1:openess4 - crab.status2:openess4 = 0"))

mod15.2 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.df, mods = ~ crab.status:preytype.3 -1, 
                   weighted = TRUE, level = 95, method = "REML")

linearHypothesis(mod15.2, c("crab.status1:preytype.3infauna - crab.status2:preytype.3infauna = 0"))
linearHypothesis(mod15.2, c("crab.status1:preytype.3mobile invertebrates - crab.status2:preytype.3mobile invertebrates = 0"))
linearHypothesis(mod15.2, c("crab.status1:preytype.3primary producer - crab.status2:preytype.3primary producer = 0"))
linearHypothesis(mod15.2, c("crab.status1:preytype.3sessile invertebrates - crab.status2:preytype.3sessile invertebrates = 0"))
linearHypothesis(mod15.2, c("crab.status1:preytype.3vertebrates - crab.status2:preytype.3vertebrates = 0"))

mod15.3 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.df, mods = ~ crab.status:latitude.pos -1, 
                   weighted = TRUE, level = 95, method = "REML")

linearHypothesis(mod15.3, c("crab.status1:latitude.pos - crab.status2:latitude.pos = 0"))

#To look at the change per degree of latitude you need to use the qm model where you haven't forced the intercept.
#Also how you determine if the effects are statisctically different from zero (Qm)
#The predict fcn can be used to get effect sizes for a range of latitudes.
mod15.3qm <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.df, mods = ~ crab.status:latitude.pos, 
                     weighted = TRUE, level = 95, method = "REML") 

predict.rma(mod15.3qm, newmods = c(0,0), addx = TRUE)  #effect size when both crab and latitude = 0, equivalent to intercept 
predict.rma(mod15.3qm, newmods = c(10,0), addx = TRUE) #effect size when invasive:latitude = 10 
predict.rma(mod15.3qm, newmods = c(20,0), addx = TRUE)
predict.rma(mod15.3qm, newmods = c(30,0), addx = TRUE)
predict.rma(mod15.3qm, newmods = c(40,0), addx = TRUE)
predict.rma(mod15.3qm, newmods = c(50,0), addx = TRUE)

predict.rma(mod15.3qm, newmods = c(0,0), addx = TRUE)  #effect size when both crab and latitude = 0, equivalent to intercept 
predict.rma(mod15.3qm, newmods = c(0,10), addx = TRUE) #effect size when native:latitude = 10 
predict.rma(mod15.3qm, newmods = c(0,20), addx = TRUE)
predict.rma(mod15.3qm, newmods = c(0,30), addx = TRUE)
predict.rma(mod15.3qm, newmods = c(0,40), addx = TRUE)
predict.rma(mod15.3qm, newmods = c(0,50), addx = TRUE)

table(indirect.df$openess, indirect.df$crab.status)
table(indirect.df$preytype.3, indirect.df$crab.status)
table(indirect.df$crab.status)

##### PUBLICATION BIAS AND NON-INDEPENDENCE TESTING ####

direct.rma <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.df, weighted = TRUE, method = "REML", level = 95, slab = paste(point)) 
indirect.rma <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.df, weighted = TRUE, method = "REML", level = 95, slab = paste(point)) 

#Cumulative MA 
cumulative.direct <- cumul(direct.rma, order=order(direct.df$vi))
cumulative.indirect <- cumul(indirect.rma, order=order(indirect.df$vi))

#Forest plots
forest.cumul.rma(cumulative.direct)
forest.cumul.rma(cumulative.indirect)

#Create subsets of direct and indirect dataframes with only "large studies" as determined by eye from forest plots and repeat 
direct.ls.df <- subset(direct.df, direct.df$vi<0.70)
indirect.ls.df <- subset(indirect.df, indirect.df$vi<0.85)

direct.ls.rma <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, weighted = TRUE, method = "REML", slab = paste(point)) 
indirect.ls.rma <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, weighted = TRUE, method = "REML", slab = paste(point)) 

#Cumulative meta-analysis 
cumulative.direct.ls <- cumul(direct.ls.rma, order=order(direct.ls.df$vi))
cumulative.indirect.ls <- cumul(indirect.ls.rma, order=order(indirect.ls.df$vi))

forest.cumul.rma(cumulative.direct.ls)
forest.cumul.rma(cumulative.indirect.ls)

#Funnel and QQ plots 

funnel.rma(direct.rma)
funnel.rma(indirect.rma)

funnel.rma(direct.ls.rma)
funnel.rma(indirect.ls.rma)

qqnorm(direct.rma)
qqnorm(indirect.rma)

qqnorm(direct.ls.rma)
qqnorm(indirect.ls.rma)

#Rank correlation and Reg Tests

ranktest.rma(direct.rma)
ranktest.rma(indirect.rma)

ranktest.rma(direct.ls.rma)
ranktest.rma(indirect.ls.rma)

regtest.rma(direct.rma)
regtest.rma(indirect.rma)

regtest.rma(direct.ls.rma)
regtest.rma(indirect.ls.rma)

#Trim & Fill 

trimfill(direct.rma)
trimfill(indirect.rma)

trimfill(direct.ls.rma)
trimfill(indirect.ls.rma)

#Fail Safe Numbers 

fsn(yi=adjusted_effect_size, vi=vi, data=direct.df, type="Rosenthal") #file drawer analysis (missing studies have an effect size = 0)
fsn(yi=adjusted_effect_size, vi=vi, data=direct.df, type="Orwin") #missing studies have an effect size > 0 

fsn(yi=adjusted_effect_size, vi=vi, data=indirect.df, type="Rosenthal") #file drawer analysis (missing studies have an effect size = 0)
fsn(yi=adjusted_effect_size, vi=vi, data=indirect.df, type="Orwin") #missing studies have an effect size > 0 

fsn(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, type="Rosenthal") #file drawer analysis (missing studies have an effect size = 0)
fsn(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, type="Orwin") #missing studies have an effect size > 0 

fsn(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, type="Rosenthal") #file drawer analysis (missing studies have an effect size = 0)
fsn(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, type="Orwin") #missing studies have an effect size > 0 

#Leave one out for all studies (change rev_effects to direct.df or indirect.df throughout)

study_factors <- levels(indirect.df$article)
rma_results <- data.frame()
for (study_number in study_factors) {
  all_but_one_study_rev_effects <- indirect.df[indirect.df$article!=study_number,]
  results <- rma.uni(yi=adjusted_effect_size, vi=vi, data=all_but_one_study_rev_effects, mods = NULL, weighted = TRUE, subset = NULL, method = "REML")
  new_data_frame_row <- data.frame(removed_study_id=study_number, coef(summary(results)))
  rownames(new_data_frame_row) <- c(study_number)
  rma_results <- rbind(rma_results, new_data_frame_row)
}

forest.default(x = rma_results$estimate, ci.lb = rma_results$ci.lb, ci.ub = rma_results$ci.ub)

#Leave-one-out for all experiments

l1o.d <- leave1out(direct.rma)
forest.default(x = l1o.d$estimate, ci.lb = l1o.d$ci.lb, ci.ub = l1o.d$ci.ub)

l1o.i <- leave1out(indirect.rma)
forest.default(x = l1o.i$estimate, ci.lb = l1o.i$ci.lb, ci.ub = l1o.i$ci.ub)

##### FUNNEL PLOTS - SUPPLEMENTARY FIGURE #####

# plots need "direct" "indirect" labels! 

tiff("FunnelPlots-SuppFigX.tiff", width=4.5, height=6, units="in",
     pointsize=10, res=600, compression="lzw") #exports figure to working directory with proper formatting 

par(mfrow=c(2,1), oma=c(4, 4, 2, 2), mar=c(1, 1, 1, 1), 
    tck=0.03, mgp=c(3,0.4,0), family="sans", lty = 1, col = "black") 

funnel.rma(direct.rma, cex = direct.df$inv.var/4, level = c(95),  pch = 21, bg = 240, xlim = c(-6,4), ylim = c(0, 2.0), digits = 1, 
           back = "white")
funnel.rma(indirect.rma, cex = indirect.df$inv.var/4, level = c(95),  pch = 21, bg = 240, xlim = c(-4,6), ylim = c(0, 1.5), digits = 1, 
           back = "white")

mtext("Standardized mean difference (g)", side = 1, outer = TRUE, cex = 1, line = 2.2)
mtext("Standard error", side = 2, outer = TRUE, cex = 1, line = 2.2)

dev.off()

##### CUMULATIVE FOREST PLOTS - SUPPLEMENTARY FIGURE ####

#Without point annotation (for FIGURE)
direct.df["vi.round"] <- (round(direct.df$vi, 2))
indirect.df["vi.round"] <- (round(indirect.df$vi, 2))

direct.rma <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.df, weighted = TRUE, method = "REML", level = 95, slab = paste(direct.df$vi.round))
indirect.rma <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.df, weighted = TRUE, method = "REML", level = 95, slab = paste(indirect.df$vi.round)) 

#Cumulative MA 
cumulative.direct <- cumul(direct.rma, order=order(direct.df$vi))
cumulative.indirect <- cumul(indirect.rma, order=order(indirect.df$vi))

#Create dataframes of estimates and CI's so you can use "forestplot" pkg

direct.est <- cumulative.direct$estimate
direct.cilb <- cumulative.direct$ci.lb
direct.ciub <- cumulative.direct$ci.ub
cumulative.direct.df <- cbind.data.frame(direct.est,direct.cilb,direct.ciub)
cumulative.direct.df["blank"] <- NA

indirect.est <- cumulative.indirect$estimate
indirect.cilb <- cumulative.indirect$ci.lb
indirect.ciub <- cumulative.indirect$ci.ub
cumulative.indirect.df <- cbind.data.frame(indirect.est,indirect.cilb,indirect.ciub)
cumulative.indirect.df["blank"] <- NA

library(forestplot)

tiff("CumulativeForestPlots-SuppFigX.tiff", width=6, height=6, units="in",
     pointsize=10, res=600, compression="lzw") #exports figure to working directory with proper formatting 

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
grid.text("Standardized mean difference (g)", 0.5, 0.02, gp=gpar(cex = 1.0))

pushViewport(viewport(layout.pos.col = 1))
forestplot(cumulative.direct.df$blank, cumulative.direct.df$direct.est, 
                cumulative.direct.df$direct.cilb, cumulative.direct.df$direct.ciub, 
           zero = 0, title = "Direct effects", lineheight = "auto", graphwidth = unit(7, "cm"), 
           hrzl_lines = list("340" = gpar(lwd = 0.8, col = "red", lty = 2, columns = 2)),
           ci.vertices = TRUE, boxsize = 0.25, lwd.zero = 2.5, lwd.xaxis = 1.5, cex = 1,
           col = fpColors(box = "black", lines = "darkgray", summary = "black", zero = "black", text = "black", axes = "black", hrz_lines = "black"))
popViewport()

pushViewport(viewport(layout.pos.col = 2))
forestplot(cumulative.indirect.df$blank, cumulative.indirect.df$indirect.est, 
                cumulative.indirect.df$indirect.cilb, cumulative.indirect.df$indirect.ciub,
           zero = 0, title = "Indirect effects", lineheight = "auto", graphwidth = unit(7, "cm"), graph.pos = 1,
           hrzl_lines = list("106" = gpar(lwd = 0.8, col = "red", lty = 2, columns = 1)),
           ci.vertices = TRUE, boxsize = 0.25, lwd.zero = 2.5, lwd.xaxis = 1.5, cex = 1,
           col = fpColors(box = "black", lines = "darkgray", summary = "black", zero = "black", text = "black", axes = "black", hrz_lines = "black"))
popViewport()
popViewport()

dev.off()

#### RERUN WITH LS DATASETS (DIRECT) ####

#Are indirect and direct still significantly different? 
rev_effects2 <- rbind(direct.ls.df, indirect.ls.df) 

effect.qm2 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects2,
                     mods = ~ factor(effect),
                     weighted = TRUE, level = 95, method = "REML") #without dummy variable (for Qm test)


effect.qm <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects2,
                      mods = ~ factor(effect) -1,
                      weighted = TRUE, level = 95, method = "REML") #without dummy variable (for Qm test)

#With intercept (for Qm stats)
mod6.1qm.ls <-rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, 
                   mods = ~ factor(crab.status), 
                   weighted = TRUE, level = 95, method = "REML") #DIRECT, CRAB STATUS


mod6.2qm.ls <-rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, 
                   mods = ~ factor(openess), 
                   weighted = TRUE, level = 95, method = "REML") #DIRECT, OPENESS


mod6.3qm.ls <-rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, 
                   mods = ~ factor(preytype.3), 
                   weighted = TRUE, level = 95, method = "REML") #DIRECT, PREY TYPE


#Without intercept (for Figure 2, effect sizes, and CIs)
mod6.1.ls <-rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, 
                    mods = ~ factor(crab.status) -1, 
                    weighted = TRUE, level = 95, method = "REML") #DIRECT, CRAB STATUS


mod6.2.ls <-rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, 
                    mods = ~ factor(openess) -1, 
                    weighted = TRUE, level = 95, method = "REML") #DIRECT, OPENESS


mod6.3.ls <-rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, 
                    mods = ~ factor(preytype.3) -1, 
                    weighted = TRUE, level = 95, method = "REML") #DIRECT, PREY TYPE


#Subsetting by each moderator to see if there is a difference between INVASIVE and NATIVE

#QM STATS for INVASIVE v. NATIVE (mod) for DIRECT and EXP DESIGN subsets:
mod8.1.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status), 
                  weighted = TRUE, subset = (openess=="1"), level = 95, method = "REML") 

mod8.2.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status),
                  weighted = TRUE, subset = (openess=="2"), level = 95, method = "REML") 

mod8.3.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status), 
                  weighted = TRUE, subset = (openess=="3"), level = 95, method = "REML") 

mod8.4.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status), 
                  weighted = TRUE, subset = (openess=="4"), level = 95, method = "REML") 

mod8.5.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status), 
                  weighted = TRUE, subset = (openess=="5"), level = 95, method = "REML") 

#Put -1 back in to get effect sizes and errors relative to a dummy variable 
mod9.1.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status) -1, 
                     weighted = TRUE, subset = (openess=="1"), level = 95, method = "REML") 

mod9.2.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status) -1,
                     weighted = TRUE, subset = (openess=="2"), level = 95, method = "REML") 

mod9.3.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status) -1, 
                     weighted = TRUE, subset = (openess=="3"), level = 95, method = "REML") 

mod9.4.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status) -1, 
                     weighted = TRUE, subset = (openess=="4"), level = 95, method = "REML") 

mod9.5.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status) -1, 
                     weighted = TRUE, subset = (openess=="5"), level = 95, method = "REML") 

#QM STATS for INVASIVE v. NATIVE (mod) for DIRECT and PREY TYPE subsets:

mod10.1.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status), 
                   weighted = TRUE, subset = (preytype.3=="infauna"), level = 95, method = "REML") 

mod10.2.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status), 
                   weighted = TRUE, subset = (preytype.3=="mobile invertebrates"),level = 95, method = "REML") 

mod10.3.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status), 
                   weighted = TRUE, subset = (preytype.3=="sessile invertebrates"), level = 95, method = "REML") 

mod10.4.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status), 
                   weighted = TRUE, subset = (preytype.3=="primary producer"), level = 95, method = "REML") 

mod10.5.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status), 
                   weighted = TRUE, subset = (preytype.3=="vertebrates"), level = 95, method = "REML")

# Repeat with dummy variable
mod11.1.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status) -1, 
                      weighted = TRUE, subset = (preytype.3=="infauna"), level = 95, method = "REML") 

mod11.2.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status) -1, 
                      weighted = TRUE, subset = (preytype.3=="mobile invertebrates"),level = 95, method = "REML") 

mod11.3.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status) -1, 
                      weighted = TRUE, subset = (preytype.3=="sessile invertebrates"), level = 95, method = "REML") 

mod11.4.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status) -1, 
                      weighted = TRUE, subset = (preytype.3=="primary producer"), level = 95, method = "REML") 

mod11.5.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=direct.ls.df, mods = ~ factor(crab.status) -1, 
                      weighted = TRUE, subset = (preytype.3=="vertebrates"), level = 95, method = "REML")


#Mixed effect model with interaction for effect of LATITUDE:CRAB STATUS on DIRECT effects

#NEED TO CREATE A NEW DF WITHOUT CAT.5 OPENESS EXPERIMENTS
latitude.d.ls <- subset(direct.ls.df, openess!="5")

mod12.1qm.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=latitude.d.ls, mods = ~ factor(crab.status):latitude.pos, 
                     weighted = TRUE, level = 95, method = "REML")

mod12.1.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=latitude.d.ls, mods = ~ factor(crab.status):latitude.pos -1, 
                        weighted = TRUE, level = 95, method = "REML") #without intercept for effect sizes and CIs 

#### RERUN WITH LS DATASETS (INDIRECT) ####

#With intercept (for Qm stats)
mod7.1qm.ls <-rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, 
                   mods = ~ factor(crab.status), 
                   weighted = TRUE, level = 95, method = "REML") #INDIRECT, CRAB STATUS

mod7.2qm.ls <-rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, 
                   mods = ~ factor(openess), 
                   weighted = TRUE, level = 95, method = "REML") #INDIRECT, OPENESS 

mod7.3qm.ls <-rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, 
                   mods = ~ factor(preytype.3), 
                   weighted = TRUE, level = 95, method = "REML") #INDIRECT, PREY TYPE


#Without intercept (for Figure 2, effect sizes, and CIs)
mod7.1.ls <-rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, 
                 mods = ~ factor(crab.status) - 1, 
                 weighted = TRUE, level = 95, method = "REML") #INDIRECT, CRAB STATUS

mod7.2.ls <-rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, 
                 mods = ~ factor(openess) - 1, 
                 weighted = TRUE, level = 95, method = "REML") #INDIRECT, OPENESS 

mod7.3.ls <-rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, 
                 mods = ~ factor(preytype.3) - 1, 
                 weighted = TRUE, level = 95, method = "REML") #INDIRECT, PREY TYPE


#Subsetting by each moderator to see if there is a difference between INVASIVE and NATIVE

#QM STATS for INVASIVE v. NATIVE (mod) for INDIRECT and EXP DESIGN subsets:
mod8.6.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status), 
                  weighted = TRUE, subset = (openess=="1"), level=95, method = "REML") 

mod8.7.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status), 
                  weighted = TRUE, subset = (openess=="2"),level=95, method = "REML") 

mod8.8.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status), 
                  weighted = TRUE, subset = (openess=="3"),level=95, method = "REML") 

mod8.9.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status), 
                  weighted = TRUE, subset = (openess=="4"), level=95,method = "REML") 

#Put -1 back in to get effect sizes and errors relative to a dummy variable 
mod9.6.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status) -1, 
                  weighted = TRUE, subset = (openess=="1"),level=95, method = "REML") 

mod9.7.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status) -1, 
                  weighted = TRUE, subset = (openess=="2"),level=95, method = "REML") 

mod9.8.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status) -1, 
                  weighted = TRUE, subset = (openess=="3"),level=95, method = "REML") 

mod8.9.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status) -1, 
                  weighted = TRUE, subset = (openess=="4"),level=95, method = "REML") 


#QM STATS for INVASIVE v. NATIVE (mod) for INDIRECT and PREY TYPE subsets:
mod10.6.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status), 
                   weighted = TRUE, subset = (preytype.3=="infauna"), level = 95, method = "REML") 

mod10.7.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status), 
                   weighted = TRUE, subset = ( preytype.3=="mobile invertebrates"),level = 95, method = "REML") 

mod10.8.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status), 
                   weighted = TRUE, subset = ( preytype.3=="sessile invertebrates"),level = 95, method = "REML") 

mod10.9.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status), 
                   weighted = TRUE, subset = (preytype.3=="primary producer"),level = 95, method = "REML") 

mod10.10.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status), 
                    weighted = TRUE, subset = (preytype.3=="vertebrates"),level = 95, method = "REML")

# Repeat with dummy variable
mod11.6.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status) -1, 
                   weighted = TRUE, subset = (preytype.3=="infauna"),level = 95, method = "REML") 

mod11.7.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status) -1, 
                   weighted = TRUE, subset = (preytype.3=="mobile invertebrates"),level = 95, method = "REML") 

mod11.8.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status) -1, 
                   weighted = TRUE, subset = (preytype.3=="sessile invertebrates"),level = 95, method = "REML") 

mod11.9.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status) -1, 
                   weighted = TRUE, subset = (preytype.3=="primary producer"),level = 95, method = "REML") 

mod11.10.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=indirect.ls.df, mods = ~ factor(crab.status) -1, 
                    weighted = TRUE, subset = (preytype.3=="vertebrates"),level = 95, method = "REML")


#Mixed effect model with interaction for effect of LATITUDE:CRAB STATUS on INDIRECT effects

latitude.i.ls <- subset(indirect.ls.df)

mod13.1qm.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=latitude.i.ls, mods = ~ factor(crab.status):latitude.pos, 
                     weighted = TRUE, level = 95, method = "REML")

mod13.1.ls <- rma.uni(yi=adjusted_effect_size, vi=vi, data=latitude.i.ls, mods = ~ factor(crab.status):latitude.pos -1, 
                      weighted = TRUE, level = 95, method = "REML") #without intercept for effect sizes and CIs 

#### EARLY WORK ######
#(1) SINGLE MODERATORS ONLY - NO SUBSETS 

mod1.1 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = cbind(crab.species), weighted = TRUE, subset = NULL, 
                  method = "REML", control = list(tau2.max = 5000), intercept = FALSE) #mod: crab species only, no subsets

mod1.2 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = cbind(crab.family), weighted = TRUE, subset = NULL, 
                  method = "REML")

mod1.3 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = cbind(crab.status), weighted = TRUE, subset = NULL, 
                  method = "REML")

mod1.4 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = cbind(crab.region), weighted = TRUE, subset = NULL, 
                  method = "REML")

mod1.5 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = cbind(preytype.1), weighted = TRUE, subset = NULL, 
                  method = "REML")

mod1.6 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = cbind(preytype.2), weighted = TRUE, subset = NULL, 
                  method = "REML")

mod1.8 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = cbind(prey.region), weighted = TRUE, subset = NULL, 
                  method = "REML")

mod1.9 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = cbind(exp.region), weighted = TRUE, subset = NULL, 
                  method = "REML")

mod1.10 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = cbind(lab.field), weighted = TRUE, subset = NULL, 
                   method = "REML")

mod1.11 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = cbind(openess), weighted = TRUE, subset = NULL, 
                   method = "REML")

mod1.12 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = cbind(duration), weighted = TRUE, subset = NULL, 
                   method = "REML")

mod1.13 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = cbind(dens.diff), weighted = TRUE, subset = NULL, 
                   method = "REML")

mod1.14 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = cbind(measure), weighted = TRUE, subset = NULL, 
                   method = "REML") 

mod1.15 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = cbind(effect), weighted = TRUE, subset = NULL, 
                   method = "REML")

mod1.16 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = cbind(prey.status), weighted = TRUE, subset = NULL, 
                   method = "REML")

# (2) SUB-SETTING BY FACTOR 

mod2.1 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(crab.species) - 1, 
                  weighted = TRUE, subset = NULL, method = "REML")

mod2.2 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(crab.family) - 1, 
                  weighted = TRUE, subset = NULL, method = "REML")

mod2.3 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(crab.status) - 1, 
                  weighted = TRUE, subset = NULL, method = "REML")

mod2.4 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(preytype.1) - 1, 
                  weighted = TRUE, subset = NULL, method = "REML")

mod2.5 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(preytype.2) - 1, 
                  weighted = TRUE, subset = NULL, method = "REML")

mod2.6 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(prey.status) - 1, 
                  weighted = TRUE, subset = NULL, method = "REML")

mod2.7 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(prey.region) - 1, 
                  weighted = TRUE, subset = NULL, method = "REML")

mod2.8 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(exp.region) - 1, 
                  weighted = TRUE, subset = NULL, method = "REML")

mod2.9 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(lab.field) - 1, 
                  weighted = TRUE, subset = NULL, method = "REML")

mod2.10 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(openess) - 1, 
                   weighted = TRUE, subset = NULL, method = "REML")

mod2.11 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(effect) - 1, 
                   weighted = TRUE, subset = NULL, method = "REML")

# (3) MULTIPLE SUBSETS 


#1a. Do invasive crabs have greater impact on mortality than native crabs (direct, indirect, combined)?

mod2.11 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(effect) - 1, 
                   weighted = TRUE, subset = NULL, method = "REML") #both effect types, mod = crab type

mod3.1 <-rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(effect) - 1, 
                 weighted = TRUE, subset = (crab.status=="1"), method = "REML") #invasive only, mod = effect

mod3.2 <-rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(effect) - 1, 
                 weighted = TRUE, subset = (crab.status=="2"), method = "REML") #native only, mod = effect

#2. Are certain prey types more vulnerable to predation by invasive species than they are to native species (direct only)?

#For prey type 2 (collapsed groups)
mod3.3 <-rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(preytype.2) - 1, 
                 weighted = TRUE, subset = (crab.status=="1" & effect=="1"), method = "REML") #invasive and direct effects only, mod = prey type 2

mod3.4 <-rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(preytype.2) - 1, 
                 weighted = TRUE, subset = (crab.status=="2" & effect=="1"), method = "REML") #native and direct effects only, mod = prey type 2

#For prey type 1 (detailed groups)

mod3.15 <-rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(preytype.1) - 1, 
                  weighted = TRUE, subset = (crab.status=="1" & effect=="1"), method = "REML") #invasive and direct effects only, mod = prey type 2

mod3.16 <-rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(preytype.1) - 1, 
                  weighted = TRUE, subset = (crab.status=="2" & effect=="1"), method = "REML") #native and direct effects only, mod = prey type 2

#3. Are certain invasive species more detrimental than others? 

mod3.5 <-rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(crab.species) -1,
                 weighted = TRUE, subset = (crab.status=="1"), method = "REML") #invasive only


#4. What is the effect of experiment type (openess) on the first two questions

mod2.10 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(openess) - 1, 
                   weighted = TRUE, subset = NULL, method = "REML") #all species, mod = openess

mod3.6 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(openess) - 1, 
                  weighted = TRUE, subset = (crab.status=="1"), method = "REML") #

mod3.7 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(openess) - 1, 
                  weighted = TRUE, subset = (crab.status=="2"), method = "REML")

mod3.8 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(openess) - 1, 
                  weighted = TRUE, subset = (crab.status=="1" & effect=="1"), method = "REML")

mod3.9 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(openess) - 1, 
                  weighted = TRUE, subset = (crab.status=="2" & effect=="1"), method = "REML")

mod3.10 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(preytype.2) - 1, 
                   weighted = TRUE, subset = (crab.status=="1" & effect=="1" & openess=="4"), method = "REML")

mod3.11 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(preytype.2) - 1, 
                   weighted = TRUE, subset = (crab.status=="2" & effect=="1" & openess=="4"), method = "REML")

mod3.12 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(preytype.2) - 1, 
                   weighted = TRUE, subset = (crab.status=="1" & effect=="1" & openess=="3"), method = "REML")

mod3.13 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(preytype.2) - 1, 
                   weighted = TRUE, subset = (crab.status=="2" & effect=="1" & openess=="3"), method = "REML")

#Are particular regions more suseptible to the impacts of invasive species? 

mod3.14 <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, mods = ~ factor(exp.region) - 1, 
                   weighted = TRUE, subset = (crab.status=="1"), method = "REML")


#summary(ALL MODELS)
#AIC (ALL MODELS)
#confint(ALL MODELS)

# (4 & 5) MULTIPLE MODERATOR MODELS - AIC TESTING MODELS (TABLE 1) 

#Instructions from: http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti?s[]=aic

rma.glmulti <- function(formula, data, ...) {
  rma(as.formula(paste(deparse(formula))), vi, data=data, method = "ML", ...)
}

#Direct effects only 

direct.effects <- glmulti(adjusted_effect_size ~ crab.status + openess + preytype.3 + latitude, data=rev_effects, 
                          weighted = TRUE, subset = (rev_effects$effect=="1"), level=2, fitfunction=rma.glmulti, 
                          maxsize = 1, confsetsize = 1000000,
                          crit = "aic")

#Level = 1 means its additive only (no interactions), Level = 2 includes interactions, but you have to run it for 
#ALL interactions (can't specify which ones you want). Therefore had to ignore all but the 16 models that were 
#chosen a priori and recalculate deltaAIC and weights myself. 
#Confsetsize = 6 means only 6 of best models to be included in the consensus object. If NA, all available models 
#are included.
#Does setting confsetsize mean that in model averaging only those models are used -> results change when 
#confsetsize is changed so probably, yes? If so, where to set confsetsize? 

plot(direct.effects) #IC profile with level set at 2 ICs difference 

aictable1 <- weightable(direct.effects) #AIC table with weights
aictable1

summary(direct.effects@objects[[1]]) #summary of best model #not sure how to interpret?

plot(direct.effects, type = "s")

#Indirect effects only 
#NOTE: Some redundant predictors are being dropped because the combinations don't exist in 
#the data, becuase the sample size is so much smaller. The warning comes from the meta-analysis package 
#I am using (metafor) not glmulti. See http://www.metafor-project.org/doku.php/updates for more info

indirect.effects <- glmulti(adjusted_effect_size ~ crab.status + openess + preytype.3 + latitude,
                            data=rev_effects, 
                            maxsize = 2, minsize = 1, confsetsize = 30000000,
                            weighted = TRUE, subset = (rev_effects$effect=="2"), level=2, fitfunction=rma.glmulti, 
                            crit = "aic")

plot(indirect.effects) #IC profile with level set at 2 ICs difference 

aictable2 <- weightable(indirect.effects) #AIC table with weights
aictable2

summary(indirect.effects@objects[[1]]) #summary of best model. 

plot(indirect.effects, type = "s")




#Model Averaging (to do this I needed to make the glmutli package handle rma objects again)
#http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti

setOldClass("rma.uni")
setMethod('getfit','rma.uni',function(object,...) {
  if (object$knha==FALSE) {
    cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=100000)
  } else {
    cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=object$k-object$p)
  }
})

round(coef.glmulti(direct.effects), 3) #just rounds decimal places to 3 in table output
round(coef.glmulti(indirect.effects), 3)

#COMBINED MODEL

#Purpose is to combine models used above to run a single MA with required interactions and moderators, rather than running many seperate models
#Indirect and direct would still be seperate
#Moderators included are: crab.status, prey.type3, openess, and latitude:crab.status

direct.combinedmod <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, 
                              mods = ~ openess + crab.status + preytype.3 + crab.status:latitude, 
                              weighted = TRUE, subset = (effect=="1"), method = "REML")

#Interpretting output: 
# LINE 1: Intercept is est. avg. effect (EAE) when openness = 1 for an invasive crab with infaunal prey. 
#LINES 2-5: Openess 2-5 are how much higher/lower the EAE is for each level, when crab is invasive and prey are infauna
#LINE 6: How much higher/lower the EAE is for native crabs when openess = 1 and prey are infaunal
#LINES 7-10: How much higher/lower EAE is for each prey type, when status is held at invasive and openess = 1 
#LINES 11-12: how much higher/lower the contrast between levels of crab status are for all latitudes


direct.combinedmod.woint <- rma.uni(yi=adjusted_effect_size, vi=vi, data=rev_effects, 
                                    mods = ~ openess + crab.status + preytype.3 + crab.status:latitude -1, 
                                    weighted = TRUE, subset = (effect=="1"), method = "REML")