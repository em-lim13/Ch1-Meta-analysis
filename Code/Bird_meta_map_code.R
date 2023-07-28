# This script is to reproduce the analyses from Grames et al. (2022?). The effect 
# of insect food availability on songbird reproductive success and chick body 
# condition: evidence from a systematic review and meta-analysis. Ecology Letters.

# All data needed to run the script are contained in the same directory.

# Functions --------------------------------------------------------------------

#' Calculate I2, the proportion of heterogeneity attributed to true variance in 
#' effects as opposed to sampling error
#' @param mod a model object output from metafor::rma (and related functions)
calc_i2 <- function(mod){
  W <- diag(1/mod$vi)
  X <- model.matrix(mod)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2 <- 100 * sum(mod$sigma2) / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P)))
  return(I2)
}

#' Determine the size of squares for plotting effect sizes
#' @param mod a model object output from metafor::rma (and related functions)
make_squares <- function(mod){
  wi <- weights(mod)
  psize <- wi/sum(wi, na.rm = TRUE)
  rng <- max(psize, na.rm = TRUE) - min(psize, na.rm = TRUE)
  k <- length(wi)
  if (rng <= .Machine$double.eps^0.5) {
    psize <- rep(1, k)
  }else {
    psize <- (psize - min(psize, na.rm = TRUE))/rng
    psize <- (psize * 1) + 0.5
  }
  if (all(is.na(psize))){
    psize <- rep(1, k)
  } 
  return(psize)
}

#' Generate forest plots showing effect sizes, labels, and confidence intervals
#' along with heterogeneity test results
#' @param mod a model object output from metafor::rma (and related functions)
#' @param labels a character vector of length equal to number of studies in mod with appropriate labels
plot_forest <- function(mod, labels){
  par(xpd=F)
  xvals <- seq(-5, 5, 0.5)
  upci <- max(max(mod$yi + 1.96*sqrt(mod$vi), na.rm = T), mod$ci.ub)
  lowci <- min(min((mod$yi - 1.96*sqrt(mod$vi)), na.rm=T), mod$lb)
  xlims <- c(xvals[max(which(xvals<lowci))], xvals[min(which(xvals>upci))])
  ylims <- c(1, nrow(tmpdat))
  
  metafor::forest(mod, slab=labels, 
                  cex.lab=1, xlab="Effect")
  
  text(x=xlims[1], y=length(weights(mod))+c(3), cex=1, adj=0, 
       bquote(paste("Q = ", .(formatC(mod$QE, digits=2, format="f")),
                    "; df = ", .(formatC(mod$k - mod$p, digits=2, format="f")),
                    "; p = ", .(formatC(mod$QEp, digits=2, format="f")))))

  I2 <- calc_i2(mod)
  
  text(x=xlims[1], y=length(weights(mod))+c(2), 
       adj=0, bquote(paste(I^2, " = ", .(formatC(I2, digits=2, format="f")), "; ", tau^2, " = ", .(formatC(mod$tau2, digits=2, format="f")))),
       cex=1)
}

# Workspace setup --------------------------------------------------------------

# Packages required in addition to loaded libraries:
  # metafor -- meta-analyses
  # ape -- phylogenetic tree
  # phylosignal -- plotting phylogenetic tree
  # sp -- for plotting study locations
  # rgdal -- loaded by sp
  # vegan -- for estimating accumulation of measures
  # corpcor -- matrix manipulation
  # Matrix -- positive definite matrix calculation
  # readxl -- for data loading; data are in xlsx to preserve citation formatting

library(dplyr) # for %>%
library(metaSEM) # for %&%
library(visNetwork) # optional, for interactive network plot
library(maptools) # for wrld_simpl map
library(ggplot2) # for heatmaps

# Define our color palettes
pal <- c("#464E29", "#c09129", "#ad3429", 
         "#6b140b", "#b65d35", "#656043",
         "#164057", "#9c8c73", "#3d4753")
pal2 <- c("#fedbaf", "#fa983a", 
          "#f96815", "#c05a1d", "#864c24", 
          "#833826", "#812e27", "#7f2427", 
          "#732729", "#66292b")         

# Read in data -----------------------------------------------------------------
metadata <- data.frame(readxl::read_xlsx("./study_characteristics.xlsx"))
effects <- read.csv("./effect_sizes.csv")
measures <- read.csv("./effect-measurements.csv")
traits <- read.delim("https://figshare.com/ndownloader/files/5631081", sep="\t")[1:9993,] # drop a few NA rows

# Were studies included in the meta-analytic structural equation model?
# Some articles are included in the narrative review but we could not calculate
# an effect size or the measures they used were not commonly used
metadata$included <- metadata$Article_ID%in%effects$Article_ID

# Bird species studied ---------------------------------------------------------
# Read in tree from VertLife for plotting
phylo <- ape::read.nexus("./output.nex")
tree <- phylo[[1]]
tree$tip.label<-gsub("_"," ",tree$tip.label)

counts <- table(metadata$Species.name)

phylo_count <- phylobase::phylo4d(tree,  
                               data.frame(count=as.numeric(counts[match(tree$tip.label, names(counts))])))

# Color indicates log(number) of studies of each species
# Figure 2b
par(xpd=T)
phylosignal::gridplot(
  phylo_count, show.data.axis=T,
  tree.xlim=c(-1000, 1000), 
  tree.ladderize=T, tree.open.crown=T,
  trait.bg.col="#00000000",
  cell.col = colorRampPalette(pal2)(256),
  tree.ratio=6, show.box=F,
  tree.type = "fan", tip.font=1,
  show.trait = F,
  tip.cex = .6, tree.open.angle=90)

# Figure 2a
par(xpd=F, pty="s", las=1)
families <- table(metadata$included, metadata$Family)
barplot(families[2:1,order(colSums(families))], xlim=c(0, 150), 
        horiz=T, col=paste("#164057", c("ff", "88"), sep=""),
          border="#16405700", add=F, cex.names=1)

# Note: for publication, these plots have been post-hoc combined in Inkscape

# Habitats / study locations ---------------------------------------------------
# Many studies have more than one habitat type, so can be represented more than once
for(i in 1:nrow(metadata)){
  habitats <- trimws(unlist(strsplit(metadata$Habitat.type[i], ";")))
  for(j in 1:length(habitats)){
    if(i==1 & j==1){
      habitat_data <- cbind(metadata$included[i], habitats[j])
    }else{
      habitat_data <- rbind(habitat_data, cbind(metadata$included[i], habitats[j]))
    }
  }
}

habitats <- table(habitat_data[,1], habitat_data[,2])

# We also want to order habitats by their overall classification
habitat_order <- read.csv("./habitat_groups.csv")
habitat_order <- habitat_order[!habitat_order$habitat_type=="",]
habitats <- habitats[2:1,(match(habitat_order$habitat_type, colnames(habitats)))]

#colors <- pal[(factor(habitat_order$habitat_group))]
colors <- pal[habitat_order$plot_order]

colors <- rev(colors)
colnames(habitats) <- gsub("subtropical/tropical", "(sub)tropical", colnames(habitats))

# Figure 2c
barplot(rev(colSums(habitats)), xlim=c(0, 160), axes=F,
        horiz=T, col=paste(colors, "88", sep=""), border="#00000000")
axis(1, at=seq(0, 160, 20))
barplot(rev(habitats[1,]), horiz=T, col=colors, 
        add=T, border="#00000000", axes=F,names.arg = NA)
legend("right", legend=c("forest", "artificial terrestrial", "artificial aquatic",
                         "savanna", "shrubland", "grassland", "freshwater wetland",
                         "rocky habitats", "marine wetland"),
       col=pal, pch=15, pt.cex=2, bty="n")

# Study locations --------------------------------------------------------------
# Load a base map for plotting
data(wrld_simpl)
world <- spTransform(wrld_simpl, CRS("+proj=times"))
world <- world[!world$NAME == "Antarctica", ]

# Includes haphazard coordinates based on georeferenced place names
studysites <- unique(data.frame(as.numeric(metadata$Approximate.latitude[metadata$included]), 
                    as.numeric(metadata$Approximate.longitude[metadata$included])))
names(studysites) <- c("latitude", "longitude")
sp::coordinates(studysites) <- ~ longitude + latitude

# Note, this and subsequent sp::proj4string calls will likely return a warning 
# message because of the migration from PROJ4 to PROJ6; we can safely ignore it 
# since it is simply because of losing a comment that is contained in the WKT 
# string that is not in the PROJ4 string
# See: http://rgdal.r-forge.r-project.org/articles/PROJ6_GDAL3.html
sp::proj4string(studysites) <- sp::proj4string(wrld_simpl)
studysites <- sp::spTransform(studysites, sp::proj4string(world))

unstudysites <- unique(data.frame(as.numeric(metadata$Approximate.latitude[!metadata$included]), 
                    as.numeric(metadata$Approximate.longitude[!metadata$included])))
names(unstudysites) <- c("latitude", "longitude")
sp::coordinates(unstudysites) <- ~ longitude + latitude
sp::proj4string(unstudysites) <- sp::proj4string(wrld_simpl)
unstudysites <- sp::spTransform(unstudysites, sp::proj4string(world))

# Pick a nice color for the ocean
blue <- paste("#b9d1df", sep="")

# Figure 1a
plot(world, bg=blue, border=blue, lty=0,
     col="white")
plot(world, border = blue, lty=1, lwd=.25,
     col= paste("#B4775B", "15", sep=""), add=T)
# Add study sites that were not included in the meta-analysis first
points(unstudysites, pch = 17, cex = 1, lwd=.05,
       col = paste("#B65D35", "bb", sep=""))
# And then add those that were on top of that
points(studysites, pch = 17, cex = 1, lwd=.05,
       col = paste("#164057", "aa", sep=""))
legend(x = c(-14000000), y=c(-4800000), 
       legend=c("Included in meta-analysis",
       "Not included in meta-analysis"),
       bty="n", pch=17, col=c("#164057", "#B65D35"), 
       pt.cex=2, text.col = "grey40")

# Subset to just europe for an inset plot because there are so many studies
europe <- world[world$REGION==150,]
europe <- europe[!europe$NAME%in%c("Russia", "Denmark", "Iceland"),]

b <- as(raster::extent(europe), 'SpatialPolygons')
raster::crs(b) <- raster::crs(world)

# Figure 1b
plot(b, col=blue)
plot(world, border = blue, lwd=0.001,
     col= "white", add=T)
plot(world, border = blue, lwd=0.35,
     col= paste("#164057", "15", sep=""), add=T)
points(unstudysites,
       pch = 17,
       cex = .85, lwd=.05,
       col = paste("#164057", "bb", sep=""))
points(studysites,
       pch = 17,
       cex = .85, lwd=.05,
       col = paste("#B65D35", "aa", sep=""))

# Measurements of body condition and reproductive success ----------------------

# These are all the measures, including for articles not in the meta-analysis
# but that did meet inclusion criteria for the systematic review
articles <- unique(measures$Article_ID)
repro_measures <- c(measures$Reclassified.response.measure[measures$Response.type=="reproductive success"],
                measures$Reclassified.predictor.measure[measures$Predictor.type=="reproductive success"])
bc_measures <- c(measures$Reclassified.response.measure[measures$Response.type=="body condition"],
                measures$Reclassified.predictor.measure[measures$Predictor.type=="body condition"])

# Set up placeholder matrices for article x repro and article x body condition
repro <- unique(repro_measures)
bc <- unique(bc_measures)
repro.dat <- array(dim=c(length(articles), length(repro)))
bc.dat <- array(dim=c(length(articles), length(bc)))

for(i in 1:length(articles)){
  tmp <- measures[measures$Article_ID==articles[i],]
  repro.dat[i,] <- as.numeric(repro %in% append(tmp$Reclassified.predictor.measure, tmp$Reclassified.response.measure))
  bc.dat[i,] <- as.numeric(bc %in%  append(tmp$Reclassified.predictor.measure, tmp$Reclassified.response.measure))
}

# Estimate accumulation of new measures
est.repro <- vegan::poolaccum(repro.dat, permutations = 2000)
est.bc <- vegan::poolaccum(bc.dat, permutations = 2000)

# Figure S4
plot(apply(est.repro$jack1, 1, mean), type="n", ylim=c(0,80), xlim=c(0,180), axes=F,
     ylab="Extrapolated total count of measures",
     xlab="Number of included articles")
axis(1, at=seq(0, 180, 20)); axis(2, at=seq(0, 80, 10))
ci <- apply(est.repro$jack1, 1, quantile, c(0.025, 0.975))
newx <- seq(1, length(articles)-2, 1)
polygon(x=c(newx, rev(newx)), y=c(ci[1,], rev(ci[2,])),
        col=paste("#AD3429", "66", sep=""), border = F)
lines(apply(est.repro$jack1, 1, mean), lwd=4, col="#AD3429")

ci <- apply(est.bc$jack1, 1, quantile, c(0.025, 0.975))
newx <- seq(1, length(articles)-2, 1)
polygon(x=c(newx, rev(newx)), y=c(ci[1,], rev(ci[2,])),
        col=paste("#164057", "66", sep=""), border = F)
lines(apply(est.bc$jack1, 1, mean), lwd=4, col="#164057")
legend("topleft",
       legend=c("Mean", "Reproductive success", "Body condition"),
       lwd=c(2, NA, NA, NA, NA), pch=c(NA, 15, 15, 15),
       col=c("black", c("#AD3429", "#164057")), bty="n",
       pt.cex = c(NA, 2,2,2))


# Variable selection -----------------------------------------------------------

# Take all the measurements used in any study with a calculable effect
all_variables <- unique(append(
  append(measures$Reclassified.predictor.measure[measures$Predictor.type=="food"], 
                               measures$Reclassified.predictor.measure),
                               measures$Reclassified.response.measure))
# Move the common food measures to the start to make the zero counting easier
all_variables[1:5] <- all_variables[c(1,2,5,4,3)]

# If body condition measures are uncorrected, they aren't really condition
all_variables <- all_variables[!grepl("uncorrected", all_variables)]

# Set up placeholder matrix
overlap <- array(dim=c(length(all_variables), length(all_variables)))
rownames(overlap) <- colnames(overlap) <- all_variables

# Count number of studies for each pairwise combination
for(i in 1:nrow(overlap)){
  for(j in 1:ncol(overlap)){
    overlap[i,j] <- sum(measures$Reclassified.predictor.measure==rownames(overlap)[i] & measures$Reclassified.response.measure==colnames(overlap)[j], na.rm=T) +
      sum(measures$Reclassified.response.measure==rownames(overlap)[i] & measures$Reclassified.predictor.measure==colnames(overlap)[j], na.rm=T)
  }
}

diag(overlap) <- 1

zeros <- apply(overlap,2,function(x){sum(x==0)})
stopping_criteria <- sum(zeros[-c(1:4)]>0) # 1:4 are the main food measures

# Drop variables with the most missing pairwise effects
while(stopping_criteria>5){
  fewest_pairs <- which.max(zeros[-c(1:4)])+4
  overlap <- overlap[-fewest_pairs, -fewest_pairs]
  zeros <- apply(overlap,2,function(x){sum(x==0)})
  stopping_criteria <- sum(zeros[-c(1:4)]>0)
}

rownames(overlap)

# Missing data -----------------------------------------------------------------

# Even though we are now only missing a few effects, we still need a complete
# matrix for the meta-analytic structural equation model so we add in an
# artificial study with sample size of 2 and effect of 0.001 (i.e. essentially
# no effect whatsoever that could affect the model)

variables <- unique(append(effects$Predictor.measure, effects$Response.measure))
lookups <- as.matrix(table(variables, variables))
lookups[upper.tri(lookups, diag=F)] <- seq(1, length(lookups[upper.tri(lookups, diag=F)]))
diag(lookups) <- 0
Cell <- c()

# We need each pairwise combination of variables to be its own cell
for(i in 1:nrow(effects)){
  Cell[i] <- max(c(lookups[colnames(lookups)==effects$Response.measure[i], rownames(lookups)==effects$Predictor.measure[i]],
                   lookups[rownames(lookups)==effects$Predictor.measure[i], colnames(lookups)==effects$Response.measure[i]],
                   lookups[colnames(lookups)==effects$Response.measure[i], colnames(lookups)==effects$Predictor.measure[i]],
                   lookups[colnames(lookups)==effects$Predictor.measure[i], colnames(lookups)==effects$Response.measure[i]]))
}

effects$Cell <- Cell

# We are missing some pairwise correlations and we need a value for every
# pairwise correlation to fit the model, so we replace missing values with
# very, very tiny effects that have a very small sample size
missingcorrs <- which(!(1:max(lookups) %in% Cell))
missing_vars <- array(dim=c(length(missingcorrs), 2))
for(i in 1:length(missingcorrs)){
  missing_vars[i,] <- c(names(which(rowSums(lookups==missingcorrs[i])>0)), 
                        names(which(colSums(lookups==missingcorrs[i])>0)))
  effects <- data.frame(rbind(effects, c(Article_ID=9999,
                                         Study_ID=9999,
                                         ES_ID=paste("miss", missingcorrs[i], sep=""), 
                                         Habitat.type=NA,
                                         Family=NA,
                                         Species.name=NA,
                                         Response.type=NA,
                                         Response.measure=missing_vars[i,1],
                                         Predictor.type=NA,
                                         Predictor.measure=missing_vars[i,2],
                                         r=0.001, # very small effect
                                         N=2, # with very little weight
                                         Cell=missingcorrs[i])))
}

# Data prep --------------------------------------------------------------------

effects$diet <- traits$Diet.Inv[match(effects$Species.name, traits$Scientific)]
effects$func.guild <- traits$Diet.5Cat[match(effects$Species.name, traits$Scientific)]
effects$Latitude <- abs(as.numeric(metadata$Approximate.latitude[match(effects$Article_ID, metadata$Article_ID)]))

# Because some habitats were included infrequently, we reclassify
# the main set of habitats into open or semi-open (farmland, grassland, savanna,
# etc.), forest (all types, tropical, temperate, etc.), and wetland 
effects$Habitat.type[grep("arable", effects$Habitat.type)] <- "open"
effects$Habitat.type[grep("pasture", effects$Habitat.type)] <- "open"
effects$Habitat.type[grep("rural gardens", effects$Habitat.type)] <- "open"
effects$Habitat.type[grep("grassland", effects$Habitat.type)] <- "open"
effects$Habitat.type[grep("savanna", effects$Habitat.type)] <- "open"
effects$Habitat.type[grep("shrub", effects$Habitat.type)] <- "open"
effects$Habitat.type[grep("other", effects$Habitat.type)] <- "open" # sand quarry
effects$Habitat.type[grep("forest", effects$Habitat.type)] <- "forest"
effects$Habitat.type[grep("bogs", effects$Habitat.type)] <- "wetland"
effects$Habitat.type[grep("lagoon", effects$Habitat.type)] <- "wetland"
effects$Habitat.type[grep("rivers", effects$Habitat.type)] <- "wetland"

# And a dummy variable for if the habitat is forest (most common type)
effects$forest <- as.numeric(effects$Habitat.type=="forest") -.5

# Clean up data structure
effects$Cell <- factor(effects$Cell)
effects$Study_ID <- factor(effects$Study_ID)
effects$ES <- as.numeric(effects$r)
effects$inv_n <- 1/as.numeric(effects$N)
effects$Latitude <- as.numeric(effects$Latitude)
effects$diet <- as.numeric(effects$diet)

# Add dummy coded matrix indicating which constructs are included 
effects$Predictor.measure <- factor(append(effects$Predictor.measure, effects$Response.measure))[1:nrow(effects)]
effects$Response.measure <- factor(append(effects$Response.measure, as.character(effects$Predictor.measure)))[1:nrow(effects)]
pred.vars <- table(effects$ES_ID, effects$Predictor.measure)
resp.vars <- table(effects$ES_ID, effects$Response.measure)
dummy.mat <- as.data.frame.matrix(pred.vars + resp.vars)
effects <- cbind(effects, dummy.mat)
colnames(effects) <- gsub(" ", "_", colnames(effects))

# Fill in missing values 
effects$Latitude[is.na(effects$Latitude)] <- mean(effects$Latitude, na.rm=T)
effects$Habitat.type[is.na(effects$Habitat.type)] <- "forest"
effects$diet[is.na(effects$diet)] <- mean(effects$diet, na.rm=T)
effects$func.guild[is.na(effects$func.guild)] <- FALSE

effects$insectivore <- effects$func.guild=="Invertebrate"

# Publication bias -------------------------------------------------------------

# We only include the main effects in the publication bias plot, since there
# is no reason to expect there is bias in how e.g. correlation between number
# of fledglings and hatching success is reported

shades <- paste(paste(pal[7], c("05", "20", "60"), sep=""), "#ffffffff")
shades <- unlist(strsplit(shades, " "))

# Note: the warning message about dropping NA studies are the ones that
# we introduced as temporary placeholders for missing pairwise correlations so
# we can ignore the warning message in this case
rma.food <- metafor::rma(effects$ES[effects$Predictor.type=="food"], 
                         effects$inv_n[effects$Predictor.type=="food"])
# Figure S5
metafor::funnel.rma(rma.food,
                    shade = shades[1:5],
                    level = c(0.05, 0.05, 0.01, 0.01, 0.001),
                    cex.lab=1.25, col=paste(pal[7], "cc", sep=""),
                    ylim = c(0, .5), back = "white", legend=F)
pb <- metafor::regtest(rma.food$yi, rma.food$vi, model="lm")
legend("topleft", legend=c(paste("z =", round(pb$zval, 3)), paste("p =", round(pb$pval, 3))), bty="n", cex=1.5)
legend("topright", legend=c("0.05 < p < 1.00",
                            "0.01 < p < 0.05",
                            ".001 < p < 0.10", "Studies"),
       pch=c(15, 15, 15, 19), col=c(shades[c(1,3,5)], pal[7]), pt.cex=2,
       bty="n", cex=1.2)

# Pooled correlation matrix ----------------------------------------------------
# These are the unadjusted effects to get a sense of how the data are structured
Step1 <- metafor::rma.mv(yi=ES, V=inv_n,
                        data=effects,
                        random=list(~1|Study_ID,~1|ES_ID),
                        method="ML", mods=~Cell-1)

summary(Step1)
calc_i2(Step1) # I2 for unadjusted

# then we save the asymptotic variance-covariance matrix
Acov_UnadjES <- Step1$vb

# now we create a symmetrical correlation matrix
Matrix_UnadjES <- corpcor::vec2sm(Step1$b, diag = F)
diag(Matrix_UnadjES) <- 1

colnames(Matrix_UnadjES) <- variables
rownames(Matrix_UnadjES) <- variables

# Must be positive-definite for MASEM
Matrix_UnadjES <- as.matrix(Matrix::nearPD(Matrix_UnadjES, keepDiag = T)$mat)

# Meta-regressions -------------------------------------------------------------

# For the meta-regressions, we selected number of fledglings, which was the 
# outcome measure with the most studies for the two most commonly used food 
# measures: biomass/abundance and supplemental feeding

# We use the binary classification of forest here because it was the most
# frequently studied habitat type. Other covariates are latitude, diet (proportion
# of the diet that is invertebrates), and a binary classification if a species
# is classified as an insectivore

## fledglings ~ biomass/abundance --------------------------------------------
# Note: throughout this section there will be many warning messages about 'type' 
# being obsolete because of converting between factors/characters/numbers, but, 
# we just want to plot a blank canvas set up with the  appropriate axes so we 
# can ignore this warning

par(pty="s", mfrow=c(2,3))
tmpdat <- effects[effects$Response.measure=="number of fledglings" & effects$Predictor.measure=="biomass_abundance",]
tmpdat <- tmpdat[!colSums(apply(tmpdat[,c(15, 17, 18)], 1, is.na)>0),]
tmpdat <- tmpdat[order(tmpdat$ES),]
rma1 <- metafor::rma.mv(ES, inv_n, 
                        mods = ~ Latitude * forest + diet + insectivore, 
                        random=list(~1|Study_ID/ES_ID), 
                        data=tmpdat)
summary(rma1)

# Diet effect at mean latitude, non-insectivorous, and in forests

# Figure 4a
metafor::regplot(rma1, mod="diet", bg=c("#D8645A", "#1C362B")[factor(tmpdat$forest)], 
                 pch=21, type="n", pred=F, shade=F, main="Fledglings ~ Biomass/Abundance", xlab=
                   "Proportion of Diet Insectivorous")

newx <- seq(0, 110, 1)
dat0 <- as.matrix(cbind(diet=newx, 
                        forest=rep(-.5, length(newx)),
                        'Latitude:forest'=mean(tmpdat$Latitude)*-.5, 
                        Latitude=rep(mean(tmpdat$Latitude), length(newx)),
                        insectivore=rep(FALSE, length(newx))))
dat1 <- as.matrix(cbind(diet=newx, forest=rep(0.5, length(newx)), 
                        'Latitude:forest'=mean(tmpdat$Latitude)*0.5, 
                        Latitude=rep(mean(tmpdat$Latitude), length(newx)),
                        insectivore=rep(FALSE, length(newx))))
pred0 <- predict(rma1, newmods = dat0)
pred1 <- predict(rma1, newmods = dat1)

polygon(x=c(newx, rev(newx)), y=c(pred0$ci.lb, rev(pred0$ci.ub)), col="#D8645A88")
abline(lm(pred0$pred ~ newx), col="#D8645A", lwd=2)

polygon(x=c(newx, rev(newx)), y=c(pred1$ci.lb, rev(pred1$ci.ub)), col="#1C362B88")
abline(lm(pred1$pred ~ newx), col="#1C362B88", lwd=2)
points(tmpdat$ES ~ tmpdat$diet, cex=make_squares(rma1)*2, pch=21, 
       bg=c("#D8645A", "#1C362B")[factor(tmpdat$forest)])

# Latitude effect at mean diet, non-insectivorous, and in forests
# Figure 4b
metafor::regplot(rma1, mod="Latitude", bg=c("#D8645A", "#1C362B")[factor(tmpdat$forest)], pch=21, type="n", pred=F, shade=F, main="Fledglings ~ Biomass/Abundance", xlab=
                   "Latitude (Absolute)")

newx <- seq(0, 70, 1)
dat0 <- as.matrix(cbind(Latitude=newx, 
                        forest=rep(-.5, length(newx)),
                        'Latitude:forest'=newx*-.5, 
                        diet=rep(mean(tmpdat$diet), length(newx)),
                        insectivore=rep(FALSE, length(newx))))
dat1 <- as.matrix(cbind(Latitude=newx, 
                        forest=rep(0.5, length(newx)),
                        'Latitude:forest'=newx*0.5, 
                        diet=rep(mean(tmpdat$diet), length(newx)),
                        insectivore=rep(FALSE, length(newx))))
pred0 <- predict(rma1, newmods = dat0)
pred1 <- predict(rma1, newmods = dat1)

polygon(x=c(newx, rev(newx)), y=c(pred0$ci.lb, rev(pred0$ci.ub)), col="#D8645A88")
abline(lm(pred0$pred ~ newx), col="#D8645A", lwd=2)

polygon(x=c(newx, rev(newx)), y=c(pred1$ci.lb, rev(pred1$ci.ub)), col="#1C362B88")
abline(lm(pred1$pred ~ newx), col="#1C362B", lwd=2)
points(tmpdat$ES ~ tmpdat$Latitude, cex=make_squares(rma1)*2, pch=21, bg=c("#D8645A", "#1C362B")[factor(tmpdat$forest)])

# Figure 4c
plot_forest(rma1, tmpdat$Study_ID)

## fledglings ~ supplemental -------------------------------------------------

# Figure 4d
tmpdat <- effects[effects$Response.measure=="number of fledglings" & effects$Predictor.measure=="supplemental",]
tmpdat <- tmpdat[!colSums(apply(tmpdat[,c(15, 17, 18)], 1, is.na)>0),]
tmpdat <- tmpdat[order(tmpdat$ES),]

rma2 <- metafor::rma.mv(ES, inv_n, 
                        mods = ~ Latitude * forest+diet + insectivore, 
                        random=list(~1|Study_ID/ES_ID), 
                        data=tmpdat)
summary(rma2)

# Diet effect at mean latitude, non-insectivorous, and in forests
metafor::regplot(rma2, mod="diet", bg=c("#D8645A", "#1C362B")[factor(tmpdat$forest)], pch=21, type="n", pred=F, shade=F, main="Fledglings ~ Supplemental (Mealworms)", xlab=
                   "Proportion of Diet Insectivorous")

newx <- seq(0, 110, 1)
dat0 <- as.matrix(cbind(diet=newx, 
                        forest=rep(-.5, length(newx)),
                        'Latitude:forest'=mean(tmpdat$Latitude)*-.5, 
                        Latitude=rep(mean(tmpdat$Latitude), length(newx)),
                        insectivore=rep(FALSE, length(newx))))
dat1 <- as.matrix(cbind(diet=newx, forest=rep(0.5, length(newx)), 
                        'Latitude:forest'=mean(tmpdat$Latitude)*0.5, 
                        Latitude=rep(mean(tmpdat$Latitude), length(newx)),
                        insectivore=rep(FALSE, length(newx))))
pred0 <- predict(rma2, newmods = dat0)
pred1 <- predict(rma2, newmods = dat1)

polygon(x=c(newx, rev(newx)), y=c(pred0$ci.lb, rev(pred0$ci.ub)), col="#D8645A88")
abline(lm(pred0$pred ~ newx), col="#D8645A", lwd=2)

polygon(x=c(newx, rev(newx)), y=c(pred1$ci.lb, rev(pred1$ci.ub)), col="#1C362B88")
abline(lm(pred1$pred ~ newx), col="#1C362B88", lwd=2)

points(tmpdat$ES ~ tmpdat$diet, cex=make_squares(rma2)*2, pch=21, bg=c("#D8645A", "#1C362B")[factor(tmpdat$forest)])
legend("bottomleft", legend=c("Forest", "Other habitats"),
       pt.bg=c("#1C362B", "#D8645A"), pch=21, pt.cex=2, bty="n")


# Latitude effect at mean diet, non-insectivorous, and in forests
# Figure 4e
metafor::regplot(rma2, mod="Latitude", bg=c("#D8645A", "#1C362B")[factor(tmpdat$forest)], pch=21, type="n", pred=F, shade=F, main="Fledglings ~ Supplemental (Mealworms)", xlab=
                   "Latitude (Absolute)")

newx <- seq(0, 70, 1)
dat0 <- as.matrix(cbind(Latitude=newx, 
                        forest=rep(-.5, length(newx)),
                        'Latitude:forest'=newx*-.5, 
                        diet=rep(mean(tmpdat$diet), length(newx)),
                        insectivore=rep(FALSE, length(newx))))
dat1 <- as.matrix(cbind(Latitude=newx, 
                        forest=rep(0.5, length(newx)),
                        'Latitude:forest'=newx*0.5, 
                        diet=rep(mean(tmpdat$diet), length(newx)),
                        insectivore=rep(FALSE, length(newx))))
pred0 <- predict(rma2, newmods = dat0)
pred1 <- predict(rma2, newmods = dat1)

polygon(x=c(newx, rev(newx)), y=c(pred0$ci.lb, rev(pred0$ci.ub)), col="#D8645A88")
abline(lm(pred0$pred ~ newx), col="#D8645A", lwd=2)

polygon(x=c(newx, rev(newx)), y=c(pred1$ci.lb, rev(pred1$ci.ub)), col="#1C362B88")
abline(lm(pred1$pred ~ newx), col="#1C362B", lwd=2)
points(tmpdat$ES ~ tmpdat$Latitude, cex=make_squares(rma2)*2, pch=21, bg=c("#D8645A", "#1C362B")[factor(tmpdat$forest)])

# Figure 4f
plot_forest(rma2, tmpdat$Study_ID)

par(mfrow=c(1,1))
# Note: these plots have been post-hoc combined in Inkscape to be more readable,
# For example increasing the size of the axes and increasing text size in the 
# forest plots and rearranging to a vertical structure to match page layouts


# Insectivore status has no effect, but does it explain heterogeneity?
rma3 <- metafor::rma.mv(ES, inv_n, 
                        mods = ~ Latitude * forest+diet, 
                        random=list(~1|Study_ID/ES_ID), 
                        data=tmpdat)
summary(rma3) 
calc_i2(rma3)

# Heterogeneity in effects -----------------------------------------------------

# We already know there is heterogeneity and some of it can be explained by the
# covariates from the meta-regressions, but want to look at it for all the 'main'
# effects, and also the other effects used to partition variance when fitting 
# the structural model. First, we look at just bird responses to food as predictor
food.dat <- effects[effects$Predictor.type=="food",]

# Study 9999 is the non-existent study with N=2 for the missing correlations
food.dat <- food.dat[!is.na(food.dat$Study_ID) & !food.dat$Study_ID=="9999",]
mod1 <- metafor::rma.mv(yi=ES, V=inv_n,
                data=food.dat,
                random=list(~1|Study_ID/ES_ID),
                method="ML",
                mods=~ Habitat.type*Latitude + diet)
mod1
calc_i2(mod1)
# Suggests we should retain the interaction effect between habitat and latitude

# Now for all effects, including between outcome measures
# Again, drop the non-existent study for checking moderators
tmp.effects <- effects[!effects$Study_ID=="9999",]
mod2 <- metafor::rma.mv(yi=ES, V=inv_n,
                        data=tmp.effects,
                        random=list(~1|Study_ID/ES_ID),
                        method="ML",
                        mods=~ Habitat.type*Latitude + diet)
mod2
calc_i2(mod2)

# Diet does not come out as significant in these models, but is in the same direction
# as the meta-regressions and we like it as a moderator for biological reasons

# What if we also throw in functional guild as suggested by reviewer 2?
mod3 <- metafor::rma.mv(yi=ES, V=inv_n,
                        data=tmp.effects,
                        random=list(~1|Study_ID/ES_ID),
                        method="ML",
                        mods=~ Habitat.type*Latitude + diet + insectivore)
mod3
calc_i2(mod3)
# It is not significant, but, it does explain some heterogenetiy so we'll keep it

# Covariate adjusted effects ---------------------------------------------------

# Now we add in moderator variables for the selected covariates, and for which
# sets of variables are represented by a correlation
Step2 <- metafor::rma.mv(yi=ES, V=inv_n,
                         data=effects,
                         random=list(~1|Study_ID/ES_ID),
                         method="ML",
                         mods=~ Habitat.type*Latitude + diet + insectivore +
                           biomass_abundance +
                           chick_body_mass_corrected_by_age +
                           chick_body_mass_residuals +
                           chick_tarsus_corrected_by_age +
                           chick_wing_length_corrected_by_age +
                           fledging_success_rate +
                           frass +
                           hatching_success_rate +
                           nest_success +
                           number_of_fledglings +
                           number_of_nestlings +
                           supplemental)
Step2

predmean <- colMeans(metafor::model.matrix.rma(Step2))[c(1:7,20:22)] %*% Step2$b[c(1:7, 20:22),1]
cell.adj <- (metafor::model.matrix.rma(Step2))[,-c(1:7,20:22)] %*% Step2$b[-c(1:7, 20:22),1]

effects$residuals <- metafor::residuals.rma(Step2)

effects$Covadj_ES <- effects$residuals + as.numeric(predmean)

# How do the adjusted covariates map back to the unadjusted ones?s
plot(effects$Covadj_ES ~ effects$ES)

# Meta-analytic structural equation model --------------------------------------

# We now re-estimate the pooled correlation matrix using the adjusted values
Step3 <-metafor::rma.mv(yi=Covadj_ES, V=inv_n,
                        data=effects,
                        random=list(~1|Study_ID/ES_ID),
                        method="ML",
                        mods=~ Cell-1)
summary(Step3)

#Save the asymptotic variance-covariance matrix
Acov_AdjES <- Step3$vb
#Create the correlation matrix
Matrix_AdjES <- corpcor::vec2sm(Step3$b, diag = F)
diag(Matrix_AdjES) <- 1

colnames(Matrix_AdjES) <- variables
rownames(Matrix_AdjES) <- variables

Matrix_AdjES <- as.matrix(Matrix::nearPD(Matrix_AdjES, keepDiag = T)$mat)

n <- sum(as.numeric(effects$N)) 
variables <- gsub(" ", "_", variables) # bad coding practice to overwrite this

# Specify our model structure using lavaan syntax
model1 <- "# measurements
           food =~ biomass_abundance + supplemental + suppression + frass
           reproduction =~ number_of_nestlings + number_of_fledglings + fledging_success_rate + hatching_success_rate + nest_success
           condition =~  chick_body_mass_corrected_by_age + chick_wing_length_corrected_by_age + chick_tarsus_corrected_by_age + chick_body_mass_residuals
           # regressions
           reproduction ~ food
           condition ~ food
           # correlations
           reproduction ~~ condition
           food ~~1*food"

# Check that the model matches the relationships we think we specified
plot(
  model1,
  layout = "tree",
  style = "mx",
  edge.label.cex = .75,
  edge.color = "black",
  label.color = "black"
)

# Set up the matrices we need
RAM1 <- metaSEM::lavaan2RAM(
  model1,
  obs.variables = variables,
  A.notation = "on",
  S.notation = "with",
)

Smatrix <- metaSEM::as.mxMatrix(RAM1$S, name="S")
Amatrix <- metaSEM::as.mxMatrix(RAM1$A, name="S")
Fmatrix <- metaSEM::as.mxMatrix(RAM1$F, name="S")

dv <- apply(Amatrix$free, 1, any)
diag(Smatrix$free)[dv] <- TRUE # allow to vary

# Here is the version with unadjusted effect sizes
fit1 <- metaSEM::wls(
  Matrix_UnadjES,
  Acov_UnadjES,
  n = n, 
  Smatrix = Smatrix, Fmatrix = Fmatrix, Amatrix = Amatrix,
  cor.analysis = TRUE,
  intervals = "z",
  diag.constraints = T
)
summary(fit1)

# And the version corrected by habitat and response type
# Qualitatively, they are not too different, but in the adjusted version 
# we have narrower confidence intervals because some of the noise has
# been removed prior to fitting the model
fit2 <- metaSEM::wls(
  Matrix_AdjES,
  Acov_AdjES,
  n = n, 
  Smatrix = Smatrix, Fmatrix = Fmatrix, Amatrix = Amatrix,
  cor.analysis = TRUE,
  intervals = "z",
  diag.constraints = T
)
summary(fit2)

# Interactive plot of the model ------------------------------------------------

# Note: for the manuscript, we have reproduced the SEM plot in PowerPoint to
# improve readability and resolution of the figure, thus, the aesthetics do
# not match the interactive plot, however, the results are the same

# Figure 6

sem.plot <- metaSEM::meta2semPlot(fit2, labels="labels")

nodedf <- data.frame(id=append(variables, c("reproduction", "food", "condition")),
                     label=gsub("_", " ", append(variables, c("reproduction", "food", "condition"))),
                     group=1,
                     font.size=10,
                     value=c(rep(7, 13), rep(20, 3)))

edgedf <- data.frame(label=sem.plot@Pars$label,
                     from=sem.plot@Pars$lhs,
                     to=sem.plot@Pars$rhs,
                     est=sem.plot@Pars$est)

edgedf$label[edgedf$label=="" & edgedf$from!=edgedf$to] <- paste(edgedf$from[edgedf$label=="" & edgedf$from!=edgedf$to],
                                                                 edgedf$to[edgedf$label=="" & edgedf$from!=edgedf$to],sep="with")

coefs <- summary(fit2)$coefficient
coefs$label <- rownames(coefs)
edgedf <- dplyr::left_join(edgedf, coefs, by="label")                     
edgedf$lbound[edgedf$lbound==edgedf$ubound] <- NA
edgedf$ubound[edgedf$lbound==edgedf$ubound] <- NA

edges <- edgedf
nodes <- nodedf
nodes$font.size <- 25
edges$value <- edges$est
edges$font.size <- 25
edges$label <- paste(round(edges$est, 2), "\n", "(",
                     round(edges$lbound, 2), ", ",
                     round(edges$ubound, 2),")", sep="")
edges$label[is.na(edges$lbound)] <- round(edges$est, 2)[is.na(edges$lbound)]

nodes$label<- gsub(" corrected by age", "\ncorrected by age", nodes$label)

visNetwork(nodes, edges, width = "100%") %>%
  visIgraphLayout() %>%
  visNodes(
    scaling = list(
      min = min(nodes$value)+10,
      max = max(nodes$value)+10, 
      label=list(enabled="true")
    ),
    shadow = list(enabled = F)
  ) %>%
  visEdges(arrows =list(to = list(enabled = TRUE, scaleFactor = 1)),    
           scaling = list(
             min = min(edges$value),
             max = max(edges$value), 
             label=list(enabled="true")
           ),
           shadow = FALSE,
           color = list(color = "#06A9B2", highlight = "#06A9B2")
  ) %>%
  visOptions(
    highlightNearest = list(enabled = T, degree = 1, hover = T)
  ) %>%
  visGroups(groupname = "1", color = list(background = "grey",
                                          border = "#06A9B2",
                                          highlight = list(background="grey", border="#06A9B2"),
                                          hover=list(background="grey", border="#06A9B2"))) %>%
  visLayout(randomSeed = 11)


# Various plots ----------------------------------------------------------------
  
  ## Heatmap of unadjusted pooled correlations -----
  # Figure 3 above the diagonal
  hm <- lookups
  betas <- as.numeric(Step1$b)
  names(betas) <- gsub("Cell", "", rownames(Step1$b))
  betas <- betas[order(as.numeric(names(betas)))]
  betas[missingcorrs] <- NA
  
  hm[upper.tri(hm)] <- betas
  tm <- t(hm)[,-nrow(hm)]
  hm[lower.tri(hm)] <- tm[lower.tri(tm)]
  rownames(hm) <- colnames(hm) <- paste(c(101,107,106,108,109,113,102,110,114,112,111,103,104),
                                        rownames(hm), sep="_")
  hm <- hm[order(rownames(hm), decreasing = T), order(colnames(hm), decreasing = T)]
  hm[lower.tri(hm)] <- NA
  
  tmp <- as.data.frame(hm)
  tmp$Freq[tmp$Freq==0] <- NA
  
  ggplot(tmp, aes(variables, variables.1, fill = Freq)) +
    geom_tile() +
    scale_fill_gradient2(high = "#164057",
                         mid="#FBF7EE",
                         low = "#6B140B",
                         na.value = "white",
                         limits=c(-1, 1)) + # fill with color
    theme(
      line = element_blank(),
      # remove the background, tickmarks, etc
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 9
      ),
      axis.text.y = element_text(size = 9),
      axis.title = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size = 10)
    ) +
    coord_equal()+
    geom_text(aes(variables, variables.1, label = round(Freq,2)), color = "black", size = 4)
  
  # Note: we can ignore this warning message since we intentionally gave
  # it a bunch of blanks for the lower triangle
  
  # Heatmap with covariate adjusted correlations -----
  # Figure S7 above the diagonal
  hm <- lookups
  betas <- as.numeric(Step3$b)
  names(betas) <- gsub("Cell", "", rownames(Step3$b))
  betas <- betas[order(as.numeric(names(betas)))]
  betas[missingcorrs] <- NA
  
  hm[upper.tri(hm)] <- betas
  tm <- t(hm)[,-nrow(hm)]
  hm[lower.tri(hm)] <- tm[lower.tri(tm)]
  rownames(hm) <- colnames(hm) <- paste(c(101,107,106,108,109,113,102,110,114,112,111,103,104),
                                        rownames(hm), sep="_")
  hm <- hm[order(rownames(hm), decreasing = T), order(colnames(hm), decreasing = T)]
  hm[lower.tri(hm)] <- NA
  
  tmp <- as.data.frame(hm)
  tmp$Freq[tmp$Freq==0] <- NA
  
  ggplot(tmp, aes(variables, variables.1, fill = Freq)) +
    geom_tile() +
    scale_fill_gradient2(high = "#164057",
                         mid="#FBF7EE",
                         low = "#6B140B",
                         na.value = "white",
                         limits=c(-1, 1)) + # fill with color
    theme(
      line = element_blank(),
      # remove the background, tickmarks, etc
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 9
      ),
      axis.text.y = element_text(size = 9),
      axis.title = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size = 10)
    ) +
    coord_equal()+
    geom_text(aes(variables, variables.1, label = round(Freq,2)), color = "black", size = 4)
  
  # Note: we can ignore this warning message since we intentionally gave
  # it a bunch of blanks for the lower triangle
  
  # Heatmap of the number of studies underlying those correlations -----
  # Figure 3 and Figure S7 below the diagonal
  counts <- lookups
  counts[!is.na(counts)] <- NA
  
  realstudies <- effects[!effects$Article_ID=="9999",]

  for(i in 1:nrow(counts)){
    for(j in 1:ncol(counts)){
      counts[i,j] <- sum(realstudies$Predictor.measure==rownames(counts)[i] & realstudies$Response.measure==colnames(counts)[j], na.rm=T) +
        sum(realstudies$Response.measure==rownames(counts)[i] & realstudies$Predictor.measure==colnames(counts)[j], na.rm=T)
    }
  }
  
  hmN <- counts
  
  # To get things in order of food, body condition, then reproductive success
  rownames(hmN) <- colnames(hmN) <- paste(c(101,107,106,108,109,113,
                                            102,110,114,112,111,103,104),
                                          rownames(hmN), sep="_")
  hmN <- hmN[order(rownames(hmN), decreasing = T), order(colnames(hmN), decreasing = T)]
  hmN[upper.tri(hmN)] <- NA
  
  tmp <- as.data.frame(hmN)
  tmp$Freq[tmp$Freq==0] <- NA
  
  ggplot(tmp, aes(variables, variables.1, fill = Freq)) +
    geom_tile() +
    scale_fill_gradient(high = "#b8820b",
                        low = "#e0ebbc",
                        na.value = "white",
                        limits=c(0, 50)) + # fill with color
    theme(
      line = element_blank(),
      # remove the background, tickmarks, etc
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 9
      ),
      axis.text.y = element_text(size = 9),
      axis.title = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size = 10)
    ) +
    coord_equal()+
    geom_text(aes(variables, variables.1, label = round(Freq,2)), color = "black", size = 4)
  
  # Note: we can ignore this warning message since we intentionally gave
  # it a bunch of blanks for the upper triangle
  
  # Supplemental figures with RMA for each combination -----
  # Figures S8-S39 
  
  subgroupdat <- effects[!effects$Study_ID=="9999",]
  predictors <- c("biomass_abundance", "suppression", "supplemental", "frass")
  variables <- unique(append(effects$Predictor.measure, effects$Response.measure))
  responses <- variables[!variables %in% predictors]
  responses <- responses[!is.na(responses)]

    samp.size <- mean.mat <- ci.u.mat <- ci.l.mat <- array(dim=c(length(responses), length(predictors)))
  
  dir.create("SI-figures")
  setwd("./SI-figures")
  png(width=1200, height=600, units="px")
  layout(mat=matrix(c(2,1,2,4,2,3), nrow=2, ncol=3),
         heights=c(1, 2),
         widths=c(2,2,2))
  par(pty="s", xpd=T, las=1)
  n.sim <- 100
  for(i in 1:length(responses)){
    for(j in 1:length(predictors)){
      
      rm(rma1) 
      tmpdat <- subgroupdat[which(subgroupdat$Response.measure==responses[i] & subgroupdat$Predictor.measure==predictors[j]),]
      tmpdat <- tmpdat[order(tmpdat$ES, decreasing = F),]
      if(nrow(tmpdat)>1){
        
        rma1 <- metafor::rma.mv(yi=Covadj_ES, V=inv_n,
                                data=tmpdat,
                                mods=
                                random=list(~1|Study_ID,~1|ES_ID),
                                method="ML",
                                mods=~ Cell-1)
        
        xvals <- seq(-5, 5, 0.5)
        upci <- max(max(rma1$yi + 1.96*sqrt(rma1$vi), na.rm = T), rma1$ci.ub)
        lowci <- min(min((rma1$yi - 1.96*sqrt(rma1$vi)), na.rm=T), rma1$lb)
        xlims <- c(xvals[max(which(xvals<lowci))], xvals[min(which(xvals>upci))])
        ylims <- c(1, nrow(tmpdat))
        
        wi <- weights(rma1)
        psize <- wi/sum(wi, na.rm = TRUE)
        rng <- max(psize, na.rm = TRUE) - min(psize, na.rm = TRUE)
        k <- length(wi)
        if (rng <= .Machine$double.eps^0.5) {
          psize <- rep(1, k)
        }else {
          psize <- (psize - min(psize, na.rm = TRUE))/rng
          psize <- (psize * 1) + 0.5
        }
        if (all(is.na(psize))){
          psize <- rep(1, k)
        } 
        
        par(xpd=F)
        
        if(nrow(tmpdat)>20){newcex <- 1.1}else{newcex <- 1.5}
        if(nrow(tmpdat)>40){newcex <- .9}
        
        metafor::forest(rma1, slab=paste(tmpdat$Study_ID, sep=""), 
                        cex.lab=1, xlab="Effect", cex=newcex)
        
        text(x=xlims[1], y=length(wi)+c(3), cex=1.5, adj=0, 
             bquote(paste("Q = ", .(formatC(rma1$QE, digits=2, format="f")),
                          "; df = ", .(formatC(rma1$k - rma1$p, digits=2, format="f")),
                          "; p = ", .(formatC(rma1$QEp, digits=2, format="f")))))
        
        I2 <- calc_i2(rma1)
        
        text(x=xlims[1], y=length(wi)+c(2), 
             adj=0, bquote(paste(I^2, " = ", .(formatC(I2, digits=2, format="f")), "; ", tau^2, " = ", .(formatC(rma1$tau2, digits=2, format="f")))),
             cex=1.5)
        
        par(xpd=T, mar=c(5,5,0,5))
        plot(0,0, type="n", axes=F, ylab="", xlab="")
        text(x=0, y=-1, paste(responses[i], 
                              gsub("_", " ", predictors[j]), sep=" ~ "), 
             cex=2)
        
        mean.mat[i,j] <- as.numeric(rma1$beta)
        ci.l.mat[i,j] <- as.numeric(rma1$ci.lb)
        ci.u.mat[i,j] <- as.numeric(rma1$ci.ub)
        samp.size[i,j] <- as.numeric(rma1$k.all)
        if(nrow(tmpdat)>9 & rma1$QEp<0.05){

          # # If heterogeneity, then check cumulative
          meandiff <- array(dim=c(n.sim, nrow(tmpdat)))
          for(n in 1:n.sim){

            tmpdat <- tmpdat[sample(1:nrow(tmpdat)),]

            loci <- est <- upci <- c()

            for(z in 2:nrow(tmpdat)){
              tryCatch({rmatmp <- metafor::rma(yi=tmpdat$Covadj_ES[1:z], vi=tmpdat$inv_n[1:z],
                                               control=list(stepadj=1, maxiter=1000))},error=function(e){
                                               })
              est[z] <- rmatmp$b
              upci[z] <- rmatmp$ci.ub
              loci[z] <- rmatmp$ci.lb
            }
            est <- est[-1]

            for(m in 2:length(est)){
              meandiff[n,m] <- abs(diff(c(est[m],est[m-1])))
            }
          }

          yl <- c(0,max(meandiff, na.rm=T)+.05)
          xl <- c(1, ceiling(ncol(meandiff)/10)*10)

          plot(apply(meandiff, 2, mean, na.rm=T), type="n", lwd=2,
               ylab="Difference in cumulative mean estimate", xlab="Number of studies added",
               axes=F, ylim=yl, xlim=xl, cex.lab=1.5)
          axis(1)
          axis(2, at = seq(min(yl), max(yl), 0.02))
          ci <- apply(meandiff, 2, quantile, c(0.025, 0.975), na.rm=T)
          newx <- seq(1, nrow(tmpdat), 1)[-which(is.na(ci[1,]))]
          polygon(x=c(newx, rev(newx)), y=c(ci[1,-which(is.na(ci[1,]))], rev(ci[2,-which(is.na(ci[2,]))])),
                  col=paste(pal[7], "44", sep=""), border = F)
          lines(apply(meandiff, 2, mean, na.rm=T), lwd=4, col=pal[7])
          abline(h=0.05, lty=2, lwd=2)
        }else{
          plot(0,0, type="n", axes=F, ylab="", xlab="")
        }
        par(xpd=F, las=0)
        metafor::funnel.rma(rma1,
                            shade = shades[1:5],
                            cex=2.5, lwd=2.5, 
                            cex.lab=1.5, 
                            cex.axis=1.5, 
                            level = c(0.05, 0.05, 0.01, 0.01, 0.001),
                            col=paste(pal[7], "cc", sep=""),
                            ylim = c(0, .5), back = "white", 
                            legend=F)
        par(xpd=T)
        #title(ylab="Standard error", line=4, cex.lab=1.5)
        pb <- summary(lm(as.numeric(residuals(rma1)) ~ rma1$vi))$coefficients
        if(any(grep("rma", rownames(pb)))){
          if(grep("rma", rownames(pb))){
            legend("topleft", legend=c(paste("t =", round(pb['rma1$vi','t value'], 3)), paste("p =", round(pb['rma1$vi','Pr(>|t|)'], 3))), bty="n", cex=1.5)
          }
        }
      }
    }
  }
  dev.off()
  
  # Summary of supplemental figures -----
  # Figure 5
  par(mfrow=c(1,1), xpd=F)
  wi <- as.numeric(t(samp.size))
  psize <- wi/sum(wi, na.rm = TRUE)
  rng <- max(psize, na.rm = TRUE) - min(psize, na.rm = TRUE)
  k <- length(wi)
  if (rng <= .Machine$double.eps^0.5) {
    psize <- rep(1, k)
  }else {
    psize <- (psize - min(psize, na.rm = TRUE))/rng
    psize <- (psize * 1) + 0.5
  }
  if (all(is.na(psize))){
    psize <- rep(1, k)
  } 
  
  newpal <- pal[c(1,2,3,7)]
  
  summ.dat <- data.frame(cbind(mean=as.numeric(t(mean.mat)), 
                               ci.l=as.numeric(t(ci.l.mat)), ci.u=as.numeric(t(ci.u.mat)),
                               N=as.numeric(t(samp.size))))
  summ.dat$types <- rep(predictors, length(responses))
  summ.dat$response <- rep(responses, each=length(predictors))
  
  
  par(las=1)
  plot(summ.dat$mean, pch=15, col=newpal[factor(summ.dat$types)], cex=psize*2, ylim=c(-.6, .8),
       axes=F, ylab="Effect")
  axis(2)
  arrows(x0=seq(1, nrow(summ.dat)), y0=summ.dat$ci.l, y1=summ.dat$ci.u,
         length=0.1, angle=90, code=3, col=newpal[factor(summ.dat$types)], lwd=2)
  abline(v=seq(0.5, nrow(summ.dat), 4), lty=2, col="grey60", lwd=0.7)
  abline(h=0, lty=2)
  legend("bottomright", legend=c(gsub("_", "/", predictors)), col=newpal, pch=15, bty="n", pt.cex=2)
  text(seq(2, nrow(summ.dat), 4), y= -.5, responses, srt=90)

  # We made this more readable in Inkscape and changed the blue to a 
  # lighter shade to make it easier to distinguish from green
  
  # Meta-regression covariates -----
  # Figure S6
  par(las=1)
  plot(Step2$b[-c(8:19)], pch=15, 
       col=(pal[c(4,7)])[factor(Step2$b>0)], cex=2, ylim=c(-1, 2.5),
       axes=F, ylab="Effect")
  axis(2)
  arrows(x0=seq(1, 10), y0=Step2$ci.lb[-c(8:19)], 
         y1=Step2$ci.ub[-c(8:19)],
         length=0.1, angle=90, code=3, col=(pal[c(4,7)])[factor(Step2$b>0)], lwd=2)
  abline(h=0, lty=2)
  legend("bottomright",
         legend=c("Positive", "Negative"),
         col=(pal[c(7,4)]), pch=15, bty="n", pt.cex=2)
  text(seq(1, 11, 1), y= -.5, c("Intercept",
                                "Open habitat",
                                "Urban",
                                "Wetland",
                                "Latitude",
                                "Proportion insectivorous",
                                "Diet guild",
                                "Latitude*Open",
                                "Latitude*Urban",
                                "Latitude*Wetland"), srt=90)
  text(seq(1, 11, 1)[Step2$ci.ub[-c(8:19)]*Step2$ci.lb[-c(8:19)]>0], 
       y=Step2$ci.ub[-c(8:19)][Step2$ci.ub[-c(8:19)]*Step2$ci.lb[-c(8:19)]>0]+.1,
       c("*"), cex=3)
  
  # fin