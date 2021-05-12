# Seedling growth gradients interact with homogeneous disturbance regimes
# to explain tree cover discontinuities in savannas
# Code author: R. Holdo (rholdo@uga.edu)

library(dplyr)
library(ggplot2)
library(nlme)
library(lme4)
library(reshape2)
library(grid)
library(gridExtra)
library(cowplot)
library(sp)
library(spdep)
library(spatialEco)

# The following assumes you are using RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
df <- read.csv("Seedling_final.csv")
# Note: 'Subplot' refers to an individual seedling
# Extract species info to reuse later
df.species <- subset(df, Per == 'Jun 2019')
df.species <- df.species[, c(1,3)]

# What proportion of individual stems were not found at some point in the study?
100 - mean(df$Found, na.rm = TRUE) * 100
# Only 0.507 %
# How many were lost for good at some point?
# This does not count those missing only in the very last survey (N = 2)
lost <- aggregate((1 - Found) ~ Subplot, df, function(x){sum(x, na.rm = TRUE)})
names(lost)[2] <- 'Lost'
lost2plus <- subset(lost, Lost > 1)
# These were 'lost' for more than 2 surveys. How many of these were
# still missing by the last survey?
lostfinal <- subset(df, Per == "Jun 2019" & Found == 0)
lost2plus$Subplot
lostfinal$Subplot
# There is no overlap in the two groups
# There were therefore no seedlings that were never again found

# Include only cases where seedling was found
df <- subset(df, Found == 1)
df$Found <- NULL
df <- df[order(df$Subplot, df$Day), ]

# Obtain mean grass biomass across subplots
grass.ag <- aggregate(Grass ~ Subplot, df, function(x){mean(x, na.rm = TRUE)})

# For seedlings that appeared dead but came back to life,
# Recategorize as 'topkilled' instead of dead, and fill
# corresponding zeros in 'Alive' column with 1s
subplots <- unique(df$Subplot)
N <- length(subplots)
for (i in 1:N){
  sub <- subset(df, Subplot == subplots[i])
  idx <- max(which(sub$Alive == 1))
  misclass.death <- which(sub$Alive[1:idx] == 0)
  if (length(misclass.death) > 0){
    sub$Alive[misclass.death] <- 1
    sub$Topkilled[misclass.death] <- 1
  }
  if (i == 1) df.cl <- sub
  else df.cl <- rbind(df.cl, sub)
}

# Proportion of seedlings with mortality and topkill
# Make subset without Initial period
df.cl2 <- subset(df.cl, Per != "Initial") # Exclude setup survey
df.cl2$Dead <- 1 - df.cl2$Alive

# Did an individual die over the course of the survey (3 years)?
df.ag.mort <- aggregate(Dead ~ Subplot,
                        df.cl2, function(x){max(x, na.rm = TRUE)})
# Did an individual sustain damage in a given year?
df.ag.dam <- aggregate(cbind(Topkilled, Damage) ~ Subplot + Year,
                       df.cl2, function(x){max(x, na.rm = TRUE)})
# Calculate mean fire and herbivory damage in year 3
df.ag.fireherb <- aggregate(cbind(Fire, Herbivory) ~ Subplot, subset(df.cl2, Year == 2019),
                            function(x){mean(x, na.rm = TRUE)})
names(df.ag.fireherb)[2:3] <- c('F2019', 'H2019')
# Mean annual damage and topkill
df.ag <- aggregate(cbind(Topkilled, Damage) ~ Subplot, df.ag.dam, mean)
df.ag <- merge(df.ag, df.ag.fireherb, by = 'Subplot', all = TRUE)
df.ag <- merge(df.ag, df.ag.mort, by = 'Subplot', all = TRUE)

# Find maximum vertical growth across 2017, 2018, and 2019 seasons
# First, make initial height = 0 for individuals with size NA in Initial period
# These were recently topkilled at the very outset
df$Basal_cm <- ifelse(is.na(df$Basal_cm) & df$Per == 'Initial', 0, df$Basal_cm)
df$Ht_m <- ifelse(is.na(df$Ht_m) & df$Per == 'Initial', 0, df$Ht_m)
year1 <- subset(df, Year == 2017)
# Keep only seedlings that survived year 1
year1IDs <- subset(year1, Per == "May 2017" & Alive == 1)
year1 <- subset(year1, Subplot %in% year1IDs$Subplot)
year1H0 <- subset(year1, Per == "Initial")
year1H0 <- year1H0[, c(1,11)]
year1H0$Ht_m <- ifelse(is.na(year1H0$Ht_m), 0, year1H0$Ht_m)
year1Hmax <- aggregate(Ht_m ~ Subplot, year1[year1$Per != "Initial", ],
                       function(x){max(x, na.rm = TRUE)})
year1H <- merge(year1H0, year1Hmax, by = 'Subplot')
names(year1H)[2:3] <- c("H0", "Hmax")

year2 <- subset(df, Per == "May 2017" | Year == 2018)
# Keep only seedlings that survived year 2
year2IDs <- subset(year2, Per == "May 2018" & Alive == 1)
year2 <- subset(year2, Subplot %in% year2IDs$Subplot)
year2H0 <- subset(year2, Per == "May 2017")
year2H0 <- year2H0[, c(1,11)]
year2Hmax <- aggregate(Ht_m ~ Subplot, year2[year2$Per != "May 2017", ],
                       function(x){max(x, na.rm = TRUE)})
year2H <- merge(year2H0, year2Hmax, by = 'Subplot')
names(year2H)[2:3] <- c("H0", "Hmax")

year3 <- subset(df, Per == "May 2018" | Year == 2019)
# Keep only seedlings that survived year 3
year3IDs <- subset(year3, Per == "Jun 2019" & Alive == 1)
year3 <- subset(year3, Subplot %in% year3IDs$Subplot)
year3H0 <- subset(year3, Per == "May 2018")
year3H0 <- year3H0[, c(1,11)]
year3Hmax <- aggregate(Ht_m ~ Subplot, year3[year3$Per != "May 2018", ],
                       function(x){max(x, na.rm = TRUE)})
year3H <- merge(year3H0, year3Hmax, by = 'Subplot')
names(year3H)[2:3] <- c("H0", "Hmax")

# Merge together
growth <- merge(year1H, year2H, by = 'Subplot', all = TRUE)
growth <- merge(growth, year3H, by = 'Subplot', all = TRUE)
growth$deltaH1 <- growth$Hmax.x - growth$H0.x
growth$deltaH2 <- growth$Hmax.y - growth$Hmax.x
growth$deltaH3 <- growth$Hmax - growth$Hmax.y
growth <- growth[, -c(4,6)]
names(growth)[2:5] <- c('H0', 'H1', 'H2', 'H3')
growth.l <- melt(growth, id.vars = "Subplot", measure.vars = 6:8,
                 variable.name = "Year", value.name = "deltaH")
growth.l$Year <- ifelse(growth.l$Year == "deltaH1", 2017, 
                        ifelse(growth.l$Year == "deltaH2", 2018, 2019))

# Combine growth data with damage data by year
grdam <- merge(growth.l, df.ag.dam, by = c('Subplot', 'Year'), all = TRUE)
grundam <- subset(grdam, Damage == 0 & Topkilled == 0)
grundam.ag <- aggregate(deltaH ~ Subplot, grundam, 
                        function(x){mean(x, na.rm = TRUE)})
names(grundam.ag)[2] <- 'deltaH.undam'

# Combine aggregated disturbance dataset with aggregated grass biomass dataset 
df.ag <- merge(df.ag, grass.ag, by = 'Subplot', all = TRUE)
# Now combine with growth data
growth$deltaH <- rowMeans(growth[,6:8], na.rm = TRUE)
growth.all <- growth[, c(1,9)]
df.ag <- merge(df.ag, growth.all, by = 'Subplot', all = TRUE)
# Now mean growth across 3 years has been added
# Add mean growth for periods without damage
df.ag <- merge(df.ag, grundam.ag, by = 'Subplot', all = TRUE)
# Reinsert species variable
df.ag <- merge(df.ag, df.species)

# Figure for calibration of light sensor data
ldrcal <- read.csv("Light_cal.csv")
ldrreg <- lm(log(PAR) ~ log(R), ldrcal)
figA1 <- ggplot(ldrcal) + geom_point(aes(x = log(R), y = log(PAR)), 
                                     shape = 21, size = 6, stroke = 2) +
  labs(x = expression(paste('log R (', Omega, ')')), 
       y = expression(paste('log PAR (', mu, 'mol m'^'-2', ' s'^'-1',')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(2)),
        axis.title.x = element_text(vjust = 0, size = rel(2)),
        axis.text = element_text(size = rel(2)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)))
figA1 <- figA1 + geom_abline(intercept = coef(ldrreg)[1], slope = coef(ldrreg)[2], lwd = 1.2) +
  annotate("text", x = 6.5, y = 6, label = "R ^ 2 == 0.99", parse = TRUE, size = 8)

# Obtain TC data, recalculate distances, combine with transect data
tc <- read.csv("Transect_TC_data.csv")
tc <- subset(tc, Type == "Seedling" | Type == "Sensor")
# Estimate basal area, tree density and seedling density  along transects from belt transect data
# Calculate distances along transect using PCA
# Recover distance along and from transect for ground tree survey
tree <- read.csv("Transect_Tree_Density_Data.csv")
# Obtain mean crown area
tree$Crown <- pi * tree$Can1 * tree$Can2 / 4
Cr.mean <- mean(tree$Crown, na.rm = TRUE)
# Calculate basal and canopy area
tree <- tree[, -(7:10)]
tree$Basal2 <- ifelse(is.na(tree$Basal2), 0, tree$Basal2)
tree$Basal3 <- ifelse(is.na(tree$Basal3), 0, tree$Basal3)
tree$Bas <- (tree$Basal1 ^ 2 + tree$Basal2 ^ 2 + tree$Basal3 ^ 2) * pi / 4
tree <- tree[, -c(1, 7:11)]
names(tree)[2] <- "ID"
# Import seedling data from Deus' survey as well
seed <- read.csv("Transect_Seedling_Density_Data.csv")
seed <- seed[, c(2:5,11)]

# Join the seedling location and tree/seedling survey datasets
tc2 <- tc[, c(2,1,3,4,5)]
tree2 <- tree[, 1:4]
tree2$Type <- "Tree.tran" # Belt transect tree
tree2 <- tree2[, c(1,2,5,3,4)]
seed2 <- seed[, 1:4]
seed2$Type <- "Seedl.tran" # Belt transect seedling
seed2 <- seed2[, c(1,2,5,3,4)]
joint <- rbind(tc2, tree2)
joint <- rbind(joint, seed2)
tran <- unique(joint$Site)
for (i in 1:length(tran)){
  sub <- subset(joint, Site == tran[i])
  # Use seedling/sensor locations as reference for transect line
  # The transect line is found via PCA on XY data
  # Most parsimonious way to join survey seedlings and belt transect data
  sub.seedl <- subset(sub, Type == "Seedling" | Type == "Sensor")
  Xmin <- min(sub.seedl$X)
  Ymin <- min(sub.seedl$Y)
  sub.seedl$Xadj <- sub.seedl$X - Xmin
  sub.seedl$Yadj <- sub.seedl$Y - Ymin
  # Calculate best-fit transect line with PCA
  pca <- princomp(sub.seedl[, 6:7])
  pc1min <- min(pca$scores[,1])
  # Adjusted PC1 is distance along transect
  sub.seedl$Distpca <- pca$scores[,1] - pc1min
  # PC2 is distance from the transect (Lateral)
  sub.seedl$Lateral <- pca$scores[,2]
  # Apply model to tree locations to map trees onto transect
  sub.tran <- subset(sub, Type != "Seedling" & Type != "Sensor")
  sub.tran$Xadj <- sub.tran$X - Xmin
  sub.tran$Yadj <- sub.tran$Y - Ymin
  pred <- predict(pca, sub.tran[, 6:7])
  sub.tran$Distpca <- pred[,1] - pc1min
  sub.tran$Lateral <- pred[,2]
  # Reassemble
  sub <- rbind(sub.seedl, sub.tran)
  if (i == 1) {
    newjoint <- sub
  }
  else {
    newjoint <- rbind(newjoint, sub)
  } 
}
# Problems with Mbuzi Mawe - remove anything that is more than 100 m
# from original transect start - start of belt transect was apparently misplaced
newjoint <- subset(newjoint, Distpca > -100)
# Use segmented regression to identify tree cover breakpoint
newjoint2 <- subset(newjoint, Type == "Seedling" | Type == "Sensor")
tcid <- tc[,c(1,6)]
tcid <- merge(tcid, newjoint2, by = "ID")
# Flip distance direction in Makoma, MM, Simiyu and Togoro
# so that grassland-woodland transition happens as distance increases
distmax <- aggregate(Distpca ~ Site, newjoint2, max)
names(distmax)[2] <- "Dmax"
tcid <- merge(tcid, distmax)
tcid$Distpca <- ifelse(tcid$Site == "Makoma" | tcid$Site == "Mbuzi Mawe" | tcid$Site == "Simiyu" | tcid$Site == "Togoro",
                       tcid$Dmax - tcid$Distpca, tcid$Distpca)
tcid$Dmax <- NULL
# Find breakpoints (best-fit distance giving a transition from grassland to woodland)
brkdf <- data.frame(tran)
names(brkdf) <- "Site"
brkdf$Brkdist <- 0
for (i in 1:length(tran)){
  sub <- subset(tcid, Site == tran[i])
  sub <- sub[order(sub$Distpca), ]
  sub$ssq <- numeric(nrow(sub))
  for (j in 1:nrow(sub)){
    brk <- sub$Distpca[j]
    val <- ifelse(sub$Distpca <= brk, mean(sub$TC30[1:j]), mean(sub$TC30[(j + 1):nrow(sub)]))
    sub$ssq[j] <- sum((val - sub$TC30) ^ 2)
  }
  brkdf$Brkdist[i] <- sub$Distpca[which.min(sub$ssq)]
}
# Add Brkdist to tc
tcid <- merge(tcid, brkdf)
tcid$Habitat <- ifelse(tcid$Distpca <= tcid$Brkdist, "Open", "Woody")

# Transect plot for appendix
sitecode <- c("FS", "IK", "MA", "MM", "SI", "SO", "TA", "TO")
tcid$Transect <- factor(sitecode[as.numeric(as.factor(tcid$Site))])
figB1 <- ggplot(tcid) + geom_point(aes(x = Distpca, y = TC30, fill = Habitat), 
                                   shape = 21, size = 4) +
  facet_wrap(~ Transect) +
  labs(x = 'Distance (m)', 
       y = expression(paste('TC'[30]))) +
  scale_fill_hue(l=40) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(vjust = 0, size = rel(1.2)),
        axis.text = element_text(size = rel(1.2)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)))

tc <- tcid
tc <- tc[, -c(7,8,10,11,13)]
# Now apply breakpoint to transect data
transect <- subset(newjoint, Type == "Tree.tran" | Type == "Seedl.tran")
transect <- merge(transect, brkdf)
transect <- merge(transect, distmax)
transect$Distpca <- ifelse(transect$Site == "Makoma" | transect$Site == "Mbuzi Mawe" | transect$Site == "Simiyu" | transect$Site == "Togoro",
                           transect$Dmax - transect$Distpca, transect$Distpca)
transect$Habitat <- ifelse(transect$Distpca <= transect$Brkdist, "Open", "Woody")
transect <- transect[, -c(4,5,6,7,11)]
tree <- merge(tree, subset(transect, Type == "Tree.tran"), by = c("Site", "ID"))
seed <- merge(seed, subset(transect, Type == "Seedl.tran"), by = c("Site", "ID"))
# Remove temporary dataframes
rm(joint, newjoint, newjoint2, brkdf, tree2, seed2, sub, sub.tran, tc2)

# Calculate tree basal area and seedling density in woody and open habitats
# Get tree basal area data along transects
basdist <- aggregate(cbind(Distpca, Brkdist) ~ Site, tree, max)
bassum <- aggregate(Bas ~ Site + Habitat, tree, sum)
bassum <- merge(bassum, basdist)
# Get max tree height by habitat
basht <- aggregate(Ht ~ Site + Habitat, tree, max)
bassum <- merge(bassum, basht)
bassum$Area <- ifelse(bassum$Habitat == 'Open', bassum$Brkdist * 20,
                      (bassum$Distpca - bassum$Brkdist) * 20)
bassum$BA <- bassum$Bas / bassum$Area # Basal area in m2/ha
bassum <- bassum[, -c(3,4,5,7)]

# Compare tree cover with basal area
tchab <- aggregate(TC30 ~ Site + Habitat, tc, mean)
tchab <- merge(tchab, bassum)

# Seedling density per habitat type
seeddist <- aggregate(Distpca ~ Site, seed, max)
# Maximum transect length is max distance from tree or seedling dataset
seeddist <- merge(seeddist, basdist, by = 'Site')
seeddist$Distpca <- ifelse(seeddist$Distpca.x > seeddist$Distpca.y, seeddist$Distpca.x, seeddist$Distpca.y)
seeddist <- seeddist[, -(2:3)]
seedsum <- aggregate(ID ~ Site + Habitat, subset(seed, Seedl_Respr == 'S'), function(x){length(x)})
seedsum <- merge(seedsum, seeddist)
seedsum$Area <- ifelse(seedsum$Habitat == 'Open', seedsum$Brkdist * 3,
                      (seedsum$Distpca - seedsum$Brkdist) * 3)
names(seedsum)[3] <- 'N'
seedsum$Density <- seedsum$N / seedsum$Area * 1e4 # Density in seedlings / ha

# Plot for appendix figures
# Tree basal area and seedling density across habitats
# Fill in 0's for seedsum - integrate with bassum
seedsum2 <- seedsum[, c(1,2,7)]
bassum <- merge(bassum, seedsum2, all = TRUE)
bassum$Density <- ifelse(is.na(bassum$Density), 0, bassum$Density)
bassum$BA <- ifelse(is.na(bassum$BA), 0, bassum$BA)
fig2a <- ggplot(bassum, aes(x = Habitat, y = BA)) +
  geom_line(aes(group = Site), lwd=1.2) + 
  labs(x = 'Habitat', 
       y = expression(paste('Basal area (m'^2,' ha'^-1,')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1.2)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm")) +
  annotate("text", x = 1.2, y = 6, label = "P = 0.008", parse = FALSE, size = 4)

fig2b <- ggplot(bassum, aes(x = Habitat, y = Density)) +
  geom_line(aes(group = Site), lwd=1.2) + 
  labs(x = 'Habitat', 
       y = expression(paste('Seedlings ha'^-1))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1.2)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm")) +
  annotate("text", x = 2, y = 55, label = "P = 0.81", parse = FALSE, size = 4)

fig2c <- ggplot(bassum, aes(x = Habitat, y = Ht)) +
  geom_line(aes(group = Site), lwd = 1.2) + 
  ylim(0, 12) +
  labs(x = 'Habitat', 
       y = 'Max. height (m)') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1.2)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm")) +
  annotate("text", x = 1, y = 11.5, label = "P = 0.09", parse = FALSE, size = 4)

fig2 <- plot_grid(fig2a, fig2b, fig2c, labels = "AUTO", align = "v", label_size = 20,
                   nrow = 1, hjust = 0, vjust = 1, scale = 0.9)

# Stats
bassum <- bassum[order(bassum$Habitat, bassum$Site), ]
BA.test <- t.test(bassum$BA[1:8], bassum$BA[9:16], alternative = "less",
                    paired = TRUE)
Density.test <- t.test(bassum$Density[1:8], bassum$Density[9:16], alternative = "less",
                  paired = TRUE)
Height.test <- t.test(bassum$Ht[1:8], bassum$Ht[9:16], alternative = "less",
                       paired = TRUE)

# Rough calculation of crown area to be expected from current seedlings/gullivers
# moving into the adult class
Seed.mean <- mean(seedsum$Density)
Cr.area.proj <- Seed.mean * Cr.mean # Mean crown projected area (no overlap)
Cr.area.proj / 1e4
# Mean projected crown area for all escaped seedlings is about 0.034

# TC30 as a function of basal area
figC1 <- ggplot(tchab) + geom_point(aes(x = BA, y = TC30, fill = Habitat), 
                                    shape = 21, size = 6) +
  labs(x = expression(paste('Basal area (m'^2,' ha'^-1,')')), 
       y = expression(paste('TC'[30]))) +
  scale_fill_hue(l=40) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(2)),
        axis.title.x = element_text(vjust = 0, size = rel(2)),
        axis.text = element_text(size = rel(2)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)))
lpred <- summary(lm(TC30 ~ BA, tchab))
figC1 <- figC1 + geom_abline(slope = coef(lpred)[2], intercept = coef(lpred)[1], lwd = 1.2) +
  annotate("text", x = 4, y = 0.2, label = "P = 0.0016", parse = FALSE, size = 8)

# Stats
TC30.BA.mod <- summary(lme(TC30 ~ BA, data = tchab, random = ~ 1 | Site))
anova(TC30.BA.mod)
summary(lm(predict(TC30.BA.mod) ~ tchab$TC30)) # R-squared value
summary(lm(TC30 ~ BA, data = tchab))

# Test also for a relationship between tree basal area and seedling density
tchab <- merge(tchab, seedsum, all = TRUE)
Den.BA.mod <- summary(lme(Density ~ BA, data = tchab,
                          na.action = 'na.exclude', random = ~ 1 | Site))
anova(Den.BA.mod)

#--------------------------------------------------------------------------
# Combine demography data with tc data
df.ag <- merge(df.ag, tc, by.x = 'Subplot', by.y = 'ID')
# Make binary version of fire/herbivory variables
df.ag$F2019.bin <- ifelse(df.ag$F2019 > 0, 1, 0)
df.ag$H2019.bin <- ifelse(df.ag$H2019 > 0, 1, 0)

# Test for relationship between tree cover and damage/mortality/growth
# Use distance as autoregressive covariate
# Problem with spatial data: two identical spatial coordinates - jitter
duplicated(df.ag[, 14:15])
df.ag$X[398] <- df.ag$X[398] + 1
df.ag$X[694] <- df.ag$X[694] + 1

# Get summary data for growth, topkill etc.
df.ag.tr <- aggregate(cbind(Topkilled, Damage, F2019.bin, H2019.bin, Dead, Grass, deltaH, deltaH.undam) ~
                         Site, df.ag, function(x){round(mean(x), 3)})
df.ag.mn <- aggregate(cbind(Topkilled, Damage, F2019.bin, H2019.bin, Dead, Grass, deltaH, deltaH.undam) ~
                        1, df.ag.tr, function(x){round(mean(x), 3)})
df.ag.se <- aggregate(cbind(Topkilled, Damage, F2019.bin, H2019.bin, Dead, Grass, deltaH, deltaH.undam) ~
                        1, df.ag.tr, function(x){round(sd(x) / sqrt(8), 3)})
df.ag.min <- aggregate(cbind(Topkilled, Damage, F2019.bin, H2019.bin, Dead, Grass, deltaH, deltaH.undam) ~
                        1, df.ag.tr, function(x){round(min(x), 3)})
df.ag.max <- aggregate(cbind(Topkilled, Damage, F2019.bin, H2019.bin, Dead, Grass, deltaH, deltaH.undam) ~
                         1, df.ag.tr, function(x){round(max(x), 3)})
df.ag.sum <- rbind(df.ag.mn, df.ag.se, df.ag.min, df.ag.max)
df.ag.sum$Metric <- c("Mean", "SE", "Min", "Max")
df.ag.sum <- df.ag.sum[, c(9,7,8,2,3,4,1,5,6)]
df.ag.sum$Grass <- round(df.ag.sum$Grass)
write.csv(df.ag.sum, "Table_S1.csv", row.names = FALSE)

# Obtain mean damage to fire or herbivore damaged seedlings
df.ag.sub.fi <- subset(df.ag, F2019 > 0)
mean(df.ag.sub.fi$F2019)
sd(df.ag.sub.fi$F2019)
min(df.ag.sub.fi$F2019)
max(df.ag.sub.fi$F2019)
df.ag.sub.he <- subset(df.ag, H2019 > 0)
mean(df.ag.sub.he$H2019)
sd(df.ag.sub.he$H2019)
min(df.ag.sub.he$H2019)
max(df.ag.sub.he$H2019)
rm(df.ag.sub.fi, df.ag.sub.he)

# Count number of height reversals per seedling across surveys
df <- df[order(df$Subplot, df$Day), ]
subplots <- unique(df$Subplot)
for (i in 1:length(subplots)){
  sub <- subset(df, Subplot == subplots[i])
  hts <- sub$Ht_m
  hts <- hts[!is.na(hts)]
  htsdf <- data.frame(subplots[i])
  names(htsdf) <- 'ID'
  htsdf$Reversal <- mean((hts[2:length(hts)] - hts[1:(length(hts) - 1)]) < 0)
  if (i == 1) htsdf.all <- htsdf
  else htsdf.all <- rbind(htsdf.all, htsdf)
}
mean(htsdf.all$Reversal)
sd(htsdf.all$Reversal)

df$Tran <- factor(sitecode[as.numeric(as.factor(df$Transect))])
# Get mean day per survey period
perday <- aggregate(Day ~ Per + Tran, df, mean)
df.ag.ht <- aggregate(Ht_m ~ Per + Tran, df, mean)
df.ag.ht <- merge(df.ag.ht, perday)

# Find expected height distribution without disturbance and compare with actual
hts <- df[, c(6,1,11)]
hinit <- subset(hts, Per == 'Initial')
hinit$Per <- NULL
hfinal <- subset(hts, Per == 'Jun 2019')
hfinal$Per <- NULL
htgr <- df.ag[, c(1,9)]
hinit <- merge(hinit, hfinal, by = 'Subplot')
names(hinit)[2:3] <- c('H0', 'Hf')
hinit <- merge(hinit, htgr)
hinit$Hf.pred <- hinit$H0 + hinit$deltaH.undam * 3

fig3 <- ggplot(NULL, aes(x = Day, y = Ht_m)) +
  geom_line(data = df, aes(group = Subplot)) + facet_wrap(~ Tran) +
  geom_line(data = df.ag.ht, col = 'red', lwd = 1.2) +
  labs(x = 'Day', 
       y = 'Height (m)') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"))

# For Gaussian models, fit spatial lme models with site random effects
# Exponential autocorrelation models
fit.deltaH.exp <- lme(deltaH ~ TC30, data = df.ag,
            random = ~ 1 | Site,
            corr = corSpatial(form = ~ X + Y, type ="exponential", nugget = F), 
            method = "ML", na.action = 'na.exclude')
fit.deltaH.undam.exp <- lme(deltaH.undam ~ TC30, data = df.ag,
            random = ~ 1 | Site,
            corr = corSpatial(form = ~ X + Y, type ="exponential", nugget = F), 
            method = "ML", na.action = 'na.exclude')
fit.Grass.exp <- lme(Grass ~ TC30, data = df.ag,
            random = ~ 1 | Site,
            corr = corSpatial(form = ~ X + Y, type ="exponential", nugget = F), 
            method = "ML", na.action = 'na.exclude')
fit.Topkill.exp <- lme(Topkilled ~ TC30, data = df.ag,
            random = ~ 1 | Site,
            corr = corSpatial(form = ~ X + Y, type ="exponential", nugget = F), 
            method = "ML", na.action = 'na.exclude')
fit.Damage.exp <- lme(Damage ~ TC30, data = df.ag,
            random = ~ 1 | Site,
            corr = corSpatial(form = ~ X + Y, type ="exponential", nugget = F), 
            method = "ML", na.action = 'na.exclude')
# Spherical autocorrelation models
fit.deltaH.sph <- lme(deltaH ~ TC30, data = df.ag,
            random = ~ 1 | Site,
            corr = corSpatial(form = ~ X + Y, type ="spherical", nugget = F), 
            method = "ML", na.action = 'na.exclude')
fit.deltaH.undam.sph <- lme(deltaH.undam ~ TC30, data = df.ag,
            random = ~ 1 | Site,
            corr = corSpatial(form = ~ X + Y, type ="spherical", nugget = F), 
            method = "ML", na.action = 'na.exclude')
fit.Grass.sph <- lme(Grass ~ TC30, data = df.ag,
            random = ~ 1 | Site,
            corr = corSpatial(form = ~ X + Y, type ="spherical", nugget = F), 
            method = "ML", na.action = 'na.exclude')
fit.Topkill.sph <- lme(Topkilled ~ TC30, data = df.ag,
            random = ~ 1 | Site,
            corr = corSpatial(form = ~ X + Y, type ="spherical", nugget = F), 
            method = "ML", na.action = 'na.exclude')
fit.Damage.sph <- lme(Damage ~ TC30, data = df.ag,
            random = ~ 1 | Site,
            corr = corSpatial(form = ~ X + Y, type ="spherical", nugget = F), 
            method = "ML", na.action = 'na.exclude')
AIC(fit.deltaH.exp, fit.deltaH.sph)
AIC(fit.deltaH.undam.exp, fit.deltaH.undam.sph)
AIC(fit.Grass.exp, fit.Grass.sph)
AIC(fit.Topkill.exp, fit.Topkill.sph)
AIC(fit.Damage.exp, fit.Damage.sph)

# Little difference - exponential models generally fit better
summary(fit.deltaH.exp)
summary(fit.deltaH.undam.exp)
summary(fit.Grass.exp)
summary(fit.Topkill.exp)
summary(fit.Damage.exp)

# Re-examine growth patterns for more common species
df.acator <- subset(df.ag, Species == 'ACATOR')
df.comafr <- subset(df.ag, Species == 'COMAFR')
fit.deltaH.undam.acator <- lme(deltaH.undam ~ TC30, data = df.acator,
                            random = ~ 1 | Site,
                            corr = corSpatial(form = ~ X + Y, type ="exponential", nugget = F), 
                            method = "ML", na.action = 'na.exclude')
fit.deltaH.undam.comafr <- lme(deltaH.undam ~ TC30, data = df.comafr,
                               random = ~ 1 | Site,
                               corr = corSpatial(form = ~ X + Y, type ="exponential", nugget = F), 
                               method = "ML", na.action = 'na.exclude')
# Results hold for the two most common species

shapes <- c(0, 1, 2, 3, 4, 15, 16, 17)
df.ag$Transect <- factor(sitecode[as.numeric(as.factor(df.ag$Site))])
# Plots
fig4a <- ggplot(df.ag) + geom_point(aes(x = TC30, y = deltaH.undam, shape = Transect), 
                                    size = 2) +
  labs(x = expression(paste('TC'[30])), 
       y = expression(paste('Growth (m y'^-1,')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(0.8)),
        axis.title.x = element_text(vjust = 0, size = rel(0.8)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8))) +
  scale_shape_manual(labels = sitecode, values = shapes) +
  guides(shape = guide_legend(override.aes = list(size=2)), size = "none")
df.ag$Pred <- fitted(fit.deltaH.undam.exp)
lpred <- summary(lm(Pred ~ TC30, df.ag))
fig4a <- fig4a + geom_abline(slope = coef(lpred)[2], intercept = coef(lpred)[1], 
                             col = 'red', size = 1.5) +
  annotate("text", x = 0.4, y = -0.25, label = "P = 0.0002", parse = FALSE, size = 4)

fig4b <- ggplot(df.ag) + geom_point(aes(x = TC30, y = Grass, shape = Transect), 
                                    size = 2) +
  labs(x = expression(paste('TC'[30])), 
       y = expression(paste('Grass biomass (g m'^'-2',')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(0.8)),
        axis.title.x = element_text(vjust = 0, size = rel(0.8)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8))) +
  scale_shape_manual(labels = sitecode, values = shapes) +
  guides(shape = guide_legend(override.aes = list(size=2)), size = "none")
df.ag$Pred <- fitted(fit.Grass.exp)
lpred <- summary(lm(Pred ~ TC30, df.ag))
fig4b <- fig4b + geom_abline(slope = coef(lpred)[2], intercept = coef(lpred)[1], 
                             col = 'red', size = 1.5) +
  annotate("text", x = 0.45, y = 1000, label = "P = 0.22", parse = FALSE, size = 4)

fig4c <- ggplot(df.ag) + geom_point(aes(x = TC30, y = Topkilled, shape = Transect), 
                                    size = 2) +
  labs(x = expression(paste('TC'[30])), 
       y = expression(paste('Topkill (y'^'-1',')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(0.8)),
        axis.title.x = element_text(vjust = 0, size = rel(0.8)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8))) +
  scale_shape_manual(labels = sitecode, values = shapes) +
  guides(shape = guide_legend(override.aes = list(size=2)), size = "none")
df.ag$Pred <- fitted(fit.Topkill.exp)
lpred <- summary(lm(Pred ~ TC30, df.ag))
fig4c <- fig4c + geom_abline(slope = coef(lpred)[2], intercept = coef(lpred)[1], 
                             col = 'red', size = 1.5) +
  annotate("text", x = 0.4, y = 0.9, label = "P = 0.60", parse = FALSE, size = 4)

fig4d <- ggplot(df.ag) + geom_point(aes(x = TC30, y = Damage, shape = Transect), 
                                    size = 2) +
  labs(x = expression(paste('TC'[30])), 
       y = expression(paste('Damage (y'^'-1',')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(0.8)),
        axis.title.x = element_text(vjust = 0, size = rel(0.8)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8))) +
  scale_shape_manual(labels = sitecode, values = shapes) +
  guides(shape = guide_legend(override.aes = list(size=2)), size = "none")
df.ag$Pred <- fitted(fit.Damage.exp)
lpred <- summary(lm(Pred ~ TC30, df.ag))
fig4d <- fig4d + geom_abline(slope = coef(lpred)[2], intercept = coef(lpred)[1], 
                             col = 'red', size = 1.5) +
  annotate("text", x = 0.45, y = 0.9, label = "P = 0.76", parse = FALSE, size = 4)

prow <- plot_grid(
  fig4a + theme(legend.position="none"),
  fig4b + theme(legend.position="none"),
  fig4c + theme(legend.position="none"),
  fig4d + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C", "D"),
  hjust = -1,
  nrow = 2
)
legend <- get_legend(
  # Create some space to the left of the legend
  fig4a + theme(legend.box.margin = margin(0, 0, 0, 12))
)
fig4 <- plot_grid(prow, legend, rel_widths = c(3, .4), scale = 0.9)

# For logistic models, fit glms and autologistic models and obtain slopes by site
Site <- unique(df.ag$Site)
regs <- data.frame(Site)
Sites <- Site
regs$Fi.slope <- NA
regs$Fi.slope.se <- NA
regs$Fi.slope.auto <- NA
regs$Fi.slope.se.auto <- NA
regs$Fi.slope.pval.auto <- NA
regs$He.slope <- NA
regs$He.slope.se <- NA
regs$He.slope.auto <- NA
regs$He.slope.se.auto <- NA
regs$He.slope.pval.auto <- NA
regs$Mo.slope <- NA
regs$Mo.slope.se <- NA
regs$Mo.slope.auto <- NA
regs$Mo.slope.se.auto <- NA
regs$Mo.slope.pval.auto <- NA
df.ag$Dist2 <- rnorm(nrow(df.ag), df.ag$Distpca, 0.01)
for (i in 1:length(Site)){
  sub <- subset(df.ag, Site == Sites[i])
  sub <- subset(sub, !is.na(F2019))
  # Conduct both logistic and autologistic regressions
  coordinates(sub) <- ~ X + Y
  if (sum(sub$F2019.bin) > 0){
    # Fire logistic
    mod.Fi <- glm(F2019.bin ~ TC30, sub, family = binomial)
    regs$Fi.slope[i] <- coef(summary(mod.Fi))["TC30", "Estimate"]
    regs$Fi.slope.se[i] <- coef(summary(mod.Fi))["TC30", "Std. Error"]
    sub$Pred.Fi <- fitted(mod.Fi)
    # Fire autologistic
    mod.Fi.auto <- logistic.regression(sub, y = 'F2019.bin', x = 'TC30', autologistic=TRUE, 
                                       coords = coordinates(sub))
    regs$Fi.slope.auto[i] <- mod.Fi.auto$coefTable[2, "Coef"]
    regs$Fi.slope.se.auto[i] <- mod.Fi.auto$coefTable[2, "StdError"]
    regs$Fi.slope.pval.auto[i] <- mod.Fi.auto$coefTable[2, "Prob"]
  } else {
    sub$Pred.Fi <- NA
  }
  if (sum(sub$H2019.bin) > 0){
    # Herbivory logistic
    mod.He <- glm(H2019.bin ~ TC30, sub, family = binomial)
    regs$He.slope[i] <- coef(summary(mod.He))["TC30", "Estimate"]
    regs$He.slope.se[i] <- coef(summary(mod.He))["TC30", "Std. Error"]
    sub$Pred.He <- fitted(mod.He)
    # Herbivory autologistic
    mod.He.auto <- logistic.regression(sub, y = 'H2019.bin', x = 'TC30', autologistic=TRUE, 
                                       coords = coordinates(sub))
    regs$He.slope.auto[i] <- mod.He.auto$coefTable[2, "Coef"]
    regs$He.slope.se.auto[i] <- mod.He.auto$coefTable[2, "StdError"]
    regs$He.slope.pval.auto[i] <- mod.He.auto$coefTable[2, "Prob"]
  } else {
    sub$Pred.He <- NA
  }
  if (sum(sub$Dead) > 0){
    # Mortality logistic
    mod.Mo <- glm(Dead ~ TC30, sub, family = binomial)
    regs$Mo.slope[i] <- coef(summary(mod.Mo))["TC30", "Estimate"]
    regs$Mo.slope.se[i] <- coef(summary(mod.Mo))["TC30", "Std. Error"]
    sub$Pred.Mo <- fitted(mod.Mo)
    # Mortality autologistic
    mod.Mo.auto <- logistic.regression(sub, y = 'Dead', x = 'TC30', autologistic=TRUE, 
                                       coords = coordinates(sub))
    regs$Mo.slope.auto[i] <- mod.Mo.auto$coefTable[2, "Coef"]
    regs$Mo.slope.se.auto[i] <- mod.Mo.auto$coefTable[2, "StdError"]
    regs$Mo.slope.pval.auto[i] <- mod.Mo.auto$coefTable[2, "Prob"]
  } else {
    sub$Pred.Mo <- NA
  }
  if (i == 1) df.ag2 <- sub
  else df.ag2 <- rbind(df.ag2, sub)
}
df.ag2$Dist2 <- NULL
summary(lm(Fi.slope ~ 1, regs))
summary(lm(Fi.slope.auto ~ 1, regs))
summary(lm(He.slope ~ 1, regs))
summary(lm(He.slope.auto ~ 1, regs))
summary(lm(Mo.slope ~ 1, regs))
summary(lm(Mo.slope.auto ~ 1, regs))

tableE1 <- regs
tableE1[, 2:16] <- round(tableE1[, 2:16], 3)
write.csv(tableE1, "Table_E1.csv", row.names = FALSE)

# Mortality declines with tree cover
# Is mortality related to damage?
# Simplest analysis: break down dead and alive seedlings by transect
# and damage proportion - use paired t-test
df.mort <- aggregate(Damage ~ Dead + Site, df.ag, mean)
df.mort <- df.mort[order(df.mort$Dead, df.mort$Site), ]
Mort.test <- t.test(df.mort$Damage[1:8], df.mort$Damage[9:16], alternative = "less",
                      paired = TRUE)
# Plot
df.mort$Dead <- ifelse(df.mort$Dead == 0, 'Alive', 'Dead')
fig6 <- ggplot(df.mort, aes(x = Dead, y = Damage)) +
  geom_line(aes(group = Site), lwd = 2) + 
  labs(x = '', 
       y = expression(paste('Damage (y'^'-1',')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(2)),
        axis.title.x = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm")) +
  annotate("text", x = 1.5, y = 0.9, label = "P = 0.0002", parse = FALSE, size = 8)

# Fit GLMM models for fire, herbivory and mortality
fire.glm <- glm(F2019.bin ~ TC30, df.ag2, na.action = "na.exclude", family = binomial)
herb.glm <- glm(H2019.bin ~ TC30, df.ag2, na.action = "na.exclude", family = binomial)
mort.glm <- glm(Dead ~ TC30, df.ag2, na.action = "na.exclude", family = binomial)
# These models ignore fine-scale spatial autocorrelation and site
# random effects and are therefore approximations
# They serve to show aggregate trends across transects

# Plots
df.ag2$Transect <- factor(sitecode[as.numeric(as.factor(df.ag2$Site))])
df.ag2 <- as.data.frame(df.ag2)
# Add overall predictions to dataframe
df.ag2$Pred.Fi.all <- fitted(fire.glm)
df.ag2$Pred.He.all <- fitted(herb.glm)
df.ag2$Pred.Mo.all <- fitted(mort.glm)
fig5a <- ggplot(data = df.ag2) + geom_line(aes(x = TC30, y = Pred.Fi, color = Transect), 
                                    size = 1.5) +
  geom_line(aes(x = TC30, y = Pred.Fi.all), size = 1.5) +
  labs(x = expression(paste('TC'[30])), 
       y = expression(paste('2019 fire damage (y'^'-1',')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(0.8)),
        axis.title.x = element_text(vjust = 0, size = rel(0.8)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8))) +
  scale_shape_manual(labels = sitecode, values = shapes) +
  guides(color = guide_legend(override.aes = list(size = 1.5)), size = "none") +
  annotate("text", x = 0.15, y = 0.95, label = "P = 0.37", parse = FALSE, size = 4)

fig5b <- ggplot(data = df.ag2) + geom_line(aes(x = TC30, y = Pred.He, color = Transect), 
                                     size = 1.5) +
  geom_line(aes(x = TC30, y = Pred.He.all), size = 1.5) +
  labs(x = expression(paste('TC'[30])), 
       y = expression(paste('2019 herbivore damage (y'^'-1',')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(0.8)),
        axis.title.x = element_text(vjust = 0, size = rel(0.8)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8))) +
  scale_shape_manual(labels = sitecode, values = shapes) +
  guides(shape = guide_legend(override.aes = list(size = 1.5)), size = "none") +
  annotate("text", x = 0.35, y = 0.95, label = "P = 0.32", parse = FALSE, size = 4)

fig5c <- ggplot(data = df.ag2) + geom_line(aes(x = TC30, y = Pred.Mo, color = Transect), 
                                     size = 1.5) +
  geom_line(aes(x = TC30, y = Pred.Mo.all), size = 1.5) +
  labs(x = expression(paste('TC'[30])), 
       y = expression(paste('2017-19 mortality'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(0.8)),
        axis.title.x = element_text(vjust = 0, size = rel(0.8)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8))) +
  scale_shape_manual(labels = sitecode, values = shapes) +
  guides(shape = guide_legend(override.aes = list(size = 1.5)), size = "none") +
  annotate("text", x = 0.35, y = 0.4, label = "P = 0.017", parse = FALSE, size = 4)

prow <- plot_grid(
  fig5a + theme(legend.position="none"),
  fig5b + theme(legend.position="none"),
  fig5c + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C"),
  hjust = -1,
  nrow = 1
)
legend <- get_legend(
  # create some space to the left of the legend
  fig5a + theme(legend.box.margin = margin(0, 0, 0, 12))
)
fig5 <- plot_grid(prow, legend, rel_widths = c(3, .4))

# Examine microclimate data in relation to tree cover
micro <- read.csv('Daily_microclimate_data.csv')
micro.ag <- aggregate(cbind(VWCmean, Tmean, Tmax, Tmin, PARmean, PARmax) ~ ID, micro, mean)
micro.ag <- merge(micro.ag, tc, by = 'ID')

# Site random effects - intercept
mod.PARmean <- lme(PARmean ~ TC30, data = micro.ag, control = list(opt = "optim"), 
                 random = ~ 1 | Site)
mod.PARmax <- lme(PARmax ~ TC30, data = micro.ag, control = list(opt = "optim"), 
                   random = ~ 1 | Site)
mod.VWCmean <- lme(VWCmean ~ TC30, data = micro.ag, control = list(opt = "optim"), 
                   random = ~ 1 | Site)
mod.Tmean <- lme(Tmean ~ TC30, data = micro.ag, control = list(opt = "optim"), 
                   random = ~ 1 | Site)
mod.Tmax <- lme(Tmax ~ TC30, data = micro.ag, control = list(opt = "optim"), 
                 random = ~ 1 | Site)
mod.Tmin <- lme(Tmin ~ TC30, data = micro.ag, control = list(opt = "optim"), 
                 random = ~ 1 | Site)
anova(mod.PARmean)
anova(mod.PARmax)
anova(mod.VWCmean)
anova(mod.Tmean)
anova(mod.Tmax)
anova(mod.Tmin)

# Figs
micro.ag$Transect <- factor(sitecode[as.numeric(as.factor(micro.ag$Site))])
fig7a <- ggplot(micro.ag) + geom_point(aes(x = TC30, y = PARmean, shape = Transect), size = 2) +
  labs(x = expression(paste('TC'[30])), 
       y = expression(paste('Mean PAR (', mu, 'mol m'^-2,' s'^-1,')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(0.8)),
        axis.title.x = element_text(vjust = 0, size = rel(0.8)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8))) +
  scale_shape_manual(labels = sitecode, values = shapes) +
  guides(shape = guide_legend(override.aes = list(size=2)), size = "none")
micro.ag$Pred <- fitted(mod.PARmean)
lpred <- summary(lm(Pred ~ TC30, micro.ag))
fig7a <- fig7a + geom_abline(slope = coef(lpred)[2], intercept = coef(lpred)[1], 
                             col = 'red', size = 1.5) +
  annotate("text", x = 0.3, y = 900, label = "P = 0.039", parse = FALSE, size = 4)

fig7b <- ggplot(micro.ag) + geom_point(aes(x = TC30, y = PARmax, shape = Transect), size = 2) +
  labs(x = expression(paste('TC'[30])), 
       y = expression(paste('Max. PAR (', mu, 'mol m'^-2,' s'^-1,')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(0.8)),
        axis.title.x = element_text(vjust = 0, size = rel(0.8)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8))) +
  scale_shape_manual(labels = sitecode, values = shapes) +
  guides(shape = guide_legend(override.aes = list(size=2)), size = "none")
micro.ag$Pred <- fitted(mod.PARmax)
lpred <- summary(lm(Pred ~ TC30, micro.ag))
fig7b <- fig7b + geom_abline(slope = coef(lpred)[2], intercept = coef(lpred)[1], 
                             col = 'red', size = 1.5) +
  annotate("text", x = 0.3, y = 1800, label = "P = 0.071", parse = FALSE, size = 4)

fig7c <- ggplot(micro.ag) + geom_point(aes(x = TC30, y = VWCmean, shape = Transect), size = 2) +
  labs(x = expression(paste('TC'[30])), 
       y = expression(paste('Mean VWC (cm'^3,' cm'^-3,')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(0.8)),
        axis.title.x = element_text(vjust = 0, size = rel(0.8)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8))) +
  scale_shape_manual(labels = sitecode, values = shapes) +
  guides(shape = guide_legend(override.aes = list(size=2)), size = "none")
micro.ag$Pred <- fitted(mod.VWCmean)
lpred <- summary(lm(Pred ~ TC30, micro.ag))
fig7c <- fig7c + geom_abline(slope = coef(lpred)[2], intercept = coef(lpred)[1], 
                             col = 'red', size = 1.5) +
  annotate("text", x = 0.3, y = 40, label = "P = 0.0069", parse = FALSE, size = 4)

fig7d <- ggplot(micro.ag) + geom_point(aes(x = TC30, y = Tmean, shape = Transect), size = 2) +
  labs(x = expression(paste('TC'[30])), 
       y = expression(paste('Mean T (', degree, 'C)'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(0.8)),
        axis.title.x = element_text(vjust = 0, size = rel(0.8)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8))) +
  scale_shape_manual(labels = sitecode, values = shapes) +
  guides(shape = guide_legend(override.aes = list(size=2)), size = "none")
micro.ag$Pred <- fitted(mod.Tmean)
lpred <- summary(lm(Pred ~ TC30, micro.ag))
fig7d <- fig7d + geom_abline(slope = coef(lpred)[2], intercept = coef(lpred)[1], 
                             col = 'red', size = 1.5) +
  annotate("text", x = 0.3, y = 45, label = "P = 0.58", parse = FALSE, size = 4)


fig7e <- ggplot(micro.ag) + geom_point(aes(x = TC30, y = Tmin, shape = Transect), size = 2) +
  labs(x = expression(paste('TC'[30])), 
       y = expression(paste('Min. T (', degree, 'C)'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(0.8)),
        axis.title.x = element_text(vjust = 0, size = rel(0.8)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8))) +
  scale_shape_manual(labels = sitecode, values = shapes) +
  guides(shape = guide_legend(override.aes = list(size=2)), size = "none")
micro.ag$Pred <- fitted(mod.Tmin)
lpred <- summary(lm(Pred ~ TC30, micro.ag))
fig7e <- fig7e + geom_abline(slope = coef(lpred)[2], intercept = coef(lpred)[1], 
                             col = 'red', size = 1.5) +
  annotate("text", x = 0.3, y = 35, label = "P = 0.25", parse = FALSE, size = 4)

fig7f <- ggplot(micro.ag) + geom_point(aes(x = TC30, y = Tmax, shape = Transect), size = 2) +
  labs(x = expression(paste('TC'[30])), 
       y = expression(paste('Max. T (', degree, 'C)'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(0.8)),
        axis.title.x = element_text(vjust = 0, size = rel(0.8)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8))) +
  scale_shape_manual(labels = sitecode, values = shapes) +
  guides(shape = guide_legend(override.aes = list(size=2)), size = "none")
micro.ag$Pred <- fitted(mod.Tmax)
lpred <- summary(lm(Pred ~ TC30, micro.ag))
fig7f <- fig7f + geom_abline(slope = coef(lpred)[2], intercept = coef(lpred)[1], 
                             col = 'red', size = 1.5) +
  annotate("text", x = 0.3, y = 62, label = "P = 0.86", parse = FALSE, size = 4)

prow <- plot_grid(
  fig7a + theme(legend.position="none"),
  fig7b + theme(legend.position="none"),
  fig7c + theme(legend.position="none"),
  fig7d + theme(legend.position="none"),
  fig7e + theme(legend.position="none"),
  fig7f + theme(legend.position="none"),
  align = 'vh',
  labels = 'AUTO',
  hjust = -1,
  nrow = 2
)
legend <- get_legend(
  # create some space to the left of the legend
  fig7a + theme(legend.box.margin = margin(0, 0, 0, 12))
)
fig7 <- plot_grid(prow, legend, rel_widths = c(3, .4))

# Simulation model
# First use the original dataset
# Deterministic growth and stochastic disturbance
# Recalculate parameters using linear models
# Growth and mortality are related to TC
# Bootstrapped parameter estimates
df.sim <- df.ag[, c(2,3,6,8,9,12,18,19)]
ITER <- 100 # Bootstrapped iterations
# Bootstrapped coefficients
# Growth and mortality
gr <- matrix(nrow = ITER, ncol = 2)
gr.u <- matrix(nrow = ITER, ncol = 2)
mo <- matrix(nrow = ITER, ncol = 3)
gr.alt <- numeric(ITER) # Alternative parameterization where growth is unrelated to TC
gr.u.alt <- numeric(ITER)
mo.alt <- numeric(ITER)
# Damage and topkill
he <- numeric(ITER)
fi <- numeric(ITER)
tk <- numeric(ITER)
he.alt <- matrix(nrow = ITER, ncol = 2)
fi.alt <- matrix(nrow = ITER, ncol = 2)
tk.alt <- matrix(nrow = ITER, ncol = 2)
for (i in 1:ITER){
  df.boot <- df.sim[sample(nrow(df.sim), replace = TRUE), ]
  # Growth and mortality vary with TC30
  sim.gr <- lm(deltaH ~ TC30, df.boot, na.action = "na.omit")
  sim.gr.undam <- lm(deltaH.undam ~ TC30, df.boot, na.action = "na.omit")
  sim.mreg <- glm(Dead ~ TC30 + Damage, df.boot, na.action = "na.omit", family = binomial)
  gr[i, ] <- coef(sim.gr)
  gr.u[i, ] <- coef(sim.gr.undam)
  mo[i, ] <- coef(sim.mreg)
  # Disturbance does not
  he[i] <- mean(df.boot$H2019.bin, na.rm = TRUE)
  fi[i] <- mean(df.boot$F2019.bin, na.rm = TRUE)
  tk[i] <- mean(df.boot$Topkilled)
  # Now allow disturbance/topkill to vary with TC30, but not growth
  # Growth and mortality
  gr.alt[i] <- mean(df.boot$deltaH, na.rm = TRUE)
  gr.u.alt[i] <- mean(df.boot$deltaH.undam, na.rm = TRUE)
  mo.alt[i] <- mean(df.boot$Dead, na.rm = TRUE)
  # Disturbance
  hreg <- glm(H2019.bin ~ TC30, df.boot, family = binomial)
  he.alt[i, ] <- coef(hreg)
  freg <- glm(F2019.bin ~ TC30, df.boot, family = binomial)
  fi.alt[i, ] <- coef(freg)
  treg <- lm(Topkilled ~ TC30, df.boot)
  tk.alt[i, ] <- coef(treg)
}

N <- 1e4 # 1000 individuals per 10 tree cover bins
f <- c(0.1, 0.5, 1, 1.5) # Proportional disturbance reduction
# Calculate bootstrapped escape proportions
# Scenario 1: Growth varies across habitats, disturbance does not
# This is analogous to the default scenario supported by the analysis
for (s in 1:4){
  for (i in 1:ITER){
    H <- numeric(N) 
    sim <- data.frame(H)
    sim$TC <- rep(seq(from = 0, to = 0.45, by = 0.05), each = 1000)
    sim$Dam <- 0
    sim$Esc <- 0
    # Cycle through 500 time steps
    for (t in 1:1000){
      sim$H <- ifelse(sim$Dam == 0, sim$H + gr.u[i,1] + gr.u[i,2] * sim$TC,
                      sim$H + gr[i,1] + gr[i,2] * sim$TC)
      sim$Dam <- rbinom(N, size = 1, prob = he[i] * f[s]) + 
        rbinom(N, size = 1, prob = fi[i] * f[s])
      sim$Dam <- ifelse(sim$Dam == 2, 1, sim$Dam)
      topkilled <- rbinom(N, size = 1, prob = tk[i] * f[s])
      sim$H <- ifelse(topkilled == 1 & sim$Esc == 0, 0, sim$H)
      # Probability of mortality (depends on TC and damage)
      lp.m <- mo[i,1] + mo[i,2] * sim$TC + mo[i,3] * sim$Dam # Linear predictor for inverse logit
      prob.m.3y <- exp(lp.m) / (exp(lp.m) + 1) # Three-year probability of mortality - annualize
      prob.m.1y <- 1 - (1 - prob.m.3y) ^ (1/3)
      Dead <- rbinom(N, size = 1, prob = prob.m.1y)
      # Assume dead trees are removed and replaced with new seedlings on height 0
      sim$H <- ifelse(Dead == 1 & sim$Esc == 0, 0, sim$H)
      sim$Esc <- ifelse(sim$H > 2, 1, 0)
    }
    # Calculate proportion escaped for each TC class
    sim.ag <- aggregate(Esc ~ TC, sim, mean)
    sim.ag$f <- f[s]
    sim.ag$ITER <- i
    if (s == 1 & i == 1) sim1 <- sim.ag
    else sim1 <- rbind(sim1, sim.ag)
    cat(s, i, '\n')
  }
}
sim1$Scen <- "Scenario 1"

# Scenario 2: Growth is spatial, all else nonspatial
for (s in 1:4){
  for (i in 1:ITER){
    H <- numeric(N) 
    sim <- data.frame(H)
    sim$TC <- rep(seq(from = 0, to = 0.45, by = 0.05), each = 1000)
    sim$Dam <- 0
    sim$Esc <- 0
    # Cycle through 100 time steps
    for (t in 1:1000){
      sim$H <- ifelse(sim$Dam == 0, sim$H + gr.u[i,1] + gr.u[i,2] * sim$TC,
                      sim$H + gr[i,1] + gr[i,2] * sim$TC)
      sim$Dam <- rbinom(N, size = 1, prob = he[i] * f[s]) + 
        rbinom(N, size = 1, prob = fi[i] * f[s])
      sim$Dam <- ifelse(sim$Dam == 2, 1, sim$Dam)
      topkilled <- rbinom(N, size = 1, prob = tk[i] * f[s])
      sim$H <- ifelse(topkilled == 1 & sim$Esc == 0, 0, sim$H)
      # Probability of mortality is fixed
      prob.m.3y <- mo.alt[i] # Three-year probability of mortality - annualize
      prob.m.1y <- 1 - (1 - prob.m.3y) ^ (1/3)
      Dead <- rbinom(N, size = 1, prob = prob.m.1y)
      # Assume dead trees are removed and replaced with new seedlings on height 0
      sim$H <- ifelse(Dead == 1 & sim$Esc == 0, 0, sim$H)
      sim$Esc <- ifelse(sim$H > 2, 1, 0)
    }
    # Calculate proportion escaped for each TC class
    sim.ag <- aggregate(Esc ~ TC, sim, mean)
    sim.ag$f <- f[s]
    sim.ag$ITER <- i
    if (s == 1 & i == 1) sim2 <- sim.ag
    else sim2 <- rbind(sim2, sim.ag)
    cat(s, i, '\n')
  }
}
sim2$Scen <- "Scenario 2"

# Scenario 3: Growth is non-spatial, all else is spatial
for (s in 1:4){
  for (i in 1:ITER){
    H <- numeric(N) 
    sim <- data.frame(H)
    sim$TC <- rep(seq(from = 0, to = 0.45, by = 0.05), each = 1000)
    sim$Dam <- 0
    sim$Esc <- 0
    # Cycle through 100 time steps
    for (t in 1:100){
      sim$H <- ifelse(sim$Dam == 0, sim$H + gr.u.alt[i],
                      sim$H + gr.alt[i])
      # Probability of herbivory
      lp.h <- he.alt[i,1] + he.alt[i,2] * sim$TC # Linear predictor for inverse logit
      prob.h <- exp(lp.h) / (exp(lp.h) + 1) * f[s]
      # Probability of fire
      lp.f <- fi.alt[i,1] + fi.alt[i,2] * sim$TC # Linear predictor for inverse logit
      prob.f <- exp(lp.f) / (exp(lp.f) + 1) * f[s]
      # Probability of damage
      sim$Dam <- rbinom(N, size = 1, prob.h) + rbinom(N, size = 1, prob.f)
      sim$Dam <- ifelse(sim$Dam == 2, 1, sim$Dam)
      # Probability pf topkill
      prob.t <- tk.alt[i,1] + tk.alt[i,2] * sim$TC
      topkilled <- rbinom(N, size = 1, prob = prob.t * f[s])
      sim$H <- ifelse(topkilled == 1 & sim$Esc == 0, 0, sim$H)
      # Probability of mortality (depends on TC and damage)
      lp.m <- mo[i,1] + mo[i,2] * sim$TC + mo[i,3] * sim$Dam # Linear predictor for inverse logit
      prob.m.3y <- exp(lp.m) / (exp(lp.m) + 1) # Three-year probability of mortality - annualize
      prob.m.1y <- 1 - (1 - prob.m.3y) ^ (1/3)
      Dead <- rbinom(N, size = 1, prob = prob.m.1y)
      # Assume dead trees are removed and replaced with new seedlings on height 0
      sim$H <- ifelse(Dead == 1 & sim$Esc == 0, 0, sim$H)
      sim$Esc <- ifelse(sim$H > 2, 1, 0)
    }
    # Calculate proportion escaped for each TC class
    sim.ag <- aggregate(Esc ~ TC, sim, mean)
    sim.ag$f <- f[s]
    sim.ag$ITER <- i
    if (s == 1 & i == 1) sim3 <- sim.ag
    else sim3 <- rbind(sim3, sim.ag)
    cat(s, i, '\n')
  }
}
sim3$Scen <- "Scenario 3"
sim.all <- rbind(sim1, sim2)
sim.all <- rbind(sim.all, sim3)
sim.all.ag <- aggregate(Esc ~ Scen + f + TC, sim.all, mean)
sim.all.ag$ITER <- 0
sim.all.ag <- sim.all.ag[ ,c(1:3,5,4)]
sim.all <- sim.all[, c(5,3,1,4,2)]
sim.all <- rbind(sim.all.ag, sim.all)

# Save simulation results
write.csv(sim.all, "Simulation_results_final.csv", row.names = FALSE)
sim.all <- read.csv("Simulation_results_final.csv")
sim.all$f <- ifelse(sim.all$f == 0.1, 'f = 0.1',
                    ifelse(sim.all$f == 1, 'f = 1.0', 
                           ifelse(sim.all$f == 1.5, 'F = 1.5', 'f = 0.5')))
sim.all$f <- as.factor(sim.all$f)

fig8 <- ggplot(NULL, aes(x = TC, y = Esc)) + 
  geom_line(data = subset(sim.all, ITER != "0"), aes(group = ITER)) +
  geom_line(data = subset(sim.all, ITER == "0"), col = 'red', lwd = 1.2) + 
  facet_grid(rows = vars(f), cols = vars(Scen)) +
labs(x = expression(paste('TC'[30])), 
     y = 'Escape probability') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(1.5)),
        axis.title.x = element_text(vjust = 0, size = rel(1.5)),
        axis.text = element_text(size = rel(1.2)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"))

# Plot the discontinuity in tree cover across the breakpoint
# Obtain the mean tree cover plus/minus 100 m from each breakpoint
# For each transect, recalculate distance in relation to the breakpoint
tcid$Dist <- tcid$Distpca - tcid$Brkdist
tcid <- tcid[, c(1,3,13,14)]
tcid <- subset(tcid, Dist >= -150 & Dist <= 150)
# Obtain a mean TC30 across all transects
tcid$Dbin30 <- cut(tcid$Dist, breaks = seq(from = -150, to = 150, by = 30),
                  labels = seq(from = -135, to = 135, by = 30))
tcid.ag <- aggregate(TC30 ~ Dbin30, tcid, mean)
tcid.ag.se <- aggregate(TC30 ~ Dbin30, tcid, sd)
names(tcid.ag.se)[2] <- 'TC30.se'
tcid.ag.se$TC30.se <- tcid.ag.se$TC30.se / sqrt(8)
tcid.ag <- merge(tcid.ag, tcid.ag.se, by = 'Dbin30')

# Plot actual TC along transition
fig9 <- ggplot(tcid.ag, aes(x = as.numeric(as.character(Dbin30)), y = TC30)) +
  geom_line(lwd = 2) +
  geom_errorbar(aes(ymin = TC30 - TC30.se,
                    ymax = TC30 + TC30.se),
                width = 5, lwd = 1) +
  labs(x = 'Distance (m)', 
       y = expression(paste('TC'[30]))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(2)),
        axis.title.x = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"))