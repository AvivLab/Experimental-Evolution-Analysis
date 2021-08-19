# Code written on August 2021 by Kenny Ye and Maryl Lambros
# to test for significance of relative growth rate results, 
# assign strains by evolutionary strategy (generalist, 
# specialist, or not significant) based on relative growth
# rate results, and test for significance of evolutionary
# strategy by environmental fluctuation regime

setwd("C:/Users/mmlam/Desktop/BergmanLabRotation/EcoliExperimentalEvolutionPaper/frontiersInMicrobiology/Code")
growth <- read.csv('../Data/growth_strains_new.csv')

growth$Treatment <- factor(growth$Treatment, levels = c("Media", "Random", "Slow", "Fast"))

### different Linage in the same day need to be differentiated
growth$Experiment <- paste(growth$Linage, growth$Experiment, sep = "_")

growth$group <- paste(growth$Strain, growth$Experiment, sep = ":")

library(ggplot2)

### estimate the standard deviation of the residuals at 15
lm.15 <- lm(Max_deriv ~ group, data = growth, subset = Temperature == 15)
s.15 <- summary(lm.15)$sigma

### estimate the standard deviation of the residuals at 43
lm.43 <- lm(Max_deriv ~ group, data = growth, subset = Temperature == 43)
s.43 <- summary(lm.43)$sigma

### variance equalized with power transmation at 1/6,
lm.15a <- lm((Max_deriv)^(1/6) ~ group, data = growth, subset = Temperature == 15)
lm.43a <- lm((Max_deriv)^(1/6) ~ group, data = growth, subset = Temperature == 43)

avg <- aggregate(Max_deriv~group,growth,mean)
avg$reps <- aggregate(Max_deriv~group,growth,length)[,2]
avg <- cbind(avg,growth[match(avg[,"group"],growth$group),c("Treatment","Strain","Experiment","Temperature")])
avg$geometric.mean <- aggregate(Max_deriv~group,growth,function(x) exp(mean(log(x))))[,2]
names(avg)[2] <- "mean.max_deriv"

reference <- subset(avg,Treatment=="Media")
tests <- subset(avg,Treatment!="Media") 


tests <- cbind(tests,reference[match(tests$Exp,reference$Exp),c("mean.max_deriv","reps","geometric.mean")])
### means
names(tests)[c(2,9)] <- c("mu1","mu0")

### number of replicates
names(tests)[c(3,10)] <- c("n1","n0")

### geometric means
names(tests)[c(8,11)] <- c("mu1a","mu0a")

### drop the first column 
tests <- tests[,-1]

### ratio between means
tests$ratio <- tests$mu1/tests$mu0

### ratio between geometric means
tests$ratio.a <- tests$mu1a/tests$mu0a

### compare the ratios from mean and geometric means; they are extremely similar so use the means
plot(tests$ratio,tests$ratio.a);abline(0,1)

### test statistics from pooled variance
### (mu1-mu0)/(sqrt(1/n1+1/n0)*s)

tests$s <- s.15
tests$s[tests$Temperature==43] <- s.43
tests$Zstat <- (tests$mu1-tests$mu0)/(sqrt(1/tests$n1+1/tests$n0)*tests$s)

### p-val from normal approximation
tests$pval <- (1-pnorm(abs(tests$Zstat)))*2

### re-arrange columns and export

tests <- tests[,c(3,4,5,6,1,2,8,9,11,7,10,12,13:15)]
#SAVE
# write.csv(tests,"../Data/tests.csv",row.names=F)

# 
plot(tests$pval[tests$Temperature == 15], tests$ratio[tests$Temperature == 15]);abline(v=0.1, col="red")
# minimum relative growth rates greater than 1.0 that is significant for 15
min15 = min(tests$ratio[tests$pval < 0.07 & tests$Temperature == 15 & tests$ratio>=1.0])
# maximum relative growth rates greater than 1.0 that is NOT significant for 15
max15 = max(tests$ratio[tests$pval >= 0.07 & tests$Temperature == 15])
mean(min15,max15)

# minimum relative growth rates greater than 1.0 that is significant for 43
min43 = min(tests$ratio[tests$pval < 0.07 & tests$Temperature == 43 & tests$ratio>=1.0])
# maximum relative growth rates greater than 1.0 that is NOT significant for 43
max43 = max(tests$ratio[tests$pval >= 0.07 & tests$Temperature == 43])
mean(min43,max43)

temp15 <- subset(tests,Temperature==15)
temp43 <- subset(tests,Temperature==43)

### check that strain matches

sum(temp15$Strain==temp43$Strain)

ratios <- cbind(temp15[,c("Strain","Treatment","ratio")],temp43[,"ratio"])
names(ratios)[3:4] <- c("at15dgr","at43dgr")

# Add Lineage to ratios data frame:
ratios$Lineage <- ifelse(grepl("606", ratios$Strain), "606", "607")

# Assign cut off that is significant 
sig43cutoff = 1.3
sig15cutoff = 1.15
library("ggrepel")
dev.new()
#pdf("RelativeGrowthScatterPlotFor43Vs15_August2021.pdf")
ggplot(ratios,aes(x=at15dgr,y=at43dgr,col=Treatment,shape=Lineage,size=0.7))+geom_point(stroke=1.4)+
  geom_text_repel(aes(label=as.factor(Strain)),size=5)+
  geom_hline(yintercept=c(sig43cutoff),linetype="dotted")+
  geom_vline(xintercept = c(1/sig15cutoff,sig15cutoff),linetype="dotted")+xlim(0.5, 2)+ylim(1,2)+
  geom_vline(xintercept = c(1),linetype="solid")+geom_hline(yintercept = c(1), linetype="solid")+
  scale_color_manual(values=c('orange','#5BBCD6',"mediumaquamarine"))+
  scale_shape_manual(values=c(4, 16))+xlab("Relative Growth at 15C")+ylab("Relative Growth at 43C")
#dev.off()

# Add to ratios dataframe which strains are generalists vs specialist based on results:
ratios$Strategy = ifelse(ratios$at15dgr >= sig15cutoff & ratios$at43dgr >= sig43cutoff, "Generalist", 
                         ifelse(ratios$at15dgr >= sig15cutoff & ratios$at43dgr < sig43cutoff, "Specialist15", 
                         ifelse(ratios$at15dgr < sig15cutoff & ratios$at43dgr >= sig43cutoff, "Specialist43", "Other")))
#write.csv(ratios,"../Data/August2021_updatedStrategyTable.csv")
  
# Fisher Exact Test for Generalists significance for random vs periodic strains:
fisher.test(cbind(c(7,9),c(0,8)))
# significant at p-val = 0.054
# Fisher Exact Test for Specialists significance for random vs periodic strains:
fisher.test(cbind(c(6,10),c(2,6)))
# not significant at p-val = 0.667
# Fisher exact test for 607 generalists vs 606 generalists:
fisher.test(cbind(c(5,3),c(2,6)))
# Fisher exact test for Random 607 specialists to Random 606 specialists
fisher.test(cbind(c(2,2),c(0,4)))
# for 607 vs 606 strains with relative growth increase at 15C:
fisher.test((cbind(c(6,6),c(4,8))))

### test for variance
library(pander)

pander(var.test(log(at15dgr)~I(Treatment=="Random"),data=ratios))

#Remove S606-2 and S607-2 that grow less at 15C than ancestor significantly for variance testing
removedRatios <- ratios[ratios$Strain!="S606-2" & ratios$Strain!="S607-2",]
pander(var.test(log(at15dgr)~I(Treatment=="Random"),data=removedRatios))

