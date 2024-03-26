# DESCRIPTION:  This Script takes the finalized mutation effect data and creates all main text and supplementary figures
# Other than S02 B -> this is created in the script that processes raw barcode counts to fitness effects


# PACKAGES USE 
library(ggplot2)
library(ggpubr)
library('MASS')
library(viridis)
library(dplyr)
library( tidyr )
library(stats4) 
library(tidyverse)
library(lmtest)
library(moments)
library(purrr)
library(broom)
library(ggridges)
library(ggeffects)
library(microseq) # use to get sequences out from fasta files on computer
library(stringdist) # calc string distances 
library(R.utils) # unzip fastq.gz files
library(bioseq)
library(readxl)


## ALL FILES NEEDED in same folder 
# df_s_ests_wGR_myGRest.csv
# Johnson_DFEmeanDat.csv
# KRE33_pres_abs.csv
# Mut_full_info_Use.csv
# Confirmation Experiment 
# allBCs_Matches 
# transferODtable 


### UPDATE WORKING DIRECTORY TO FOLDER IN WHICH ALL DATA FILES ARE LOCATED 
setwd("~/Desktop/")
fileSave = '~/Desktop/' # UPDATE TO SUB-FOLDER WANT TO SAVE DATA TO (IF DESIRED)
                        # Must CREATE folder first 



neutralMutIDs = c(102,51,99,91,6)
minNumStrains = 5 # min num dat points to do a regression on 
df_s_ests_wGR = read.csv(  paste0('df_s_ests_wGR_myGRest.csv'))
df_s_ests_wGR = df_s_ests_wGR[,-1]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####



######## Figure 1 : Beneficial and Deleterious Mutations ##### 

# Prep data: only keep s-ests for mutations with sufficient data
df_mutsCanFitLine = df_s_ests_wGR %>%
  group_by(env, Mut_ID) %>%
  mutate(totStrainsIn = length(Strain))
df_mutsCanFitLine = subset(df_mutsCanFitLine, totStrainsIn>=minNumStrains & !(Mut_ID %in% neutralMutIDs) & Mut_ID != 28)
# remove neutral muts because want to focus on non-putative neutral, and remove mut 28 because doesnt have enough data


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

## ~~ Plot: Prop Ben and Delt Muts vs GR ####

prop_ben_delt = df_mutsCanFitLine %>%
  group_by(Strain, env, Mean_GR) %>%
  summarize(PropBen = sum(mutCall == 'Beneficial')/length(mutCall),
            PropDelt = sum(mutCall == 'Deleterious')/length(mutCall))

prop_ben_delt_g = gather(prop_ben_delt, key = 'type', value = 'Proportion',PropBen:PropDelt )


plotPben_Delt = function(envID) 
{
  colorDF = data.frame(env = c("30SC3" ,"30SC5", "30SC7" ,"37SC3", "37SC5" ,"37SC7" ,"YPD"  ), color = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'))
  
  colorChoice = subset(colorDF, env == envID)$color
  
  data = subset(prop_ben_delt_g, env == envID)
  data$type = factor(data$type, levels = c('PropBen', 'PropDelt'))
  return(  ggplot(data,aes(x = Mean_GR,y = Proportion, color = type, fill = type, linetype = type))+
             geom_point(aes(shape = type), alpha = 0.6, size  = 1.2)+
             geom_smooth(method = 'lm', formula = 'y~x', lineend = 'round', se = F)+
             scale_color_manual(values = c(colorChoice,colorChoice))+
             scale_fill_manual(values = c(colorChoice, 'white'))+
             scale_shape_manual(values = c(24,25))+
             scale_linetype_manual(values = c('solid', 'dashed'))+
             theme_classic() +
             xlab('Growth Rate (1/hr)')+
             ylab('Proportion in DFE')+
             #ggtitle(paste0(envID))+
             scale_y_continuous(expand = c(0,0), limits = c(-0.2,0.68), breaks = c(0, 0.3, 0.6))+
             scale_x_continuous(  breaks = c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4))+
             coord_cartesian( ylim=c(0,0.68)) +
             theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
                   axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'),
                   #axis.ticks.length=unit(.04, "cm"),
                   axis.title = element_blank(),
                   axis.text.x = element_text(size=10,color = 'black', family = 'Helvetica'),
                   axis.text.y = element_blank(),
                   plot.title = element_text(size=10,color = 'black', family = 'Helvetica'),
                   legend.position = 'none' ,
                   plot.margin = unit(c(0.5,0.5,0.1,0), 'lines') ) )
  
}



allPben_delt_plots = map(unique(prop_ben_delt_g$env), ~plotPben_Delt( .x))

allPben_delt_plots_a = ggarrange(plotlist = allPben_delt_plots, nrow = 2, ncol = 3)

#ggsave(paste0(fileSave,'allPben_delt_plots_a.pdf' ),allPben_delt_plots_a, width = 11.2, height = 7.55, unit = 'cm')






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure 2 : Stats of Slopes and Intercepts  ####


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
## ~~ Analysis: Fit variable slopes model ####

# fit eq (1) for each mutation (e.g, get 1 slope and 1 int per env, save Rsq/pvalues)
perMut_6Slope6int_Rsq_p = df_mutsCanFitLine %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(fit = map(data, ~ lm(avgS ~ Mean_GR * env, .)), glance = map(fit, ~ glance(.x))) %>%
  unnest(glance)
## this gives rsq for full model (eq (1))
perMut_6Slope6int_fits = perMut_6Slope6int_Rsq_p[, c('Mut_ID', 'fit' )]
names(perMut_6Slope6int_fits)[2] = paste0(names(perMut_6Slope6int_fits)[2], '_6s6i')
perMut_6Slope6int_Rsq_p = perMut_6Slope6int_Rsq_p[, c('Mut_ID', 'r.squared' ,'adj.r.squared', 'p.value', 'logLik')]
names(perMut_6Slope6int_Rsq_p)[2:5] = paste0(names(perMut_6Slope6int_Rsq_p)[2:5], '_6s6i')


# also fits eq (1) for each mutation, but saves Rsq/p for the regression in each environment (or the slopes/ints)
perMut_Env_6Slope6int_Rsq_p = df_mutsCanFitLine %>%
  nest(data = -c(Mut_ID ,env)) %>% # nest into Mut_ID, env
  mutate(fit = map(data, ~ lm(avgS ~ Mean_GR , .)), glance = map(fit, ~ glance(.x))) %>%
  unnest(glance) 

perMut_Env_6Slope6int_slopes_ints = df_mutsCanFitLine %>%
  nest(data = -c(Mut_ID ,env)) %>% # nest into Mut_ID, env
  mutate(fit = map(data, ~ lm(avgS ~ Mean_GR , .)), tidy = map(fit, ~ tidy(.x))) %>%
  unnest(tidy) 
# this gives all slopes and ints 
sub = perMut_Env_6Slope6int_slopes_ints[, c('env', 'Mut_ID', 'term', 'estimate')]
perMut_Env_6Slope6int_slopes_ints_spread = spread(sub, term, estimate)
names(perMut_Env_6Slope6int_slopes_ints_spread) = c('env', 'Mut_ID', 'intercept','slope')
remove(sub)



## ~~ Plot: Histogram of Slopes and Intercepts ####

# TO make y-axis fraction, just specify bindiwth and multiple y value by binwidth 
hist_mutSlopes = ggplot(data = perMut_Env_6Slope6int_slopes_ints_spread, aes(x = slope,color = env)) + 
  stat_binline(binwidth = 0.12, aes(color = env,y = after_stat(density)*0.12),geom="step", position = position_dodge(width = 0.02),  linewidth = 0.7, alpha = 1 )+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  geom_vline(xintercept = 0, linetype = 'solid', color = 'grey', linewidth = 0.5, lineend = 'round')+
  theme_classic()+
  xlab('Slope')+
  ylab('Fraction')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 12, family = 'Helvetica'),
        axis.text = element_text(color = 'black', size = 10, family = 'Helvetica'),
        legend.position = 'none')
#ggsave(paste0(fileSave,'hist_mutSlopes.pdf' ),hist_mutSlopes, width = 7, height = 5.4, unit = 'cm')




hist_mutInts = ggplot(data = perMut_Env_6Slope6int_slopes_ints_spread, aes(x = intercept,color = env)) + 
  stat_binline(binwidth = 0.03, aes(color = env,y = after_stat(density)*0.03),geom="step", position = position_dodge(width = 0.02),  linewidth = 0.7, alpha = 1 )+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  geom_vline(xintercept = 0, linetype = 'solid', color = 'grey', linewidth = 0.55, lineend = 'round')+
  theme_classic()+
  xlab('Intercept')+
  ylab('Fraction')+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 12, family = 'Helvetica'),
        axis.text = element_text(color = 'black', size = 10, family = 'Helvetica'),
        legend.position = 'none')

#ggsave(paste0(fileSave,'hist_mutInts.pdf' ),hist_mutInts, width = 7, height = 5.4, unit = 'cm')



# save parameters of slope and int distributions in each env
# allSlopeIntDistParams = perMut_Env_6Slope6int_slopes_ints_spread %>% 
#   group_by(env) %>%
#   summarize(SlopeDistMean = mean(slope, na.rm = T), SlopeDistSD = sd(slope, na.rm = T), , SlopeDistSkew = skewness(slope),
#             intDistMean = mean(intercept, na.rm = T), intDistSD = sd(intercept, na.rm = T), , intDistSkew = skewness(intercept))

#write.csv(allSlopeIntDistParams, 'allSlopeIntDistParams.csv')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

## ~~ Plot:  Slope Int Corr ####
slopeIntCorr = ggplot(perMut_Env_6Slope6int_slopes_ints_spread, aes(x = slope, y = intercept, color = env, group = env))+
  geom_vline(xintercept = 0, color = 'grey', linewidth = 0.3, linetype = 'solid')+
  geom_hline(yintercept = 0, color = 'grey', linewidth = 0.3, linetype = 'solid')+
  geom_point( alpha = 0.5, size = 2, shape = 16)+
  geom_smooth(method = 'lm', se = F, lineend = 'round')+
  xlab('Slope') +  ylab('Intercept')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title = element_text(color = 'black', size = 12, family = 'Helvetica'),
        axis.text =  element_text(color = 'black', size = 10, family = 'Helvetica'),
        legend.text = element_text(color = 'black', size = 10, family = 'Helvetica'),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)

#ggsave(paste0(fileSave,'slopeIntCorr.pdf' ),slopeIntCorr, width = 5.3, height = 5.4, unit = 'cm')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 






## ~~ Numbers: for Fig 1:2 ####

# Number of significant regressions (per Mut)
adjP = p.adjust(perMut_6Slope6int_Rsq_p$p.value_6s6i, method = 'BH')
sum(adjP < 0.05)/length(adjP) # 88/94 (94%)

# IQR for Rsq when fit  eq 1
quantile(perMut_6Slope6int_Rsq_p$r.squared_6s6i, prob = c(0.25,0.5,0.75))  # gives lower and upper quants with median in middle

# Number of significant regressions (per mut, environment)
perMut_Env_6Slope6int_Rsq_p$adjP = p.adjust(perMut_Env_6Slope6int_Rsq_p$p.value ,  method = 'BH')
sum(adjP < 0.05)/length(adjP) # 205/545 (38%)
sub_sigPerEnvMut = subset(perMut_Env_6Slope6int_Rsq_p, adjP<0.05)

# Percent of negative slopes (of sig ones)
sub = inner_join(perMut_Env_6Slope6int_slopes_ints_spread, sub_sigPerEnvMut, by = c("Mut_ID", 'env')) # get just the signficant regression info
sum(sub$slope < 0)/length(sub$slope) # 98% (200/205)

# Percent of positive intercepts (of sig ones)
sum(sub$intercept > 0)/length(sub$intercept) # 96% (196/205)


# Slope Intercept correlation stats  
perMut_Env_6Slope6int_slopes_ints_spread %>%
  nest(data = -env) %>%
  mutate(fit = map(data, ~ lm(intercept ~ slope , .)), glance = map(fit, ~ glance(.x))) %>%
  unnest(glance) 


# Prop Ben and Delt vs GR stats
reg_pBen = prop_ben_delt %>%
  nest(data = -env) %>%
  mutate(fit = map(data, ~ lm(PropBen ~ Mean_GR, .)), glance = map(fit, ~ glance(.x))) %>%
  unnest(glance)

reg_pDelt = prop_ben_delt %>%
  nest(data = -env) %>%
  mutate(fit = map(data, ~ lm(PropDelt ~ Mean_GR, .)), glance = map(fit, ~ glance(.x))) %>%
  unnest(glance)

reg_pBen[,c('env', 'adj.r.squared', 'p.value')]
reg_pDelt[,c('env', 'adj.r.squared', 'p.value')]



# For IQR around variable slope model 
quantile(perMut_6Slope6int_Rsq_p$r.squared_6s6i, prob = c(0.25,0.5,0.75))  # gives lower and upper quants with median in middle





# Times each mutation is beneficial and deleterious 
sub = subset(df_s_ests_wGR, !(Mut_ID %in% c(neutralMutIDs, 28)))
t = sub %>% group_by(Mut_ID) %>% summarize(timesNeutral = sum(mutCall == 'Neutral'), timesBen = sum(mutCall == 'Beneficial'), timesDelt = sum(mutCall == 'Deleterious'))

sum(t$timesBen>0)/length(t$timesBen) # 88/94 (94%)
sum(t$timesBen>0 & t$timesDelt>0)/length(t$timesBen)  # 86/94  (91%)

sum(t$timesBen == 0 & t$timesDelt == 0) # no mut is neutral in all cases 

# Range of proportion beneficial and delterious mutations
t2 = sub %>% 
  group_by(Strain, env) %>%
  summarize(pBen = sum(mutCall == 'Beneficial')/length(mutCall), 
            pDelt = sum(mutCall == 'Deleterious')/length(mutCall))
range(t2$pBen)
range(t2$pDelt)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####



# Figure 3 : Microscopic Epistasis Regressions - all Mutations ####

## first get mutation name info 
Mut_full_info_Use = read.csv( 'Mut_full_info_Use.csv')
Mut_full_info_Use = Mut_full_info_Use[, !names(Mut_full_info_Use) %in% 'X']
df_mutsCanFitLine$Mut_ID = as.numeric(df_mutsCanFitLine$Mut_ID)
df_mutsCanFitLine = inner_join(df_mutsCanFitLine,Mut_full_info_Use, by = 'Mut_ID' )

pSigDiffBins = data.frame(min = c(-1,0,0.2,0.4), max = c(0,0.2,0.4,0.6), bin = c(1,2,3,4))
yaxisScales = data.frame(mMin = c(1,57,79),mMax =  c(56,78,94), scaleMin = c(-0.09,-0.09,-0.09), scaleMax = c(0.09,0.09,0.09))
  # here, jsut using same scales for y-axis for all plots 
plotExample_includepairwiseSlopeDiff = function(rowID) 
{
  m = as.numeric(joinedRanks[rowID,'Mut_ID']) # get the mutation ID assocaited with the row of interest
  
  ## custom y axis scale for plot
  yaxisScales$isIn = rowID >= yaxisScales$mMin & rowID <= yaxisScales$mMax
  ymin = yaxisScales$scaleMin[yaxisScales$isIn]
  ymax = yaxisScales$scaleMax[yaxisScales$isIn]
  
  fitted = subset(perMut_6Slope6int_fits, Mut_ID == m)$fit_6s6i[[1]] # get the variable slope model regression
  data = subset(df_mutsCanFitLine, Mut_ID == m)
  data$predS = fitted$fitted.values # note this only works because data here is subsetted from eaxacr same DF in exact same way as when do the LM
  data$env = factor(data$env, levels = c('30SC3', '30SC5', '30SC7', '37SC3', '37SC5', '37SC7'))
  mutName = subset(Mut_full_info_Use, Mut_ID ==m)$Gene.Use
  
  colorScale_pSigDiff = rev(c('black', 'black', 'black', 'black'))
  
  propSigDiff = subset(allPairwiseSlopeDiffs, Mut_ID == m)$PropSigDiff
  pSigDiffBins$isIn = propSigDiff> pSigDiffBins$min & propSigDiff<=pSigDiffBins$max
  
  colorChoice = colorScale_pSigDiff[pSigDiffBins$isIn]
  
  
  return(  ggplot(data,aes(x = Mean_GR))+
             geom_hline(yintercept = 0, color = 'grey', linetype = 'solid', linewidth = 0.5/2)+
             # geom_point(aes(y = avgS, color = env), alpha = 0.2, size  = 0.4)+
             geom_line(aes(y = predS, color = env),   linewidth = 0.5)+
             scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = F)+
             theme_void() +
             xlab('Growth Rate (1/hr)')+
             ylab('Selection Coefficient (1/hr)')+
             # ggtitle(paste0(mutName))+
             annotate(geom="text", x=-Inf, y=-Inf, label=paste0(mutName),
                      color=colorChoice, size = 1.945, vjust = -0.6, hjust = -0.1)+
             #scale_y_continuous( breaks = c(-0.12,-0.09, -0.06,-0.03, -0.01, 0.00, 0.01,0.03, 0.06, 0.09, 0.12))+
             scale_y_continuous(limits = c(ymin,ymax), breaks = c( ymin, 0.00, ymax))+
             scale_x_continuous(limits = c(0,0.4), breaks = c(0.1,0.3))+
             theme(axis.line =element_line(color = 'black',linewidth=0.5) ,
                   axis.ticks  =element_line(color = 'black',linewidth=0.5),
                   axis.ticks.length=unit(.04, "cm"),
                   axis.title =element_blank(),
                   axis.text = element_blank(),
                   plot.title = element_text(size=5 ,vjust = -7, hjust = 0.15),
                   legend.position = 'none' ,
                   plot.margin = unit(c(0,0.1,0.25,0), 'lines') ) )
  
}

## Calculate the number of pairwise differences in slope per mutation 
# https://strengejacke.github.io/ggeffects/articles/introduction_comparisons.html#hypothesis-testing-for-slopes-of-numeric-predictors

getPropSigDifSlopes = function(data, pCutoff = 0.05)
{
  sub = data #subset(df_mutsCanFitLine, Mut_ID == m)
  sub$Mean_GR = as.numeric(sub$Mean_GR)
  m.interaction <- lm(avgS ~ Mean_GR*env, data =sub)
  
  
  # slopes_withCI = hypothesis_test(m.interaction, c("Mean_GR", "env"), test = NULL) # FOR 71, all but 37sc5 are sig
  slopedifftest = hypothesis_test(m.interaction, c("Mean_GR", "env"))
  slopedifftest$p.value_adj =   p.adjust(slopedifftest$p.value, method = 'BH') # slopedifftest$p.value
  numSigDiffs = sum(slopedifftest$p.value_adj <= pCutoff)
  numTests = length(slopedifftest$p.value_adj <= pCutoff)
  
  out = as.data.frame(slopedifftest)
  sigEnvpairs = paste(subset(out, p.value_adj<=pCutoff)$env, collapse=':')
  alltestedEnvPairs = paste(out$env, collapse=':')
  
  return(data.frame(PropSigDiff = numSigDiffs/numTests, sigEnvpairs = sigEnvpairs, alltestedEnvPairs =alltestedEnvPairs ))
  
}

allPairwiseSlopeDiffs = df_mutsCanFitLine %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(PropSigDiffs_out = map(data, ~getPropSigDifSlopes(.)))%>%
  unnest(PropSigDiffs_out)

allPairwiseSlopeDiffs = allPairwiseSlopeDiffs[,c('Mut_ID', 'PropSigDiff', 'sigEnvpairs','alltestedEnvPairs' )]


#


## Get single slope fit per mutation ##
perMut_1Slope6int_slope = df_mutsCanFitLine %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(fit = map(data, ~ lm(avgS ~ Mean_GR + env, .)), tidy = map(fit, ~ tidy(.x))) %>%
  unnest(tidy)

perMut_1Slope6int_slope = subset(perMut_1Slope6int_slope, term == 'Mean_GR')# only need slopes, not ints
perMut_1Slope6int_slope = perMut_1Slope6int_slope[, c('Mut_ID', 'estimate')]
perMut_1Slope6int_slope$rank_bySlope = rank(perMut_1Slope6int_slope$estimate, ties.method = 'first')# just make the early mutID be 1st rnak when are same, since dont want multiple at same place



# Make all main plots
joinedRanks = inner_join(allPairwiseSlopeDiffs,perMut_1Slope6int_slope, by = 'Mut_ID')
joinedRanks$rank_bySlope = factor(joinedRanks$rank_bySlope )
joinedRanks = joinedRanks[(order(joinedRanks$rank_bySlope )),] # put in slope order, so when bin by propSigdiff, will be in slope order
joinedRanks$rowID = 1:nrow(joinedRanks)
allPlots = map(joinedRanks$rowID, ~plotExample_includepairwiseSlopeDiff( .x))
#allPlots_slopeOrder = ggarrange(plotlist = allPlots, nrow = 12, ncol = 8)


## START : Make the insets - half tile plot ###
envs = unique(df_mutsCanFitLine$env)
env1 = numeric(15)
env2 = numeric(15)
counter = 0
for(e1i in 1:(length(envs)-1))
{
  e1 = envs[e1i]
  for(e2i in (e1i+1):length(envs))
  {
    counter = counter + 1
    e2 = envs[e2i]
    env1[counter] = e1
    env2[counter] = e2
  }
}

allEnvPairs_df =  data.frame(env1,env2)
allEnvPairs_df$envPair =  paste0(allEnvPairs_df$env1,'-', allEnvPairs_df$env2)
allEnvPairs_df$altenvPair =  paste0(allEnvPairs_df$env2,'-', allEnvPairs_df$env1)

# flip right 
allEnvPairs_df = allEnvPairs_df[rev(order((allEnvPairs_df$env1))), ]
allEnvPairs_df$env1 = factor(allEnvPairs_df$env1, levels = unique(allEnvPairs_df$env1) )


makeInsetTile = function(rowID)
{
  pairs = joinedRanks$sigEnvpairs[rowID]
  indivSig = str_split(pairs, ":")[[1]]
  allTested = str_split( joinedRanks$alltestedEnvPairs[rowID], ":")[[1]]
  allEnvPairs_df$sigThisMut = allEnvPairs_df$envPair %in% indivSig | allEnvPairs_df$altenvPair %in% indivSig
  allEnvPairs_df$testedThisMut = allEnvPairs_df$envPair %in% allTested | allEnvPairs_df$altenvPair %in% allTested
  allEnvPairs_df$testedThisMut = factor(allEnvPairs_df$testedThisMut, levels = c(T, F))
  
  # try filling ones that are not sig but were tested with grey
  allEnvPairs_df$sigThisMut[allEnvPairs_df$sigThisMut == F & allEnvPairs_df$testedThisMut == T]  = 'realFalse'
  allEnvPairs_df$sigThisMut = factor(allEnvPairs_df$sigThisMut, levels = c(F,T, 'realFalse'))
  allEnvPairs_df$annotate = ''
  allEnvPairs_df$annotate[allEnvPairs_df$sigThisMut == F] = 'X'
  
  plot_tile = ggplot(allEnvPairs_df, aes(x = env1, y = env2, fill = sigThisMut))+
    geom_tile( color = '#948982', linewidth = 0.2)+ # 
    scale_fill_manual(values = c('grey', 'black','white'), drop = F)+
   #scale_color_manual(values = c('black', 'transparent'), drop = F)+
    #geom_text(aes(label=annotate), color = 'black', size = 0.6)+
    scale_x_discrete(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0))+
    theme_classic()+
    theme_void()+
    theme(axis.text = element_blank(), axis.title = element_blank(), 
          legend.position = 'none' , axis.line = element_blank(), axis.ticks = element_blank(),
          plot.margin = unit(c(0.02,0.02,0,0), 'lines') )
  
  return(plot_tile)
}


# use joined ranks so in same order as all Main plots
allInsets = map(1:nrow(joinedRanks), ~makeInsetTile(.x))

## END : Make the insets - half tile plot ###


## put main plots and insets together 
mainPLot_plusInset = function(i)
{
  p1 = allPlots[[i]]
  p2 = allInsets[[i]]
  
  ## custom y axis scale for plot
  yaxisScales$isIn = i >= yaxisScales$mMin & i <= yaxisScales$mMax
  ymax = yaxisScales$scaleMax[yaxisScales$isIn]
  
  return(p1 + annotation_custom(
    ggplotGrob(p2), 
    xmin = 0.28, xmax = 0.4, ymin = 0.02, ymax = ymax+0.01 # 0.1
  ))
}

allMutPLots_w_Inset = map(1:94, ~mainPLot_plusInset(.x))
allPlots_slopeOrder_wInset = ggarrange(plotlist = allMutPLots_w_Inset, nrow = 12, ncol = 8)

## to make it so can go around legend 
sub1 = allMutPLots_w_Inset[1:10]
sub2 = allMutPLots_w_Inset[11:82]
sub3 = allMutPLots_w_Inset[83:94]

allPlots_slopeOrder_wInset_1 = ggarrange(plotlist = sub1, nrow = 2, ncol = 5)
allPlots_slopeOrder_wInset_2 = ggarrange(plotlist = sub2, nrow = 9, ncol = 8)
allPlots_slopeOrder_wInset_3 = ggarrange(plotlist = sub3, nrow = 2, ncol = 6)

# ggsave(paste0(fileSave,'allPlots_slopeOrder_wInset_1a.pdf' ),allPlots_slopeOrder_wInset_1, width = 16.5/8*5, height = 16.5/12*2, unit = 'cm')
# ggsave(paste0(fileSave,'allPlots_slopeOrder_wInset_2a.pdf' ),allPlots_slopeOrder_wInset_2, width = 16.5/8*8, height = 16.5/12*9, unit = 'cm')
# ggsave(paste0(fileSave,'allPlots_slopeOrder_wInset_3a.pdf' ),allPlots_slopeOrder_wInset_3, width = 16.5/8*6, height = 16.5/12*2, unit = 'cm')



## ~~ Plot : Barplot % Sig Diff Slopes ####
cBin1 = sum(allPairwiseSlopeDiffs$PropSigDiff == 0)
cBin2 = sum(allPairwiseSlopeDiffs$PropSigDiff > 0 & allPairwiseSlopeDiffs$PropSigDiff <= 0.2)
cBin3 = sum(allPairwiseSlopeDiffs$PropSigDiff > 0.2 & allPairwiseSlopeDiffs$PropSigDiff <= 0.4)
cBin4 = sum(allPairwiseSlopeDiffs$PropSigDiff > 0.4 & allPairwiseSlopeDiffs$PropSigDiff <= 0.6)
dfPsigDifCol = data.frame(bin = c(1,2,3,4), count = c(cBin1, cBin2,cBin3,cBin4))

histPsigDiff = ggplot(dfPsigDifCol, aes(x = as.factor(bin), y = count))+
  geom_col(fill = 'black', color = '#948982', linewidth = 0.5)+
  #scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0), limits = c(0,60), breaks = c(0,30,60))+
  ylab('Count')+
  theme_classic()+
  theme(axis.title.y = element_text(color = 'black', size = 10),axis.text = element_text(color = 'black',size = 8), plot.margin = unit(c(0,0,0,0), 'lines'),
        axis.title.x = element_blank())
#ggsave(paste0(fileSave,'histPsigDiff.pdf' ),histPsigDiff, width = 5.8, height = 2.4, unit = 'cm')





## ~~ Numbers - Figure 3 ####

#  Proportion of signficiantly different slopes 
sum(allPairwiseSlopeDiffs$PropSigDiff == 0 )/length(allPairwiseSlopeDiffs$PropSigDiff) # with BH correction, 56% of mutations show no sig differences in slopes 
mean(allPairwiseSlopeDiffs$PropSigDiff  )
max(allPairwiseSlopeDiffs$PropSigDiff )


getNumEnvs = function(input)
{
  return(as.numeric(length(strsplit(input, ':')[[1]])))
}

allPairwiseSlopeDiffs$numtestedEnvs = as.numeric(map(allPairwiseSlopeDiffs$alltestedEnvPairs, ~getNumEnvs(.x)[[1]]))
allPairwiseSlopeDiffs$numSigEnvs = as.numeric(map(allPairwiseSlopeDiffs$sigEnvpairs, ~getNumEnvs(.x)[[1]]))
allPairwiseSlopeDiffs$numNonSigEnvs = allPairwiseSlopeDiffs$numtestedEnvs - allPairwiseSlopeDiffs$numSigEnvs

totComps = sum(allPairwiseSlopeDiffs$numtestedEnvs)
totNS = sum(allPairwiseSlopeDiffs$numNonSigEnvs)

totNS/totComps # 86% (1153/1333)




# ~~~ Looking at only signficant slope muts (those that are sig in variable slopes model)
getTestCases =  perMut_Env_6Slope6int_Rsq_p[, c('env', 'Mut_ID','p.value' )]

df_muts_wSigSlopes = inner_join(df_mutsCanFitLine, getTestCases, by = c('Mut_ID', 'env'))
df_muts_wSigSlopes = subset(df_muts_wSigSlopes, p.value < 0.05)
  
df_muts_wSigSlopes_andEnoughEnvs = df_muts_wSigSlopes %>%
  group_by(Mut_ID) %>%
  mutate(numEnv = length(unique(env)))

df_muts_wSigSlopes_andEnoughEnvs = subset(df_muts_wSigSlopes_andEnoughEnvs, numEnv>=2)

allPairwiseSlopeDiffs_sigSlopes = df_muts_wSigSlopes_andEnoughEnvs %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(PropSigDiffs_out = map(data, ~getPropSigDifSlopes(.)))%>%
  unnest(PropSigDiffs_out)

allPairwiseSlopeDiffs_sigSlopes = allPairwiseSlopeDiffs_sigSlopes[,c('Mut_ID', 'PropSigDiff', 'sigEnvpairs','alltestedEnvPairs' )]

allPairwiseSlopeDiffs_sigSlopes$numtestedEnvs = as.numeric(map(allPairwiseSlopeDiffs_sigSlopes$alltestedEnvPairs, ~getNumEnvs(.x)[[1]]))
allPairwiseSlopeDiffs_sigSlopes$numSigEnvs = as.numeric(map(allPairwiseSlopeDiffs_sigSlopes$sigEnvpairs, ~getNumEnvs(.x)[[1]]))
allPairwiseSlopeDiffs_sigSlopes$numNonSigEnvs = allPairwiseSlopeDiffs_sigSlopes$numtestedEnvs - allPairwiseSlopeDiffs_sigSlopes$numSigEnvs

totComps_sub = sum(allPairwiseSlopeDiffs_sigSlopes$numtestedEnvs)
totNS_sub = sum(allPairwiseSlopeDiffs_sigSlopes$numNonSigEnvs)

totNS_sub/totComps_sub 
# 82.4% (263/319) , when do 3+ envs (>2 envs)
# 82.6% (285/345), when do 2+ envs




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure 4 : Macroscopic Epistasis (DFEs) ####

# ~~ Analysis : Get F0 values ####

# ~~ first, need to get slope and intercepts for the model out, isn't super straightforward from the dyplr tidy output because 
# ~~ each mut in diff numbers of envs, so write for self

Mut_ID = c()
env = c()
slope =c()
intercept =c()
GR = c()
Rsq  = c()
pval = c()
mean_int = c()
resVar = c()
allResiduals = c()
mcounter = 0
for(M in unique(df_mutsCanFitLine$Mut_ID))
{
  mcounter = mcounter + 1
  print(mcounter)
  subMut = subset(df_mutsCanFitLine, Mut_ID == M)
  
  if(dim(subMut)[1] > minNumStrains & length(unique(subMut$env))>1)
  {
    
    
    # ordinary least squares 
    lmout = (lm(avgS ~ Mean_GR + env, data = subMut, na.action = na.exclude)) 
    allCoeff  = lmout$coefficients
    
    
    # fit1 = lmout
    # test = subMut
    # test$predS = fit1$fitted.values
    # # plot with equal slopes, this only really works for mut with all envs measured
    # plot =  ggplot(test,aes(x = Mean_GR))+
    #   geom_hline(yintercept = 0, color = 'grey', linetype = 'dashed', linewidth = 3)+
    #   geom_point(aes(y = avgS, color = env), alpha = 0.5, size  = 3)+
    #   geom_line(aes(y = predS, color = env), linewidth = 3)+
    #   scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'))+
    #   geom_point(x = 0, y = allCoeff[['(Intercept)']], color ='#00026E', size = 7 )+
    #   geom_point(x = 0, y = allCoeff[['(Intercept)']]+allCoeff[['env30SC5']], color ='#295DCC', size = 7 )+
    #   geom_point(x = 0, y = allCoeff[['(Intercept)']]+allCoeff[['env30SC7']], color ='#90B6DD', size = 7 )+
    #   geom_point(x = 0, y = allCoeff[['(Intercept)']]+allCoeff[['env37SC3']], color ='#C21700', size = 7 )+
    #   geom_point(x = 0, y = allCoeff[['(Intercept)']]+allCoeff[['env37SC5']], color ='#FE8D26', size = 7 )+
    #   geom_point(x = 0, y = allCoeff[['(Intercept)']]+allCoeff[['env37SC7']], color ='#FBCF96', size = 7 )+
    #   xlim(0,0.4)+
    #   theme_classic() +
    #   xlab('Growth Rate')+
    #   ylab('Selection Coefficient')+
    #   ggtitle(paste0('Mutation ID: ',M))+
    #   theme(axis.line =element_line(color = 'black',linewidth=1.3) ,axis.ticks  =element_line(color = 'black',linewidth=1.3),
    #         axis.title =element_blank(),
    #         axis.text = element_text(color = 'black', size = 15),
    #         plot.title = element_text(size=20),
    #         legend.position = 'none' )
    
    
    
    
    
    for(envChoose in unique(subMut$env))
    {
      
      
      idName = paste0('env', envChoose)
      
      # calc all real ints
      test = allCoeff[-2] # remove GR term
      test[2:length(test)] = test[2:length(test)]+test[1]
      
      if(idName %in% names(allCoeff) ) # if is one of the actual ones
      {
        Mut_ID = c(Mut_ID, M)
        env = c(env, envChoose)
        slope = c(slope, allCoeff[['Mean_GR']] )
        intercept = c(intercept, allCoeff[['(Intercept)']]+allCoeff[[idName]] )
        Rsq  = c(Rsq, summary(lmout)$r.squared)
        pval = c(pval, summary(lmout)$coefficients[1,4])
        mean_int = c(mean_int, mean(test))
        #resVar = c(resVar, var(summary(lmout)$residuals[subMut$env == envChoose], na.rm = T)) # juet get resisduals for this ENV 
        resVar = c(resVar, var(summary(lmout)$residuals, na.rm = T)) #  get resisduals for all ENVs
        
        allResiduals = c(allResiduals, as.numeric(summary(lmout)$residuals))
      }else{ # is the first one, 
        Mut_ID = c(Mut_ID, M)
        env = c(env, envChoose)
        slope = c(slope, allCoeff[['Mean_GR']] )
        intercept = c(intercept, allCoeff[['(Intercept)']])
        Rsq  = c(Rsq, summary(lmout)$r.squared)
        pval = c(pval, summary(lmout)$coefficients[1,4])
        mean_int = c(mean_int, mean(test))
       # resVar = c(resVar,  var(summary(lmout)$residuals[subMut$env == envChoose], na.rm = T)) # juet get resisduals for this ENV 
         resVar = c(resVar, var(summary(lmout)$residuals, na.rm = T)) # juet get resisduals for this ENV 
        
         allResiduals = c(allResiduals, as.numeric(summary(lmout)$residuals))
      }
      
      
      
      
    }
    
    
    
  }else{
    Mut_ID = c(Mut_ID, M)
    env = c(env, envChoose)
    slope = c(slope, NA)
    intercept = c(intercept, NA)
    Rsq  = c(Rsq, NA)
    pval = c(pval, NA)
    mean_int = c(mean_int, NA)
    resVar = c(resVar,  NA)
   
  }
  
}


df_Mutslopes_1slope_6int = data.frame(Mut_ID ,
                                      env ,
                                      slope ,
                                      intercept,
                                      Rsq  ,
                                      pval ,
                                      mean_int ,
                                      resVar)
df_Mutslopes_constantSlopes = df_Mutslopes_1slope_6int


df_F0_1 = df_Mutslopes_1slope_6int %>%
  nest(data = -c(env)) %>% # nest into Mut_ID, env
  mutate(fit = map(data, ~ lm(intercept ~ 0+slope , .)), tidy = map(fit, ~ tidy(.x))) %>% # enforced intercept at 0,0
  unnest(tidy) 
df_F0_1 = df_F0_1[, c('env', 'term' ,  'estimate', 'std.error')]
df_F0 = spread(df_F0_1, term, estimate)
names(df_F0) = c('env',  'stdErr','F0')
df_F0$F0 = abs(df_F0$F0)


#write.csv(df_F0, paste0(fileSave, 'df_F0.csv'))



# ~~ Analysis : calculate adjusted GR for all strains in all envs ####

df_mutsCanFitLine_wAdjGR = inner_join(df_mutsCanFitLine, df_F0, by = 'env')
df_mutsCanFitLine_wAdjGR$adjGR = df_mutsCanFitLine_wAdjGR$Mean_GR - df_mutsCanFitLine_wAdjGR$F0




# ~~ Plot : All Regressions each env, eta = 0 and otherwise ####
plotAllRegressions = function(environment, eta_to_0 = F, title = T, useAdjGR  = F)
{
  xlim = c(0,0.8)
  subDF = subset(df_Mutslopes_1slope_6int, env == environment)
  F0 = subset(df_F0, env == environment)$F0
  xIntVal = F0
  if(eta_to_0 == T & useAdjGR == F) # if using eta to 0, change intercept to that pred by -1*slope*F0
  {
    subDF$intercept = -1*subDF$slope * F0
  }
  
  if(useAdjGR == T & eta_to_0 == F) 
  {
    subDF$intercept =  subDF$intercept + 1*subDF$slope * F0 # shift all ints by mean int pred
    xlim = c(-0.3,0.3)
    xIntVal = 0
  }
  if(eta_to_0 == T & useAdjGR == T)
  {
    subDF$intercept = 0 
    xlim = c(-0.3,0.3)
    xIntVal = 0
  }
  
  testerdata = data.frame(x = c(0.037,F0,0.38), y = 0.1)
  titleUse = paste0(environment)[title] # gives title if T, blank title if F
  
  plotout = ggplot()+
    geom_abline(data = subDF, aes(slope = slope, intercept = intercept, group = Mut_ID, color = as.factor(sign(slope))), linewidth = 0.14, alpha = 0.9,lineend='round')+
    geom_hline(yintercept = 0, color = 'grey', linewidth = 1, linetype = 'solid')+
    geom_vline(xintercept = xIntVal, color = 'grey', linewidth = 0.5, linetype = 'solid')+
    scale_color_manual( values = c('black', 'black'))+
    scale_x_continuous(limits = xlim)+
    scale_y_continuous(limits = c(-0.4,0.2))+
    xlab('Growth Rate (1/hr)')+ ylab('Selection Coefficient')+
    #geom_point(data = testerdata, aes(x = x ), y = -Inf,fill = 'grey',color = 'white', size = 8, shape = c(22,23,24))+
    ggtitle(paste0(titleUse))+
    theme_classic()+
    theme(axis.line = element_line(linewidth = 0.5, color = 'black',lineend='round'),
          axis.ticks = element_line(linewidth = 0.5, color = 'black',lineend='round'),
          axis.title = element_blank(),
          axis.text =  element_text(color = 'black', size = 8),
          legend.text = element_text(color = 'black', size = 8),
          legend.key.size = unit(0.2, 'in'),
          legend.position = 'none') #text(color = 'black', size = 10,angle = 90)
  
  
  return(plotout)
}



SingleEnv_ex = plotAllRegressions( '30SC7', useAdjGR = T, eta_to_0 = T,title = F)
#ggsave(paste0(fileSave,'SingleEnv_ex_allReg.pdf' ),SingleEnv_ex, width = 4.9, height = 4.33+0.287+0.358, unit = 'cm')





# ~~ Plot : 3 adj GR bins, plot compiled DFEs each env ####
numBins = 3

## Make ADJ GR bins 
binSeq_adjGR = seq(-0.16, 0.16, length.out = numBins+1 )# slightly over range of actual adj GRs, but allows nice #s
adjGR_binDF = data.frame(minBin = binSeq_adjGR[1:(length(binSeq_adjGR)-1)], maxBin =  binSeq_adjGR[2:(length(binSeq_adjGR))])
adjGR_binDF$binID = 1:nrow(adjGR_binDF)
names(adjGR_binDF ) = paste0(names(adjGR_binDF), '_adjGR')

## Make ABS GR bins
binSeq_absGR = seq(0,0.4, length.out = numBins+1 )
absGR_binDF = data.frame(minBin = binSeq_absGR[1:(length(binSeq_absGR)-1)], maxBin =  binSeq_absGR[2:(length(binSeq_absGR))])
absGR_binDF$binID = 1:nrow(absGR_binDF)
names(absGR_binDF ) = paste0(names(absGR_binDF), '_absGR')

## add to df_mutsCanFitLine_wAdjGR
df_mutsCanFitLine_wAdjGR$bin_absGR = NA
df_mutsCanFitLine_wAdjGR$bin_adjGR = NA
for(i in 1:max(absGR_binDF$binID_absGR))
{
  minBin_abs = absGR_binDF$minBin_absGR[i]
  maxBin_abs = absGR_binDF$maxBin_absGR[i]
  
  minBin_adj = adjGR_binDF$minBin_adjGR[i]
  maxBin_adj = adjGR_binDF$maxBin_adjGR[i]
  
  
  df_mutsCanFitLine_wAdjGR$bin_absGR[df_mutsCanFitLine_wAdjGR$Mean_GR >= minBin_abs & df_mutsCanFitLine_wAdjGR$Mean_GR < maxBin_abs] = i
  df_mutsCanFitLine_wAdjGR$bin_adjGR[df_mutsCanFitLine_wAdjGR$adjGR >= minBin_adj & df_mutsCanFitLine_wAdjGR$adjGR < maxBin_adj] = i
  
}


# if put all envs together
plot_compiledDFE_adjGR_allenv = function( binChoose, binwidthChoose = 0.03, adjGR = T)
{
  if(adjGR) # use adjGR
  {
    data = subset(df_mutsCanFitLine_wAdjGR, bin_adjGR == binChoose)
    numStrainsperenv = data %>% group_by(env) %>% summarize(numStrains = length(unique(Strain)))
    envsChoose = subset(numStrainsperenv, numStrains>=3)$env
    data = subset(data, env%in%envsChoose)
    data$env = factor(data$env, levels = c("30SC3" ,"30SC5" ,"30SC7", "37SC3" ,"37SC5" ,"37SC7"))
    
  }else{
    data = subset(df_mutsCanFitLine_wAdjGR, bin_absGR == binChoose)
    numStrainsperenv = data %>% group_by(env) %>% summarize(numStrains = length(unique(Strain)))
    envsChoose = subset(numStrainsperenv, numStrains>=3)$env
    data = subset(data, env%in%envsChoose)
    data$env = factor(data$env, levels = c("30SC3" ,"30SC5" ,"30SC7", "37SC3" ,"37SC5" ,"37SC7"))
    
  }
  
  
  
  # use stat_binline instead of staat_bin,, this actually just outlines histogram instead of creating shifted bin count
  adjGR_compDFEs = ggplot(data, aes(x = avgS))+
  # geom_histogram(bins = 17, aes(y = after_stat(density)), fill = '#74706e' , color = '#74706e')+
   # stat_binline(bins = 17, aes(color = env,group = env, y = after_stat(density)),geom="step", position = position_dodge(width = 0.002),  linewidth = 0.7, alpha = 1 )+
    geom_histogram(binwidth = binwidthChoose, aes(y = after_stat(density)*binwidthChoose), fill = 'grey' , color = 'grey')+
    stat_binline(binwidth = binwidthChoose, aes(color = env,group = env, y = after_stat(density)*binwidthChoose),geom="step", position = position_dodge(width = 0.002),  linewidth = 0.7, alpha = 1 )+
    geom_vline(xintercept = 0, color = 'grey', linetype = 'solid', linewidth = 0.7,lineend='round')+
    scale_x_continuous(limits = c(-0.2, 0.2), breaks = c(-0.15, 0, 0.15))+
    scale_y_continuous(expand = c(0,0), limits = c(0,0.8), breaks = c(0,0.3,0.6))+
    scale_color_manual(values = c("#00026E", "#295DCC", "#90B6DD", "#C21700" ,"#FE8D26", "#FBCF96"), drop = F)+
    #ggtitle(paste0('Adjusted GR Bin, ', e2))+
    theme_classic()+
    theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
          axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'), 
          axis.title = element_blank(),
          axis.text.y = element_text(color = 'black', family = 'Helvetica', size = 8),
          axis.text.x = element_blank(),
          legend.position = 'none',
          plot.margin = unit(c(0.2,0.1,0,0), 'lines'))
  
  
  
  return(adjGR_compDFEs)
}

allBinDFEs = ggarrange(plotlist = map(unique(adjGR_binDF$binID_adjGR), ~plot_compiledDFE_adjGR_allenv( .x)), nrow = 1, ncol = length(unique(adjGR_binDF$binID_adjGR)))
plotlist = map(unique(adjGR_binDF$binID_adjGR), ~plot_compiledDFE_adjGR_allenv( .x))
 # ggsave(paste0(fileSave,'compiledDFE_bin1.pdf' ),plotlist[[1]], width = 4.8-0.56, height = 2.35, unit = 'cm')
 # ggsave(paste0(fileSave,'compiledDFE_bin2.pdf' ),plotlist[[2]], width = 4.8-0.56, height = 2.35, unit = 'cm')
 # ggsave(paste0(fileSave,'compiledDFE_bin3.pdf' ),plotlist[[3]], width = 4.8-0.56, height = 2.35, unit = 'cm')




# get dist of GRs for each of the bins
plot_distGRs_eachBin = function( binChoose)
{
  
  
  ## to highlight the area of GRs in the full hist
 minBin =  subset(adjGR_binDF, binID_adjGR ==binChoose)$minBin_adjGR
 maxBin =  subset(adjGR_binDF, binID_adjGR ==binChoose)$maxBin_adjGR
 
  
  data = df_mutsCanFitLine_wAdjGR # subset(df_mutsCanFitLine_wAdjGR, bin_adjGR == binChoose)
  numStrainsperenv = data %>% group_by(env) %>% summarize(numStrains = length(unique(Strain)))
  envsChoose = subset(numStrainsperenv, numStrains>=3)$env
  data = subset(data, env%in%envsChoose)
  data$env = factor(data$env, levels = c("30SC3" ,"30SC5" ,"30SC7", "37SC3" ,"37SC5" ,"37SC7"))
  
  data$colorBy = 'outRange'
  data$colorBy[data$bin_adjGR == binChoose] = 'inRange'
  # 
  # adjGR_compDFEs = ggplot(data, aes(x = adjGR))+
  #   #geom_histogram(fill = 'transparent', color = 'black', aes(y = after_stat(density)), bins = 20)+
  # #  geom_rect(xmin=minBin, xmax=maxBin, ymin=0, ymax=6, fill='grey', color="grey", alpha=0.5)+
  #   geom_vline(xintercept = minBin, color = 'grey', linetype = 'solid', linewidth = 0.5)+
  #   geom_vline(xintercept = maxBin, color = 'grey', linetype = 'solid', linewidth = 0.5)+
  #  stat_binline(aes(y = after_stat(density)),geom="step",bins=10, linewidth = 0.2, alpha = 1)+
  #  # geom_histogram(  aes(y = after_stat(density), fill = colorBy), color = 'black', bins = 10)+
  #   geom_vline(xintercept = 0, color = 'grey', linetype = 'solid', linewidth = 0.5)+
  #   scale_x_continuous(limits = c(-0.2,0.2), breaks = c(-0.18, 0, 0.18))+
  #   scale_y_continuous(expand = c(0,0), breaks = c(0,3,6), limits = c(0,6))+
  #   #scale_fill_manual(values = c("white", "darkgrey"), drop = F)+
  #   theme_classic()+
  #   theme(axis.line = element_line(color = 'black',linewidth=0.2, lineend = 'round') ,
  #         axis.ticks  =element_line(color = 'black',linewidth=0.2, lineend = 'round'), 
  #         axis.title = element_blank(),
  #         axis.text = element_blank(),
  #         legend.position = 'none',
  #         plot.margin = unit(c(0.1,0.1,0.1,0), 'lines'),
  #         axis.ticks.length=unit(.05, "cm"))
  # 
  
  numEntries = length(data$Strain)
  binWidthChoose = 0.0355
  adjGR_compDFEs =  ggplot(data, aes(x = adjGR))+
    geom_histogram(binwidth = binWidthChoose, aes(y = after_stat(count)/numEntries, fill = colorBy), color = 'transparent')+
    stat_binline(binwidth = binWidthChoose, aes(y = after_stat(count)/numEntries),geom="step",linewidth = 0.2, alpha = 1, color = 'black')+
    geom_vline(xintercept = 0, color = 'grey', linetype = 'solid', linewidth = 0.5)+
    scale_fill_manual(values = c('#58595B', 'white'))+
    scale_y_continuous(expand = c(0,0), limits = c(0,0.25), breaks = c(0,0.2) )+
    scale_x_continuous( limits = c(-0.2,0.2), breaks = c(-0.2,0,0.2) )+
    ylab('Prop.')+
    xlab('Adj. GR')+
    theme_classic()+
    theme(axis.line = element_line(color = 'black',linewidth=0.2, lineend = 'round') ,
          axis.ticks  =element_line(color = 'black',linewidth=0.2, lineend = 'round'), 
          axis.title.x = element_text(size = 6, color = 'black', family = 'Helvetica'),
          axis.title.y = element_blank(),
          axis.text =element_text(size = 6, color = 'black', family = 'Helvetica'),
          legend.position = 'none',
          plot.margin = unit(c(0.1,0.1,0.1,0), 'lines'),
          axis.ticks.length=unit(.05, "cm"))
  
   
  
  
  return(adjGR_compDFEs)
}

grs_bin1 = plot_distGRs_eachBin(1)
grs_bin2 = plot_distGRs_eachBin(2)
grs_bin3 = plot_distGRs_eachBin(3)
#ggsave(paste0(fileSave,'grs_bin1.pdf' ),grs_bin1, width = 1.65, height = 1.2, unit = 'cm')
#ggsave(paste0(fileSave,'grs_bin2.pdf' ),grs_bin2,  width = 1.65, height = 1.2, unit = 'cm')
#ggsave(paste0(fileSave,'grs_bin3.pdf' ),grs_bin3, width = 1.65, height = 1.2, unit = 'cm')




# ~~ Analysis: Predicted DFEs ####
df_F0_2 = df_Mutslopes_1slope_6int %>%
  nest(data = -c(env)) %>% # nest into Mut_ID, env
  mutate(fit = map(data, ~ lm(intercept ~ 0+slope , .)), augment = map(fit, ~ augment(.x))) %>% # enforced intercept at 0,0
  unnest(augment) 

# Parameters needed
Mean_slopeDist = mean(df_Mutslopes_1slope_6int$slope) # can also take avg of avg in each env, yields -0.1507 (almost same)
Var_slopeDist = var(df_Mutslopes_1slope_6int$slope)
Skew_slopeDist = skewness(df_Mutslopes_1slope_6int$slope)

lm_resVar_vs_b = summary(lm(resVar ~0+slope, data =df_Mutslopes_1slope_6int )) # enforce 0,0
alpha_resVar_vs_b = -1*lm_resVar_vs_b$coefficients[1]


# for diff type- where slope indep res var
avgResVar = 0.0002291042 # mean(df_Mutslopes_1slope_6int$resVar)
# 0.0002291042 is the var of ALlresiduals saved
etaVar = var(df_F0_2$.resid) # e.g., sigma^2 pivot


dfe_parampred = function(GR, environment)
{
  
  F_0_e = df_F0[match(environment,df_F0$env ),]$F0 #subset(df_F0, env == environment)$F0
  # this syntax allows for multiple envs inputted 
  # have to either input 1 or numGR entries though
 
  # for  slope-dependent res var
  DFEmean = (GR - F_0_e) * Mean_slopeDist 
  DFEvar = -alpha_resVar_vs_b*Mean_slopeDist + etaVar + (GR - F_0_e)^2 * Var_slopeDist  
  DFEskew = (Skew_slopeDist *  (GR - F_0_e)^3 + -3*alpha_resVar_vs_b*(GR - F_0_e)/sqrt(Var_slopeDist))/(( (GR - F_0_e)^2 + (-alpha_resVar_vs_b*Mean_slopeDist +etaVar)/Var_slopeDist)^(3/2)) 
  
  
  return(data.frame(DFEmean, DFEvar, DFEskew))
  
}


# ~~ Analysis :  get DFE stats ####
oldDFEdat = df_mutsCanFitLine_wAdjGR %>%
  group_by(Strain, env, Mean_GR, adjGR, Std_err, F0) %>%
  summarize(DFEmean = mean(avgS, na.rm= T), DFEvar = var(avgS, na.rm = T), DFEskew = skewness(avgS, na.rm = T))

oldDFEdat$predDFEmean = dfe_parampred(as.numeric(oldDFEdat$Mean_GR), oldDFEdat$env)$DFEmean
oldDFEdat$predDFEvar = dfe_parampred(as.numeric(oldDFEdat$Mean_GR), oldDFEdat$env)$DFEvar
oldDFEdat$predDFEskew = dfe_parampred(as.numeric(oldDFEdat$Mean_GR), oldDFEdat$env)$DFEskew




# ~~ Analysis : Get Error in all params from Bootstrapping ####
numMutsSamp = 70
numIts = 300
set.seed(NULL)

sampleNmuts_DFEstat = function(data, numMutsSamp)
{
  indecies = 1:nrow(data)
  samps = sample(indecies, numMutsSamp, replace = T)
  samp_s = (data[samps,])$avgS
  DFEmean = mean(samp_s)
  DFEvar = var(samp_s)
  DFEskew = skewness(samp_s)

    return(data.frame(DFEmean,DFEvar,DFEskew) )
  
}

bootError = function(data) # data is subset of df_mutsCanFitLine_wAdjGR for 1 env, 1 strain
{
  # For testing, 
  # data = subset(df_mutsCanFitLine_wAdjGR, Strain == 'LK1-A02'& env == '30SC3')
  if(nrow(data) > 5)
  {
    allDFEmean = numeric(numIts)
    allDFEvar = numeric(numIts)
    allDFEskew = numeric(numIts)
    for(i in 1:numIts)
    {
      out = sampleNmuts_DFEstat(data, numMutsSamp)
      allDFEmean[i] = out$DFEmean
      allDFEvar[i] = out$DFEvar
      allDFEskew[i] = out$DFEskew
    }
    
    return(data.frame(DFEmean_SD  = sd(allDFEmean),DFEvar_SD  = sd(allDFEvar),DFEskew_SD  = sd(allDFEskew)   ))
  }else{
    return(data.frame(DFEmean_SD  =NA,DFEvar_SD  = NA,DFEskew_SD  = NA  ))
  }
  
  
}

allDFEbootSD = df_mutsCanFitLine_wAdjGR %>%
  nest(data = -c(Strain,env)) %>%
  mutate(fit = map(data, ~ bootError(.)) )%>%
  unnest(fit) 

allDFEbootSD = allDFEbootSD[, c('Strain', 'env', 'DFEmean_SD', 'DFEvar_SD', 'DFEskew_SD')]
#write.csv(allDFEbootSD,paste0( fileSave, 'allDFEbootSD_300Its.csv'))
#allDFEbootSD = read.csv(paste0( fileSave, 'allDFEbootSD_300Its.csv'))


oldDFEdat = full_join(oldDFEdat,allDFEbootSD, by = c('Strain', 'env') )
oldDFEdat$DFEmean_min = oldDFEdat$DFEmean - oldDFEdat$DFEmean_SD
oldDFEdat$DFEmean_max = oldDFEdat$DFEmean + oldDFEdat$DFEmean_SD
oldDFEdat$DFEvar_min = oldDFEdat$DFEvar - oldDFEdat$DFEvar_SD
oldDFEdat$DFEvar_max = oldDFEdat$DFEvar + oldDFEdat$DFEvar_SD
oldDFEdat$DFEskew_min = oldDFEdat$DFEskew - oldDFEdat$DFEskew_SD
oldDFEdat$DFEskew_max = oldDFEdat$DFEskew + oldDFEdat$DFEskew_SD

oldDFEdat = oldDFEdat[, !(names(oldDFEdat) %in% c('DFEmean_SD', 'DFEvar_SD', 'DFEskew_SD'))]

## ~~ END: Get Error in all params from Bootstrapping


## ~~~ Add Johnson data 
JohnsonDat = read.csv(paste0( 'JohnsonDFEmeanDat.csv'))
sub_JohnsonDat = subset(JohnsonDat, Strain %in% unique(oldDFEdat$Strain))
sub_JohnsonDat$F0 = 0.52024718 # estimated previously in excel sheet 
df_F0 = rbind(df_F0, data.frame(env = 'YPD', stdErr = NA, F0 = 0.52024718))
sub_JohnsonDat$adjGR = sub_JohnsonDat$Mean_GR - sub_JohnsonDat$F0

sub_JohnsonDat$predDFEmean = dfe_parampred(as.numeric(sub_JohnsonDat$Mean_GR), sub_JohnsonDat$env)$DFEmean
sub_JohnsonDat$predDFEvar = dfe_parampred(as.numeric(sub_JohnsonDat$Mean_GR), sub_JohnsonDat$env)$DFEvar
sub_JohnsonDat$predDFEskew = dfe_parampred(as.numeric(sub_JohnsonDat$Mean_GR), sub_JohnsonDat$env)$DFEskew


oldDFEdat = rbind(oldDFEdat,sub_JohnsonDat )



# ~~ Plot : All DFE moments ####

DFEMean_plot_rawGR = ggplot(subset(oldDFEdat), aes(x = Mean_GR,color = env))+
  geom_linerange(aes(ymin = DFEmean_min, ymax = DFEmean_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = Mean_GR - Std_err, xmax = Mean_GR + Std_err, y = DFEmean), alpha = 0.5,lineend='round')+
  geom_point(aes(y = DFEmean ), color = 'white', shape = 16)+
  geom_point(aes(y = DFEmean ), alpha = 0.5, shape = 16)+
  # geom_smooth(method = 'lm', se = F)+
  geom_line(aes(y = predDFEmean, color = env, group = env), linewidth = 1, lineend='round')+
  xlab('Background Growth Rate (1/h)')+
  ylab('DFE Mean')+
  #ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  #scale_x_continuous(limits = c(0,0.4), breaks = c(0,0.2,0.4))+
  scale_x_continuous( breaks = c(0,0.2,0.4,0.6,0.8))+
  scale_y_continuous( breaks = c(-0.03,0,0.03))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_text(size = 18, color = 'black'),
        axis.text = element_text(size = 12, color = 'black'),
        legend.position = 'none' )
#ggsave(paste0(fileSave, 'DFEMean_plot_rawGR.pdf'), DFEMean_plot_rawGR, width = 4.29, height = 4.35 )


DFEVar_plot_rawGR = ggplot(subset(oldDFEdat), aes(x = Mean_GR ,color = env))+
  geom_linerange(aes(ymin = DFEvar_min, ymax = DFEvar_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = Mean_GR - Std_err, xmax = Mean_GR + Std_err, y = DFEvar), alpha = 0.5,lineend='round')+
  geom_point(aes(y = DFEvar ), color = 'white', shape = 16)+
  geom_point(aes( y = DFEvar), alpha = 0.5, shape = 16)+
  # geom_smooth(method = 'lm', se = F)+
  geom_line(aes(y = predDFEvar, color = env, group = env), linewidth = 1,lineend='round')+
  xlab('Background Growth Rate (1/h)')+
  ylab('DFE Variance')+
  # ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  #scale_x_continuous(limits = c(0,0.4), breaks = c(0,0.2,0.4))+
  scale_x_continuous( breaks = c(0,0.2,0.4,0.6,0.8))+
  scale_y_continuous( breaks = c( 0,0.001,0.002))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_text(size = 18, color = 'black'),
        axis.text = element_text(size = 12, color = 'black'),
        legend.position = 'none' )


#ggsave(paste0(fileSave, 'DFEVar_plot_rawGR.pdf'), DFEVar_plot_rawGR, width = 4.29, height = 4.35 )



DFESkew_plot_rawGR = ggplot(subset(oldDFEdat), aes(x = Mean_GR,color = env))+
  geom_linerange(aes(ymin = DFEskew_min, ymax = DFEskew_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = Mean_GR - Std_err, xmax = Mean_GR + Std_err, y = DFEskew), alpha = 0.5,lineend='round')+
  geom_point(aes(y = DFEskew ), color = 'white', shape = 16)+
  geom_point(aes( y = DFEskew ), alpha = 0.5, shape = 16)+
  # geom_smooth(method = 'lm', se = F)+
  geom_line(aes(y = predDFEskew, color = env, group = env), linewidth = 1,lineend='round')+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Skewness')+
  #ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  #scale_x_continuous(limits = c(0,0.4), breaks = c(0,0.2,0.4))+
  scale_x_continuous( breaks = c(0,0.2,0.4,0.6,0.8))+
  scale_y_continuous( breaks = c(-3, 0,3))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_text(size = 18, color = 'black'),
        axis.text = element_text(size = 12, color = 'black'),
        legend.position = 'none' )

#ggsave(paste0(fileSave, 'DFESkew_plot_rawGR.pdf'), DFESkew_plot_rawGR, width = 4.29, height = 4.35 )


## points to add to adjGr plot
adjGR_binDF$midpoint = (adjGR_binDF$minBin_adjGR + adjGR_binDF$maxBin_adjGR)/2
dfextrapoints = data.frame(x = c(adjGR_binDF$midpoint))
predfestatextrapt = dfe_parampred(dfextrapoints$x + df_F0$F0[df_F0$env == '30SC5'], '30SC5' ) # just turn adj gr back to raw for 1 env (will work for any env)
predfestatextrapt$x = dfextrapoints$x

DFEMean_plot_adjGR = ggplot(oldDFEdat, aes(x = adjGR,color = env))+
  geom_linerange(aes(ymin = DFEmean_min, ymax = DFEmean_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEmean), alpha = 0.5,lineend='round')+
  geom_point(aes(y = DFEmean ), color = 'white', shape = 16)+
  geom_point(aes(y = DFEmean ), alpha = 0.5, shape = 16)+
  # geom_smooth(method = 'lm', se = F)+
  geom_line(aes(y = predDFEmean), color = 'black', linewidth = 1,lineend='round')+
  geom_point(data = predfestatextrapt, aes(x = x, y = DFEmean),  color = 'black', fill = 'white', shape = c(22,23,24), size = 3, stroke = 1.3)+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Mean')+
  #ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(-0.03,0,0.03))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  = element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )




DFEVar_plot_adjGR = ggplot(oldDFEdat, aes(x = adjGR ,color = env))+
  geom_linerange(aes(ymin = DFEvar_min, ymax = DFEvar_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEvar), alpha = 0.5,lineend='round')+
  geom_point(aes( y = DFEvar ), color = 'white', shape = 16)+# underlay so dont see cross of error bar in point
  geom_point(aes( y = DFEvar), alpha = 0.5, shape = 16)+
  # geom_smooth(method = 'lm', se = F)+
  geom_line(aes(y = predDFEvar), color = 'black', linewidth = 1,lineend='round')+
  geom_point(data = predfestatextrapt, aes(x = x, y = DFEvar), color = 'black', fill = 'white', shape = c(22,23,24), size = 3, stroke = 1.3)+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Variance')+
  # ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(0,0.001,0.002))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )




DFESkew_plot_adjGR = ggplot(oldDFEdat, aes(x = adjGR,color = env))+
  geom_linerange(aes(ymin = DFEskew_min, ymax = DFEskew_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEskew), alpha = 0.5,lineend='round')+
  geom_point(aes( y = DFEskew ), color = 'white', shape = 16)+# underlay so dont see cross of error bar in point
  geom_point(aes( y = DFEskew ), alpha = 0.5, shape = 16)+
  # geom_smooth(method = 'lm', se = F)+
  geom_line(aes(y = predDFEskew), color = 'black', linewidth = 1,lineend='round')+
  geom_point(data = predfestatextrapt, aes(x = x, y = DFEskew), color = 'black', fill = 'white', shape = c(22,23,24), size = 3, stroke = 1.3)+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Skewness')+
  #ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(-3,0,3))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )



## local rep
x_adj = 0.25--0.2
x_raw = 0.8

width_Mean_GR = 7.1
adjGR_width = width_Mean_GR*(x_adj/x_raw) + 0.2

# 
# ggsave(paste0(fileSave,'DFEMean_plot_adjGR.pdf' ),DFEMean_plot_adjGR, width = adjGR_width , height =4.37-0.4, unit = 'cm')
# ggsave(paste0(fileSave,'DFEvar_plot_adjGR.pdf' ),DFEVar_plot_adjGR, width = adjGR_width, height = 4.37-0.4, unit = 'cm')
# ggsave(paste0(fileSave,'DFEskew_plot_adjGR.pdf' ),DFESkew_plot_adjGR, width = adjGR_width, height = 4.37-0.4, unit = 'cm')
# 
# ggsave(paste0(fileSave,'DFEMean_plot_rawGR.pdf' ),DFEMean_plot_rawGR, width = width_Mean_GR, height = 4.37-0.4, unit = 'cm')
# ggsave(paste0(fileSave,'DFEvar_plot_rawGR.pdf' ),DFEVar_plot_rawGR, width = width_Mean_GR, height = 4.37-0.4, unit = 'cm')
# ggsave(paste0(fileSave,'DFEskew_plot_rawGR.pdf' ),DFESkew_plot_rawGR, width = width_Mean_GR, height = 4.37-0.4, unit = 'cm')
# 
# 
# 


# ~~ Numbers : Figure 4 ####

# mean and sd of eta dist in all envs 
df_F0_2 %>%
  group_by(env)%>%
  summarize(mean = mean(.resid),sd = sd(.resid))



# Get Rsq values
summary(lm(DFEmean~ predDFEmean, data = subset(oldDFEdat)))
summary(lm(DFEvar~ predDFEvar, data = subset(oldDFEdat)))
summary(lm(DFEskew~ predDFEskew, data = subset(oldDFEdat)))




# ggarrange(plotlist = list(DFEMean_plot_rawGR ,DFEVar_plot_rawGR,DFESkew_plot_rawGR, DFEMean_plot_adjGR ,DFEVar_plot_adjGR,DFESkew_plot_adjGR), nrow = 2, ncol = 3)

# 
# ggplot(df_Mutslopes_1slope_6int, aes(x = resVar, fill = env))+
#   geom_histogram(alpha = 0.1, position = 'identity')+
#   theme_classic()



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
##        S U P P L E M E N T ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
growthRate_data = df_mutsCanFitLine_wAdjGR %>%
  group_by(Strain, env) %>%
  summarize(Mean_GR = mean(Mean_GR), Std_err = mean(Std_err), adjGR = mean(adjGR))
colorDF = data.frame(env = c("30SC3" ,"30SC5", "30SC7" ,"37SC3", "37SC5" ,"37SC7" ,"YPD"  ), color = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'))

# checked that var is NA for all of them, so means there is only 1 val per group

# Figure S1 : Growth Rate Histograms ####
# ~~ Growth Rate Histogram, per Env ####
perEnvGRhist = function(envChoose)
{
  df = subset(growthRate_data, env == envChoose)
  col = subset(colorDF, env == envChoose)$color
  grHist = ggplot(data = df, aes(x = Mean_GR,color = env)) + 
    #stat_bin(geom="step",bins=14, position = position_dodge(width = 0.02), linewidth = 0.5)+
    stat_binline(bins = 15, aes(color = env,y = after_stat(density)),geom="step", position = position_dodge(width = 0.02),  linewidth = 0.7, alpha = 1 )+
    scale_color_manual(values = col)+
    scale_fill_manual(values =col)+
    scale_x_continuous(limits = c(0,0.4), breaks = c(0,0.2,0.4))+
    scale_y_continuous(limits = c(0,15), expand = c(0,0), breaks = c(0,7,14))+
    #geom_vline(xintercept = 0, linetype = 'dashed', color = 'black', linewidth = 1)+
    theme_classic()+
    xlab('Background Growth Rate (1/h)')+
    ylab('Density')+
    theme(axis.line =element_line(color = 'black',linewidth=0.5) ,
          axis.ticks  =element_line(color = 'black',linewidth=0.5), 
          axis.title = element_blank(),
          axis.text.x = element_text(color = 'black', size = 10, family = 'Helvetica'),
          axis.text.y = element_blank(),
          legend.position = 'none')
}

allGRhist = map(unique(growthRate_data$env), ~perEnvGRhist(.x))


allGrhist_a = ggarrange(plotlist = allGRhist, nrow = 1, ncol = 6)
#ggsave(paste0(fileSave, 'allGrhist_a.pdf'), allGrhist_a, width = 17, height = 5 , unit = 'cm')




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure S2 : Data Quality Checks ####

# ~~ Plot: Replicate Correlation (BAckgrounds and Mutations) ####
  ## In Script processing individual BC s-ests to final s-ests 


# ~~ Plot: Confirmation Experiment ####

## DATA NEEDED
allBCs_Matches = read.csv("allBarcodeMatches.csv")
transferODtable = read.csv("TransferODTable.csv")
#df_s_ests_wGR = read.csv(  paste0('df_s_ests_wGR_myGRest.csv'))


## Analysis
transferODtable$Strain = substr(transferODtable$Sample, 4,10)
transferODtable$Rep = substr(transferODtable$Sample, 1,2)

mutInfoDF = subset(transferODtable, (Rep %in% c('1-', '2-')))

transferODtable = subset(transferODtable, !(Rep %in% c('1-', '2-'))) # get rid of rows with the libary IDs 
transferODtable$env = NA
transferODtable$env[transferODtable$ID %in% 1:60] = 'e1'
transferODtable$env[transferODtable$ID %in% 61:102] = 'e2'
names(transferODtable)[names(transferODtable) == 'ID'] = 'samp'
transferODtable$samp = as.factor(transferODtable$samp)


# ~~ Add sample and Mut Info to data frame
names(transferODtable)[names(transferODtable) == 'samp'] = 'sampID'
transferODtable$sampID = as.numeric(transferODtable$samp)

df_MappedReads_use = inner_join(allBCs_Matches, transferODtable, by = 'sampID')
mutInfoDF = mutInfoDF[, c('ID', 'Strain', 'Rep')]
names(mutInfoDF) = c('allMatchedSample_name','Mut_ID', 'Rep' )
mutInfoDF$Rep = substr(mutInfoDF$Rep,1,1)
mutInfoDF$allMatchedSample = paste0('sample_', mutInfoDF$allMatchedSample_name)
mutInfoDF$Rep = paste0('R', mutInfoDF$Rep)
mutInfoDF = mutInfoDF[, c('allMatchedSample', 'Mut_ID', 'Rep')]


df_MappedReads_withMut = inner_join(df_MappedReads_use, mutInfoDF, by = c('allMatchedSample','Rep'))
# when join by matched sample AND rep, are filtering out all reads from wrong rep
# filters out ~27% of rows


# ~~ Subset to reads interested in
sub_Reads = subset(df_MappedReads_withMut, Strain %in% c('LK6-A05', 'LK2-D07')) # only using these because are only two that correctly map (that have good reference pool sequenced)
sub_Reads = sub_Reads[, c( "allDistances"  , 'allMatchedBC',   "Time",
                           'Strain', 'Rep', 'env', 'Mut_ID')]

names(sub_Reads) = c('Distance', 'Barcode', "Time", 'Strain', 'Rep', 'env', 'Mut_ID')

# ~~ Use Clustering Algorithm to get final BC_IDs
uniqueBCs = unique(sub_Reads$Barcode)
trimmedBC = substr(uniqueBCs, 1, 28) # have to trim so all are 28 bps
trimmed_uniqueBC  = trimmedBC[nchar(trimmedBC) == 28]
uniqueBCs_1 = uniqueBCs[nchar(trimmedBC) == 28] # save these bc are bcs which actually amtch with data frame 

clust = seq_cluster(dna(trimmed_uniqueBC), threshold = 0.1, method = "complete")

df_clust = data.frame(uniqueBCs_1, clust)
names(df_clust) = c('Barcode', 'BC_ID')

BCcounts_thisSamp = inner_join(sub_Reads, df_clust, by = 'Barcode')


# ~~ Get final counts for BCs
neutralMutIDs = c(102,51,99,91,6)

BCcounts_thisSamp_raw = BCcounts_thisSamp %>%
  group_by(Time, Strain, Rep, env, Mut_ID, BC_ID) %>%
  summarize(count = length(Strain), .groups = 'keep')

BCcounts_thisSamp_raw = BCcounts_thisSamp_raw %>%
  group_by(Time, Strain, Rep, env) %>%
  mutate(totCount_thisTime_Strain_rep_env = sum(count))

BCcounts_thisSamp_raw$freq = BCcounts_thisSamp_raw$count / BCcounts_thisSamp_raw$totCount_thisTime_Strain_rep_env
BCcounts_thisSamp_raw$isNeutral = BCcounts_thisSamp_raw$Mut_ID %in% neutralMutIDs


# ~~ Subset data for min initial freq, calc median neutral 
BCcounts_thisSamp_raw = subset(BCcounts_thisSamp_raw, Time != 'T-1')
BCcounts_thisSamp_raw$timeNum = substr(BCcounts_thisSamp_raw$Time,2,2)

BCcounts_thisSamp_raw = BCcounts_thisSamp_raw %>%
  group_by(Strain, BC_ID, env, Rep,Mut_ID) %>%
  mutate(initialFreq = freq[which.min(timeNum)])

sub = subset(BCcounts_thisSamp_raw,initialFreq > 1e-04 )

sub = sub %>%
  group_by(Strain, env, Rep, Time) %>%
  mutate(medianNeutral = median(count[isNeutral == T]), na.rm = T)

sub$relCount = sub$count / sub$medianNeutral


# ~ Put data into time intervals
myBC_counts_T0 = subset(sub, Time == 'T0')
myBC_counts_T1 = subset(sub, Time == 'T1')
myBC_counts_T2 = subset(sub, Time == 'T2')
myBC_counts_T3 = subset(sub, Time == 'T3')
myBC_counts_T4 = subset(sub, Time == 'T4')

myBC_counts_TimeInt_1 = full_join(myBC_counts_T0, myBC_counts_T1, by = c('env','Strain', 'Mut_ID', 'Rep',  'BC_ID','isNeutral'))
myBC_counts_TimeInt_1$timeInt = rep(1, times = length(myBC_counts_TimeInt_1$env))
myBC_counts_TimeInt_2 = full_join(myBC_counts_T1, myBC_counts_T2, by =  c('env','Strain', 'Mut_ID', 'Rep',  'BC_ID','isNeutral'))
myBC_counts_TimeInt_2$timeInt = rep(2, times = length(myBC_counts_TimeInt_2$env))
myBC_counts_TimeInt_3 = full_join(myBC_counts_T2, myBC_counts_T3, by =  c('env','Strain', 'Mut_ID', 'Rep',  'BC_ID','isNeutral'))
myBC_counts_TimeInt_3$timeInt = rep(3, times = length(myBC_counts_TimeInt_3$env))
myBC_counts_TimeInt_4 = full_join(myBC_counts_T3, myBC_counts_T4,by =  c('env','Strain', 'Mut_ID', 'Rep',  'BC_ID','isNeutral'))
myBC_counts_TimeInt_4$timeInt = rep(4, times = length(myBC_counts_TimeInt_4$env))


# ~ Put all together
my_bc_counts_full_TimeIntervals = rbind(myBC_counts_TimeInt_1,myBC_counts_TimeInt_2,myBC_counts_TimeInt_3,myBC_counts_TimeInt_4)


# ~ Get S-ests for all BCs
# ~~~ assuming 12 hour transfers, dont see her say different anywhere
my_bc_counts_full_TimeIntervals$s_est = (1/12)*log(  my_bc_counts_full_TimeIntervals$relCount.y  /my_bc_counts_full_TimeIntervals$relCount.x  )


# ~ Compare to my s-ests  from Main Exp 
df_s_est_repsPooled = my_bc_counts_full_TimeIntervals %>%
  group_by( env, Strain, Mut_ID, isNeutral ) %>%
  summarize(avgS = mean(s_est, na.rm = T), varS = var(s_est, na.rm = T), seS = sd(s_est, na.rm = T)/sqrt(sum(!is.na(s_est))))


df_s_est_repsPooled$env[df_s_est_repsPooled$env == 'e1'] = '30SC5'
df_s_est_repsPooled$env[df_s_est_repsPooled$env == 'e2'] = '37SC7'

sub_df_s_ests_wGR = subset(df_s_ests_wGR,  Strain %in% c( "LK1-C09" ,"LK1-H02", "LK2-D07", "LK5-C04" ,"LK6-A05") & Mut_ID %in% c( 6 ,  10 , 66 , 71 , 99,  117 ,127) & env %in% c('30SC5', '37SC7'))

sub_df_s_ests_wGR$Mut_ID = as.factor(sub_df_s_ests_wGR$Mut_ID)
df_s_est_repsPooled$Mut_ID = as.factor(df_s_est_repsPooled$Mut_ID)

joined = full_join(sub_df_s_ests_wGR, df_s_est_repsPooled, by = c('Strain', 'env', 'Mut_ID'))
joined = subset(joined, !is.na(avgS.x) & !is.na(avgS.y))


confExp = ggplot(joined, aes(x = avgS.x, y = avgS.y, color = env,fill = env, shape = factor(Mut_ID)))+
  geom_linerange(aes(xmin = avgS.x -  seS.x, xmax = avgS.x +  seS.x), alpha = 0.9)+
  geom_linerange(aes(ymin = avgS.y -  seS.y, ymax = avgS.y +  seS.y), alpha = 0.9)+
  geom_point(size = 2.1, alpha = 1, color = 'white')+
  geom_point(size = 2.1, alpha = 0.75)+
  scale_color_manual(values = c('#295DCC', '#FBCF96'), drop = FALSE)+
  scale_alpha_manual(values = c(0.3,0.9))+
  scale_fill_manual(values = c('#295DCC', '#FBCF96'), drop = FALSE)+
  scale_shape_manual(values = c(21,22,23,24,25,4))+
  xlab('Main Experiment')+
  ylab('Confirmation Experiment')+
  geom_abline(slope = 1, intercept = 0, color = 'grey', linewidth = 0.5)+
  theme_classic()+
  scale_x_continuous(breaks = c(-0.1,0,0.1))+
  scale_y_continuous(breaks = c(-0.1,0,0.1))+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black'),
        axis.title = element_text(color = 'black', size = 11),
        axis.text =  element_text(color = 'black', size = 11),
        legend.text = element_text(color = 'black', size = 11),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)

summary(lm(avgS.y ~ avgS.x, data = joined))
# Rsq = 0.73 when I do exlude NA in getting avgS (which is the right way to do it)

#ggsave(paste0(fileSave, 'confExp.pdf'), confExp, width = 2.8/1.35,height =3.2/1.2)




## Get barcode association file from BCcounts_thisSamp
BCassociation = BCcounts_thisSamp %>%
  group_by(Barcode, Strain, Rep, env, Mut_ID) %>%
  summarise(numTimes = length(Distance))

BCassociation = BCassociation[, names(BCassociation) != c('numTimes')]
BCassociation$env[BCassociation$env == 'e1'] = '30SC5'
BCassociation$env[BCassociation$env == 'e2'] = '37SC7'

#write.csv(BCassociation, 'BCassociation.csv')



# ~~ Plot: Correlation Correction ####
threshNumMuts_forCor = 10 
df_s_ests_wGR_sub = subset(df_mutsCanFitLine_wAdjGR ,  !is.na(avgS))
rawCor = c()
correctedCor = c()
testCor = c()
Mut_ID = c()
env = c()


sigma_sq_bg =   mean(df_s_ests_wGR_sub$Std_err^2)#
sigma_sq_mut =   mean(df_s_ests_wGR_sub$seS^2)#

for(e in unique(df_s_ests_wGR_sub$env))
{
  subEnv = subset(df_s_ests_wGR_sub, env == e)
  print(e)
  for(m in unique(subEnv$Mut_ID))
  {
    subEnv_Mut = subset(subEnv, Mut_ID == m )
    
    if(nrow(subEnv_Mut)>threshNumMuts_forCor)
    {
      
      
      cov_e_F_s = cov( subEnv_Mut$Mean_GR,subEnv_Mut$avgS)
      
      pre_var_F = var(subEnv_Mut$Mean_GR)
      pre_var_s = var(subEnv_Mut$avgS)
      
      postVar_s = pre_var_F - sigma_sq_bg
      postVar_F = pre_var_s - sigma_sq_bg -  sigma_sq_mut
      
      
      if(postVar_s<0 | postVar_F < 0)#post_corrected_var_dif_obs_exp<0 | post_corrected_var_predFit<0)
      {
        rawCor = c(rawCor,  cor( subEnv_Mut$Mean_GR,subEnv_Mut$avgS))
        correctedCor = c(correctedCor, 'corrected_varBelow0')
        Mut_ID = c(Mut_ID, m)
        env = c(env, e)
      }else{
        
        post_corrected_cor = (cov_e_F_s + sigma_sq_bg)/(sqrt((postVar_s)*(postVar_F)))
        
        
        rawCor = c(rawCor,  cor( subEnv_Mut$Mean_GR,subEnv_Mut$avgS))
        correctedCor = c(correctedCor, post_corrected_cor)
        Mut_ID = c(Mut_ID, m)
        env = c(env, e)
      }
      
      
      
    }else{
      
      rawCor = c(rawCor, 'belowPointThresh')
      correctedCor = c(correctedCor, 'belowPointThresh')
      Mut_ID = c(Mut_ID, m)
      env = c(env, e)
    }
    
    
  }
  
}


df_corCorrection = data.frame(env, Mut_ID, rawCor, correctedCor)

## look at the actual and corrected corrs

sub_df_corCorrection = subset(df_corCorrection,!(rawCor %in% 'belowPointThresh') &  !(correctedCor %in% c('belowPointThresh', 'corrected_varBelow0'))) # remove the character values 
sub_df_corCorrection$rawCor = as.numeric(sub_df_corCorrection$rawCor)
sub_df_corCorrection$correctedCor = as.numeric(sub_df_corCorrection$correctedCor)

CorrCorrelation_Plot = ggplot(sub_df_corCorrection, aes(x = rawCor, y = correctedCor, color = env))+
  geom_point(alpha = 0.5,shape =16)+
  xlab('Raw Correlation')+ylab('Corrected Correlation')+
  scale_x_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
  geom_abline(slope = 1, intercept = 0, color = 'grey', linewidth = 0.5)+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5) ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5), 
        axis.title = element_text(color = 'black', size = 11, family = 'Helvetica'),
        axis.text = element_text(color = 'black', size = 11, family = 'Helvetica'),
        legend.position = 'none')
    # ~~ remove all un-real scale corr (e.g., above or below scale 1. Retain in corr calc though)

#ggsave(paste0(fileSave,'CorrCorrelation_Plot.pdf' ),CorrCorrelation_Plot, width = 4.9359, height = 6.4953+0.2758, unit = 'cm')

summary(lm(correctedCor~rawCor,data = sub_df_corCorrection))
## Stats
length(df_corCorrection$correctedCor) # 545 cases (~94*6), of those
sum(df_corCorrection$correctedCor == 'belowPointThresh')/length(df_corCorrection$correctedCor) # 6% (34/545) have less than 10 points for doing a correlation
sum(df_corCorrection$correctedCor == 'corrected_varBelow0')/length(df_corCorrection$correctedCor) # 37% (200/545) have corrected variance below 0 --> suggesting no real correlation . Even lower than percent with N.S regressions, so results are even better with this, consistent with plo showing corrected corr is stronger magnitude 
22/length(df_corCorrection$correctedCor) # 4% have corrected corr outside -1,1





# ~~ Std_err, per Env ####
perEnvStdErrhist = function(envChoose)
{
  # using growthrate data from proces
  df = subset(growthRate_data, env == envChoose)
  col = subset(colorDF, env == envChoose)$color
  grHist = ggplot(data = df, aes(x = log10(Std_err),color = env)) + 
    #stat_bin(geom="step",bins=14, position = position_dodge(width = 0.02), linewidth = 0.5)+
    stat_binline(bins = 20, aes(color = env, y = after_stat(density)),geom="step", position = position_dodge(width = 0.02),  linewidth = 0.7, alpha = 1 )+
    scale_color_manual(values = col)+
    scale_fill_manual(values =col)+
    scale_x_continuous(limits = c(-4.2,-0.0002), breaks = c(-4,-2,0))+
    scale_y_continuous(limits = c(0,3.3), expand = c(0,0), breaks = c(0,3))+
    # scale_x_continuous(limits = c(-0.007,0.06), breaks = c(0,0.05))+
    #scale_y_continuous(limits = c(0,250), expand = c(0,0), breaks = c(0,100,200))+
    #geom_vline(xintercept = 0, linetype = 'dashed', color = 'black', linewidth = 1)+
    theme_classic()+
    xlab('Background Growth Rate (1/h)')+
    ylab('Density')+
    theme(axis.line =element_line(color = 'black',linewidth=0.5) ,
          axis.ticks  =element_line(color = 'black',linewidth=0.5), 
          axis.title = element_blank(),
          axis.text.x = element_text(color = 'black', size = 10, family = 'Helvetica'),
          axis.text.y = element_blank(),
          legend.position = 'none')
}

allSTDEhist = map(unique(growthRate_data$env), ~perEnvStdErrhist(.x))


allSTDEhist_a = ggarrange(plotlist = allSTDEhist, nrow = 1, ncol = 6)

#ggsave(paste0(fileSave, 'allSTDEhist_a.pdf'), allSTDEhist_a, width = 17, height = 3, unit = 'cm')


# ~~ Std_err, per Env ####

perEnvMutSEhist = function(envChoose)
{
  df = subset(df_mutsCanFitLine_wAdjGR, env == envChoose)
  col = subset(colorDF, env == envChoose)$color
  grHist = ggplot(data = df, aes(x = log10(seS),color = env)) + 
    #stat_bin(geom="step",bins=14, position = position_dodge(width = 0.02), linewidth = 0.5)+
    stat_binline(bins = 20, aes(color = env, y = after_stat(density)),geom="step", position = position_dodge(width = 0.02),  linewidth = 0.7, alpha = 1 )+
    scale_color_manual(values = col)+
    scale_fill_manual(values =col)+
    scale_x_continuous(limits = c(-4.2,-0.0002), breaks = c(-4,-2,0))+
    scale_y_continuous(limits = c(0,3.3), expand = c(0,0), breaks = c(0,3))+
    #  scale_x_continuous(limits = c(-0.007,0.06), breaks = c(0,0.05))+
    #scale_y_continuous(limits = c(0,250), expand = c(0,0), breaks = c(0,100,200))+
    #geom_vline(xintercept = 0, linetype = 'dashed', color = 'black', linewidth = 1)+
    theme_classic()+
    xlab('Background Growth Rate (1/h)')+
    ylab('Density')+
    theme(axis.line =element_line(color = 'black',linewidth=0.5) ,
          axis.ticks  =element_line(color = 'black',linewidth=0.5), 
          axis.title = element_blank(),
          axis.text.x = element_text(color = 'black', size = 10, family = 'Helvetica'),
          axis.text.y = element_blank(),
          legend.position = 'none')
}

allSTDE_mutshist = map(unique(growthRate_data$env), ~perEnvMutSEhist(.x))


allSTDE_mutshist_a = ggarrange(plotlist = allSTDE_mutshist, nrow = 1, ncol = 6)

#ggsave(paste0(fileSave, 'allSTDE_mutshist_a.pdf'), allSTDE_mutshist_a, width = 17, height = 3 , unit = 'cm')


# 
# ggplot(df_s_ests_wGR, aes(x = Mean_GR, y = Std_err,color = env, group = env))+
#   geom_point(alpha = 0.5, shape= 16)+
#   scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
#   geom_smooth(method = 'lm', se = F)+
#   scale_y_continuous(limits = c(0,0.025))+
#   theme_classic()

# df_s_ests_wGR %>%
#   nest(data = -c(env)) %>%
#   mutate(lmout = map(data, ~lm(Std_err~Mean_GR, data = .)), tidy = map(lmout, ~glance(.))) %>%
#   unnest(tidy)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure S3 : GxG, GxE interactions, %var in sign####

# ~~ Plot: GxG Effect Re-assortment ####

dfUse =  subset(df_mutsCanFitLine_wAdjGR ,  !is.na(avgS))

countMutSigns_GxG = dfUse %>%
  group_by(Mut_ID,env)%>%
  summarize(numStrainsBen_in = sum(mutCall == 'Beneficial' ),
            numStrainsDelt_in = sum(mutCall == 'Deleterious' ), 
            numStrainsNeutral_in = sum(mutCall == 'Neutral' ), 
            numStrainsIn = length(mutCall))

test_GxG = countMutSigns_GxG %>%
  group_by(env) %>%
  summarize(percent_alwaysBen = sum(numStrainsBen_in == numStrainsIn)/length(numStrainsBen_in),
            percent_alwaysDelt = sum(numStrainsDelt_in == numStrainsIn)/length(numStrainsBen_in),
            percent_alwaysNeutral = sum(numStrainsNeutral_in == numStrainsIn)/length(numStrainsBen_in),
            percent_ben_and_delt = sum(numStrainsBen_in != 0 & numStrainsDelt_in != 0 )/length(numStrainsBen_in)
            , nMut = length(numStrainsBen_in))


test_GxG$other = 1 - test_GxG$percent_alwaysBen -  test_GxG$percent_alwaysDelt -   test_GxG$percent_alwaysNeutral - test_GxG$percent_ben_and_delt

test_GxG = test_GxG[, c('env', 'percent_alwaysBen', 'percent_alwaysDelt','percent_alwaysNeutral', 'percent_ben_and_delt', 'other')]

test_GxG_g = gather(test_GxG, 'test', 'proportion',percent_alwaysBen:other )

test_GxG_g$test = factor(test_GxG_g$test, levels = c('percent_alwaysBen', 'percent_alwaysDelt','percent_alwaysNeutral','other', 'percent_ben_and_delt'))

propMutTypes_gxG = ggplot(test_GxG_g, aes(x = env, y = proportion, fill = test))+
  geom_bar(position = 'stack',stat="identity")+
  # scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  scale_fill_manual(values = c('#1a9b67','#9b1a4e','#6DB3F5','grey', '#e3ac5b'), drop = FALSE)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.5,1.0))+
  xlab('Environment')+
  ylab('Proportion of Mutations')+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line =element_line(color = 'black',linewidth =0.5, lineend = 'round') ,axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'),
        axis.text.y = element_text(color = 'black', size = 10) ,
        axis.text.x =  element_blank(),
        legend.position = 'none' )

range(test_GxG$percent_ben_and_delt)


#ggsave(paste0(fileSave ,'propMutTypes_gxG.pdf'), plot  = propMutTypes_gxG, width = 6.35, height =5, unit = 'cm')



# ~~ Plot: GxE Effect Re-assortment ####

countMutSigns_GxE = dfUse %>%
  group_by(Mut_ID,Strain)%>%
  summarize(numEnvsBen_in = sum(mutCall == 'Beneficial' ),
            numEnvsDelt_in = sum(mutCall == 'Deleterious' ), 
            numEnvsNeutral_in = sum(mutCall == 'Neutral' ), 
            numEnvsIn = length(mutCall))


test_GxE = countMutSigns_GxE %>%
  group_by(Strain) %>%
  summarize(percent_alwaysBen = sum(numEnvsBen_in == numEnvsIn)/length(numEnvsBen_in),
            percent_alwaysDelt = sum(numEnvsDelt_in == numEnvsIn)/length(numEnvsBen_in),
            percent_alwaysNeutral = sum(numEnvsNeutral_in == numEnvsIn)/length(numEnvsBen_in),
            percent_ben_and_delt = sum(numEnvsBen_in != 0 & numEnvsDelt_in != 0 )/length(numEnvsBen_in), 
            nmut = length(numEnvsBen_in))

range(test_GxE$percent_ben_and_delt)


test_GxE$other = 1 - test_GxE$percent_alwaysBen -  test_GxE$percent_alwaysDelt -   test_GxE$percent_alwaysNeutral - test_GxE$percent_ben_and_delt

test_GxE = test_GxE[, c('Strain', 'percent_alwaysBen', 'percent_alwaysDelt','percent_alwaysNeutral', 'percent_ben_and_delt', 'other')]

test_GxE_g = gather(test_GxE, 'test', 'proportion',percent_alwaysBen:other )

test_GxE_g$test = factor(test_GxE_g$test, levels = c('other','percent_alwaysBen', 'percent_alwaysDelt','percent_alwaysNeutral', 'percent_ben_and_delt'))


test_GxE_g = test_GxE_g %>%
  group_by(Strain) %>%
  mutate(totProp_BenDelt = sum(proportion[test == 'percent_ben_and_delt']))

test_GxE_g = test_GxE_g[rev(order(test_GxE_g$totProp_BenDelt )), ]

test_GxE_g$test = factor(test_GxE_g$test, levels = c('percent_alwaysBen', 'percent_alwaysDelt','percent_alwaysNeutral','other', 'percent_ben_and_delt'))

test_GxE_g$Strain = factor(test_GxE_g$Strain, levels = unique(test_GxE_g$Strain))


propMutTypes_GxE = ggplot(test_GxE_g, aes(x = Strain, y = proportion, fill = test))+
  geom_bar(position = 'stack',stat="identity")+
  # scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  scale_fill_manual(values = c('#1a9b67','#9b1a4e','#6DB3F5','grey', '#e3ac5b'), drop = FALSE)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0),breaks = c(0.0,0.5,1.0))+
  xlab('Strain')+
  ylab('Proportion of Mutations')+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line =element_line(color = 'black',linewidth =0.5, lineend = 'round') ,axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'),
        axis.text.y = element_blank(),
        axis.text.x =  element_blank(),
        legend.position = 'none' )

#ggsave(paste0(fileSave ,'propMutTypes_GxE.pdf'), plot  = propMutTypes_GxE, width = 11.5, height = 5, unit = 'cm')




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####


# Figure S4 : Best Microscopic Epistasis Model, Differences in Distributions of Slopes and Intercepts ####

# ~~ Analysis: Fitting Alternative Models of Global Epistasis ####

## ~~~ Constant slopes (1 intercept per env, same slope all envs)
perMut_1Slope6int_Rsq_p = df_mutsCanFitLine %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(fit = map(data, ~ lm(avgS ~ Mean_GR + env, .)), glance = map(fit, ~ glance(.x))) %>%
  unnest(glance)
## this gives rsq for full model 

perMut_1Slope6int_fits =perMut_1Slope6int_Rsq_p[, c('Mut_ID', 'fit' )]
names(perMut_1Slope6int_fits)[2] = paste0(names(perMut_1Slope6int_fits)[2], '_1s6i')


perMut_1Slope6int_Rsq_p = perMut_1Slope6int_Rsq_p[, c('Mut_ID', 'r.squared' ,'adj.r.squared', 'p.value', 'logLik')]
names(perMut_1Slope6int_Rsq_p)[2:5] = paste0(names(perMut_1Slope6int_Rsq_p)[2:5], '_1s6i')


## ~~~ 0 slope 1 int per env
perMut_0Slope6int_Rsq_p = df_mutsCanFitLine %>%
  nest(data = -Mut_ID) %>% # nest with Mut_ID (e.g nest all cols other than mutID into grouper mutID), alone this makes tibble with cols = Mut_ID and data, which stores all data of that mutID
  mutate(fit = map(data, ~ lm(avgS ~  env, .)), 
         glance = map(fit, ~ glance(.x))) %>%
  unnest(glance)

perMut_0Slope6int_fits =perMut_0Slope6int_Rsq_p[, c('Mut_ID', 'fit' )]
names(perMut_0Slope6int_fits)[2] = paste0(names(perMut_0Slope6int_fits)[2], '_0s6i')


perMut_0Slope6int_Rsq_p = perMut_0Slope6int_Rsq_p[, c('Mut_ID', 'r.squared' ,'adj.r.squared', 'p.value', 'logLik')]
names(perMut_0Slope6int_Rsq_p)[2:5] = paste0(names(perMut_0Slope6int_Rsq_p)[2:5], '_0s6i')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

# ~~ Analysis : Best Model Choice ####


## ~~~ Gather all fits into a list 
df_list <- list(perMut_0Slope6int_fits, perMut_1Slope6int_fits, perMut_6Slope6int_fits)
allFits_EpiModels = df_list %>% purrr::reduce(full_join, by='Mut_ID')

likelihood_ratio_test_0s6i_v_1s6i <- function(test,x) {
  sub = test[x,]
  return(lrtest(sub$fit_0s6i[[1]], sub$fit_1s6i[[1]])$'Pr(>Chisq)'[2])
}

likelihood_ratio_test_1s6i_v_6s6i <- function(test,x) {
  sub = test[x,]
  return(lrtest(sub$fit_1s6i[[1]], sub$fit_6s6i[[1]])$'Pr(>Chisq)'[2])
}

# Apply the function to each row of the tibble
allFits_EpiModels$rowNum = 1:length(allFits_EpiModels$Mut_ID)
resultLRT_Test <- rowwise(allFits_EpiModels) %>%
  mutate(LRT_p_0s6i_v1s6i =  likelihood_ratio_test_0s6i_v_1s6i(allFits_EpiModels,rowNum),
         LRT_p_1s6i_v6s6i =  likelihood_ratio_test_1s6i_v_6s6i(allFits_EpiModels,rowNum))

resultLRT_Test$padj_0s6i_v1s6i = p.adjust(resultLRT_Test$LRT_p_0s6i_v1s6i, method = 'BH')
resultLRT_Test$padj_1s6i_v6s6i = p.adjust(resultLRT_Test$LRT_p_1s6i_v6s6i, method = 'BH')


# Identify the best model for each mutation 
pCutoff = 0.05
resultLRT_Test$bestModel = 'none'
resultLRT_Test$bestModel[resultLRT_Test$padj_0s6i_v1s6i>pCutoff] = '1s6i' # take 0s6i model as special case of 1 slope 6 intercept model (it is just 1, size 0 slope)
resultLRT_Test$bestModel[resultLRT_Test$bestModel == 'none' & resultLRT_Test$padj_1s6i_v6s6i>pCutoff] = '1s6i'
resultLRT_Test$bestModel[resultLRT_Test$bestModel == 'none' & resultLRT_Test$padj_1s6i_v6s6i<pCutoff] = '6s6i'


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

# ~~ Plot : Best Model ####
barPLot_bestModel = ggplot( resultLRT_Test, aes(x = bestModel))+
  geom_bar(aes(y= after_stat(count)/sum(after_stat(count))*100),fill = 'white', color = 'black',  linewidth = 1)+
  scale_x_discrete(drop = FALSE)+
  scale_y_continuous(breaks = c(0,30,60))+
  theme_classic() +
  xlab('Test')+
  ylab('% Best Model')+
  theme(axis.line =element_line(color = 'black',linewidth=1) ,
        axis.ticks  =element_line(color = 'black',linewidth=1), 
        axis.title = element_text(color = 'black', size = 12, family = 'Helvetica'),
        axis.text = element_text(color = 'black', size = 10, family = 'Helvetica'))

#ggsave(paste0(fileSave,'barPLot_bestModel.pdf' ),barPLot_bestModel, width = 8.4, height = 5.1-1.5, unit = 'cm')


# Proportion best fit each model : note that if include condition for numEnvs > 2, then 6s6i is slightly more 
sum(resultLRT_Test$bestModel == '1s6i')/length(resultLRT_Test$bestModel) # 60% (56/94)
sum(resultLRT_Test$bestModel == '6s6i')/length(resultLRT_Test$bestModel) # 40% (38/94)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

# ~~ Plot : Violin of Adjusted Rsq  ####
# ~~~ Gather all Rsq and Pvals into a list 
df_list <- list(perMut_0Slope6int_Rsq_p, perMut_1Slope6int_Rsq_p, perMut_6Slope6int_Rsq_p)
allRsq_Pval_EpiModels = df_list %>% purrr::reduce(full_join, by='Mut_ID')

justRsq_EpiModels = allRsq_Pval_EpiModels[,c('Mut_ID', "adj.r.squared_1s6i", "adj.r.squared_6s6i")]
names(justRsq_EpiModels) = c('MutID', '1s6i', '6s6i')
justRsq_EpiModels_g = gather(justRsq_EpiModels, key = 'Model', value = 'Rsq', '1s6i':'6s6i')

deltaRsqPlot = ggplot(data =  justRsq_EpiModels_g, aes(x = Model, y =Rsq)) +
  geom_point( position = 'jitter', alpha = 0.1, shape = 16)+
  #geom_boxplot( colour = 'black', alpha = 0.4, linewidth = 0.4, fill = 'grey') +
  #geom_violin(fill = 'transparent', linewidth = 1)+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.7,
               colour = "black")+
  #scale_fill_manual(values = c('palegreen4','#debe92','#92b2de' ))+
  theme_classic()+
  ylab('Adjusted Rsq')+
  xlab('Model')+
  theme(axis.line =element_line(color = 'black',linewidth=1) ,
        axis.ticks  =element_line(color = 'black',linewidth=1), 
        axis.title = element_text(color = 'black', size = 12),
        axis.text = element_text(color = 'black', size = 10) ,
        legend.position = 'none')

#ggsave(paste0(fileSave,'deltaRsqPlot.pdf' ),deltaRsqPlot, width = 4.2, height = 5.1, unit = 'cm')





# ~~ Plot: Slope-Intercept Correlation, Constant Slope model ####

slopeIntCorr_1s_6int = ggplot(df_Mutslopes_1slope_6int, aes(x = slope, y = intercept, color = env, group = env))+
  geom_vline(xintercept = 0, color = 'black', linewidth = 0.75, linetype = 'dashed')+
  geom_hline(yintercept = 0, color = 'black', linewidth = 0.75, linetype = 'dashed')+
  geom_point( alpha = 0.5, size = 2, shape = 16)+
  geom_smooth(method = 'lm', se = F)+
  xlab('Slope') +  ylab('Intercept')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  # scale_x_continuous(breaks = c(-0.2,0,0.2))+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 1, color = 'black'),
        axis.ticks = element_line(linewidth = 1, color = 'black'),
        axis.title = element_text(color = 'black', size = 12, family = 'Helvetica'),
        axis.text =  element_text(color = 'black', size = 10, family = 'Helvetica'),
        legend.text = element_text(color = 'black', size = 10, family = 'Helvetica'),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)

#ggsave(paste0(fileSave,'slopeIntCorr_1s_6int.pdf' ),slopeIntCorr_1s_6int, width = 5, height = 5.1, unit = 'cm')





# ~~ Analysis: Stat Differences in Slope and Intercept Distributions ####

df_allMutSlopes_6a_6b = perMut_Env_6Slope6int_slopes_ints_spread
envs = unique(df_allMutSlopes_6a_6b$env)

env1 = c()
env2 = c()
slopePval_KS = c()
intPval_KS = c()
slopePval_t = c()
intPval_t = c()
slopePval_F = c()
intPval_F = c()

for(e1i in 1:(length(envs)-1))
{
  e1 = envs[e1i]
  sub_e1 = subset(df_allMutSlopes_6a_6b, env == e1)
  for(e2i in (e1i+1):length(envs))
  {
    e2 = envs[e2i]
    
    sub_e2 = subset(df_allMutSlopes_6a_6b, env == e2)
    
    # join to only test for same set of muts 
    joined = inner_join(sub_e1,sub_e2, by = c('Mut_ID'))
    
    slopeTest_ks = ks.test(joined$slope.x, joined$slope.y) # ftest = var.test, t.test for mean diff, ks test for overall similarlity 
    intTest_ks = ks.test(joined$intercept.x, joined$intercept.y)
    
    slopeTest_t = t.test(joined$slope.x, joined$slope.y) # ftest = var.test, t.test for mean diff, ks test for overall similarlity 
    intTest_t = t.test(joined$intercept.x, joined$intercept.y)
    
    slopeTest_f = var.test(joined$slope.x, joined$slope.y) # ftest = var.test, t.test for mean diff, ks test for overall similarlity 
    intTest_f = var.test(joined$intercept.x, joined$intercept.y)
    
    
    
    env1 = c(env1, e1)
    env2 = c(env2, e2)
    
    slopePval_KS = c(slopePval_KS, slopeTest_ks$p.value )
    intPval_KS = c(intPval_KS, intTest_ks$p.value)
    slopePval_t = c(slopePval_t,slopeTest_t$p.value)
    intPval_t = c(intPval_t,intTest_t$p.value)
    slopePval_F = c(slopePval_F,slopeTest_f$p.value)
    intPval_F = c(intPval_F,intTest_f$p.value)
    
  }
}


dfSlopeInt_ksTEst = data.frame(
  env1 ,
  env2 ,
  slopePval_KS ,
  intPval_KS ,
  slopePval_t ,
  intPval_t ,
  slopePval_F ,
  intPval_F)

dfSlopeInt_ksTEst$intPval_ks_adj = p.adjust(dfSlopeInt_ksTEst$intPval_KS , method = 'BH')
dfSlopeInt_ksTEst$slopePval_ks_adj = p.adjust(dfSlopeInt_ksTEst$slopePval_KS, method = 'BH')

dfSlopeInt_ksTEst$intPval_t_adj = p.adjust(dfSlopeInt_ksTEst$intPval_t , method = 'BH')
dfSlopeInt_ksTEst$slopePval_t_adj = p.adjust(dfSlopeInt_ksTEst$slopePval_t, method = 'BH')

dfSlopeInt_ksTEst$intPval_f_adj = p.adjust(dfSlopeInt_ksTEst$intPval_F , method = 'BH')
dfSlopeInt_ksTEst$slopePval_f_adj = p.adjust(dfSlopeInt_ksTEst$slopePval_F, method = 'BH')


dfSlopeInt_ksTEst$sigInt_KS =  dfSlopeInt_ksTEst$intPval_ks_adj<0.05
dfSlopeInt_ksTEst$sigSlope_KS = dfSlopeInt_ksTEst$slopePval_ks_adj<0.05 

dfSlopeInt_ksTEst$sigInt_t =  dfSlopeInt_ksTEst$intPval_t_adj<0.05
dfSlopeInt_ksTEst$sigSlope_t = dfSlopeInt_ksTEst$slopePval_t_adj<0.05 


dfSlopeInt_ksTEst$sigInt_f =  dfSlopeInt_ksTEst$intPval_f_adj<0.05
dfSlopeInt_ksTEst$sigSlope_f = dfSlopeInt_ksTEst$slopePval_f_adj<0.05 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

# ~~ Plot : Stat Differences in Slope and Int Distributions ####

p1_ks = ggplot(dfSlopeInt_ksTEst, aes(x = env1, y = env2, fill = sigSlope_KS))+
  geom_tile(color = 'white', linewidth = 0.5)+
  scale_fill_manual(values = c('grey','black'))+
  geom_text(aes(label=  round( slopePval_ks_adj, digits = 2) ),color = 'white', size = 2)+
  # ggtitle('Slope - KS-test')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5) ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        legend.position = 'none')


p2_ks = ggplot(dfSlopeInt_ksTEst, aes(x = env1, y = env2, fill = sigInt_KS))+
  geom_tile(color = 'white', linewidth = 0.5)+
  scale_fill_manual(values = c('grey','black'))+
  geom_text(aes(label=  round(intPval_ks_adj, digits = 2) ),color = 'white', size = 2)+
  #  ggtitle('Intercept - KS-test')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5) ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        legend.position = 'none')


p1_t = ggplot(dfSlopeInt_ksTEst, aes(x = env1, y = env2, fill = sigSlope_t))+
  geom_tile(color = 'white', linewidth = 0.5)+
  scale_fill_manual(values = c('grey','black'))+
  geom_text(aes(label=  round( slopePval_t_adj, digits = 2) ),color = 'white', size = 2)+
  # ggtitle('Slope - t-test')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5) ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        legend.position = 'none')


p2_t = ggplot(dfSlopeInt_ksTEst, aes(x = env1, y = env2, fill = sigInt_t))+
  geom_tile(color = 'white', linewidth = 0.5)+
  scale_fill_manual(values = c('grey','black'))+
  geom_text(aes(label=  round(intPval_t_adj, digits = 2) ),color = 'white', size = 2)+
  # ggtitle('Intercept - t-test')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5) ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        legend.position = 'none')

p1_f = ggplot(dfSlopeInt_ksTEst, aes(x = env1, y = env2, fill = sigSlope_f))+
  geom_tile(color = 'white', linewidth = 0.5)+
  scale_fill_manual(values = c('grey','black'))+
  geom_text(aes(label=  round( slopePval_f_adj, digits = 2) ),color = 'white', size = 2)+
  #ggtitle('Slope - F-test')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5) ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        legend.position = 'none')


p2_f = ggplot(dfSlopeInt_ksTEst, aes(x = env1, y = env2, fill = sigInt_f))+
  geom_tile(color = 'white', linewidth = 0.5)+
  scale_fill_manual(values = c('grey','black'))+
  geom_text(aes(label=  round(intPval_f_adj, digits = 2) ),color = 'white', size = 2)+
  #ggtitle('Intercept - F-test')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5) ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        legend.position = 'none')

alltests_slope_int_pairwise = ggarrange(plotlist = list(p1_ks, p1_t,p1_f,p2_ks, p2_t,p2_f), nrow = 2, ncol = 3)
# ggsave(paste0(fileSave, 'alltests_slope_int_pairwise.pdf'), alltests_slope_int_pairwise, width = 14.8,height =6, unit = 'cm' )


# ~~~ Plot of residuals vs gr ####
lm_resVar_vs_b = summary(lm(resVar ~0+slope, data =df_Mutslopes_1slope_6int )) # enforce 0,0
alpha_resVar_vs_b = -1*lm_resVar_vs_b$coefficients[1]

residualVariance_vs_slope = ggplot(df_Mutslopes_1slope_6int, aes(x = slope, y = resVar))+
  geom_point(shape = 16, alpha = 0.5)+
  geom_smooth(method = 'lm', formula = 'y~0+x', se = F, color = 'grey')+
  xlab('Slope')+
  ylab('Residual Variance')+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title = element_text(size=10, color = 'black', family = 'Helvetica'),
        axis.text =  element_text(size=8, color = 'black', family = 'Helvetica'),
        legend.position = 'none') 

#ggsave(paste0(fileSave, 'residualVariance_vs_slope.pdf'), residualVariance_vs_slope, width = 4.2, height = 6, unit = 'cm')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure S5 : All Mutation Microscopic Epistasis (w Points) ####


### Plot: all model best fits  with POINTS ###
yaxisScales = data.frame(mMin = c(1,57,81),mMax =  c(56,80,94), scaleMin = c(-0.12,-0.09,-0.05), scaleMax = c(0.12,0.09, 0.05))
# have y-axis scale vary by rows 
plotExample_wPoints = function(rowID) 
{
  m = as.numeric(joinedRanks[rowID,'Mut_ID'])
  
  ## custom y axis scale for plot
  yaxisScales$isIn = rowID >= yaxisScales$mMin & rowID <= yaxisScales$mMax
  ymin = yaxisScales$scaleMin[yaxisScales$isIn]
  ymax = yaxisScales$scaleMax[yaxisScales$isIn]
  
  fitted = subset(perMut_6Slope6int_fits, Mut_ID == m)$fit_6s6i[[1]] # get the variable slope model regression
  data = subset(df_mutsCanFitLine, Mut_ID == m)
  data$predS = fitted$fitted.values # note this only works because data here is subsetted from eaxacr same DF in exact same way as when do the LM
  data$env = factor(data$env, levels = c('30SC3', '30SC5', '30SC7', '37SC3', '37SC5', '37SC7'))
  mutName = subset(Mut_full_info_Use, Mut_ID ==m)$Gene.Use
  
  
  #  colorScale_pSigDiff = rev(c('#f8db73', '#f7bae1', '#d42bd6', '#680185'))
  
  #  colorScale_pSigDiff = rev(c('#f8db73', '#d42bd6', '#680185', 'black'))
  colorScale_pSigDiff = rev(c('black', 'black', 'black', 'black'))
  
  propSigDiff = subset(allPairwiseSlopeDiffs, Mut_ID == m)$PropSigDiff
  pSigDiffBins$isIn = propSigDiff> pSigDiffBins$min & propSigDiff<=pSigDiffBins$max
  
  colorChoice = colorScale_pSigDiff[pSigDiffBins$isIn]
  
  
  return(  ggplot(data,aes(x = Mean_GR))+
             geom_hline(yintercept = 0, color = 'grey', linetype = 'solid', linewidth = 0.5/2)+
             geom_point(aes(y = avgS, color = env), alpha = 0.2, size  = 0.4, shape = 16)+
             geom_line(aes(y = predS, color = env),   linewidth = 0.5)+
             scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = F)+
             theme_void() +
             xlab('Growth Rate (1/hr)')+
             ylab('Selection Coefficient (1/hr)')+
             # ggtitle(paste0(mutName))+
             annotate(geom="text", x=-Inf, y=-Inf, label=paste0(mutName),
                      color=colorChoice, size = 1.945, vjust = -0.6, hjust = -0.1)+
             #scale_y_continuous( breaks = c(-0.12,-0.09, -0.06,-0.03, -0.01, 0.00, 0.01,0.03, 0.06, 0.09, 0.12))+
             scale_y_continuous(limits = c(ymin,ymax), breaks = c( ymin, 0.00, ymax))+
             scale_x_continuous(limits = c(0,0.4), breaks = c(0.1,0.3))+
             theme(axis.line =element_line(color = 'black',linewidth=0.5) ,
                   axis.ticks  =element_line(color = 'black',linewidth=0.5),
                   axis.ticks.length=unit(.04, "cm"),
                   axis.title =element_blank(),
                   axis.text = element_blank(),
                   plot.title = element_text(size=5 ,vjust = -7, hjust = 0.15),
                   legend.position = 'none' ,
                   plot.margin = unit(c(0,0.1,0.25,0), 'lines') ) )
  
}

allPlots = map(joinedRanks$rowID, ~plotExample_wPoints( .x)) 
## muts 1-56, use -0.09:0.09
## muts 57:80 use -0.05:0.05
## muts 81:94 use - 0.03 : 0.03
# arrange only by slope
allPlots_slopeOrder_wPoints = ggarrange(plotlist = allPlots, nrow = 12, ncol = 8)

#ggsave(paste0(fileSave,'allPlots_slopeOrder_wPoints.pdf' ),allPlots_slopeOrder_wPoints, width = 16.5, height = 16.5, unit = 'cm')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure S6 : Growth Rate Re-assortment ####

# ~~ Plot : Correlation of GRs for all Environment Pairs ####
e1Vec = c()
e2Vec = c()
grCorr = c()

allE = unique(growthRate_data$env)
for(e1i in 1:(length(allE)-1))
{
  e1 = allE[e1i]
  sub1 = subset(growthRate_data, env == e1)
  for(e2i in (e1i+1):length(allE))
  {
    e2 = allE[e2i]
    sub2 = subset(growthRate_data, env == e2)
    joined = inner_join(sub1, sub2, by = c('Strain'))
    
    e1Vec = c(e1Vec, e1)
    e2Vec = c(e2Vec, e2)
    grCorr = c(grCorr, cor(joined$Mean_GR.x, joined$Mean_GR.y))
    
  }
}
dfGrCorr = data.frame(e1Vec, e2Vec, grCorr)

allEnvPairs = expand.grid(allE,allE)
names(allEnvPairs) = c('e1Vec', 'e2Vec')

dfGrCorr_wNA = full_join(dfGrCorr, allEnvPairs, by = c('e1Vec', 'e2Vec'))
dfGrCorr_wNA = subset(dfGrCorr_wNA, e1Vec != e2Vec & e1Vec != '37SC7' & e2Vec != '30SC3')

grCorr_tile = ggplot(dfGrCorr_wNA, aes(x = e1Vec, y = e2Vec, fill = grCorr))+
  geom_tile(color = 'white', linewidth = 0.5)+
  geom_text(aes(label=round(grCorr, digits = 2)), color = 'black', size = 3)+
  scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white', midpoint = 0, na.value = 'grey')+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 12, family = 'Helvetica'),
        axis.text = element_text(color = 'black', size = 10, family = 'Helvetica'),
        legend.position = 'right')

#ggsave(paste0(fileSave ,'grCorr_tile.pdf'), plot  = grCorr_tile, width = 10, height = 7, unit = 'cm')







# ~~ Plot : Growth Rate Re-assortment ####
##~~ add median GR
sumGR = growthRate_data %>%
  group_by(env) %>%
  mutate(medianGR = median(Mean_GR), medianAdjGR = median(adjGR))

sumGR = sumGR[order(sumGR$medianGR), ]
sumGR$env = factor(sumGR$env, levels = unique(sumGR$env))


## stats for changing side
sumGR$topSErange = sumGR$Mean_GR + sumGR$Std_err
sumGR$bottomSErange = sumGR$Mean_GR - sumGR$Std_err

sumGR$belowMedian = sumGR$topSErange < sumGR$medianGR
sumGR$aboveMedian = sumGR$bottomSErange >= sumGR$medianGR

sumGR = sumGR %>%
  group_by(Strain) %>%
  mutate(changeSideofMedian = sum(belowMedian)>0 & sum(aboveMedian)>0)

sumGR$relGR = sumGR$Mean_GR - sumGR$medianGR
sumGR$env = factor(sumGR$env, levels = c('30SC3', '30SC5', '30SC7', '37SC3', '37SC5', '37SC7'))

relGR_reassort = ggplot(sumGR, aes(x = env))+
  geom_line(aes( y = relGR,group = Strain, color = changeSideofMedian),alpha = 0.5)+
  geom_point(aes( y = relGR, group = Strain, color = changeSideofMedian),alpha = 0.5, shape = 16)+
  scale_color_manual(values = c('black', '#993254'))+#beb1b5
  geom_hline(color = 'grey', linewidth = 0.5, yintercept = 0 )+
  xlab("Environment")+
  ylab('Background Relative Growth Rate (1/h)')+
  theme_classic()+
  theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
        axis.title = element_text(color = 'black', size = 12, family = 'Helvetica'),
        axis.text = element_text(color = 'black', size = 10, family = 'Helvetica'),
        legend.position = 'none')

#ggsave(paste0(fileSave ,'relGR_reassort.pdf'), plot  = relGR_reassort, width = 8, height = 7, unit = 'cm')
sum(sumGR$changeSideofMedian)/length(sumGR$changeSideofMedian)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure S7 : Mutation Effect Re-assortment ####

## first, determine order so can have ones with different shared y axes 
minrange = 0
maxrange = 0
whichAxis = c()
Schoose_all = c()
count = 0
for(Schoose in unique(df_s_ests_wGR$Strain))
{
  count = count + 1
  sub = subset(df_s_ests_wGR, Strain == Schoose)
  
  sub = sub %>%
    group_by(env) %>%
    mutate(medianS = median(avgS))
  
  sub = sub[order(sub$medianS), ]
  sub$env = factor(sub$env, levels = unique(sub$env))
  
  
  ## stats for changing side
  sub$topSErange = sub$avgS + sub$seS
  sub$bottomSErange = sub$Mean_GR - sub$medianS
  
  
  sub$belowMedian = sub$topSErange < sub$medianS
  sub$aboveMedian = sub$bottomSErange >sub$medianS
  
  
  sub = sub %>%
    group_by(Mut_ID) %>%
    mutate(changeSideofMedian = sum(belowMedian)>0 & sum(aboveMedian)>0)
  
  sub$relS = sub$avgS - sub$medianS
  
  minrels = min(sub$relS)
  maxrels = max(sub$relS)
  
  if(minrels > -0.12 & maxrels < 0.12)
  {
    whichAxis = c(whichAxis, '0.06')
  }else{
    whichAxis = c(whichAxis, '0.12')
  }
  
  
  Schoose_all = c(Schoose_all,Schoose )
  
  
}

df_order = data.frame(Schoose_all, whichAxis)
df_order = df_order[order(df_order$whichAxis),]
df_order$Schoose_all = factor(df_order$Schoose_all, levels = unique(df_order$Schoose_all))

## now, make plots
minrange = 0
maxrange = 0
plist_strainMutrearrange = list()
plist_strainMutrearrange_relS = list()
count = 0
for(Schoose in unique(df_order$Schoose))
{
  count = count + 1
  sub = subset(df_s_ests_wGR, Strain == Schoose)
  
  sub = sub %>%
    group_by(env) %>%
    mutate(medianS = median(avgS))
  
  sub = sub[order(sub$medianS), ]
  sub$env = factor(sub$env, levels = unique(sub$env))
  
  
  ## stats for changing side
  sub$topSErange = sub$avgS + sub$seS
  sub$bottomSErange = sub$Mean_GR - sub$medianS
  
  
  sub$belowMedian = sub$topSErange < sub$medianS
  sub$aboveMedian = sub$bottomSErange >sub$medianS
  
  
  sub = sub %>%
    group_by(Mut_ID) %>%
    mutate(changeSideofMedian = sum(belowMedian)>0 & sum(aboveMedian)>0)
  
  propChangeSide =  sum(sub$changeSideofMedian)/length(sub$changeSideofMedian)

  plist_strainMutrearrange[[count]] = ggplot(sub, aes(x = env))+
    geom_line(aes( y = avgS,group = Mut_ID, color = changeSideofMedian),alpha = 0.5, lineend = 'round')+
    geom_point(aes( y = avgS), color = 'white' , alpha = 1, size = 0.5, shape = 16)+ # add underlay for points 
    geom_point(aes( y = avgS, color = changeSideofMedian),alpha = 0.5, size = 0.5, shape = 16)+
    scale_color_manual(values = c('#beb1b5', '#993254'))+
    geom_segment(aes(y = medianS, xend=after_stat(x)+0.15 ,yend=after_stat(y)),linewidth = 2, color = 'black', lineend = 'round')+
    geom_segment(aes(y = medianS, xend=after_stat(x)-0.15 ,yend=after_stat(y)),linewidth = 2, color = 'black', lineend = 'round')+
    xlab("Environment")+
    ylab('Fitness Effect (1/h)')+
    ggtitle(paste0(Schoose))+
    theme_classic()+
    theme(axis.line =element_line(color = 'black',linewidth=0.5, lineend = 'round') ,
          axis.ticks  =element_line(color = 'black',linewidth=0.5, lineend = 'round'), 
          axis.title = element_blank(),
          axis.text.x = element_text(color = 'black', size = 9, family = 'Helvetica', angle = 45, hjust = 0.9),
          axis.text.y = element_text(color = 'black', size = 9, family = 'Helvetica'),
          plot.title = element_text(color = 'black', size = 9, family = 'Helvetica'),
          legend.position = 'none')
  
  
  sub$relS = sub$avgS - sub$medianS
  sub$env = factor(sub$env, levels = c('30SC3', '30SC5', '30SC7', '37SC3', '37SC5', '37SC7'))
  
  minrange = min(minrange, min(sub$relS))
  maxrange = max(maxrange, max(sub$relS))

  
  plist_strainMutrearrange_relS[[count]] = ggplot(sub, aes(x = env))+
    geom_line(aes( y = relS,group = Mut_ID, color = changeSideofMedian),alpha = 0.5, linewidth = 0.2, lineend = 'round')+
    geom_point(aes( y = relS, group = Mut_ID), color = 'white' , alpha = 1, size = 0.35, shape = 16)+ # add underlay for points 
    geom_point(aes( y = relS, group = Mut_ID, color = changeSideofMedian, fill =changeSideofMedian ),alpha = 0.5, size = 0.4, shape = 16)+
    scale_x_discrete(drop = F)+
    scale_y_continuous(limits = c(-0.18, 0.15), breaks = c(-0.12, 0 , 0.12))+
    #scale_y_continuous( breaks = axisUse_breaks)+
    scale_color_manual(values = c('black', '#993254'))+
    geom_hline(yintercept = 0 , color = 'grey', linewidth = 0.5)+
    annotate(geom="text", x=-Inf, y=-Inf, label=paste0(round(propChangeSide, 2)*100),
             color='#993254', size = 1.945, vjust = -0.6, hjust = -0.1)+
    annotate(geom="text", x=-Inf, y=-Inf, label=paste0(', ', round(1-propChangeSide, 2)*100, '%'),
             color='black', size = 1.945, vjust = -0.6, hjust = -0.5)+
    xlab("Environment")+
    ylab('Relative Fitness Effect (1/h)')+
    ggtitle(paste0(Schoose))+#,subtitle =  paste0(round(propChangeSide, digits = 2)))+
    theme_classic()+
    theme(axis.line =element_line(color = 'black',linewidth=0.3, lineend = 'round') ,
          axis.ticks  =element_line(color = 'black',linewidth=0.3, lineend = 'round'), 
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(color = 'black',vjust = -4, hjust = 0.1, size = 9, family = 'Helvetica'),
          #plot.subtitle = element_text(color = 'black', size = 7, family = 'Helvetica'),
          legend.position = 'none',
          plot.margin = unit(c(0,0.1,0.25,0), 'lines') )
  
  
  
  
  
  
  
}

allStrains_Mutreassort = ggarrange(plotlist = plist_strainMutrearrange_relS, nrow = 7, ncol = 6)
#ggsave(paste0(fileSave ,'allStrains_Mutreassort.pdf'), plot  = allStrains_Mutreassort, width = 16.5, height = 15, unit = 'cm')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure S8 : All mutation regressions, Eta distributions ####

# ~~ Plot : All mutations regressions, each environment ####
allPlots_eta_to_0 = map(unique(df_mutsCanFitLine_wAdjGR$env), ~plotAllRegressions( .x, eta_to_0 = T, useAdjGR = T))
allPlots_realInt = map(unique(df_mutsCanFitLine_wAdjGR$env), ~plotAllRegressions( .x, eta_to_0 = F, useAdjGR = T))

allEnv_alLReg_eta0 = ggarrange(plotlist = allPlots_eta_to_0, nrow = 2, ncol = 3)
allEnv_alLReg_etaReal = ggarrange(plotlist = allPlots_realInt, nrow = 2, ncol = 3)
#ggsave(paste0(fileSave,'allEnv_alLReg_etaReal.pdf' ),allEnv_alLReg_etaReal, width = 16, height = 12, unit = 'cm')



# ~~ Plot : F0 vs Mean GR ####

subgr = df_mutsCanFitLine %>%
  group_by(Strain, env) %>%
  summarize(Mean_GR = mean(Mean_GR),)

perEnvMeanGR = subgr %>%
  group_by(env) %>%
  summarize(Mean_GR_thisEnv = mean(Mean_GR))


df_F0_wMeanGR = inner_join(df_F0,perEnvMeanGR, by = 'env' )
f0_vs_meanGR = ggplot(df_F0_wMeanGR, aes(x =Mean_GR_thisEnv, y = F0, color = env))+
  geom_abline(slope =1, intercept = 0 , color = 'grey', linewidth = 0.5, lineend = 'round')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'), drop = FALSE)+
  geom_point(size = 2, shape = 16)+
  xlab('Average GR')+
  ylab('Pivot GR')+
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
        axis.title = element_text(color = 'black', size = 10, family = 'Helvetica'),
        axis.text =  element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.text = element_text(color = 'black', size = 8, family = 'Helvetica'),
        legend.key.size = unit(0.2, 'in'),
        legend.position = 'none') #text(color = 'black', size = 10,angle = 90)

#ggsave(paste0(fileSave,'f0_vs_meanGR.pdf' ), f0_vs_meanGR, width = 6, height = 7.4, unit = 'cm')



# ~~ Analysis : get all eta (e.g residuals of slope int corr) ####
df_F0_2 = df_Mutslopes_1slope_6int %>%
  nest(data = -c(env)) %>% # nest into Mut_ID, env
  mutate(fit = map(data, ~ lm(intercept ~ 0+slope , .)), augment = map(fit, ~ augment(.x))) %>% # enforced intercept at 0,0
  unnest(augment) 


# ~~ Plot : all Eta Distributions ####
plotEtaDist = function(envChoose)
{
  sub = subset(df_F0_2, env == envChoose)
  
  plot = ggplot(sub, aes(x = .resid))+
    geom_histogram(bins =25, alpha =1, position = 'identity', aes(y = after_stat(density)), fill = '#595857')+
    geom_vline(xintercept = 0, color = 'grey', linewidth = 0.5)+
    stat_function(fun = dnorm, args = list(mean = mean(sub$.resid), sd = sd(sub$.resid)), linewidth = 1)+
    scale_y_continuous(limits = c(0,45), breaks = c(0,20,40))+
    scale_x_continuous(limits = c(-0.062,0.062), breaks = c(-0.05,0,0.05))+
    # ggtitle(paste0(envChoose))+
    theme_classic()+
    theme(axis.line = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
          axis.ticks = element_line(linewidth = 0.5, color = 'black', lineend = 'round'),
          axis.title = element_blank(), 
          axis.text =  element_blank(),
          legend.text = element_text(color = 'black', size = 8, family = 'Helvetica'),
          legend.key.size = unit(0.2, 'in'),
          legend.position = 'none') #text(color = 'black', size = 10,angle = 90)
  
  
  return(plot)
}

allEtaDist = map(unique(df_F0_2$env), ~plotEtaDist(.x))
allEtaDist_a = ggarrange(plotlist = allEtaDist, nrow = 2, ncol = 3)
#ggsave(paste0(fileSave, 'allEtaDist.pdf'), allEtaDist_a, width = 16.5, height = 5, unit = 'cm')






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure S9 : Differences in Matched DFEs ####

# for each strain, get its matched, closest in GR strain in a different env

getClosestStrain_DeltaDFE_adjGR = function(StrainChoice, envChoice, compEnvSet = c('30SC3', '30SC5', '30SC7', '37SC3', '37SC5', '37SC7'))
{
  focalStrainInfo = subset(oldDFEdat, Strain == StrainChoice & env == envChoice )
  focalStrainGR = focalStrainInfo$adjGR
  
  compEnvSet = compEnvSet[compEnvSet!= envChoice] # make sure cant compare to same env
  allPossibleComps = subset(oldDFEdat, env %in% compEnvSet & Strain !=StrainChoice ) # or to itself in diff env
  
  if(nrow(allPossibleComps)<1)
  {
    returnDF = data.frame(deltaGR =  NA,
                          deltaMean = NA, 
                          deltaVar = NA,
                          deltaSkew =NA,
                          KS_p_unjoined = NA,
                          KS_p_joined =NA)
    
    return(returnDF)
  }
  
  allPossibleComps$deltaGR = abs(allPossibleComps$adjGR - focalStrainGR )
  
  compChoice =   allPossibleComps[which.min(allPossibleComps$deltaGR), ]
  
  
  
  ### ~~~ for doing KS test ~~~ ##
  dfe1 = subset(df_s_ests_wGR, Strain == StrainChoice & env == envChoice)
  dfe2 = subset(df_s_ests_wGR, Strain == compChoice$Strain & env == compChoice$env)
  
 ks_unjoined =  ks.test(dfe1$avgS, dfe2$avgS)
  
 joined = inner_join(dfe1,dfe2, by = 'Mut_ID')
 
 ks_joined =  ks.test(joined$avgS.x, joined$avgS.y)
 
  ### ~~~~~~~~~~~~~~~~~~~~~~~~~ ## 
  
  
  returnDF = data.frame(deltaGR =  compChoice$deltaGR,
                        deltaMean = abs(compChoice$DFEmean - focalStrainInfo$DFEmean), 
                        deltaVar = abs(compChoice$DFEvar - focalStrainInfo$DFEvar),
                        deltaSkew = abs(compChoice$DFEskew - focalStrainInfo$DFEskew),
                        KS_p_unjoined = ks_unjoined$statistic,
                        KS_p_joined = ks_joined$statistic)
  
  
  
  
  
  return(returnDF)
}



getClosestStrain_DeltaDFE_rawGR = function(StrainChoice, envChoice, compEnvSet = c('30SC3', '30SC5', '30SC7', '37SC3', '37SC5', '37SC7'))
{
  focalStrainInfo = subset(oldDFEdat, Strain == StrainChoice & env == envChoice)
  focalStrainGR = focalStrainInfo$Mean_GR
  
  compEnvSet = compEnvSet[compEnvSet!= envChoice] # make sure cant compare to same env
  allPossibleComps = subset(oldDFEdat, env %in% compEnvSet & Strain !=StrainChoice ) # or to itself in diff env
  
  if(nrow(allPossibleComps)<1)
  {
    returnDF = data.frame(deltaGR =  NA,
                          deltaMean = NA, 
                          deltaVar = NA,
                          deltaSkew =NA,
                          KS_p_unjoined = NA,
                          KS_p_joined =NA)
    
    return(returnDF)
  }
  
  allPossibleComps$deltaGR = abs(allPossibleComps$Mean_GR - focalStrainGR )
  
  compChoice =   allPossibleComps[which.min(allPossibleComps$deltaGR), ]
 
   ### ~~~ for doing KS test ~~~ ##
  dfe1 = subset(df_s_ests_wGR, Strain == StrainChoice & env == envChoice)
  dfe2 = subset(df_s_ests_wGR, Strain == compChoice$Strain & env == compChoice$env)
  
  ks_unjoined =  ks.test(dfe1$avgS, dfe2$avgS)
  
  joined = inner_join(dfe1,dfe2, by = 'Mut_ID')
  
  ks_joined =  ks.test(joined$avgS.x, joined$avgS.y)
  
  ### ~~~~~~~~~~~~~~~~~~~~~~~~~ ## 
  
  
  returnDF = data.frame(deltaGR =  compChoice$deltaGR,
                        deltaMean = abs(compChoice$DFEmean - focalStrainInfo$DFEmean), 
                        deltaVar = abs(compChoice$DFEvar - focalStrainInfo$DFEvar),
                        deltaSkew = abs(compChoice$DFEskew - focalStrainInfo$DFEskew),
                        KS_p_unjoined = ks_unjoined$statistic,
                        KS_p_joined = ks_joined$statistic)
  
  return(returnDF)
}


getClosestStrain_DeltaDFE_sameStrainDiffEnv = function(StrainChoice, envChoice, compEnvSet = c('30SC3', '30SC5', '30SC7', '37SC3', '37SC5', '37SC7'))
{
  focalStrainInfo = subset(oldDFEdat, Strain == StrainChoice & env == envChoice)
  focalStrainGR = focalStrainInfo$Mean_GR
  
  compEnvSet = compEnvSet[compEnvSet!= envChoice] # make sure cant compare to same env
  allPossibleComps = subset(oldDFEdat, env %in% compEnvSet & Strain ==StrainChoice ) # ocomp to itself
  
  if(nrow(allPossibleComps)<1)
  {
    returnDF = data.frame(deltaGR =  NA,
                          deltaMean = NA, 
                          deltaVar = NA,
                          deltaSkew =NA,
                          KS_p_unjoined = NA,
                          KS_p_joined =NA)
    
    return(returnDF)
  }else{
    allPossibleComps$deltaGR = abs(allPossibleComps$Mean_GR - focalStrainGR )
    
    compChoice =   allPossibleComps[which.min(allPossibleComps$deltaGR), ]
    ### ~~~ for doing KS test ~~~ ##
    dfe1 = subset(df_s_ests_wGR, Strain == StrainChoice & env == envChoice)
    dfe2 = subset(df_s_ests_wGR, Strain == compChoice$Strain & env == compChoice$env)
    
    ks_unjoined =  ks.test(dfe1$avgS, dfe2$avgS)
    
    joined = inner_join(dfe1,dfe2, by = 'Mut_ID')
    
    ks_joined =  ks.test(joined$avgS.x, joined$avgS.y)
    
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~ ## 
    
    
    returnDF = data.frame(deltaGR =  compChoice$deltaGR,
                          deltaMean = abs(compChoice$DFEmean - focalStrainInfo$DFEmean), 
                          deltaVar = abs(compChoice$DFEvar - focalStrainInfo$DFEvar),
                          deltaSkew = abs(compChoice$DFEskew - focalStrainInfo$DFEskew),
                          KS_p_unjoined = ks_unjoined$statistic,
                          KS_p_joined = ks_joined$statistic)
    
    return(returnDF)
  }
 
}



### ~~~ ALL ENVS ##
dfeDatComps_0 = data.frame(Strain = oldDFEdat$Strain, env = oldDFEdat$env)
dfeDatComps_0 = subset(dfeDatComps_0, env != 'YPD')

dfeDatComps_adjGR = dfeDatComps_0 %>%
  rowwise() %>%
  mutate(output = list(getClosestStrain_DeltaDFE_adjGR(Strain,env)))%>%
  unnest_wider(output)
  
dfeDatComps_adjGR$type = 'adjGR'

dfeDatComps_rawGR = dfeDatComps_0 %>%
  rowwise() %>%
  mutate(output = list(getClosestStrain_DeltaDFE_rawGR(Strain,env)))%>%
  unnest_wider(output)

dfeDatComps_rawGR$type = 'rawGR'


dfeDatComps_sameStrain = dfeDatComps_0 %>%
  rowwise() %>%
  mutate(output = list(getClosestStrain_DeltaDFE_sameStrainDiffEnv(Strain,env)))%>%
  unnest_wider(output)

dfeDatComps_sameStrain$type = 'sameStrainDiffEnv'

joined_dfeComp = rbind(dfeDatComps_adjGR,dfeDatComps_rawGR,dfeDatComps_sameStrain )


## what if imposed gr diff cuttoff

# joined_dfeComp = rbind(subset(dfeDatComps_adjGR, deltaGR<0.025),
#                        subset(dfeDatComps_rawGR, deltaGR<0.025)
#                        ,dfeDatComps_sameStrain )

## ~~~~ end gr diff cuttoff 

# First, check that have very similar distributions of delta GR
# ggplot(joined_dfeComp, aes(x = deltaGR, color = type))+
#   stat_binline(bins = 15, aes(y = after_stat(density)),geom="step", linewidth = 0.7, alpha = 1 )+
#   theme_classic()
## same strian diff env has much wider dist of delta GRs, as expected
## adj and raw are very similr though


# Now, look at stats of DFE sim 
meandf = joined_dfeComp %>% 
  group_by(type) %>% 
  summarize(meanVal_u = mean(deltaMean), 
            meanSD = mean(deltaVar), 
            meanSkew = mean(deltaSkew),
            meanKS = mean(KS_p_unjoined))



pMean = ggplot(joined_dfeComp, aes(x = deltaMean, color = type))+
  stat_binline(bins = 11, aes(y = after_stat(count), color = type),geom="step", linewidth = 1, alpha = 0.6 )+
 # geom_histogram(bins = 11, position = 'identity', aes(fill = type),color = 'transparent', alpha = 0.5)+
  geom_point(data = meandf, aes(x = meanVal_u, y = 0,color = type), shape = 17, size = 3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.position = 'none', 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

pVar = ggplot(joined_dfeComp, aes(x = deltaVar, color = type))+
  stat_binline(bins = 11, aes(y = after_stat(count), color = type),geom="step", linewidth = 1, alpha = 0.6 )+
 # geom_histogram(bins = 11, position = 'identity', aes(fill = type),color = 'transparent', alpha = 0.5)+
  geom_point(data = meandf, aes(x = meanSD, y = 0,color = type), shape = 17, size = 3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.position = 'none', 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

pSkew = ggplot(joined_dfeComp, aes(x = deltaSkew, color = type))+
  stat_binline(bins = 11, aes(y = after_stat(count), color = type),geom="step", linewidth = 1, alpha = 0.6 )+
 # geom_histogram(bins = 11, position = 'identity', aes(fill = type),color = 'transparent', alpha = 0.5)+
  geom_point(data = meandf, aes(x = meanSkew, y = 0,color = type), shape = 17, size = 3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.position = 'none', 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())


pKS = ggplot(joined_dfeComp, aes(x = KS_p_unjoined, color = type))+
  stat_binline(bins = 11, aes(y = after_stat(count), color = type),geom="step", linewidth = 1, alpha = 0.6 )+
  # geom_histogram(bins = 11, position = 'identity', aes(fill = type),color = 'transparent', alpha = 0.5)+
  geom_point(data = meandf, aes(x = meanKS, y = 0,color = type), shape = 17, size = 3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.position = 'none', 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

comp_DFEs_meanGR_vs_adjGR_allenv = ggarrange(plotlist = list(pMean,pVar,pSkew,pKS), nrow = 1, ncol = 4)

#ggsave(paste0(fileSave,'comp_DFEs_meanGR_vs_adjGR_allEnvs.pdf' ),comp_DFEs_meanGR_vs_adjGR_allenv, width = 17, height = 5, unit = 'cm')





### ~~~ Just most different envs ##
dfeDatComps_0 = data.frame(Strain = oldDFEdat$Strain, env = oldDFEdat$env)
dfeDatComps_0 = subset(dfeDatComps_0, env %in% c('30SC5', '37SC7'))

dfeDatComps_adjGR = dfeDatComps_0 %>%
  rowwise() %>%
  mutate(output = list(getClosestStrain_DeltaDFE_adjGR(Strain,env, compEnvSet = c('30SC5', '37SC7'))))%>%
  unnest_wider(output)

dfeDatComps_adjGR$type = 'adjGR'

dfeDatComps_rawGR = dfeDatComps_0 %>%
  rowwise() %>%
  mutate(output = list(getClosestStrain_DeltaDFE_rawGR(Strain,env, compEnvSet = c('30SC5', '37SC7'))))%>%
  unnest_wider(output)

dfeDatComps_rawGR$type = 'rawGR'


dfeDatComps_sameStrain = dfeDatComps_0 %>%
  rowwise() %>%
  mutate(output = list(getClosestStrain_DeltaDFE_sameStrainDiffEnv(Strain,env, compEnvSet = c('30SC5', '37SC7'))))%>%
  unnest_wider(output)

dfeDatComps_sameStrain$type = 'sameStrainDiffEnv'

joined_dfeComp = rbind(dfeDatComps_adjGR,dfeDatComps_rawGR,dfeDatComps_sameStrain )



# Now, look at stats of DFE sim 
meandf = joined_dfeComp %>% 
  group_by(type) %>% 
  summarize(meanVal_u = mean(deltaMean, na.rm = T), 
            meanSD = mean(deltaVar, na.rm = T), 
            meanSkew = mean(deltaSkew, na.rm = T),
            meanKS = mean(KS_p_unjoined, na.rm = T))



pMean = ggplot(joined_dfeComp, aes(x = deltaMean))+
 stat_binline(bins = 11, aes(y = after_stat(count), color = type),geom="step", linewidth = 1, alpha = 0.6 )+
 # geom_histogram(bins = 11, position = 'identity', aes(fill = type),color = 'transparent', alpha = 0.5)+
  geom_point(data = meandf, aes(x = meanVal_u, y = 0,color = type), shape = 17, size = 3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.position = 'none', 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

pVar = ggplot(joined_dfeComp, aes(x = deltaVar, color = type))+
  stat_binline(bins = 11, aes(y = after_stat(count), color = type),geom="step", linewidth = 1, alpha = 0.6 )+
 # geom_histogram(bins = 11, position = 'identity', aes(fill = type),color = 'transparent', alpha = 0.5)+
  geom_point(data = meandf, aes(x = meanSD, y = 0,color = type), shape = 17, size = 3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.position = 'none', 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

pSkew = ggplot(joined_dfeComp, aes(x = deltaSkew, color = type))+
  stat_binline(bins = 11, aes(y = after_stat(count), color = type),geom="step", linewidth = 1, alpha = 0.6 )+
 # geom_histogram(bins = 11, position = 'identity', aes(fill = type),color = 'transparent', alpha = 0.5)+
  geom_point(data = meandf, aes(x = meanSkew, y = 0,color = type), shape = 17, size = 3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.position = 'none', 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())



pKS = ggplot(joined_dfeComp, aes(x = KS_p_unjoined, color = type))+
  stat_binline(bins = 11, aes(y = after_stat(count), color = type),geom="step", linewidth = 1, alpha = 0.6 )+
  # geom_histogram(bins = 11, position = 'identity', aes(fill = type),color = 'transparent', alpha = 0.5)+
  geom_point(data = meandf, aes(x = meanKS, y = 0,color = type), shape = 17, size = 3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.position = 'none', 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())


comp_DFEs_meanGR_vs_adjGR = ggarrange(plotlist = list(pMean,pVar,pSkew,pKS), nrow = 1, ncol = 4)

#ggsave(paste0(fileSave,'comp_DFEs_meanGR_vs_adjGR_envsMostDiff.pdf' ),comp_DFEs_meanGR_vs_adjGR, width = 17, height = 5, unit = 'cm')







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure S10 : Effect of Missing Measures on DFE moments ####
# ~~ Here, went through each environment, clustered mutations for whether have a measurement in each strain
# and then manually chose a set of 40 mutations measured in as many strains as possible


ENVCHOOSE = '37SC5'
# make matrix of 0 and 1 for if have measure for every mut in every strain
df_mut_strain = df_mutsCanFitLine[, c('Mut_ID', "Strain", 'env')]

df_mut_strain_env =  subset(df_mut_strain, env == ENVCHOOSE)[, c('Mut_ID', "Strain")]
df_mut_strain_env$x_val = 1


a = unique(df_mutsCanFitLine$Mut_ID)
b = unique(df_mutsCanFitLine$Strain)
allMuts_Strains = expand.grid(a,b)
names(allMuts_Strains) = c('Mut_ID', 'Strain')

joined_all_measured = full_join(df_mut_strain_env, allMuts_Strains, by = c('Mut_ID', 'Strain'))
joined_all_measured$x_val[is.na(joined_all_measured$x_val) ] = 0

joined_all_measured_spread = spread(joined_all_measured, key = Strain, value = x_val)

binary_matrix = joined_all_measured_spread[, !(names(joined_all_measured_spread) %in% 'Mut_ID')]
binary_matrix = as.matrix(binary_matrix)
rownames(binary_matrix) = joined_all_measured_spread$Mut_ID

## Plot unordered 
#image(as.matrix(binary_matrix), col = c("red","blue"))

# as.data.frame(binary_matrix) %>% mutate(Mut_ID = rownames(binary_matrix)) %>%
#   pivot_longer(-Mut_ID) %>%
#   ggplot(aes(x = Mut_ID, y = reorder(name, desc(name)), fill = as.factor(value)))+
#   geom_tile()+
#   scale_fill_manual(name = "Code", values = c("red","blue"))+
#   theme(axis.text.x = element_text(angle = 90))+
#   labs(y = "")


## cluster
row_clusters <- hclust(dist(binary_matrix), method = "complete")
col_clusters <- hclust(dist(t(binary_matrix)), method = "complete")
ordered_matrix <- binary_matrix[row_clusters$order, col_clusters$order]



# plot ordered 
#image(as.matrix(ordered_matrix), col = c("red","blue"))
orderedDF = as.data.frame(ordered_matrix) %>% mutate(Mut_ID = as.factor(rownames(ordered_matrix))) %>%
  pivot_longer(-Mut_ID) 
orderedDF$Mut_ID = factor(orderedDF$Mut_ID, levels = unique(orderedDF$Mut_ID))
orderedDF$name = factor(orderedDF$name, levels = unique(orderedDF$name))

ggplot(orderedDF, aes(x = Mut_ID, y = name, fill = as.factor(value)))+
  geom_tile()+
  scale_fill_manual(name = "Code", values = c("red","blue"))+
  theme(axis.text.x = element_text(angle = 90))+
  labs(y = "")


## need to choose visually 

## ~~~~ 30SC3 
#  ~ all strains
#  middels ~ 1/2 of mutations (3-88)
#paste(shQuote(levels(orderedDF$Mut_ID)[28:67]), collapse=", ") # get 40 mutations
# c('3', '121', '59', '123', '115', '105', '67', '68', '1', '85', '103', '84', '55', '127', '124', '114', '110', '108', '104', '101', '92', '81', '80', '72', '66', '58', '44', '35', '17', '15', '14', '10', '2', '7', '18', '46', '11', '62', '29', '88')

SubStrainEnv_30SC3 = subset(df_mutsCanFitLine_wAdjGR, env == '30SC3' & Mut_ID %in% c('3', '121', '59', '123', '115', '105', '67', '68', '1', '85', '103', '84', '55', '127', '124', '114', '110', '108', '104', '101', '92', '81', '80', '72', '66', '58', '44', '35', '17', '15', '14', '10', '2', '7', '18', '46', '11', '62', '29', '88'))
# this is rough subset to mostly 1's, but now get rid of the strains without all muts measured
SubStrainEnv_30SC3 = SubStrainEnv_30SC3 %>% 
  group_by(Strain) %>%
  mutate(numMutsMeasured = length(unique(Mut_ID)))

SubStrainEnv_30SC3 = subset(SubStrainEnv_30SC3, numMutsMeasured == 40) # 40 is the total put in above
length(unique(SubStrainEnv_30SC3$Strain)) # get 30 strains
#write.csv(SubStrainEnv_30SC3, paste0(fileSave, 'SubStrainEnv_30SC3.csv'))


## ~~~~ 30SC5 
#  ~ all strains
#  middle ~ 1/2 of mutations (1-62)
#paste(shQuote(levels(orderedDF$Mut_ID)[25:64]), collapse=", ") # get 40 mutations
# c('18', '127', '124', '123', '121', '115', '114', '110', '108', '105', '104', '92', '88', '82', '81', '80', '72', '70', '68', '66', '62', '58', '55', '46', '44', '35', '29', '27', '17', '15', '14', '11', '10', '7', '2', '3', '97', '21', '22', '69')

SubStrainEnv_30SC5 = subset(df_mutsCanFitLine_wAdjGR, env == '30SC5' & Mut_ID %in%  c('18', '127', '124', '123', '121', '115', '114', '110', '108', '105', '104', '92', '88', '82', '81', '80', '72', '70', '68', '66', '62', '58', '55', '46', '44', '35', '29', '27', '17', '15', '14', '11', '10', '7', '2', '3', '97', '21', '22', '69'))

# this is rough subset to mostly 1's, but now get rid of the strains without all muts measured
SubStrainEnv_30SC5 = SubStrainEnv_30SC5 %>% 
  group_by(Strain) %>%
  mutate(numMutsMeasured = length(unique(Mut_ID)))

SubStrainEnv_30SC5 = subset(SubStrainEnv_30SC5, numMutsMeasured == 40) # 40 is the total put in above
length(unique(SubStrainEnv_30SC5$Strain)) # get 36 strains
#write.csv(SubStrainEnv_30SC5, paste0(fileSave, 'SubStrainEnv_30SC5.csv'))


## ~~~~ 30SC7 
#  ~ all strains
#  top ~ 1/3 of mutations (1-62)
#paste(shQuote(levels(orderedDF$Mut_ID)[43:82]), collapse=", ") # get 40 mutations
# c('127', '124', '123', '121', '115', '114', '110', '108', '105', '104', '103', '101', '92', '88', '85', '82', '81', '80', '72', '70', '66', '62', '58', '55', '46', '44', '35', '29', '27', '18', '17', '14', '11', '10', '2', '7', '69', '3', '15', '22')

SubStrainEnv_30SC7 = subset(df_mutsCanFitLine_wAdjGR, env == '30SC7' & Mut_ID %in%  c('127', '124', '123', '121', '115', '114', '110', '108', '105', '104', '103', '101', '92', '88', '85', '82', '81', '80', '72', '70', '66', '62', '58', '55', '46', '44', '35', '29', '27', '18', '17', '14', '11', '10', '2', '7', '69', '3', '15', '22'))
# this is rough subset to mostly 1's, but now get rid of the strains without all muts measured
SubStrainEnv_30SC7 = SubStrainEnv_30SC7 %>% 
  group_by(Strain) %>%
  mutate(numMutsMeasured = length(unique(Mut_ID)))

SubStrainEnv_30SC7 = subset(SubStrainEnv_30SC7, numMutsMeasured == 40) # 40 is the total put in above
length(unique(SubStrainEnv_30SC7$Strain)) # get 40 strains
#write.csv(SubStrainEnv_30SC7, paste0(fileSave, 'SubStrainEnv_30SC7.csv'))


## ~~~~ 37SC3 
#  ~ all strains
#  top ~ 1/2 of mutations (1-62)
#paste(shQuote(levels(orderedDF$Mut_ID)[23:62]), collapse=", ") # get 40 mutations
# c('27', '85', '105', '11', '26', '82', '121', '115', '68', '88', '92', '122', '17', '55', '84', '110', '15', '80', '1', '101', '83', '46', '18', '127', '124', '114', '108', '81', '72', '62', '44', '35', '29', '2', '10', '7', '104', '123', '50', '3')

SubStrainEnv_37SC3 = subset(df_mutsCanFitLine_wAdjGR, env == '37SC3' & Mut_ID %in%  c('27', '85', '105', '11', '26', '82', '121', '115', '68', '88', '92', '122', '17', '55', '84', '110', '15', '80', '1', '101', '83', '46', '18', '127', '124', '114', '108', '81', '72', '62', '44', '35', '29', '2', '10', '7', '104', '123', '50', '3'))
# this is rough subset to mostly 1's, but now get rid of the strains without all muts measured
SubStrainEnv_37SC3 = SubStrainEnv_37SC3 %>% 
  group_by(Strain) %>%
  mutate(numMutsMeasured = length(unique(Mut_ID)))

SubStrainEnv_37SC3 = subset(SubStrainEnv_37SC3, numMutsMeasured == 40) # 40 is the total put in above
length(unique(SubStrainEnv_37SC3$Strain)) # get 23 strains
#write.csv(SubStrainEnv_37SC3, paste0(fileSave, 'SubStrainEnv_37SC3.csv'))


## ~~~~ 37SC5 
#  bottom ~ 1/2 of mutations (1-62)
#paste(shQuote(levels(orderedDF$Mut_ID)[8:47]), collapse=", ") # get 40 mutations
# c('92', '82', '105', '84', '26', '55', '67', '123', '122', '83', '50', '58', '127', '27', '115', '18', '121', '1', '104', '101', '7', '85', '72', '81', '62', '80', '46', '110', '68', '15', '124', '114', '108', '44', '35', '29', '17', '2', '10', '60')


SubStrainEnv_37SC5 = subset(df_mutsCanFitLine_wAdjGR, env == '37SC5' & Mut_ID %in%  c('92', '82', '105', '84', '26', '55', '67', '123', '122', '83', '50', '58', '127', '27', '115', '18', '121', '1', '104', '101', '7', '85', '72', '81', '62', '80', '46', '110', '68', '15', '124', '114', '108', '44', '35', '29', '17', '2', '10', '60'))

# this is rough subset to mostly 1's, but now get rid of the strains without all muts measured
SubStrainEnv_37SC5 = SubStrainEnv_37SC5 %>% 
  group_by(Strain) %>%
  mutate(numMutsMeasured = length(unique(Mut_ID)))

SubStrainEnv_37SC5 = subset(SubStrainEnv_37SC5, numMutsMeasured == 40) # 40 is the total put in above
length(unique(SubStrainEnv_37SC5$Strain)) # get 16 strains
#write.csv(SubStrainEnv_37SC5, paste0(fileSave, 'SubStrainEnv_37SC5.csv'))




### ~~~~~ 37SC7 
# bottom ~half of strains + first ~1/2 of mutations
# c('LK4-B12', 'LK2-E08', 'LK2-D05', 'LK4-B01', 'LK5-F08', 'LK5-H12', 'LK3-C04', 'LK2-B07', 'LK3-B08', 'LK3-C12', 'LK5-D09', 'LK2-D07', 'LK6-A05', 'LK1-C09', 'LK4-H11', 'LK2-A10', 'LK5-G01', 'LK1-B05', 'LK5-G03', 'LK2-A06', 'LK3-H11')
# c('82', '88', '121', '123', '84', '58', '92', '50', '27', '29', '3', '55', '1', '11', '83', '105', '62', '72', '18', '104', '101', '114', '35', '124', '108', '81', '44', '17', '15', '10', '2', '7', '26', '127', '70', '80', '85', '46', '122', '103')

#paste(shQuote(levels(orderedDF$name)[1:21]), collapse=", ")
#paste(shQuote(levels(orderedDF$Mut_ID)[1:50]), collapse=", ")

SubStrainEnv_37SC7 = subset(df_mutsCanFitLine_wAdjGR,env == '37SC7' &  Mut_ID %in% c('82', '88', '121', '123', '84', '58', '92', '50', '27', '29', '3', '55', '1', '11', '83', '105', '62', '72', '18', '104', '101', '114', '35', '124', '108', '81', '44', '17', '15', '10', '2', '7', '26', '127', '70', '80', '85', '46', '122', '103'))
# this is rough subset to mostly 1's, but now get rid of the strains without all muts measured
SubStrainEnv_37SC7 = SubStrainEnv_37SC7 %>% 
  group_by(Strain) %>%
  mutate(numMutsMeasured = length(unique(Mut_ID)))
SubStrainEnv_37SC7 = subset(SubStrainEnv_37SC7, numMutsMeasured == 40) # 40 is the total put in above
length(unique(SubStrainEnv_37SC7$Strain)) # get 16 strains
#write.csv(SubStrainEnv_37SC7, paste0(fileSave, 'SubStrainEnv_37SC7.csv'))




### ~~~~~~ PUT ALL TOGETHER

allSubDat = rbind(SubStrainEnv_30SC3,SubStrainEnv_30SC5,SubStrainEnv_30SC7,SubStrainEnv_37SC3,SubStrainEnv_37SC5,SubStrainEnv_37SC7)


dfeQualCheck_DF = allSubDat %>%
  group_by(Strain, env, Mean_GR, adjGR, Std_err, F0) %>%
  summarize(DFEmean = mean(avgS, na.rm= T), DFEvar = var(avgS, na.rm = T), DFEskew = skewness(avgS, na.rm = T))

dfeQualCheck_DF$predDFEmean = dfe_parampred(as.numeric(dfeQualCheck_DF$Mean_GR), dfeQualCheck_DF$env)$DFEmean
dfeQualCheck_DF$predDFEvar = dfe_parampred(as.numeric(dfeQualCheck_DF$Mean_GR), dfeQualCheck_DF$env)$DFEvar
dfeQualCheck_DF$predDFEskew = dfe_parampred(as.numeric(dfeQualCheck_DF$Mean_GR), dfeQualCheck_DF$env)$DFEskew


dfeMean_qualCheck = ggplot(dfeQualCheck_DF, aes(x = adjGR,color = env))+
  #geom_linerange(aes(ymin = DFEmean_min, ymax = DFEmean_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEmean), alpha = 0.5,lineend='round')+
  geom_point(aes(y = DFEmean ),shape = 16, color = 'white')+
  geom_point(aes(y = DFEmean ),shape = 16, alpha = 0.5)+
  # geom_smooth(method = 'lm', se = F)+
  geom_line(aes(y = predDFEmean), color = 'black', linewidth = 1,lineend='round')+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Mean')+
  #ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(-0.03,0,0.03))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  = element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )




dfeVar_qualCheck = ggplot(dfeQualCheck_DF, aes(x = adjGR ,color = env))+
  #  geom_linerange(aes(ymin = DFEvar_min, ymax = DFEvar_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEvar), alpha = 0.5,lineend='round')+
  geom_point(aes( y = DFEvar ),shape = 16, color = 'white')+# underlay so dont see cross of error bar in point
  geom_point(aes( y = DFEvar), shape = 16,alpha = 0.5)+
  # geom_smooth(method = 'lm', se = F)+
  geom_line(aes(y = predDFEvar), color = 'black', linewidth = 1,lineend='round')+
  #geom_point(data = predfestatextrapt, aes(x = x, y = DFEvar), color = 'black', fill = 'white', shape = c(22,23,24), size = 3, stroke = 1.3)+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Variance')+
  # ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(0,0.001,0.002))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )




dfeSkew_qualCheck = ggplot(dfeQualCheck_DF, aes(x = adjGR,color = env))+
  #geom_linerange(aes(ymin = DFEskew_min, ymax = DFEskew_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEskew), alpha = 0.5,lineend='round')+
  geom_point(aes( y = DFEskew ), shape = 16,color = 'white')+# underlay so dont see cross of error bar in point
  geom_point(aes( y = DFEskew ),shape = 16, alpha = 0.5)+
  # geom_smooth(method = 'lm', se = F)+
  geom_line(aes(y = predDFEskew), color = 'black', linewidth = 1,lineend='round')+
  #geom_point(data = predfestatextrapt, aes(x = x, y = DFEskew), color = 'black', fill = 'white', shape = c(22,23,24), size = 3, stroke = 1.3)+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Skewness')+
  #ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(-3,0,3))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )


# ggsave(paste0(fileSave,'dfeMean_qualCheck.pdf' ),dfeMean_qualCheck, width = 3 , height =3, unit = 'cm')
# ggsave(paste0(fileSave,'dfeVar_qualCheck.pdf' ),dfeVar_qualCheck, width = 3, height = 3, unit = 'cm')
# ggsave(paste0(fileSave,'dfeSkew_qualCheck.pdf' ),dfeSkew_qualCheck, width = 3, height = 3, unit = 'cm')




# ~~~ Number of Measurements per Strain cor with GR #####
checkNumMutsMeasured = df_mutsCanFitLine_wAdjGR %>%
  group_by(Strain, env, Mean_GR, adjGR) %>%
  summarize(numMutsMeasured = length(unique(Mut_ID)))


plotNumMutPerStrain = ggplot(checkNumMutsMeasured, aes(x = Mean_GR, y =numMutsMeasured, color = env ))+
  geom_point(shape = 16, alpha = 0.5)+
  geom_smooth(method = 'lm', se = F)+
  xlab("Background GR")+
  ylab('Number of Measured Mutations')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        legend.position = 'none' )
# ggsave(paste0(fileSave,'plotNumMutPerStrain.pdf' ),plotNumMutPerStrain, width = 3, height = 3, unit = 'cm')



checkNumMutsMeasured %>%
  nest(data = -(env))%>%
  mutate(lmout = map(data, ~lm(numMutsMeasured~Mean_GR, .)), glance = map(lmout, ~glance(.))) %>%
  unnest(glance)

# ~~~ Hist of number of mutations per strain ####
dfe_dat_test = df_mutsCanFitLine_wAdjGR %>%
  group_by(Strain, env, Mean_GR, adjGR, Std_err) %>%
  summarize(DFEmean = mean(avgS, na.rm= T), DFEvar = var(avgS, na.rm = T), DFEskew = skewness(avgS, na.rm = T), nPoints = sum(!is.na(avgS)))


numMutsPerStrain = ggplot(dfe_dat_test, aes(x = nPoints))+
  stat_binline(bins = 9, fill = 'transparent', aes(y = after_stat(count)))+
  theme_classic()+
  geom_vline(xintercept = median(dfe_dat_test$nPoints), linewidth = 0.5, color = 'grey')+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )

#ggsave(paste0(fileSave,'numMutsPerStrain.pdf' ),numMutsPerStrain, width = 3, height = 3, unit = 'cm')


quantile(dfe_dat_test$nPoints, prob = c(0.25,0.5,0.75))  # gives lower and upper quants with median in middle


# ~~~ Reduced Strain Set ####
#remove all strains with < 60 muts 


dfe_dat_test = subset(dfe_dat_test, nPoints >= 60)

dfe_dat_test$predDFEmean = dfe_parampred(as.numeric(dfe_dat_test$Mean_GR), dfe_dat_test$env)$DFEmean
dfe_dat_test$predDFEvar = dfe_parampred(as.numeric(dfe_dat_test$Mean_GR), dfe_dat_test$env)$DFEvar
dfe_dat_test$predDFEskew = dfe_parampred(as.numeric(dfe_dat_test$Mean_GR), dfe_dat_test$env)$DFEskew



p1 =ggplot(dfe_dat_test, aes(x = adjGR,color = env))+
  #geom_linerange(aes(ymin = DFEmean_min, ymax = DFEmean_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEmean), alpha = 0.5,lineend='round')+
  geom_point(aes(y = DFEmean ),shape = 16, color = 'white')+
  geom_point(aes(y = DFEmean ),shape = 16, alpha = 0.5)+
  # geom_smooth(method = 'lm', se = F)+
  geom_line(aes(y = predDFEmean), color = 'black', linewidth = 1,lineend='round')+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Mean')+
  #ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(-0.03,0,0.03))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  = element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )



p2 = ggplot(dfe_dat_test, aes(x = adjGR ,color = env))+
  #  geom_linerange(aes(ymin = DFEvar_min, ymax = DFEvar_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEvar), alpha = 0.5,lineend='round')+
  geom_point(aes( y = DFEvar ), shape = 16,color = 'white')+# underlay so dont see cross of error bar in point
  geom_point(aes( y = DFEvar), shape = 16,alpha = 0.5)+
  # geom_smooth(method = 'lm', se = F)+
  geom_line(aes(y = predDFEvar), color = 'black', linewidth = 1,lineend='round')+
  #geom_point(data = predfestatextrapt, aes(x = x, y = DFEvar), color = 'black', fill = 'white', shape = c(22,23,24), size = 3, stroke = 1.3)+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Variance')+
  # ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(0,0.001,0.002))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )


p3 =ggplot(dfe_dat_test, aes(x = adjGR,color = env))+
  #geom_linerange(aes(ymin = DFEskew_min, ymax = DFEskew_max), alpha = 0.5,lineend='round')+
  geom_linerange(aes(xmin = adjGR - Std_err, xmax = adjGR + Std_err, y = DFEskew), alpha = 0.5,lineend='round')+
  geom_point(aes( y = DFEskew ),shape = 16, color = 'white')+# underlay so dont see cross of error bar in point
  geom_point(aes( y = DFEskew ),shape = 16, alpha = 0.5)+
  # geom_smooth(method = 'lm', se = F)+
  geom_line(aes(y = predDFEskew), color = 'black', linewidth = 1,lineend='round')+
  #geom_point(data = predfestatextrapt, aes(x = x, y = DFEskew), color = 'black', fill = 'white', shape = c(22,23,24), size = 3, stroke = 1.3)+
  xlab('Adjusted Background Growth Rate (1/h)')+
  ylab('DFE Skewness')+
  #ggtitle('With all mutations measured')+
  scale_color_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'),drop = F)+
  scale_fill_manual(values = c('#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96', '#6a9c67'), drop = F)+
  scale_x_continuous(limits = c(-0.2,0.25), breaks = c(-0.2,0.0,0.2))+
  scale_y_continuous( breaks = c(-3,0,3))+
  theme_classic()+
  theme(axis.line = element_line(color = 'black',linewidth=0.5,lineend='round') ,
        axis.ticks  =element_line(color = 'black',linewidth=0.5,lineend='round'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none' )

#ggarrange(plotlist = list(p1,p2,p3), nrow =1, ncol = 3)

# ggsave(paste0(fileSave,'dfeMean_qualCheck_reducedStrains.pdf' ),p1, width = 3 , height =3, unit = 'cm')
# ggsave(paste0(fileSave,'dfeVar_qualCheck_reducedStrains.pdf' ),p2, width = 3, height = 3, unit = 'cm')
# ggsave(paste0(fileSave,'dfeSkew_qualCheck_reducedStrains.pdf' ),p3, width = 3, height = 3, unit = 'cm')







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# Figure S11 : QTL analysis ####

kre33_data =  read.csv(paste0( 'KRE33_pres_abs.csv'))
# has 0,1 for all strains (by vs rm allele)
# for the 4 multi-hit qtl sites seen both in miulos and Jerison paper

testDat = inner_join(oldDFEdat, kre33_data, by = 'Strain')

## ~~ on Strain GRs ####
qtls = c('KRE33_chr14_376315','chr14_470303','chr12_646707','chr15_154799')

allQTL_percVarInGR_exp = data.frame()
for(q in qtls)
{
  percExp_GR_byQTL= c()
  env = c()
  for(e in unique(testDat$env))
  {
    data = subset(testDat,env == e)
    data_2 = data[, names(data) %in% c('env', 'Mean_GR', q)]
    names(data_2)[3] = 'QTL'
    
    out = summary(lm(Mean_GR ~ QTL  , data = data_2 ) )
    percExp_GR_byQTL = c(percExp_GR_byQTL, out$adj.r.squared)
    env = c(env,e)
  }
  sub1 = data.frame(env, percExp = percExp_GR_byQTL, qtl = as.character(q))
  allQTL_percVarInGR_exp = rbind(allQTL_percVarInGR_exp, sub1)
}

allQTL_percVarInGR_exp$percExp[allQTL_percVarInGR_exp$percExp < 0 ] = 0 # just set to 0 if negative

allQTL_percVarInGR_exp = allQTL_percVarInGR_exp %>%
  group_by(qtl) %>%
  mutate(avgPexp = mean(percExp))

allQTL_percVarInGR_exp = allQTL_percVarInGR_exp[rev(order(allQTL_percVarInGR_exp$avgPexp)), ]
allQTL_percVarInGR_exp$qtl = factor(allQTL_percVarInGR_exp$qtl, levels = unique(allQTL_percVarInGR_exp$qtl))


## ~~ For each Mutation ####
testDat_muts = inner_join(df_mutsCanFitLine_wAdjGR, kre33_data, by = 'Strain')

## look at all qtls added 
getp_val_allQTLs_aboveGR = function(data)
{
  # KRE33_chr14_376315
  # chr12_646707
  # chr14_470303
  # chr15_154799
  
  if(length(unique(data$KRE33_chr14_376315)) < 2 |length(unique(data$chr12_646707)) < 2 |length(unique(data$chr14_470303)) < 2  |length(unique(data$chr15_154799)) < 2 | length(data$Strain)<0)
  {
    return(data.frame(sig_GR = NA, P_explained_GR = NA, P_exp_allQTLS_together = NA))
  }else{
    lmout = (lm(avgS ~ Mean_GR + env +  as.factor(KRE33_chr14_376315)+ as.factor(chr12_646707)+ as.factor(chr14_470303)+ as.factor(chr15_154799), data = data)) 
    af = anova(lmout)
    afss <- af$"Sum Sq"
    percExpOut = (cbind(af,PctExp=afss/sum(afss)*100)) 
    
     ## only add signficant ones 
    return(data.frame(sig_GR = percExpOut$`Pr(>F)`[1], 
                       P_explained_GR = percExpOut$PctExp[1]*(percExpOut$`Pr(>F)`[1]<0.05),  
                      P_explained_Env = (percExpOut$PctExp[rownames(percExpOut) == 'env'])*(percExpOut$`Pr(>F)`[rownames(percExpOut) == 'env']<0.05),
                      P_exp_allQTLS_together = sum(percExpOut$PctExp[!(rownames(percExpOut) %in% c('Mean_GR', 'env', 'Residuals'))] * (percExpOut$`Pr(>F)`[!(rownames(percExpOut) %in% c('Mean_GR', 'env', 'Residuals'))]<0.05 ))  ))
    
    
  }
  
}

allQTls_additive_sig_perMut_env_aboveGR = testDat_muts %>% 
  nest(data = -c(Mut_ID)) %>%
  mutate(sig_df = map(data, ~getp_val_allQTLs_aboveGR(.))[[1]]) %>% 
  unnest(sig_df)

allQTls_additive_sig_perMut_env_aboveGR$p_adj_GR = p.adjust(allQTls_additive_sig_perMut_env_aboveGR$sig_GR, method = 'BH')
allQTls_additive_sig_perMut_env_aboveGR = allQTls_additive_sig_perMut_env_aboveGR[rev(order(allQTls_additive_sig_perMut_env_aboveGR$P_exp_allQTLS_together)), ]
allQTls_additive_sig_perMut_env_aboveGR$Mut_ID = factor(allQTls_additive_sig_perMut_env_aboveGR$Mut_ID, levels = unique(allQTls_additive_sig_perMut_env_aboveGR$Mut_ID))

sub = allQTls_additive_sig_perMut_env_aboveGR[,c('Mut_ID', 'P_explained_GR','P_explained_Env', 'P_exp_allQTLS_together')]

## put in mut names
Mut_full_info_Use$Mut_ID = as.factor(Mut_full_info_Use$Mut_ID)
sub = inner_join(sub, Mut_full_info_Use, by = 'Mut_ID')


## put Env+GR = %var explained by Eq2
sub2 = sub
sub2$P_exp_eq2 = sub$P_explained_Env + sub$P_explained_GR
sub2 = sub2[,c('Mut_ID','Gene.Use', 'P_exp_eq2' , 'P_exp_allQTLS_together')]

gathered2 = gather(sub2,key = 'Type', value = 'PctExp',  P_exp_eq2:P_exp_allQTLS_together)
gathered2 = gathered2 %>%
  group_by(Mut_ID) %>%
  mutate(totPctExp = sum(PctExp))
gathered2 = gathered2[rev(order(gathered2$totPctExp)),]

gathered2$Mut_ID = factor(gathered2$Mut_ID, levels = unique(gathered2$Mut_ID))
gathered2$Gene.Use = factor(gathered2$Gene.Use, levels = unique(gathered2$Gene.Use))

gathered2= subset(gathered2, !is.na(totPctExp))
qtlPlot = ggplot(gathered2, aes(x = Gene.Use, y = PctExp, fill = Type))+
  geom_col(color = 'white', linewidth = 0.2)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c( '#c24674', '#46c294'), drop = F)+
  theme(axis.text.x = element_text(angle = 90, size = 5, vjust = 0.5, color = 'black'), 
        axis.text.y = element_text(size = 5),
        legend.position = 'none')
#ggsave(paste0(fileSave, 'qtlPlot.pdf'), qtlPlot, width = 18, height = 6.4, unit = 'cm')


## ~~ Example Mutations ####
## 93, 95, 78

## 93
data = subset(testDat_muts, Mut_ID == 93)
lmout = (lm(avgS ~ Mean_GR + env +  as.factor(KRE33_chr14_376315)+ as.factor(chr12_646707)+ as.factor(chr14_470303)+ as.factor(chr15_154799), data = data)) 
af = anova(lmout)
afss <- af$"Sum Sq"
percExpOut = (cbind(af,PctExp=afss/sum(afss)*100))
# highest effect qtl is KRE33_chr14_376315

## 95
data = subset(testDat_muts, Mut_ID == 95)
lmout = (lm(avgS ~ Mean_GR + env +  as.factor(KRE33_chr14_376315)+ as.factor(chr12_646707)+ as.factor(chr14_470303)+ as.factor(chr15_154799), data = data)) 
af = anova(lmout)
afss <- af$"Sum Sq"
percExpOut = (cbind(af,PctExp=afss/sum(afss)*100)) 
# highest effect qtl is KRE33_chr14_376315


## 78
data = subset(testDat_muts, Mut_ID == 78)
lmout = (lm(avgS ~ Mean_GR + env +  as.factor(KRE33_chr14_376315)+ as.factor(chr12_646707)+ as.factor(chr14_470303)+ as.factor(chr15_154799), data = data)) 
af = anova(lmout)
afss <- af$"Sum Sq"
percExpOut = (cbind(af,PctExp=afss/sum(afss)*100)) # env and GR are both NS
# highest effect qtl is chr15_154799


pList = list()
mutsChoose = c(93, 95, 78)
corespQTLS = c('KRE33_chr14_376315', 'KRE33_chr14_376315', 'chr15_154799')

for(i in 1:length(mutsChoose))
{
  mChoose = mutsChoose[i]# 80
  #eChoose = sub_sigKre33[i,]$env #'37SC3'
  name = paste0(mChoose, ': ', subset(Mut_full_info_Use, Mut_ID == mChoose )$Gene.Use)
  qtlChoice = corespQTLS[i]
  
  sub = subset(testDat_muts, Mut_ID == mChoose)
  sub = sub[,names(sub) %in% c('env', 'Mut_ID', 'avgS', 'Mean_GR', qtlChoice)]
  names(sub)[5]= 'QTL'
  
  ## get predS
  best = '1s6i'
  fitted = resultLRT_Test[resultLRT_Test$Mut_ID == mChoose, paste0('fit_', best)][[1]][[1]]
  data = subset(df_mutsCanFitLine, Mut_ID == mChoose)
  data$predS = fitted$fitted.values # note this only works because data here is subsetted from eaxacr same DF in exact same way as when do the LM
  data$env = factor(data$env, levels = c('30SC3', '30SC5', '30SC7', '37SC3', '37SC5', '37SC7'))

  pList[[i]] = ggplot(sub, aes(x = Mean_GR, y = avgS))+
    geom_point(size = 1,  aes(color = env, shape =   as.factor(QTL), alpha = as.factor(QTL)))+
    # scale_color_manual(values = c('#f6d53c', '#3c5df6'))+
    geom_line(data = data, aes(x = Mean_GR, y = predS, color = env))+
    scale_color_manual(values = c( '#00026E','#295DCC','#90B6DD', '#C21700', '#FE8D26', '#FBCF96'),drop = F)+
    scale_alpha_discrete(range = c(0.34, 0.9))+
    ggtitle(paste0( name, ', ', qtlChoice))+
    geom_hline(yintercept = 0, color = 'grey', linewidth = 1.2)+
    theme_classic()+
    theme(legend.position = 'none', axis.title = element_blank(),
          axis.text = element_text(color = 'black', size = 8))
}

exampleQTLs = ggarrange(plotlist = pList, nrow = 1, ncol = 3)
#ggsave(paste0(fileSave, 'exampleQTLs.pdf'), exampleQTLs, width= 17, height = 4, unit = 'cm')




