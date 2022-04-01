library(lme4)
library(stats)

ParentPath = '/Users/firasmidani/Box_Duke/project_davidlab/LAD_LAB_Personnel/Firas/FIRAS_EXPERIMENTS/_2022_03_31_bacteroides_data_analysis'
setwd(ParentPath)

tab = read.table('./tables/pipe_data_to_R.txt',sep='\t',header=TRUE)
tab = subset(tab, Substrate_Group == "Complex Sugars" | Substrate_Group == "Simple Sugars")

# Likelihood Ratio Tests with Mixed Effects Model
m_0= lmer('SimpsonEvenness ~ log10(CorrectedMixedCount) + 1|TimePoint',data=tab,REML=FALSE)
m_1 = lmer('SimpsonEvenness ~ Substrate_Group + 1|TimePoint',data=tab,REML=FALSE)
m_full = lmer('SimpsonEvenness ~ log10(CorrectedMixedCount) + Substrate_Group + 1|TimePoint',data=tab,REML=FALSE)

anova(m_0,m_full) 
anova(m_1,m_full) 

# Logistic Regression
m_lg = glm(formula = 'SimpsonEvenness ~ log10(CorrectedMixedCount) + TimePoint', data=tab,family=gaussian(link = 'logit' ))
summary(m_lg)

#tab = read.table('./tables/max_counts_joined.txt',sep='\t',header=TRUE)
#tab = subset(tab, Substrate_Group == "Complex Sugars" | Substrate_Group == "Simple Sugars")
#tab

#lmer('log10(MaxCount) ~ Culture_Type + Substrate_Group + 1|Substrate',data=tab,REML=FALSE)
