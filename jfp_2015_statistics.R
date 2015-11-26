# filename:  jfp_2015_statistics.R
# description:  descriptions for statsitical tests performed for our 2015 JFP manuscript in R
# author:  Laura M. Carroll; lmc297@cornell.edu

# Abbreviations for specific data sets used in this study:
# pc45 = plate counts 45:  differential plating over time at 45C in logCFU/g
# pc40 = plate counts 40:  differential plating over time at 40C in logCFU/g
# rpoH45 = rpoH log copy numbers over time  at 45C
# rpoH40 = rpoH log copy numbers over time  at 40C
# dnaK45 = dnaK log copy numbers over time  at 45C
# dnaK40 = dnaK log copy numbers over time  at 40C
# clpB45 = clpB log copy numbers over time  at 45C
# clpB40 = clpB log copy numbers over time  at 40C
# ibpA45 = ibpA log copy numbers over time  at 45C
# ibpA40 = ibpA log copy numbers over time  at 40C
# htpG45 = htpG log copy numbers over time  at 45C
# htpG40 = htpG log copy numbers over time  at 40C
# degP45 = degP log copy numbers over time  at 45C
# degP40 = degP log copy numbers over time  at 40C

# Percent Injury
# Used to calculate the percent of the population that is injured relative to the total population via differential plating
# Used for pc45 and pc40
# Make sure to that XLD and TSA colony counts are in CFU/g, not log CFU/g
# To get the proportion of the population that is "healthy" at a time point:
# healthy_proportion = (XLD CFU/g)/(TSA CFU/g)
# To get the proportion of the population that is "injured" at a time point:
# injured_proportion = 1-healthy_proportion
# To convert that to a percent:
# percent_injury=100*injured_proportion

# Fold Change
# Used to convert log copy numbers to fold change values
# data = data frame with "time", "trial", and "log_copy" values
# First convert log copy numbers to copy numbers
data$copy=2^data$log_copy
# To get fold change, divide copy number at each time by the copy number of the corresponding trial at 0 min
# For example, divide copy number of [trial A, 5 min] by copy number of [trial A, 0 min]
# Divide copy number of [trial C, 240 min] by copy number of [trial C, 0 min]
# Put these values in a colunn called data$fc_uncorrected
# Now correct these values to not bias them upward; we will call this column data$fc_corrected
data$fc_corrected<-ifelse((data$fc_uncorrected<1), -1/(data$fc_uncorrected), data$fc_uncorrected)

# Bartlett's Test for homogeneity of variances
# Used to test for equal variances for pc45, pc40, rpoH45, rpoH40, dnaK45, dnaK40, clpB45, clpb40, ibpA45, ibpA40, htpG45, htpG40, degP45, degP40
# Example usage:
# data = data frame with "Time" and "LogCount" columns
# In my study, "Time" took the value 0, 5, 10, 15, 30, 60, 90, 180, or 240 (min)
# In my study, "LogCount" took log CFU/g values for differential plating data, and log copy numbers for RT-qPCR (genetic) data 
Time = as.factor(data$Time) # treat time as a factor to compare data variances at discrete points
bartlett.test(data$LogCount~Time, data = data)
# If p is less than alpha, reject null; variances time point(s) is/are unequal
# If p is greater than alpha; do no reject null; variances at all time points are equal
# In our case, the null hypothesis of equal variances could not be rejected in any of the data sets (P>0.05)
# It is also possible to use Levene's Test; I used leveneTest in the 'car' package
library(car)
leveneTest(data$LogCount~data$Time, data = data)
# leveneTest (unsurprisingly) gave higher P values
#I used bartlett.test to be more confident that, should there have been a departure from equal variances, I would find it

# One-Way Dunnett'Test for statitical significance of injury via differential plating
# Dunnett's Test corrects for multiple comparisons
# Used to test for significantly lower log CFU/g counts on selective media (XLD) at various time points compared to non-selective media (TSA) as a measure of injury
# Use for pc45, pc40
# Example Usage:
# data = data frame with "Time", "LogCount", "Trial", and "TimeCode" columns
# "Time" took the value 0, 5, 10, 15, 30, 60, 90, 180, or 240 (min)
# "LogCount" took log CFU/g values for TSA media
# "Trial": a letter corresponding to the biological replicate (biological rep 1 = trial "A", biological rep 2 = trial "B", etc.)
# "TimeCode" took the values "TSA", "XLDq", where q = 0, 5, 10, 15, 30, 60, 90, 180, or 240.  This was used to compare XLD at each time point to TSA at 0 min.
# First build a model using lme4 that accounts for the fact that samples were independent, but came from a batch culture
library(lme4)
data.lmer = lmer(LogCount ~ TimeCode + (1|Trial), data=data)
# Use multcomp to perform a one-way Dunnett's test, comparing log CFU/g values at each XLD time point to TSA log CFU/g values
library(multcomp)
data.dunnett = summary(glht(data.lmer, linfct=mcp(TimeCode = "Dunnett"), alternative = "less"))
summary(data.dunnett)
# If P < alpha, reject H0; XLD counts at that time point are significantly lower than TSA counts
# If P > alpha, do not reject H0; XLD counts at that time point are not significantly lower than TSA counts

# Two-Way Dunnett's Tests for changes in transcript level relative to 0 min
# Dunnett's test corrects for multiple comparisons
# Used to test for significantly higher/lower trancript levels relative to initial transcript levels at 0 min
# Used for rpoH45, rpoH40, dnaK45, dnaK40, clpB45, clpB40, ibpA45, ibpA40, htpG45, htpG40, degP45, degP40
# Example usage:
# data = data frame with "Time", "LogCount", and "Trial" columns
# "Time" takes value of 0, 5, 10, 15, 30, 60, 90, 180, or 240 (min)
# "LogCount" takes log copy number at that time point
# "Trial" takes character corresponding to the batch culture that each particular sample came from ("A", "B", etc.)
library(lme4)
library(multcomp)
# First build a model to account for the fact that samples were independent, but came from a batch culture
data.lmer = lmer(LogCount ~ Time + (1|Trial), data=data)
data.dunnett = summary(glht(data.lmer, linfct=mcp(TimeFactor = "Dunnett")))
summary(data.dunnett)
# If P < alpha, reject H0; transcript levels at that time point are significantly different than initial transcript levels
# If P > alpha, do not reject H0; transcript levels at that time point are not significantly different than initial transcript levels

# Copy ratio for modeling
# data = data frame with "time", "trial", and "log_copy" values
# First convert log copy numbers to copy numbers
data$copy=2^data$log_copy
# To get copy ratio, divide copy number at each time by the copy number of the corresponding trial at 0 min
# For example, divide copy number of [trial A, 5 min] by copy number of [trial A, 0 min]
# Divide copy number of [trial C, 240 min] by copy number of [trial C, 0 min]
# Put these values in a colunn called data$copy_ratio
# This allows us to set the value of A in the nonlinear model (below) to 1

# Nonlinear Model Creation
# Used to fit a nonlinear model to copy ratio numbers (described above in "Copy ratio for modeling") if genetic data meets model assumptions
# Rylander, M., Y. Feng, K. Zimmermann, and K. Diller. 2010. Measurement and mathematical modeling of thermally induced injury and heat shock protein expression kinetics in normal and cancerous prostate cells. International Journal of Hyperthermia. 26:748-764.
# Model assumptions:  1) increase in transcript level, 2) decrease to initial levels
# We were conservative with our assumptions; we mandated that: 
# 1) log copy numbers must show a statistically significant increase via the Two-Way Dunnett's test described above (P<0.05)
# 2) A decrease in log copy number within 240 min to levels not statistically different from initial levels via the Two-Way Dunnett's test described above (P>0.05)
# It is, however, possible to model average increases and decreases, even if not statistically significant
# ibpA at 45 and ibpA at 40 met these assumptions
# data = data frame with "time" and "copy_ratio" columns
# First load nls2 package
library(nls2)
# Now create a grid of values to test for parameters alpha, beta, and gamma
parameter_grid = expand.grid(list(alpha = seq(-5,5,by=0.01), beta = seq(-5,5,by=0.01), gamma = seq(0,3,by=0.01)))
# This can become computationally expensive if you don't have a good idea of initial values
model.nls2 = nls2(data$copy_ratio~exp(alpha*data$time-beta*(data$time^gamma)), data = data, start = parameter_grid, algorithm = "brute-force", nls.control(warnOnly = TRUE), trace = TRUE)
# This can take a while; you may need to adjust initial values
summary(model.nls2)

# Linear Model Creation
# Used to fit a linear model to copy ratio numbers (described above in "Copy ratio for modeling")
# data = data frame with "time" and "copy_ratio" columns
model.lm = lm(data$copy_ratio~data$time, data = data)
summary(model.lm)

# Intercept-Only Model Creation
# Used to fit a linear intercept-only model to copy ratio numbers (described above in "Copy ratio for modeling")
# data = data frame with "time" and "copy_ratio" columns
model.intercept.lm = lm(data$copy_ratio~1, data = data)
summary(model.intercept.lm)

# Corrected Akaike Information Criteria (AICc) for Nonlinear, Linear, and Intercept-Only Models
# AICc is used for model selection; the lower the AICc, the better the model
# Wagenmakers, E. J., and S. Farrell. 2004. AIC model selection using Akaike weights. Psychonomic Bulletin & Review. 11:192-196.
# Load AICcmodavg package
library(AICcmodavg)
# Now run AICc on each model
aicnls=AICc(model.nls2)
aiclm=AICc(model.lm)
aicint=AICc(model.intercept.lm)
# In this example, let's say model.nls2 has the lowest AICc value!
# To get AICc weights:
# First get deltaAICc
deltaAICc_A=aiclm-aicnls # subtract best-scoring model from others; in this example, nls2 is best
deltaAICc_B=aicint-aicnls # subtract best-scoring model from others
# Now to get weight of best model:
AICcw.nls = (exp(-0.5*0))/(exp(-0.5*0)+exp(-0.5*deltaAICc_A)+exp(-0.5*deltaAICc_B))
AICcw.nls
# To get weight of lm:
AICcw.lm = (exp(-0.5*deltaAICc_A))/(exp(-0.5*0)+exp(-0.5*deltaAICc_A)+exp(-0.5*deltaAICc_B))
AICcw.lm
# To get weight of intercept-only:
AICcw.int = (exp(-0.5*deltaAICc_B))/(exp(-0.5*0)+exp(-0.5*deltaAICc_A)+exp(-0.5*deltaAICc_B))
AICcw.int

# Bayesian Information Criteria (BIC) for Nonlinear, Linear, and Intercept-Only Models
# BIC is used for model selection; the lower the BIC, the better the model
# BIC penalizes more heavily for additional parameters
# Wagenmakers, E. J., and S. Farrell. 2004. AIC model selection using Akaike weights. Psychonomic Bulletin & Review. 11:192-196.
# To get BIC values:
BIC(model.nls2, model.intercept.lm)
BIC(model.nls2, model.lm)
# To get delta BIC values and BIC weights, substitute AICc values above with BIC values
