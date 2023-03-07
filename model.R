# importing libraries
library(ranger)
library(data.table)
library(dplyr)
library(caret)
library(performanceEstimation)
source("function.R")

# parameters
ntrees <- 1000
seed <- 0
set.seed(seed)

########### Loading training data ############
# loading parsed data & data preprocessing
preal1 <- preprocess.2("output/snv-parse-real1.txt") 
preal2p1 <- preprocess.2("output/snv-parse-real2_part1.txt")

truth6 <- fread("data/real1/real1_truth.bed") %>%
  mutate(label = 1) 
truth7 <- fread("data/real2_part1/real2_truth_chr1to5.bed") %>%
  mutate(label = 1) 

combined6 <- left_join(preal1, truth6, by = c("Chr"="V1", "START_POS_REF"="V2", "END_POS_REF"="V3")) %>%
  mutate(label = ifelse(is.na(label), 0, label)) %>%
  mutate(label = factor(label))
combined7 <- left_join(preal2p1, truth7, by = c("Chr"="V1", "START_POS_REF"="V2", "END_POS_REF"="V3")) %>%
  mutate(label = ifelse(is.na(label), 0, label)) %>%
  mutate(label = factor(label))

# overall combined dataframe (syn1-5,real1p1, real2p1)
pall <- rbind(combined6,combined7)

# save processed data
# fwrite(pall, "train.txt", sep = '\t')

########## Training random forest model ##########
# load processed data
# pall <- fread("train.txt")
# pall$label <- as.factor(pall$label)

# truth labels
y <- as.factor(m1$label)

# creating feature subsets, removing labels from training data
m1 <- subset(m1, select = c(Mutect2,Freebayes,Vardict,Varscan,FILTER_Mutect2,FILTER_Freebayes,FILTER_Vardict,FILTER_Varscan,m2_MQ,f_MQMR,vs_SSC,vs_SPV,vd_SSF,vd_MSI))

# fitting random forest model
start = Sys.time()
rf.m1 <- ranger(formula = y ~ ., data = m1, y = y, num.trees = ntrees, seed = seed, importance = 'permutation', replace = TRUE, splitrule = "hellinger", keep.inbag = TRUE, save.memory = TRUE, probability=TRUE)
print( Sys.time() - start )

# saving random forest model
save(rf.m1,file="final_model.rds")

