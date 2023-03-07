# importing libraries
library(ranger)
library(data.table)
library(dplyr)
library(caret)
library(performanceEstimation)
source("function.R")
source("parse-vcf-snv-wrapper.R")

######### generating predictions ##########
# loading random forest model
load("final_model.rds")

# parsing data
# parse.df = parse.snv(sample.dir='data/real2_part2')
# write.table(parse.df, 'output/snv-parse-real2_part2.txt', row.names = F, quote = F, sep = '\t')

# preprocessing test data (real2p2)
test <- preprocess.2("output/snv-parse-real2_part2.txt") %>%
  mutate(feature = paste(Chr, START_POS_REF, sep= ".")) 

# creating feature subset for test data (real2p2)
test.m1 <- subset(test, select = c(Mutect2,Freebayes,Vardict,Varscan,FILTER_Mutect2,FILTER_Freebayes,FILTER_Vardict,FILTER_Varscan,m2_MQ,f_MQMR,vs_SSC,vs_SPV,vd_SSF,vd_MSI))

# generating predictions for test data (real2p2)
probabilities <- predict(rf.m1, data = test.m1, type='response', verbose = TRUE)$predictions[,"1"]
output <- test[,1:3]
output$Prediction[probabilities > .5] = 1
output$Prediction[probabilities <= .5] = 0

# saving predictions (prediction for all positions)
fwrite(output,"prediction/team4_real2_part2.predictions_submission.bed", sep = '\t')

# saving predictions (positives only, following format of given truth files)
mutations.only <- output[output$Prediction == 1,1:3]
fwrite(mutations.only,"prediction/team4_real2_part2.predictions_mutations.only_submission.bed", sep = '\t')


######## f1 stats #########
# generating f1 stats for test data (real2p2)
# truth <- fread( ) # truth file for real2p2
# f1stats <- calc.F1(mutations.only, truth[,1:3])
