# importing libraries
library(ranger)
library(data.table)
library(dplyr)
source("function.R")
source("parse-vcf-snv-wrapper.R")

######### generating predictions ##########
# loading random forest model
load("final_model.RData")

# parsing data
# parse.df = parse.snv(sample.dir='data/real2_part2')
# write.table(parse.df, 'output/snv-parse-real2_part2.txt', row.names = F, quote = F, sep = '\t')

# preprocessing test data (real2p2)
test <- preprocess.2("output/snv-parse-real2_part2.txt") 

# creating feature subset for test data (real2p2)
test.m1 <- subset(test, select = c(Mutect2,Freebayes,Vardict,Varscan,FILTER_Mutect2,FILTER_Freebayes,FILTER_Vardict,FILTER_Varscan,m2_MQ,f_MQMR,vs_SSC,vs_SPV,vd_SSF,vd_MSI))

# generating predictions for test data (real2p2)
probabilities <- as.data.frame(predict(rf.m1, data = test.m1, type='response', verbose = TRUE)$predictions)
pred <- max.col(probabilities)-1
output <- test[,1:3]
output$Prediction = pred

# saving predictions (prediction for all positions)
fwrite(output,"prediction/team4_real2_part2.predictions_submission.bed", sep = '\t')

# saving predictions (positives only, following format of given truth files)
mutations.only <- output[output$Prediction == 1,1:3]
fwrite(mutations.only,"prediction/team4_real2_part2.predictions_mutations.only_submission.bed", sep = '\t')


######## f1 stats #########
# generating f1 stats for test data (real2p2)
# truth <- fread( ) # truth file for real2p2
# truth$label <- 1
# truth <- left_join(test, truth, by=c("Chr"="V1", "START_POS_REF"="V2", "END_POS_REF"="V3"))
# cols = c('label')
# truth[ , (cols) := lapply(.SD, nafill, fill=0), .SDcols = cols]
# f1stats <- calc.F1(as.data.frame(mutations.only), as.data.frame(truth[,c(1:3,23)]))
