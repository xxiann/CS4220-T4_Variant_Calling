#importing libraries
library(ranger)
library(data.table)
library(dplyr)
library(caret)

#parameters
ntrees <- 1000
seed <- 0
set.seed(seed)

#loading data
# psyn1 <- fread("C:/Users/zhich/Downloads/extracolumns-20230226T083050Z-001/extracolumns/syn1-processed.txt")
# psyn2 <- fread("C:/Users/zhich/Downloads/extracolumns-20230226T083050Z-001/extracolumns/syn2-processed.txt")
# psyn3 <- fread("C:/Users/zhich/Downloads/extracolumns-20230226T083050Z-001/extracolumns/syn3-processed.txt")
# psyn4 <- fread("C:/Users/zhich/Downloads/extracolumns-20230226T083050Z-001/extracolumns/syn4-processed.txt")
# psyn5 <- fread("C:/Users/zhich/Downloads/extracolumns-20230226T083050Z-001/extracolumns/syn5-processed.txt")
preal1 <- fread("C:/Users/zhich/Downloads/version2-20230304T174310Z-001/version2/real1-processed.txt")
# preal2p1 <- fread("C:/Users/zhich/Downloads/extracolumns-20230226T083050Z-001/extracolumns/real2_part1-processed.txt")

# truth1 <- fread("C:/Users/zhich/Downloads/syn1/syn1_truth.bed")
# truth2 <- fread("C:/Users/zhich/Downloads/syn2/syn2_truth.bed")
# truth3 <- fread("C:/Users/zhich/Downloads/syn3/syn3_truth.bed")
# truth4 <- fread("C:/Users/zhich/Downloads/syn4/syn4_truth.bed")
# truth5 <- fread("C:/Users/zhich/Downloads/syn5/syn5_truth.bed")
truth6 <- fread("C:/Users/zhich/Downloads/real1/real1_truth.bed")
# truth7 <- fread("C:/Users/zhich/Downloads/real2_part1/real2_truth_chr1to5.bed")

cols = c('label')

#data preprocessing
# psyn1$feature <- paste(psyn1$Chr, psyn1$START_POS_REF, sep= ".")
# truth1$feature <- paste(truth1$V1, truth1$V2, sep= ".")
# truth1$label <- 1
# truth1 = subset(truth1, select = -c(V1,V2,V3) )
# psyn1 = subset(psyn1, select = -c(Chr,START_POS_REF,END_POS_REF) )
# combined1 <- left_join(psyn1, truth1)
# combined1[ , (cols) := lapply(.SD, nafill, fill=0), .SDcols = cols]
# 
# psyn2$feature <- paste(psyn2$Chr, psyn2$START_POS_REF, sep= ".")
# truth2$feature <- paste(truth2$V1, truth2$V2, sep= ".")
# truth2$label <- 1
# truth2 = subset(truth2, select = -c(V1,V2,V3) )
# psyn2 = subset(psyn2, select = -c(Chr,START_POS_REF,END_POS_REF) )
# combined2 <- left_join(psyn2, truth2)
# combined2[ , (cols) := lapply(.SD, nafill, fill=0), .SDcols = cols]
# 
# psyn3$feature <- paste(psyn3$Chr, psyn3$START_POS_REF, sep= ".")
# truth3$feature <- paste(truth3$V1, truth3$V2, sep= ".")
# truth3$label <- 1
# truth3 = subset(truth3, select = -c(V1,V2,V3) )
# psyn3 = subset(psyn3, select = -c(Chr,START_POS_REF,END_POS_REF) )
# combined3 <- left_join(psyn3, truth3)
# combined3[ , (cols) := lapply(.SD, nafill, fill=0), .SDcols = cols]
# 
# psyn4$feature <- paste(psyn4$Chr, psyn4$START_POS_REF, sep= ".")
# truth4$feature <- paste(truth4$V1, truth4$V2, sep= ".")
# truth4$label <- 1
# truth4 = subset(truth4, select = -c(V1,V2,V3) )
# psyn4 = subset(psyn4, select = -c(Chr,START_POS_REF,END_POS_REF) )
# combined4 <- left_join(psyn4, truth4)
# combined4[ , (cols) := lapply(.SD, nafill, fill=0), .SDcols = cols]
# 
# psyn5$feature <- paste(psyn5$Chr, psyn5$START_POS_REF, sep= ".")
# truth5$feature <- paste(truth5$V1, truth5$V2, sep= ".")
# truth5$label <- 1
# truth5 = subset(truth5, select = -c(V1,V2,V3) )
# psyn5 = subset(psyn5, select = -c(Chr,START_POS_REF,END_POS_REF) )
# combined5 <- left_join(psyn5, truth5)
# combined5[ , (cols) := lapply(.SD, nafill, fill=0), .SDcols = cols]

preal1$feature <- paste(preal1$Chr, preal1$START_POS_REF, sep= ".")
truth6$feature <- paste(truth6$V1, truth6$V2, sep= ".")
truth6$label <- 1
truth6 = subset(truth6, select = -c(V1,V2,V3) )
preal1 = subset(preal1, select = -c(Chr,START_POS_REF,END_POS_REF) )
combined6 <- left_join(preal1, truth6)
combined6[ , (cols) := lapply(.SD, nafill, fill=0), .SDcols = cols]

# preal2p1$feature <- paste(preal2p1$Chr, preal2p1$START_POS_REF, sep= ".")
# truth7$feature <- paste(truth7$V1, truth7$V2, sep= ".")
# truth7$label <- 1
# truth7 = subset(truth7, select = -c(V1,V2,V3) )
# preal2p1 = subset(preal2p1, select = -c(Chr,START_POS_REF,END_POS_REF) )
# combined7 <- left_join(preal2p1, truth7)
# combined7[ , (cols) := lapply(.SD, nafill, fill=0), .SDcols = cols]

#overall combined dataframe (syn1-5,real1p1, real2p1)
# pall <- rbindlist(list(combined1,combined2,combined3,combined4,combined5,combined6,combined7))
fwrite(combined6, "C:/Users/zhich/Downloads/preal1.txt")

#creating feature subsets
pall <- fread("C:/Users/zhich/Downloads/preal1.txt")
m1 <- subset(pall, select = c(Mutect2,Freebayes,Vardict,Varscan,FILTER_Mutect2,FILTER_Freebayes,FILTER_Vardict,FILTER_Varscan,m2_MQ,f_MQMR,vs_SSC,vs_SPV,vd_SSF,vd_MSI))
# m2 <- subset(pall, select = c(CALL_mutect2,CALL_freebayes,CALL_vardict,CALL_varscan,vs_SSC,vs_SPV,vd_SSF,vd_MSI,FILTER,avgMQ))
# og <- subset(pall, select = -c(FILTER,avgMQ,feature,label))

#truth labels
y <- as.factor(pall$label)

#generating weights to balance sampling
w <- 1/table(y)
w <- w/sum(w)
weights <- rep(0, nrow(pall))
weights[y == 0] <- w['0']
weights[y == 1] <- w['1']
table(weights, y)

#fitting random forest model
start = Sys.time()
rf.m1 <- ranger(formula = y ~ ., data = m1, y = y, num.trees = ntrees, seed = seed, importance = 'permutation', replace = TRUE, splitrule = "hellinger", keep.inbag = TRUE, save.memory = TRUE, probability=TRUE, case.weights = weights)
print( Sys.time() - start )

#saving random forest model
save(rf.m1,file="C:/Users/zhich/Downloads/rf.m1.preal1.caseweight.RData")

#loading random forest model
load("C:/Users/zhich/Downloads/rf.m1.preal1.RData")

#preprocessing test data (real2p2)
preal2p1 <- fread("C:/Users/zhich/Downloads/version2-20230304T174310Z-001/version2/real2_part1-processed.txt")
preal2p1$feature <- paste(preal2p1$Chr, preal2p1$START_POS_REF, sep= ".")
preal2p1 = subset(preal2p1, select = -c(Chr,START_POS_REF,END_POS_REF) )

#preprocessing truth labels **only for testing on syn1-5, real1p1, real2p1 data**
treal2p1 <- fread("C:/Users/zhich/Downloads/real2_part1/real2_truth_chr1to5.bed")
treal2p1$feature <- paste(treal2p1$V1, treal2p1$V2, sep= ".")
treal2p1$label <- 1
treal2p1 = subset(treal2p1, select = -c(V1,V2,V3) )
creal2p1 <- left_join(preal2p1, treal2p1)
cols = c('label')
creal2p1[ , (cols) := lapply(.SD, nafill, fill=0), .SDcols = cols]

#creating feature subset for test data (real2p2)
real2p1m1 <- subset(preal2p1, select = c(Mutect2,Freebayes,Vardict,Varscan,FILTER_Mutect2,FILTER_Freebayes,FILTER_Vardict,FILTER_Varscan,m2_MQ,f_MQMR,vs_SSC,vs_SPV,vd_SSF,vd_MSI))

#generating predictions for test data (real2p2)
probabilities <- as.data.frame(predict(rf.m1, data = real2p1m1, type='response', verbose = TRUE)$predictions)
pred.real2p2 <- max.col(probabilities)-1

#saving predictions as .bed file
# output.real2p2 = data.frame(Chr=preal2p2$Chr,START_POS_REF=preal2p2$START_POS_REF,END_POS_REF=preal2p2$END_POS_REF,Prediction=pred.real2p2)
# fwrite(output.real2p2,"C:/Users/zhich/Downloads/real2_part2.predictions.bed", sep = '\t')

#confusion matrix **only for testing on syn1-5, real1p1, real2p1 data**
confusionMatrix(table(max.col(probabilities)-1,creal2p1$label))