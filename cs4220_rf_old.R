#code unused in final implementation of our rf model

#N/A values data processing
preprocess.df <- function(file){
  df <- fread(file)
  ## filling N/A values using prior knowledge
  # med <- median(df$vs_SSC, na.rm = T)
  # df <- df %>%
  #   mutate(m2_MQ = ifelse(is.na(m2_MQ), ifelse(!is.na(f_MQMR), f_MQMR, 0), m2_MQ),
  #          f_MQMR = ifelse(is.na(f_MQMR), ifelse(!is.na(m2_MQ), m2_MQ, 0), f_MQMR)) %>%
  #   mutate(vs_SSC = ifelse(is.na(vs_SSC), med, vs_SSC),
  #          vs_SPV = ifelse(is.na(vs_SPV), 1, vs_SPV),
  #          vd_SSF = ifelse(is.na(vd_SSF), 1, vd_SSF),
  #          vd_MSI = ifelse(is.na(vd_MSI), 0, vd_MSI),
  #          FILTER = FILTER_Mutect2 + FILTER_Freebayes + FILTER_Vardict + FILTER_Varscan,
  #          avgMQ = (m2_MQ+f_MQMR)/2) 
  ## filling N/A values using mean of column
  df <- as.data.frame(df)
  for(i in 1:ncol(df)){
    df[is.na(df[,i]), i] <- mean(df[,i], na.rm = TRUE)
  }
  return(setDT(df))
}

preprocess.2 <- function(file){
  ## filling N/A values using prior knowledge, generating columns for caller detection
  df <- read.csv(file)
  med <- median(df$vs_SSC, na.rm = T)
  ## obtaining VC
  x <- data.frame(t(data.frame(strsplit(as.character(df$REF_MFVdVs), "/"), row.names =c("Mutect2","Freebayes","Vardict","Varscan"))))
  rownames(x) <- c(1:nrow(x))
  x[!(is.na(x)|x=="NA")] = TRUE
  x[(is.na(x)|x=="NA")] = FALSE
  for(i in 1:4){x[,i] <- as.logical(x[,i])}
  df <- cbind(df,x)
  df <- df %>%
    mutate(m2_MQ = ifelse(is.na(m2_MQ), ifelse(!is.na(f_MQMR), f_MQMR, 0), m2_MQ),
           f_MQMR = ifelse(is.na(f_MQMR), ifelse(!is.na(m2_MQ), m2_MQ, 0), f_MQMR)) %>%
    mutate(vs_SSC = ifelse(is.na(vs_SSC), med, vs_SSC),
           vs_SPV = ifelse(is.na(vs_SPV), 1, vs_SPV),
           vd_SSF = ifelse(is.na(vd_SSF), 1, vd_SSF),
           vd_MSI = ifelse(is.na(vd_MSI), 0, vd_MSI),
           FILTER = FILTER_Mutect2 + FILTER_Freebayes + FILTER_Vardict + FILTER_Varscan,
           avgMQ = (m2_MQ+f_MQMR)/2) 
  return(df)
}

#fitting rf model (regular rf as opposed to probability tree growing in our final model)
rf.og <- ranger(formula = y ~ ., data = og, y = y, num.trees = ntrees, importance = 'permutation')
rf.m1 <- ranger(formula = y ~ ., data = m1, y = y, num.trees = ntrees, importance = 'permutation')
rf.m2 <- ranger(formula = y ~ ., data = m2, y = y, num.trees = ntrees, importance = 'permutation')

#predicting test data (using regular rf model)
pred.og <- predict(rf.og, data = real2p1og)
pred.m1 <- predict(rf.m1, data = real2p1m1)
pred.m2 <- predict(rf.m2, data = real2p1m2)

#confusion matrix (using regular rf model)
confusionMatrix(as.factor(creal2p1$label), pred.og$predictions)
confusionMatrix(as.factor(creal2p1$label), pred.m1$predictions)
confusionMatrix(as.factor(creal2p1$label), pred.m2$predictions)