calc.F1 = function(pred,truth) {
  
  # append 'Chr' and 'StartPos'
  predv = paste(pred[,1],pred[,2])
  truthv = paste(truth[,1],truth[,2])
  
  res = data.frame(matrix(nrow = 1,ncol = 6))
  colnames(res) = c('TP','FP','FN','Precision','Recall','F1')
  
  res$TP = sum(predv %in% truthv)
  res$FP = sum(!(predv %in% truthv))
  res$FN = sum(!(truthv %in% predv))
  
  res$Precision = res$TP/(res$TP + res$FP)
  res$Recall = res$TP/(res$TP + res$FN)
  res$F1 = (2*res$Precision*res$Recall)/(res$Precision + res$Recall)
  
  return(res)
}

preprocess.2 <- function(file){
  df <- fread(file)
  
  # obtaining medium values for each row
  med <- median(df$vs_SSC, na.rm = T)
  med.1 <- median(df$m2_MQ, na.rm = T)
  med.2 <- median(df$f_MQMR, na.rm = T)
  med.3 <- median(df$vs_SPV, na.rm = T)
  med.4 <- median(df$vd_SSF, na.rm = T)
  med.5 <- median(df$vd_MSI, na.rm = T)
  
  # obtaining detected VC
  x <- data.frame(t(data.frame(strsplit(df$REF_MFVdVs, "/"), row.names =c("Mutect2","Freebayes","Vardict","Varscan"))))
  rownames(x) <- c(1:nrow(x))
  x[!(is.na(x)|x=="NA")] = TRUE
  x[(is.na(x)|x=="NA")] = FALSE
  for(i in 1:4){x[,i] <- as.logical(x[,i])}
  df <- cbind(df,x)
  
  # NA filling
  df <- df %>%
    mutate(m2_MQ = ifelse(is.na(m2_MQ), ifelse(!is.na(f_MQMR), f_MQMR, med.1), m2_MQ),
           f_MQMR = ifelse(is.na(f_MQMR), ifelse(!is.na(m2_MQ), m2_MQ, med.2), f_MQMR)) %>%
    mutate(vs_SSC = ifelse(is.na(vs_SSC), med, vs_SSC),
           vs_SPV = ifelse(is.na(vs_SPV), ifelse(!is.na(vd_SSF), vd_SSF, med.3), vs_SPV),
           vd_SSF = ifelse(is.na(vd_SSF), ifelse(!is.na(vs_SPV), vs_SPV, med.4), vd_SSF),
           vd_MSI = ifelse(is.na(vd_MSI), med.5, vd_MSI)
           )

  return(df)
}
