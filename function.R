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

#' @description Reads in parsed files, filters based on criteria
#' @param count: No of detection by the 4 variant callers, default 2
#' @param mapq: Min. average mapping quality, default 40
#' @param p.value: p.value threshold, default 0.5
#' @param p.v: No. of features meeting p.value threshold, default 1 (out of 2)
#' @param msi: Min microsatellite instability, default 1
filtering <- function(file, count=2, mapq=40, p.value=0.05, p.v = 1, msi=1){
  df <- read.csv(file) 
  m = nrow(df)
  x <- df %>%
    mutate(rownumber = seq(1:nrow(.))) %>%
    mutate(m2_MQ = ifelse(is.na(m2_MQ), ifelse(!is.na(f_MQMR), f_MQMR, 0), m2_MQ),
           f_MQMR = ifelse(is.na(f_MQMR), ifelse(!is.na(m2_MQ), m2_MQ, 0), f_MQMR),
           vs_SPV = ifelse((vs_SPV <= p.value)&(!is.na(vs_SPV)), TRUE, FALSE),
           vd_SSF = ifelse((vd_SSF <= p.value)&(!is.na(vd_SSF)), TRUE, FALSE)) %>%
    mutate(avgMQ = (m2_MQ+f_MQMR)/2) %>%
    transmute(rownumber = rownumber,
              FILTER = FILTER_Mutect2 + FILTER_Freebayes + FILTER_Vardict + FILTER_Varscan,
              MQ = ifelse(avgMQ >= mapq, TRUE, FALSE),
              P.V = vs_SPV + vd_SSF,
              MSI = ifelse((vd_MSI > msi)&(!is.na(vd_MSI)), TRUE, FALSE)) %>%
    filter((FILTER >= count ) | (MQ | (P.V >= p.v) | MSI)) %>%
    select(rownumber)
  
  df <- df[x$rownumber,]
  print(paste("from:", m, "to:", nrow(df)))
  return(df)
}



