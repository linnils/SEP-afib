#Load libraries
library(survival)
library(haven)
library(medflex)
library(flextable)
library(officer)
library(dplyr)
library(formattable)
library(reshape2)
library(flexsurv)


rm(list = ls()) #Clear the workspace


#Read inn dataset
arb <- read_sas("filepath")

#Subset data to relevant variables
arb <- subset(arb, select = c("age4", "kjonn","edu4","bmi4", "kol4", "hypertgr",
                              "alc_freq","smoke4","hard4","af", "obsaf", 
                              "prev_mi", "dm4", "sleep4", "id"))

#For analysis of single mediator
#arb$myvar <- arb$smoke4


########Surv()##################################################################
#                          MEDIATION ANALYSIS                                  #
################################################################################


analyse <- function(arb, txt, bootstrapRounds = 1000){ ## BEGINNING OF analyse()
  
  
  #Defining maximum follow-up time
  maxFollowUpTimeTemp <- max(arb$obsaf)
  
  
  #Define main function - Natural effects model
  fit_nem <- function(arb,maxFollowUpTimeTemp)
  {   
    
    message("SURVREG:", appendLF = FALSE)
    
    #Fit mediation model
    arb$edu4temp <- arb$edu4
    
    #Analysis with multiple mediators
    fitImp <- survreg(Surv(obsaf,af) ~ factor(edu4temp) + bmi4 + 
    factor(hypertgr) + factor(smoke4) + factor(hard4) + factor(alc_freq) 
    + age4,data=arb)
    
    summary(fitImp) 
    fitImp
    
    #Do dataset expansion with counterfactual education level
    tempData1<-arb
    tempData1$edu40 <- tempData1$edu4
    tempData2<-arb
    tempData2$edu40 <- ifelse(tempData2$edu4 == 1, 2, tempData2$edu4)
    tempData2$edu40 <- ifelse(tempData2$edu4 == 2, 3, tempData2$edu40)
    tempData2$edu40 <- ifelse(tempData2$edu4 == 3, 4, tempData2$edu40)
    tempData2$edu40 <- ifelse(tempData2$edu4 == 4, 1, tempData2$edu40)
    tempData2$edu4temp <- tempData2$edu40
    tempData3<-arb
    tempData3$edu40 <- ifelse(tempData3$edu4 == 1, 3, tempData3$edu4)
    tempData3$edu40 <- ifelse(tempData3$edu4 == 2, 4, tempData3$edu40)
    tempData3$edu40 <- ifelse(tempData3$edu4 == 3, 1, tempData3$edu40)
    tempData3$edu40 <- ifelse(tempData3$edu4 == 4, 2, tempData3$edu40)
    tempData3$edu4temp <- tempData3$edu40
    tempData4<-arb
    tempData4$edu40 <- ifelse(tempData4$edu4 == 1, 4, tempData4$edu4)
    tempData4$edu40 <- ifelse(tempData4$edu4 == 2, 1, tempData4$edu40)
    tempData4$edu40 <- ifelse(tempData4$edu4 == 3, 2, tempData4$edu40)
    tempData4$edu40 <- ifelse(tempData4$edu4 == 4, 3, tempData4$edu40)
    tempData4$edu4temp <- tempData4$edu40
    
    #Prediction and imputation of counterfactual survival
    linPredTemp2 <- predict(fitImp,newdata=tempData2,type="linear")
    linPredTemp3 <- predict(fitImp,newdata=tempData3,type="linear")
    linPredTemp4 <- predict(fitImp,newdata=tempData4,type="linear")
    simTimesTemp2 <- rweibull(nrow(tempData2),
                              shape=1/fitImp$scale,exp(linPredTemp2))
    simTimesTemp3 <- rweibull(nrow(tempData3),
                              shape=1/fitImp$scale,exp(linPredTemp3))
    simTimesTemp4 <- rweibull(nrow(tempData4),
                              shape=1/fitImp$scale,exp(linPredTemp4))
    
    
    
    # Truncate times to maxFollowUpTimeTemp
    # Set Event if simulated time is less than maxFollowUpTimeTemp, 
    # otherwise non-event
    tempData2$af    <- ifelse(simTimesTemp2 < maxFollowUpTimeTemp, 
                              yes = 1            , no = 0)
    tempData2$obsaf <- ifelse(simTimesTemp2 < maxFollowUpTimeTemp, 
                              yes = simTimesTemp2, no = maxFollowUpTimeTemp)
    
    tempData3$af    <- ifelse(simTimesTemp3 < maxFollowUpTimeTemp, 
                              yes = 1            , no = 0)
    tempData3$obsaf <- ifelse(simTimesTemp3 < maxFollowUpTimeTemp, 
                              yes = simTimesTemp3, no = maxFollowUpTimeTemp)
    
    tempData4$af    <- ifelse(simTimesTemp4 < maxFollowUpTimeTemp, 
                              yes = 1            , no = 0)
    tempData4$obsaf <- ifelse(simTimesTemp4 < maxFollowUpTimeTemp, 
                              yes = simTimesTemp4, no = maxFollowUpTimeTemp)
    
    
    
    expData <- rbind(tempData1,tempData2,tempData3,tempData4)
    
    
    message("COXPH:", appendLF = FALSE)
    
    #Fit natural effects model
    fit_nem <- coxph(Surv(obsaf,af) ~ factor(edu4) + factor(edu40) + age4,
                     data=expData)
    
    
    return(summary(fit_nem)$coefficients)
  } #END OF FUNCTION fit_nem()
  
  #Get par estimates
  tempFitNEM <- fit_nem(arb,maxFollowUpTimeTemp=maxFollowUpTimeTemp)
  tempFitNEM
  
  #Running 50 simulations that simulates counterfactual follow-up time
  # Collecting all 50 results in one table, outtable
  # These are main results for IE and DE
  Nimp <- 50
  outTable <- array(NA,dim=c(dim(tempFitNEM),Nimp))
  message("Running ", Nimp, " simulations: ")
  for(j in 1:Nimp)
  {
    outTable[,,j] <- fit_nem(arb,maxFollowUpTimeTemp=maxFollowUpTimeTemp)
    message(".", appendLF = FALSE)
  }
  saveRDS(outTable, "outTable.RDS")
  outTable <- readRDS("outTable.RDS") 
  
  # Combining results from 50 to 1
  # Use only column 1 and 2 from results, i.e. beta and exp(beta)
  library(Amelia)
  temp <- mi.meld(q=outTable[,1,],se=outTable[,2,],byrow=F)
  tempOut <- tempFitNEM[,1:2]
  tempOut[,1] <- temp$q.mi
  tempOut[,2] <- temp$se.mi
  tempOut
  
  IE_0 <- tempOut[1:3,1]
  DE_0 <- tempOut[4:6,1]
  TE_0 <- IE_0+DE_0
  Q_0 <- IE_0/TE_0
  tempOut0 <- matrix(NA,nrow=12,ncol=1)
  tempOut0 <- rbind(IE_0,DE_0,TE_0,Q_0)
  tempOut0
  
  
  ### Get bootstrap SDs
  G <- bootstrapRounds
  outputObj <- array(NA,dim=c(dim(tempFitNEM),G))
  j <- 1
  while(j <= G)
  {
    message("Iteration ",j," of ",G," comp. ")
    tempData <- arb[sample(1:nrow(arb),replace=TRUE),]
    temp <- try(fit_nem(tempData,maxFollowUpTimeTemp),silent=FALSE)
    if( any(is.na(temp)) ){
      message("Apparent non-convergence - drawing new sample...")
    } else if( "try-error" %in% class(temp) ){
      message("Apparent non-convergence - drawing new sample...")
    } else {
      outputObj[,,j] <- temp
      j <- j + 1
    }
  }
  saveRDS(outputObj, "bootstrap-outputObj.rds")
  outputObj <- readRDS("bootstrap-outputObj.rds")
  
  outTable0 <- tempFitNEM[,1:2]
  outTable0[,1] <- apply(outputObj[,1,],1,mean,na.rm=T)
  outTable0[,2] <- apply(outputObj[,1,],1,sd,na.rm=T)
  colnames(outTable0) <- c("Beta","SD")
  outTable0
  print(outTable0)
  
  IE <- outputObj[1:3,1,]
  DE <- outputObj[4:6,1,]
  TE <- IE+DE
  Q <- IE/TE
  
  #Creating output table with estimates for IE, DE and TE, and SD and 95%CI
  outTable <- matrix(NA,nrow=12,ncol=4)
  outTable[1,]  <- c(mean(IE[1,]),sd(IE[1,]),quantile(IE[1,],c(.025,.975)))
  outTable[2,]  <- c(mean(IE[2,]),sd(IE[2,]),quantile(IE[2,],c(.025,.975)))
  outTable[3,]  <- c(mean(IE[3,]),sd(IE[3,]),quantile(IE[3,],c(.025,.975)))
  
  outTable[4,]  <- c(mean(DE[1,]),sd(DE[1,]),quantile(DE[1,],c(.025,.975)))
  outTable[5,]  <- c(mean(DE[2,]),sd(DE[2,]),quantile(DE[2,],c(.025,.975)))
  outTable[6,]  <- c(mean(DE[3,]),sd(DE[3,]),quantile(DE[3,],c(.025,.975)))
  
  outTable[7,]  <- c(mean(TE[1,]),sd(TE[1,]),quantile(TE[1,],c(.025,.975)))
  outTable[8,]  <- c(mean(TE[2,]),sd(TE[2,]),quantile(TE[2,],c(.025,.975)))
  outTable[9,]  <- c(mean(TE[3,]),sd(TE[3,]),quantile(TE[3,],c(.025,.975)))
  
  outTable[10,] <- c(mean( Q[1,]),sd( Q[1,]),quantile(Q[1,],c(.025,.975)))
  outTable[11,] <- c(mean( Q[2,]),sd( Q[2,]),quantile(Q[2,],c(.025,.975)))
  outTable[12,] <- c(mean( Q[3,]),sd( Q[3,]),quantile(Q[3,],c(.025,.975)))
  
  #Summarizing main results
  Table_final <- matrix(NA,nrow=12,ncol=5)
  Table_final[,2:5] <-outTable
  Table_final[1:3,1] <- IE_0
  Table_final[4:6,1] <- DE_0
  Table_final[7:9,1] <- TE_0
  Table_final[10:12,1] <- Q_0
  Table_final
  
  #Formatting table. Applying labels
  rownames(Table_final) <- c("Indirect effect","Indirect effect",
                             "Indirect effect","Direct effect","Direct effect",
                             "Direct effect","Total effect","Total effect",
                             "Total effect","Q2","Q3","Q4")
  colnames(Table_final) <- c("Beta","Beta_boots","StdDev","Low","High")
  print(Table_final)
  #saveRDS("Table_final-%s.rds")
  
  #Exponentiate to get HRs
  outHRs <- exp(Table_final[1:9, ])
  outHRs[,3] <- Table_final[1:9, "StdDev"]
  colnames(outHRs) <- c("HR","HR_boots","SD","Low","High")
  
  outHRs<-formattable(outHRs,format="f",digits=2)
  outHRs
  
  # Define the value labels for edu4
  edu4_labels <- setNames(c("Primary education", 
                            "Upper secondary education", 
                            "University/college <4 years", 
                            "University/college >=4 years"),
                          c(1, 2, 3, 4))
  
  
  # Convert the matrix to a data frame
  outHRs_df <- as.data.frame(outHRs)
  
  # Create the HR_CI column with formatted strings
  outHRs_df$HR_CI <- sprintf("%.2f (%.2f, %.2f)", 
                             outHRs_df$HR, outHRs_df$Low, outHRs_df$High)
  
  # Map the row names to edu4 values and effect types
  edu4_values <- c(2, 3, 4, 2, 3, 4, 2, 3, 4)
  effect_types <- rep(c("Indirect effect", "Direct effect", "Total effect"), 
                      each = 3)
  
  # Create the Education and Effect columns
  outHRs_df$Education <- edu4_labels[edu4_values]
  outHRs_df$Effect <- effect_types
  
  # Reshape the data frame to have separate columns for each effect type
  outHRs_wide <- dcast(outHRs_df, Education ~ Effect, value.var = "HR_CI")
  outHRs_wide <- outHRs_wide[order(match(outHRs_wide$Education, edu4_labels)), ]
  outHRs_wide <- outHRs_wide[, c("Education", "Total effect", "Direct effect", 
                                 "Indirect effect")]
  
  # Create a Word document using officer
  library(officer)
  doc <- read_docx()
  
  # Add a table to the Word document
  doc <- doc %>% 
    body_add_table(value = outHRs_wide, style = "table_template")
  
  doc <- doc |> body_add_table(value = outHRs_df, style="table_template")
  
  # Save the document
  setwd("filepath")
  print(doc, target = sprintf("mediationanalysis-%s.docx", txt))
  
} ### END OF analyse()
