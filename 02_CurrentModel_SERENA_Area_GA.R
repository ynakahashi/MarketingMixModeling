### library ##############################
library(GA)

### reg for GA ##############################
regGA <- function(dat, modelno, predVars, outdata=F){
  modelno <- paste(model_car, model_region, model_num, sep=" ")
  
  sink(file=paste(progdir, "Lsts/model_", modelno, ".lst",sep=""))
  predvar_list <- ""; keepvars <- predVars[[1]]
  
  if(length(predVars)==1){
    predvar_list <- paste(predvar_list, predVars[[length(predVars)]], sep="")
  } else {
    for(k in 1:(length(predVars)-1)) {
      predvar_list <- paste(predvar_list, predVars[[k]], " + ", sep="")
      keepvars <- c(keepvars, predVars[[k+1]])
    }
    
    predvar_list <- paste(predvar_list, predVars[[length(predVars)]], sep="")
  }
  
  model <- RunCodeRes(paste("glm(",target, "~", predvar_list, ", family=gaussian, data=dat)", sep=''))
  
  print(paste("Number of records: ", dim(dat)[1], sep='')); cat("\n")
  #corr <- correlation(dat, predVars)
  
  print(summary(model))
  print("Model processing complete"); cat("\n")
  sink(file=NULL)
  return(model)
}

### mod and ROI function ##############################
getModRoi <- function(x){
  
  indexSelected <- which(x == 1, arr.ind = T)
  predVarsSelected <- predVarsAll[indexSelected]
  predVars <- c(fixedVars, predVarsSelected)
  mod <- regGA(datLrn, model_num, predVars)
  
  if (length(which(predVarsSelected == "Profit_EX_logit")) > 0) {
    # in case "Profit_EX_logit" exists in predVarsSelected
    datOut <- calcAttri02(mod, datLrn, length(predVarsSelected) - 1)
  }
  else {
    # in case "Profit_EX_logit" is not selected
    datOut <- calcAttri02(mod, datLrn, length(predVarsSelected))
  }
  
  datAtt <- meltAttri(datOut)
  datROI <- data.frame(calcROI_02(datAtt, customer, week_start, roi_end))
  
  return(list(mod, datROI))
}


### AIC function ##############################
getFitness <- function(x){
  
  ModRoi <- getModRoi(x)
  
  mod <- ModRoi[1][[1]]
  datROI <- ModRoi[2][[1]]
  
  if ( (max(datROI$ROI, na.rm = T) > 5.0)  ||
       (min(datROI$ROI, na.rm = T) <= 0.0) ||
       (length(datROI$ROI[!is.na(datROI$ROI)]) == 0) ) {
    
    return(0.0)
  }
  else {
    
    return(1.0/mod$aic)
  }
}

### model number ##############################
model_num <- "201706_Tokyo-KanagawaEX_AG_02"

### fixed variables ##############################
fixedVars <- list(
  ## Seasonal Effect
  # "YW_201652", 
  # "YW_201701",
  "WEEK_17", "WEEK_32", "WEEK_51",
  # "YW_201345", 
  # "FRSTWKS_2014","NEWLYLAUNCH",
  # "YEAR_2012", "YEAR_2013", "YEAR_2014", "YEAR_2015","YW_201252", 
  # "YW_201253", "YW_201404_05", "YW_201610_23",
  # "Aft_YW_201624",
  # "YW_201405_06", "YW_201617",
  "Jan", "Feb", "Mar", "Apr", "May", "Jun", 
  "Jul", "Aug", "Sep", "Oct",
  # "Nov",
  "Dec",
  "Week_idx_02", "Week_idx_03", "Week_idx_04", "Week_idx_05", 
  "FiscalYearEnd", "GW", "OBON", "SW", 
  "nenmatsunenshi",
  "NEWYEAR.DAYS", "HolidayNum",
  "RushDemand_2014", "PreConsume_2014",
  # "RushDemand_2014_02",
  # "PreConsume_2014_02",
  "TaxRateUp",
  ## Socio
  # "GasPrice",
  # "AVERAGE.NIKKEI",
  # "UnemploymentRate",
  # "BigMinorChangeFlg.CPREV3",
  # "SpecialEditionFlg.CPREV3",
  ## Carlife event
  # "CarLifeInv",   
  # "CarLifeInv_02",
  "CarLife_ATANH",
  # "CompFMCFLG",
  # "FMC27",
  "Push",
  "ConservBuying"
)

### variables to be optimized ##############################
predVarsAll <- list(
  "Profit_EX_logit",
  ## FMI
  # "TVGRP_CarModel.ADSTOCK",
  "TVGRP_CarModel.LOG",
  "TVGRP_CarModel.LAG.1.LOG",
  "TVGRP_CarModel.LAG.2.LOG",
  "TVGRP_CarModel.LAG.3.LOG",
  # "TVGRP_Branding.ADSTOCK",
  "TVGRP_Branding.LOG",
  "TVGRP_Branding.LAG.1.LOG",
  "TVGRP_Branding.LAG.2.LOG",
  "TVGRP_Branding.LAG.3.LOG",
  # "RATIO.TO.COMP.ADSTOCK",
  # "Ratio.to.comp",
  # "DIGITAL_CARS_TOTAL.ADSTOCK",
  "DIGITAL_CARS_TOTAL.LOG",
  #"DIGITAL_CARS_TOTAL",
  #"DIGITAL_CARS_TOTAL.LAG.1",
  #"DIGITAL_CARS_TOTAL.LAG.2",
  #"DIGITAL_CARS_TOTAL.LAG.3",
  #"DIGITAL_CARS_TOTAL.LAG.4",
  "DIGITAL_CARS_TOTAL.LAG.1.LOG",
  "DIGITAL_CARS_TOTAL.LAG.2.LOG",
  "DIGITAL_CARS_TOTAL.LAG.3.LOG",
  "DIGITAL_CARS_TOTAL.LAG.4.LOG",
  #"Digi_Car_TV_Car_ON.LOG",
  #"Digi_Car_TV_Car_OFF.LOG",
  #"Digi_Brand_TV_Car_ON.LOG",
  #"Digi_Brand_TV_Car_OFF.LOG",
  #"DIGITAL_BRANDING_TOTAL.ADSTOCK",
  #"DIGITAL_BRANDING_TOTAL",
  "DIGITAL_BRANDING_TOTAL.LOG",
  #"DIGITAL_BRANDING_TOTAL.LAG.1",
  #"DIGITAL_BRANDING_TOTAL.LAG.2",
  #"DIGITAL_BRANDING_TOTAL.LAG.3",
  #"DIGITAL_BRANDING_TOTAL.LAG.4",
  #"DIGITAL_BRANDING_TOTAL.LAG.6",
  "DIGITAL_BRANDING_TOTAL.LAG.1.LOG",
  "DIGITAL_BRANDING_TOTAL.LAG.2.LOG",
  "DIGITAL_BRANDING_TOTAL.LAG.3.LOG",
  "DIGITAL_BRANDING_TOTAL.LAG.4.LOG",
  #"FlyerNum.LAG.D.3.LOG",
  #"FlyerNum.ADSTOCK",
  "FlyerNum.LOG",
  "FlyerNum.LAG.1.LOG",
  "FlyerNum.LAG.2.LOG",
  "FlyerNum.LAG.3.LOG",
  "FlyerNum.LAG.4.LOG",
  # "FlyerNum.LAG.8.LOG",
  #"NewsPaper.ADSTOCK",
  # "NewsPaper.LOG",
  # "NewsPaper.LAG.1.LOG",
  # "NewsPaper.LAG.2.LOG",
  # "NewsPaper.LAG.3.LOG",
  # "NewsPaper.LAG.4.LOG",
  #"NewsPaper_Branding.ADSTOCK",
  #"NewsPaper_Branding.LOG",
  # "NewsPaper_Branding.LAG.1.LOG",
  # "NewsPaper_Branding.LAG.2.LOG",
  # "NewsPaper_Branding.LAG.3.LOG",
  # "NewsPaper_Branding.LAG.4.LOG",
  # #"Magazine.ADSTOCK",
  #"Magazine.LOG",
  #"Magazine.LAG.1.LOG",
  #"Magazine.LAG.2.LOG",
  #"Magazine.LAG.3.LOG",
  #"Magazine.LAG.4.LOG",
  #"Magazine.LAG.8.LOG",
  #"Magazine_Branding.ADSTOCK",
  #"Magazine_Branding.LOG",
  #"Magazine_Branding.LAG.1.LOG",
  #"Magazine_Branding.LAG.2.LOG",
  #"Magazine_Branding.LAG.3.LOG",
  #"Magazine_Branding.LAG.4.LOG",
  #"Radio.ADSTOCK",
  #"Radio.LOG",
  #"Radio.LAG.1.LOG",
  #"Radio.LAG.2.LOG",
  #"Radio.LAG.3.LOG",
  #"Radio.LAG.4.LOG",
  #"Radio_Branding.ADSTOCK",
  #"Radio_Branding.LOG",
  #"Radio_Branding.LAG.1.LOG",
  #"Radio_Branding.LAG.2.LOG",
  #"Radio_Branding.LAG.3.LOG",
  #"Radio_Branding.LAG.4.LOG",
  #"RadioAll.LOG",
  "MASS3_Carmodel.LOG",
  "MASS3_Carmodel.LAG.1.LOG",
  "MASS3_Carmodel.LAG.2.LOG",
  "MASS3_Carmodel.LAG.3.LOG",
  "MASS3_Branding.LOG",
  "MASS3_Branding.LAG.1.LOG",
  "MASS3_Branding.LAG.2.LOG",
  "MASS3_Branding.LAG.3.LOG",
  # "VME_AchieveUP_SERENA",
  # "VME_Lease_all",
  # "VME_CAR_SX_Ship_RegiUP",
  "VME_CAR_SERENA_Ship_Regi_UP1",
  "VME_CAR_SERENA_Ship_Regi_UP2",
  "VME_CAR_SERENA_Op_WKTKUP",
  "VME_CAR_SERENA_TaxChangeUP",
  "VME_CAR_SERENA_OP_5",
  "VME_CCL_Intrade_Op_TakeoverUP_NSX",
  "VME_CCL_Intrade_Op_MileageUP_NSXJE",
  "VME_CCL_Intrade_TakeoverUP_SERENA",
  "VME_CCL_tokusenUP_NS",
  # "VME_CCL_RetentionUP_NSXLCuCaJ",
  "VME_CCL_RetentionUP_NSXLCuCaJ_01",
  "VME_CCL_RetentionUP_NSXLCuCaJ_02",
  "VME_LowInt_SERENA_NAC",
  #"VME_LowInt_SERENA_BVC",
  "VME_LowInt_SERENA_BVC_01",
  "VME_LowInt_SERENA_BVC_02"
  #"VME_LowInt_SERENA_Pres"
)

### initial variables ##############################
predVarsInit <- list(
  "TVGRP_CarModel.LOG",
  "DIGITAL_BRANDING_TOTAL.LAG.4.LOG",
  "FlyerNum.LAG.2.LOG",
  "VME_CAR_SERENA_Ship_Regi_UP1",
  "VME_LowInt_SERENA_BVC_02"
)

### initial settings for genetic algorithm
populationSize <- 20
generationNum <- 20

### elite solution (best found so far)
eliteSolution <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
                  0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
                  0, 0, 0, 0, 0, 1)

### similar solutions to eliteSolution
otherSolution <- matrix(eliteSolution, populationSize - 1, length(predVarsAll), byrow = T)
for (i in 1:(populationSize - 1)) {
  for (j in 1:length(predVarsAll)) {
    if (runif(1) < 1.0 / length(predVarsAll)) {
      if (otherSolution[i, j] == 0) otherSolution[i, j] <- 1
      else otherSolution[i, j] <- 0
    }
  }
}

### random solutions
#otherSolution <- matrix(sample(c(0, 1),
#                              length(predVarsAll) * (populationSize - 1),
#                              replace = T),
#                        populationSize - 1,
#                        length(predVarsAll))

### initial population is generated with eliteSolution and otherSolution
initPop <- rbind(eliteSolution, otherSolution)

### genetic operations and parameters are set appropriately
gaControl("binary" = list(selection = "gabin_tourSelection",
                          crossover = "gabin_uCrossover"))

### genetic algorithm running
runGA <- ga(type = "binary",
            fitness = getFitness,
            suggestions = initPop,
            pmutation = ga_pmutation,
            popSize = populationSize,
            maxiter = generationNum,
            nBits = length(predVarsAll),
            monitor = T)

### summary shows up and the data is saved
summary(runGA)
opt_indexselected <- which(runGA@solution[1,] == 1, arr.ind = T)
opt_predVarsSelected <- predVarsAll[opt_indexselected]
opt_fitness <- runGA@fitnessValue
opt_aic <- 1.0 / opt_fitness
opt_ModRoi <- getModRoi(runGA@solution[1,])
opt_result <- list(opt_predVarsSelected, opt_aic, opt_ModRoi[1][[1]], opt_ModRoi[2][[1]])
save(opt_result, file="GA_SizePopXXX_NumGenXXX.Rdata")
