# this is the unfinish summary result script 
 
# need to run 01., 02, calibrate_timevarying, 02_2 and 03 before running this script 
rm(list = ls()) 
basePath <- getwd()

# Load useful libraries

library("readr")
library("dplyr")
library("tidyr")
library("purrr")
library("gridExtra")
# Setup directories after setting working directory to source file 
# directory 

Rcode <- file.path(basePath, "03. Code")

DataFolder <- file.path(basePath, "01. DATA/model input" )
OutputFolder <- file.path(basePath, "04. Output" )
ResultsFolder <- file.path(basePath, "04. Output/Results")
projectFolder <- file.path(basePath)

project_name <- "HCVModel"

source(file.path(Rcode, "/Functions/HCV_model.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R"))


load(file.path(ResultsFolder, paste0("results_", "2022-10-28_23-00-37"), 
                         paste0("results_", "2022-10-28_23-00-37", ".rda")))
#### Section 1. population calibration ####
options(warn=-1)
#### MSM total population ####
popS <- popResults_MidYear(HCV, bestResults, Population = NULL,
                           Disease_prog = NULL , 
                           Cascade = NULL, param = paramResults, endYear = 50)%>%
  as.data.frame()




popS_range <- popResults_range(HCV, popS, Population = NULL, 
                               Disease_prog = NULL,
                               Cascade = NULL, end_Y = 50)
# real world data
male_pop <-read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"), 
                              "Malepop.csv"), header = TRUE)%>%
  as.data.frame()%>%mutate(realPop = Total.male.population*0.02,
                           time = Year - HCV$cabY + 1,
                           low = Total.male.population*0.01,
                           up = Total.male.population*0.025)%>%
  dplyr::select(realPop, time, low, up)%>%as.data.frame()
# plot 
totalPopPlot <- indicatorPlot(popS_range, ylabel = "Population size", 
                              facetPlot = NULL, 
                              calibration_Y = HCV$cabY,
                              rangeun = "t", 
                              xlimits = c(HCV$startYear, 
                                          HCV$startYear + 21, 1),
                              groupPlot = NULL, 
                              observationData = male_pop, 
                              simulateYear = HCV$simY) + 
  ggtitle("Overall population") 
# adjust the y-axis 
totalPopPlot  <- totalPopPlot + expandy(popS_range, 49)

# population by group

# notification data

MSMHIVnum <-read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"), 
                               "HIV_num.csv"), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year-2003, realPop = number,
                           low = lower,
                           up = upper)
MSMHIVnum <- MSMHIVnum%>%filter(population !="HIV+")
MSMHIVnum <- factorpop(MSMHIVnum)
MSMPop <- popResults_MidYear(HCV, bestResults,Population = HCV$popNames,
                             Disease_prog = NULL , 
                             Cascade = NULL, param = paramResults, endYear = 50)%>%
  as.data.frame() 
MSMPop[MSMPop < 0] = NaN

MSMPop_range <- popResults_range(HCV, MSMPop, Population = HCV$popNames,
                                 Disease_prog = NULL , 
                                 Cascade = NULL, end_Y = 50)

MSMPopPlot <- indicatorPlot(MSMPop_range, ylabel = "Population size by MSM", 
                            facetPlot = population,
                            calibration_Y = 2004,
                            rangeun = "y", 
                            xlimits = c(HCV$startYear, 
                                        HCV$startYear+26, 1),
                            groupPlot = NULL, 
                            observationData = MSMHIVnum, 
                            simulateYear = HCV$simY) + 
  geom_vline(xintercept = 2022 - HCV$cabY + 1, size = 1.2, linetype = "dashed")

MSMPopPlot <- MSMPopPlot + 
  facet_custom (~population, scales = "free", ncol = 2, 
                scale_overrides = 
                  list(scale_new(1,scale_y_continuous(limits = c(60000,300000))),
                    scale_new(2,scale_y_continuous(limits = c(0,6000))),
                    scale_new(3,scale_y_continuous(limits = c(0, 10000))),
                    scale_new(4,scale_y_continuous(limits =c(0,60000))))) + 
  theme_Publication_facet()

#### MSM HIV prevalence #### 
# # MSM HIV prevalence 

# HIV prevalence: real world data 
HIVPrev <- read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"), 
                              "HIVPrev.csv"), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year-2003, realPop = HIV.prevalence*100,
                           up = upper*100,
                           low = lower*100)
# subtract number of HIV infected 
HIVInfected <- popResults_MidYear(HCV, bestResults,
                                  Population = c("HIV+", "HIV+d"),
                                  Disease_prog = NULL, 
                                  Cascade = NULL, param = paramResults, 
                                  endYear = 50) 

HIVInf <- HIVInfected%>%ungroup()

# sum the number of HIV infected 
HIVInf <-HIVInf%>% dplyr::group_by(year)%>%
  summarise_at(.vars = c(colnames(popS)[-1]),sum)
# get the total number overtime
HIVTotal <- popResults_MidYear(HCV, bestResults,
                               Population = NULL,
                               Disease_prog = NULL, 
                               Cascade = NULL, param = paramResults, 
                               endYear = 50)
# calculate HIV prevalence 
HIV <- cbind(year = seq(HCV$startYear , 50-1 ,1),
             as.data.frame(HIVInf[, -1] / HIVTotal[ ,-1])*100)%>%
  tibble::as_tibble() 
# get the quantile of HIV prevalence 

HIV_range <- popResults_range(HCV, HIV, Population = NULL,
                                       Disease_prog = NULL , 
                                       Cascade = NULL, end_Y = 50)  
# Plot: HIV prevalence
HIVPrevPlot <- indicatorPlot(HIV_range, 
                             ylabel = "HIV prevlaence (%)",
                             xlimits = c(HCV$startYear, 
                                         HCV$startYear+21, 1),
                             calibration_Y = 2004,
                             rangeun = "y", 
                             groupPlot = NULL, 
                             facetPlot = NULL,
                             observationData = HIVPrev, 
                             simulateYear = HCV$simY)
# adjust the y-axis
HIVPrevPlot <- HIVPrevPlot + scale_y_continuous(limits = c(0,25))

#### MSM HIV diagnosis % #### 

# HIV diagnosis numbers: real world data
HIVdiag <- read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"), 
                              "HIVdiag.csv"), header = TRUE)%>%
  as.data.frame()%>%
  mutate(time = year - 2003, realPop = HIV.diagnosis.rate*100,
         low = lower*100,
         up = upper*100)

# subtract the number of diagnosed in the model
HIVInfected <- popResults_MidYear(HCV, bestResults,Population = c("HIV+d"),
                                  Disease_prog = NULL, 
                                  Cascade = NULL, param = paramResults,
                                  endYear = 50)
# calculate number of HIV diagnosed  
HIVDiagRate <- cbind(year = seq(HCV$startYear , 50-1 ,1),
                     as.data.frame(HIVInfected[, -c(1,2)]/ HIVInf[ ,-1])*100)%>%
  tibble::as_tibble() 
  

# get the quantile of HIV diagnosed number  

HIVDiagRate_range <- popResults_range(HCV, HIVDiagRate, Population = NULL,
                              Disease_prog = NULL , 
                              Cascade = NULL, end_Y = 50) 


HIVDiagPlot <- indicatorPlot(HIVDiagRate_range , 
                             ylabel = "HIV Diagnosis among HIV positive MSM (%)",
                             xlimits = c(HCV$startYear, 
                                         HCV$startYear+21, 1),
                             calibration_Y = 2004,
                             rangeun = "y", 
                             groupPlot = NULL, 
                             facetPlot = NULL,
                             observationData = HIVdiag, 
                             simulateYear = HCV$simY) + 
  ggtitle("HIV undiagnosis among HIV positive MSM (%)") 

HIVDiagPlot <- HIVDiagPlot + scale_y_continuous(limits = c(0, 100)) 



##### population transition ##### 
# subtract the transition number between sub-population 
bestResults$newpop_tran

#### array to dataframe 
PopTransition <- as.data.frame.table(bestResults$newpop_tran)%>%
  mutate(timestep = c(rep(seq(HCV$startYear ,30-HCV$timestep,
                              HCV$timestep) , each = HCV$npops*HCV$npops)),
         from = Var1,
         To = Var2)%>%select(-c(Var1, Var2, Var3))
# giving the average number to first time step 
impute<- PopTransition%>%filter(timestep >1 &timestep <2)%>%
  group_by(from, To)%>%
  mutate(total_pop = sum(Freq)/length(Freq))%>%
  select(total_pop)%>%ungroup()

impute <- impute[c(1:as.numeric(HCV$npops*HCV$npops)),]

PopTransition[c(1:as.numeric(HCV$npops*HCV$npops)), "Freq"] <- impute$total_pop 

# apply the same data sorting to paramResults 
paramTransition <- lapply(paramResults, function(x) x$newpop_tran)
paramTransition_data <- lapply(paramTransition, function(x){ 
  
  test <- as.data.frame.table(x)%>%
    mutate(timestep = c(rep(seq(HCV$startYear ,30-HCV$timestep,
                                HCV$timestep) , each = HCV$npops*HCV$npops)),
           from = Var1,
           To = Var2)%>%select(-c(Var1, Var2, Var3))

  impute<- test%>%filter(timestep >1 &timestep <2)%>%
    group_by(from, To)%>%
    mutate(total_pop = sum(Freq)/length(Freq))%>%
    select(total_pop)%>%ungroup()
  
  impute <- impute[c(1:as.numeric(HCV$npops*HCV$npops)),]
  
  test[c(1:as.numeric(HCV$npops*HCV$npops)), "Freq"] <- impute$total_pop
  
  return(test)
  
  }) 

# convert the column(Freq) in each simulation to columns in a dataframe 
library(data.table)
paramTransition_temp <- lapply(paramTransition_data, function(x)x$Freq)%>%
  as.data.table(matrix(unlist(.), ncol=length(.), byrow=FALSE))
colnames(paramTransition_temp) <- paste0("set",seq(1, HCV$numberSamples,1)) 
# reorder columns to align with other dataframes 
# set1, set10, set2, set3.....
# using str_sort
paramTransition_temp <-paramTransition_temp%>%
  select(str_sort(colnames(.), numeric = FALSE))




# bind transition of best fit and param 
PopTransition_all <- cbind(timestep = PopTransition$timestep, 
                           from = PopTransition$from,  
                           to = PopTransition$To, 
                           best = PopTransition$Freq, 
                          paramTransition_temp)%>%
  mutate(year = rep(rep(seq(1, 30-1,1),each = 1/HCV$timestep),
                    each = HCV$npops*HCV$npops))




## sum up the number in a year
PPTranToHIV <- list()
PPTranToPrEP <- list()
PPTranToHIVP <- list()
PPTranToHIVD <- list()

# subtract colnames 

cname <- c(colnames(PopTransition_all)[4:(4+HCV$numberSamples)]) 

# sum up the number between years and split by from status 

for (i in unique(PopTransition$from)) { 
  PPTranToHIV[[i]] <- PopTransition_all%>%filter(from == i, to == "HIV-")%>%
        group_by(year)%>%summarise_at(.vars = cname, sum)
}

for (i in unique(PopTransition$from)) { 
  PPTranToPrEP[[i]] <- PopTransition_all%>%filter(from == i, to == "HIV-PrEP")%>%
    group_by(year)%>%summarise_at(.vars = cname, sum)
}

for (i in unique(PopTransition$from)) { 
  PPTranToHIVP[[i]] <- PopTransition_all%>%filter(from == i, to == "HIV+")%>%
    group_by(year)%>%summarise_at(.vars = cname, sum)
}
for (i in unique(PopTransition$from)) { 
  PPTranToHIVD[[i]] <- PopTransition_all%>%filter(from == i, to == "HIV+d")%>%
    group_by(year)%>%summarise_at(.vars = cname, sum)
}





# PrEP drop 

# number of HIV PrEP transit to HIV- 
PrEPtoHIVN <- PPTranToHIV[["HIV-PrEP"]]

# number of people in PrEP  

MSMPop_PrEP <- popResults_MidYear(HCV, bestResults,
                                  Population = c("HIV-PrEP"),
                                  Disease_prog = NULL , 
                                  Cascade = NULL, param = paramResults,
                                  endYear = 50)%>%as.data.frame() 



PrEPDropInc <- cbind(year = seq(HCV$startYear, 50 -1 , 1),
                     as.data.frame(PrEPtoHIVN[, -1]/
                                     MSMPop_PrEP[ ,-c(1,2)]*100))%>%
  tibble::as_tibble() 



# summarise q5/q95 and median for PrEP drop out

PrEPDropInc_range <- popResults_range(HCV, PrEPDropInc, Population = NULL,
                                   Disease_prog = NULL , 
                                   Cascade = NULL, end_Y = 50)

# real world data
HIVPrEPdrop <- read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"), 
                                  "PrEP_dropout.csv"), header = TRUE)%>%
  as.data.frame()%>%
  mutate(time = year - 2003, realPop = drop.out*100,
         low = lower*100,
         up = upper*100)





PrEPdropPlot <- indicatorPlot(PrEPDropInc_range,
                              ylabel = "drop out rate(%)",
                              xlimits = c(HCV$startYear,
                                          HCV$startYear+21, 1),
                              calibration_Y = HCV$cabY,
                              rangeu = "y",
                              groupPlot = NULL,
                              observationData = HIVPrEPdrop, 
                              simulateYear = HCV$simY) +
  ggtitle("PrEP drop out rate") 

PrEPdropPlot <- PrEPdropPlot + expandy(PrEPDropInc_range, 49)


##### number of HIV acquisition among PrEP users #####
PrEPtoHIVP <- PPTranToHIVP[["HIV-PrEP"]]

# summarise q5/q95 and median 
PrEPtoHIVP_range <- popResults_range(HCV, PrEPtoHIVP, Population = NULL,
                 Disease_prog = NULL , 
                 Cascade = NULL, end_Y = 50)

PrEPHIVPlot <- indicatorPlot(PrEPtoHIVP_range ,
                              ylabel = "Number of HIV infected",
                              xlimits = c(HCV$startYear,
                                          HCV$startYear+21, 1),
                              calibration_Y = HCV$cabY,
                              rangeu = "t",
                              groupPlot = NULL,
                              observationData = NULL, 
                              simulateYear = HCV$simY) +
  ggtitle("number of HIV acquisition among PrEP users") +
  scale_y_continuous(limits =c(0, 10))


##### HIV incidence among PrEP users #####

PrEPHIVInc <- cbind(year = seq(HCV$startYear, 50 -1 , 1),
                     as.data.frame(PrEPtoHIVP[, -1]/
                                     MSMPop_PrEP[ ,-c(1,2)]*100))%>%
  tibble::as_tibble()  



PrEPHIVInc_range <- popResults_range(HCV, PrEPHIVInc, Population = NULL,
                                      Disease_prog = NULL , 
                                      Cascade = NULL, end_Y = 50)

PrEPHIVIncPlot <- indicatorPlot(PrEPHIVInc_range ,
                             ylabel = "HIV incidence among PrEP users (%)",
                             xlimits = c(HCV$startYear,
                                         HCV$startYear+21, 1),
                             calibration_Y = HCV$cabY,
                             rangeu = "y",
                             groupPlot = NULL,
                             observationData = NULL, 
                             simulateYear = HCV$simY) +
  ggtitle("HIV incidence among PrEP users") 


##### HIV incidence #####
HIVNtoHIVP <- PPTranToHIVP[["HIV-"]] 

# number of HIV negative 
MSM_HIVN <- popResults_MidYear(HCV, bestResults,
                                  Population = c("HIV-PrEP", "HIV-"),
                                  Disease_prog = NULL , 
                                  Cascade = NULL, param = paramResults,
                                  endYear = 50)%>%as.data.frame()%>%
  group_by(year)%>%summarise_at(.vars = c(colnames(popS)[-1]), sum)

HIVInc <- cbind(year = seq(HCV$startYear, 50 -1 , 1),
               as.data.frame((HIVNtoHIVP[, -1] + PrEPtoHIVP[ ,-1])/
                               MSM_HIVN[ ,-1]*1000)) 

HIVInc_range <- popResults_range(HCV, HIVInc, Population = NULL,
                                     Disease_prog = NULL , 
                                     Cascade = NULL, end_Y = 50)

# real world data 

MSMHIVInc <-read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"), 
                                "HIVInc.csv"), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year-2003, 
                           realPop = HIV.incidence*1000,
                           up = upper*1000,
                           low = lower*1000)



HIVIncPlot <- indicatorPlot(HIVInc_range ,
                                ylabel = "HIV incidence (%)",
                                xlimits = c(HCV$startYear,
                                            HCV$startYear+21, 1),
                                calibration_Y = HCV$cabY,
                                rangeu = "y",
                                groupPlot = NULL,
                                observationData = MSMHIVInc, 
                                simulateYear = HCV$simY) +
  ggtitle("Overall HIV incidence") 


HIVIncPlot <- HIVIncPlot + scale_y_continuous(limits = c(0,100))


##### number of HIV diagnosis ##### 

HIVDiagnum <- PPTranToHIVD[["HIV+"]]

# summarise q5/q95 and median 
HIVPtoHIVD_range <- popResults_range(HCV, HIVDiagnum, Population = NULL,
                                     Disease_prog = NULL , 
                                     Cascade = NULL, end_Y = 50) 

# real world data 
MSMHIVDiagnum <-read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"), 
                               "HIV_num_num.csv"), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year-2003, 
                           realPop = diagnumber,
                           up = upper,
                           low = lower)

# HIV diagnosed number  
HIVDiagnumPlot <- indicatorPlot(HIVPtoHIVD_range ,
                             ylabel = "Number of HIV diagnosed",
                             xlimits = c(HCV$startYear,
                                         HCV$startYear+21, 1),
                             calibration_Y = HCV$cabY,
                             rangeu = "y",
                             groupPlot = NULL,
                             observationData = MSMHIVDiagnum, 
                             simulateYear = HCV$simY) +
  ggtitle("Number of HIV diagnosed") + scale_y_continuous(limits =c(0, 3000))




####Section 2. HCV-related calibration #### 

##### HCV prevalence #####
# number of susceptible 
tempNOTInfected_all <- popResults_MidYear(HCV, bestResults,Population = NULL,
                                          Disease_prog = NULL, 
                                          Cascade = c("s", "cured"), 
                                          param = paramResults ,endYear = 50)%>%
  as.data.frame()%>%group_by(year)%>%
  summarise_at(.vars = c(colnames(popS)[-1]),sum) 


# number of all population 
tempTotal_all <- popResults_MidYear(HCV, bestResults,Population = NULL,
                                    Disease_prog = NULL, 
                                    Cascade = NULL, param = paramResults, 
                                    endYear = 50)



tempPrev_all <- cbind(year = seq(HCV$startYear , 50-1 ,1),
                      as.data.frame((tempTotal_all[, -1] - 
                                       tempNOTInfected_all[ ,-1])/ 
                                      tempTotal_all[ ,-1]*100))%>%
  tibble::as_tibble() 

# summarise q5/q95 and median for prevalence 
tempPrev_all_range <- popResults_range(HCV, tempPrev_all, Population = NULL,
                                       Disease_prog = NULL , 
                                       Cascade = NULL, end_Y = 50)

# HCV prevalence plot: general MSM
totalPrevPlot <- indicatorPlot(tempPrev_all_range, 
                               ylabel = "Overall HCV prevalence (%)",
                               xlimits = c(HCV$startYear, 
                                           HCV$startYear+21, 1),
                               calibration_Y = 2004,
                               rangeun = "y", 
                               groupPlot = NULL, 
                               facetPlot = NULL, 
                               simulateYear = HCV$simY) +
  ggtitle("Overall HCV prevalence")

totalPrevPlot <- totalPrevPlot + expandy(tempPrev_all_range, HCV$startYear +21)

# HCV prevalence in sub-MSM groups 

tempNOTInfected <- popResults_MidYear(HCV, bestResults,
                                      Population = HCV$popNames,
                                      Disease_prog = NULL, 
                                      Cascade = c("s", "cured"), 
                                      param = paramResults,
                                      endYear = 50)%>%
  group_by(year, population)%>%summarise_at(.vars = c(colnames(popS)[-1]), sum)



tempTotal <- popResults_MidYear(HCV, bestResults,
                                Population = HCV$popNames,
                                Disease_prog = NULL, 
                                Cascade = NULL, 
                                param = paramResults, 
                                endYear = 50) 



tempPrev <- cbind(year = rep(seq(HCV$startYear , 50-1 ,1), each = HCV$npops),
                  population = HCV$popNames,
                  
                  as.data.frame(100*(tempTotal[, -c(1,2)] - 
                                   tempNOTInfected[ ,-c(1,2)])/ 
                                  tempTotal[ ,-c(1,2)]))%>%
  tibble::as_tibble() 


# summarise q5/q95 and median for HCV prevalence in MSM subgroups

tempPrev_range <- popResults_range(HCV, tempPrev, Population = HCV$popNames,
                                   Disease_prog = NULL , 
                                   Cascade = NULL, end_Y = 50)


# population convert to factors 
# be careful about the NaN value 

tempPrev_rang <- tempPrev_range%>%FactorPop(., HCV$popNames)

c(unique(tempPrev_range$population))
# real world data 

MSMHCVPrev <-read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"), 
                                "HCVPrev.csv"), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year-2003, 
                           realPop = HCV.prevalence*100,
                           up = upper*100,
                           low = lower*100)

popPrevPlot <- indicatorPlot(tempPrev_range, 
                             ylabel = "Population HCV prevalence (%)",
                             xlimits = c(HCV$startYear, 
                                         HCV$startYear+21, 1),
                             calibration_Y = 2004,
                             rangeun = "y", 
                             groupPlot = NULL, 
                             facetPlot = population,
                             observationData = MSMHCVPrev, 
                             simulateYear = HCV$simY) +
  ggtitle("HCV prevalence by population") 

popPrevPlot <- popPrevPlot + 
  facet_custom (~population,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(tempPrev_range, HCV$startYear +21)[1],
                                                     20))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(tempPrev_range, HCV$startYear +21)[2],
                                                     limit_facet(tempPrev_range, HCV$startYear +21)[6]))),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(tempPrev_range, HCV$startYear +21)[3],
                                                     limit_facet(tempPrev_range, HCV$startYear +21)[7]))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(tempPrev_range, HCV$startYear +21)[4],
                                                     limit_facet(tempPrev_range, HCV$startYear +21)[8]
                              )))))

HCV$timestep
##### HCV incidence ##### 

# total population #

HCVInfect_all <- indicatorResults(HCV, bestResults, "newInfections", pop="all",
                           paramR = paramResults, range = "y",
                           endY = 50)

# reorder columns and only keep best and parameterset 
colset <- HCVInfect_all%>%select(., contains("Set"))%>%
  select(., str_sort(paste0("set", seq(1, HCV$numberSamples,1)),numeric = FALSE))

HCVInfect <- HCVInfect_all%>%select(year, best, colnames(colset))



HCVInc <- cbind(year = seq(HCV$startYear , 50-1 ,1),
                as.data.frame(HCVInfect[, -1] / 
                                 popS[ ,-1]*1000))%>%
  tibble::as_tibble() 


# get aggregate columns 
HCVInc_range <- popResults_range(HCV, HCVInc, Population = NULL,
                 Disease_prog = NULL , 
                 Cascade = NULL, end_Y = 50)

HCVIncPlot <- indicatorPlot(HCVInc_range, 
                            ylabel = "HCV incidence per 1000",
                            xlimits = c(HCV$startYear, 
                                        HCV$startYear+21, 1),
                            calibration_Y = 2004,
                            rangeun = "y", 
                            groupPlot = NULL, 
                            facetPlot = NULL,
                            observationData = NULL, 
                            simulateYear = HCV$simY) +
  ggtitle("HCV Incidence") +  expandy(HCVInc_range, HCV$startYear +21)



# HCV incidence by groups # 

HCVInfectPop <- indicatorResults(HCV, bestResults, "newInfections", 
                                pop=HCV$popNames,
                                paramR = paramResults, range = "y",
                                endY = 50)%>%select(year, population, best,
                                                      colnames(colset) )

HCVPopInc <- cbind(year = rep(seq(HCV$startYear , 50-1 ,1), each = HCV$npops),
                   population = HCVInfectPop$population,
                as.data.frame(HCVInfectPop[, -c(1,2)] / 
                                MSMPop[ ,-c(1,2)]*1000))%>%
  tibble::as_tibble() 


# get aggregate columns 

HCVPopInc_range <- popResults_range(HCV, HCVPopInc, Population = HCV$popNames,
                 Disease_prog = NULL , 
                 Cascade = NULL, end_Y = 50)%>%as.data.frame()



# real world data 

MSMHCVInc <-read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"), 
                                "HCVInc.csv"), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year-2003, 
                           realPop = HCV.incidence*1000,
                           up = upper*1000,
                           low = lower*1000)


HCVPopIncPlot <- indicatorPlot(HCVPopInc_range, 
                            ylabel = "HCV incidence per 1000 by MSM groups",
                            xlimits = c(HCV$startYear, 
                                        HCV$startYear+21, 1),
                            calibration_Y = 2004,
                            rangeun = "y", 
                            groupPlot = NULL, 
                            facetPlot = population,
                            observationData = MSMHCVInc, 
                            simulateYear = HCV$simY) +
  ggtitle("HCV Incidence by MSM groups") + 
  facet_custom (~population,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(HCVPopInc_range, HCV$startYear +21)[1],
                                                     20))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(HCVPopInc_range, HCV$startYear +21)[2],
                                                     limit_facet(HCVPopInc_range, HCV$startYear +21)[6]))),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(HCVPopInc_range, HCV$startYear +21)[3],
                                                     limit_facet(HCVPopInc_range, HCV$startYear +21)[7]))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(HCVPopInc_range, HCV$startYear +21)[4],
                                                     limit_facet(HCVPopInc_range, HCV$startYear +21)[8]
                                                   )))))


HCVPopInc_ran

##### HCV treatment initiation ##### 
# general #
TreatmentInit <- indicatorResults(HCV, bestResults, "newTreatment", 
                 pop="all",
                 paramR = paramResults, range = "y",
                 endY = 50)

TreatmentInitPlot <- indicatorPlot(TreatmentInit, 
                          ylabel = "Number of initiated HCV treatment each year", 
                          facetPlot = NULL,
                          calibration_Y = 2004,
                          rangeun = "y", 
                          xlimits = c(HCV$startYear, 
                                      HCV$startYear+21, 1),
                          groupPlot = NULL, 
                          observationData = NULL, 
                          simulateYear = HCV$simY) + 
  ggtitle("Overall number of HCV treatment initiated")




#### treatment initiation rate ####

colset <- TreatmentInit%>%select(., contains("Set"))%>%
  select(., str_sort(paste0("set", seq(1, HCV$numberSamples,1)),numeric = FALSE))

TreatmentInit <- TreatmentInit%>%select(year, best, colnames(colset))
temptreatment_all <- cbind(year = seq(HCV$startYear , 50-1 ,1),
                      as.data.frame(TreatmentInit[ , -1]/(tempTotal_all[, -1] - tempNOTInfected_all[ ,-1])*100))%>%
  tibble::as_tibble()


TreatmentInit_range <- popResults_range(HCV, temptreatment_all, Population = NULL,
                                 Disease_prog = NULL , 
                                 Cascade = NULL, end_Y = 50)

TreatmentInitPerPlot <- indicatorPlot(TreatmentInit_range, 
                                   ylabel = "HCV treatment initiated %", 
                                   facetPlot = NULL,
                                   calibration_Y = 2004,
                                   rangeun = "y", 
                                   xlimits = c(HCV$startYear, 
                                               HCV$startYear+21, 1),
                                   groupPlot = NULL, 
                                   observationData = NULL, 
                                   simulateYear = HCV$simY) + 
  ggtitle("HCV treatment initiated %") + expandy(TreatmentInit_range, HCV$startYear +21)

# sub-group  in indicatorResults function:q5/q95
TreatmentInitPop <- indicatorResults(HCV, bestResults, "newTreatment", 
                                  pop=NULL,
                                  paramR = paramResults, range = "y",
                                  endY = 50)%>%select(year, population, best,
                                                             colnames(colset))
TreatmentInitPop_range <- popResults_range(HCV, TreatmentInitPop , Population = HCV$popNames,
                                           Disease_prog = NULL , 
                                           Cascade = NULL, end_Y = 50)

TreatmentInitPopPlot <- indicatorPlot(TreatmentInitPop_range, 
                                   ylabel = "Number of initiated HCV treatment each year", 
                                   facetPlot = population,
                                   calibration_Y = 2004,
                                   rangeun = "y", 
                                   xlimits = c(HCV$startYear, 
                                               HCV$startYear+21, 1),
                                   groupPlot = NULL, 
                                   observationData = NULL, 
                                   simulateYear = HCV$simY) + 
  ggtitle("Number of HCV treatment initiated by groups") +
  facet_custom (~population,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(TreatmentInitPop_range, HCV$startYear +21)[1],
                                                     limit_facet(TreatmentInitPop_range, HCV$startYear +21)[5]))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(TreatmentInitPop_range, HCV$startYear +21)[2],
                                                     limit_facet(TreatmentInitPop_range, HCV$startYear +21)[6]))),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(TreatmentInitPop_range, HCV$startYear +21)[3],
                                                     limit_facet(TreatmentInitPop_range, HCV$startYear +21)[7]))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(TreatmentInitPop_range, HCV$startYear +21)[4],
                                                     limit_facet(TreatmentInitPop_range, HCV$startYear +21)[8]
                                                   )))))





# cumulative treatment initiated #

# overall population #
TreatmentInit_cumulative <-  TreatmentInit%>%
  select(year, best,colnames(colset))%>%
  mutate(across(c(best, colnames(colset)), ~cumsum(.x)))

TreatmentInit_cumulativeRange <- popResults_range(HCV, TreatmentInit_cumulative, 
                                                  Population = NULL,
                                                  Disease_prog = NULL, 
                                                  Cascade = NULL, end_Y = 50)

# real world data
HCVTreat <- read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"),
                               "HCVtreat.csv"),
                     header = TRUE)%>%as.data.frame()%>%
  mutate(realPop = number.of.treatment.initiated,time = year-2003,
         low = lower,
         up = upper)%>%
  dplyr::select(realPop, time, low, up)%>%as.data.frame()


totalCumTreatmentPlot <- indicatorPlot(TreatmentInit_cumulativeRange, 
                                       ylab("Cumulative number of HCV treatment initiated"), 
                                       facetPlot = NULL,
                                       calibration_Y = 2004,
                                       rangeun = "y", 
                                       xlimits = c(HCV$startYear, 
                                                   HCV$startYear+21,1),
                                       groupPlot = NULL, 
                                       observationData = NULL, 
                                       simulateYear = HCV$simY) + 
  ggtitle("Overall treatment initiated") + 
  expandy(TreatmentInit_cumulativeRange, HCV$startYear +21)



# MSM sub-groups #

TreatmentInitPop_cumulative <- TreatmentInitPop%>%group_by(population)%>%
  mutate(across(c(best, colnames(colset)), ~cumsum(.x))) 

TreatmentInitPop_cumulativeRange <- popResults_range(HCV, 
                                                     TreatmentInitPop_cumulative, 
                                                     Population = "all",
                                                     Disease_prog = NULL, 
                                                     Cascade = NULL, end_Y = 50)
PopCumTreatmentPlot <- indicatorPlot(TreatmentInitPop_cumulativeRange, 
                                     ylab("Cumulative number of HCV treatment initiated"), 
                                     facetPlot = population,
                                     calibration_Y = 2004,
                                     rangeun = "y", 
                                     xlimits = c(HCV$startYear, 
                                                 HCV$startYear+21, 1),
                                     groupPlot = NULL, 
                                     observationData = HCVTreat, 
                                     simulateYear = HCV$simY) + 
  ggtitle("Cumulative HCV treatment initiated by MSM group") + 
  facet_custom (~population,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(TreatmentInitPop_cumulativeRange, HCV$startYear +21)[1],
                                                     limit_facet(TreatmentInitPop_cumulativeRange, HCV$startYear +21)[5]))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(TreatmentInitPop_cumulativeRange, HCV$startYear +21)[2],
                                                     limit_facet(TreatmentInitPop_cumulativeRange, HCV$startYear +21)[6]))),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(TreatmentInitPop_cumulativeRange, HCV$startYear +21)[3],
                                                     limit_facet(TreatmentInitPop_cumulativeRange, HCV$startYear +21)[7]))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(TreatmentInitPop_cumulativeRange, HCV$startYear +21)[4],
                                                     limit_facet(TreatmentInitPop_cumulativeRange, HCV$startYear +21)[8]
                                                   )))))


##### HCV SVR ##### 

# general # 
options(warn=-1)

CureInit <- indicatorResults(HCV, bestResults, "newCured", 
                                  pop="all",
                                  paramR = paramResults, range = "y",
                                  endY = 50)

CureInitPlot <- indicatorPlot(CureInit, 
                              ylabel = "Number of Cured from HCV", 
                              facetPlot = NULL,
                              calibration_Y = 2004,
                              rangeun = "y", 
                              xlimits = c(HCV$startYear, HCV$startYear+21, 1),
                              groupPlot = NULL, 
                              observationData = NULL,
                              simulateYear = HCV$simY) + 
  ggtitle("Overall number of HCV cured")

CureInitPlot <- CureInitPlot + expandy(CureInit, HCV$startYear+21)


# sub-group 
CureInitPop <- indicatorResults(HCV, bestResults, "newCured", 
                                     pop=HCV$popNames,
                                     paramR = paramResults, range = "y",
                                     endY = 50)

CureInitPopPlot <- indicatorPlot(CureInitPop, 
                                 ylabel = "Number of people cure from HCV", 
                                 facetPlot = population,
                                      calibration_Y = 2004,
                                      rangeun = "y", 
                                      xlimits = c(HCV$startYear, 
                                                  HCV$startYear+21, 1),
                                      groupPlot = NULL, 
                                      observationData = NULL, 
                                      simulateYear = HCV$simY) + 
  ggtitle("Number of people cure from HCV by groups")


CureInitPopPlot  <- CureInitPopPlot  + 
  facet_custom (~population,scales = "free", ncol = 2,
                scale_overrides = list(
                  scale_new(1,
                            scale_y_continuous(
                              limits = c(0,limit_facet(CureInitPop, HCV$startYear+21)[1+HCV$npops]))),
                  scale_new(2,
                            scale_y_continuous(
                              limits = c(0,limit_facet(CureInitPop, HCV$startYear+21)[2+HCV$npops]))),
                  scale_new(3, scale_y_continuous(
                    limits = c(0,limit_facet(CureInitPop, HCV$startYear+21)[3+HCV$npops]))),
                  scale_new(4,scale_y_continuous(
                    limits = c(0,limit_facet(CureInitPop, HCV$startYear+21)[4+HCV$npops])))))

# cumulative number of HCV cure#

# overall population #
CureInit_cumulative <-  CureInit%>%
  select(year, best,colnames(colset))%>%
  mutate(across(c(best, colnames(colset)), ~cumsum(.x)))

CureInit_cumulativeRange <- popResults_range(HCV, CureInit_cumulative, 
                                                  Population = NULL,
                                                  Disease_prog = NULL, 
                                                  Cascade = NULL, end_Y = 50)

# real world data
HCVsvr <- read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"),
                             "HCVSVR.csv"),
                   header = TRUE)%>%as.data.frame()%>%
  mutate(realPop = number.of.SVR,time = year-2003,
         low = lower, 
         up = upper)%>%
  dplyr::select(realPop, time, low, up)%>%as.data.frame()


totalCumCurePlot <- indicatorPlot(CureInit_cumulativeRange, 
                                  ylab("Cumulative number cured from HCV"), 
                                  facetPlot = NULL,
                                  calibration_Y = 2004,
                                  rangeun = "y", 
                                  xlimits = c(HCV$startYear, HCV$startYear+21,1),
                                  groupPlot = NULL, 
                                  observationData = HCVsvr,
                                  simulateYear = HCV$simY) + 
  ggtitle("Overall cumulative number cured from HCV") + 
  scale_y_continuous(limits = c(0, 10000))

totalCumCurePlot <- totalCumCurePlot + expandy(CureInit_cumulativeRange, HCV$startYear+21)


# MSM sub-groups #

CureInitPop_cumulative <- CureInitPop%>%group_by(population)%>%
  mutate(across(c(best, colnames(colset)), ~cumsum(.x))) 

CureInitPop_cumulativeRange <- popResults_range(HCV, 
                                                     CureInitPop_cumulative, 
                                                     Population = "all",
                                                     Disease_prog = NULL, 
                                                     Cascade = NULL, end_Y = 50)
PopCumCurePlot <- indicatorPlot(CureInitPop_cumulativeRange, 
                                     ylab("Cumulative number of people cured from HCV by MSM groups"), 
                                     facetPlot = population,
                                     calibration_Y = 2004,
                                     rangeun = "y", 
                                     xlimits = c(HCV$startYear, 
                                                 HCV$startYear+21, 1),
                                     groupPlot = NULL, 
                                     observationData = NULL, 
                                     simulateYear = HCV$simY) + 
  ggtitle("Cumulative number of people cured from HCV by MSM groups") 


PopCumCurePlot <- PopCumCurePlot + 
  facet_custom (~population,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(CureInitPop_cumulativeRange, HCV$startYear +21)[1],
                                                     limit_facet(CureInitPop_cumulativeRange, HCV$startYear +21)[5]))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(CureInitPop_cumulativeRange, HCV$startYear +21)[2],
                                                     limit_facet(CureInitPop_cumulativeRange, HCV$startYear +21)[6]))),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(CureInitPop_cumulativeRange, HCV$startYear +21)[3],
                                                     limit_facet(CureInitPop_cumulativeRange, HCV$startYear +21)[7]))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(CureInitPop_cumulativeRange, HCV$startYear +21)[4],
                                                     limit_facet(CureInitPop_cumulativeRange, HCV$startYear +21)[8]
                                                   )))))

####Section 3.detail of components and flow #### 

##### Number of infections #####
###### total ######
cascade_NotS <- HCV$cascade_name[-6]
Infected_all <- popResults_MidYear(HCV, bestResults,Population = NULL,
                                   Disease_prog = NULL, 
                                   Cascade = cascade_NotS, 
                                   param = paramResults ,endYear = 50)%>%
  as.data.frame()%>%group_by(year)%>%
  summarise_at(.vars = c(colnames(popS)[-1]),sum) 


Infected_all_range <-popResults_range(HCV, Infected_all, Population = NULL,
                                      Disease_prog = NULL , 
                                      Cascade = NULL, end_Y = 50)

Infected_allPlot <- indicatorPlot(Infected_all_range, 
                                  ylabel = "Number of HCV infections",
                                  xlimits = c(HCV$startYear, 
                                              HCV$startYear+21, 1),
                                  calibration_Y = 2004,
                                  rangeun = "y", 
                                  groupPlot = NULL, 
                                  facetPlot = NULL, 
                                  simulateYear = HCV$simY) +
  ggtitle("Number of HCV infections") + 
  expandy(Infected_all_range, HCV$startYear + 21)

###### sub-pop ######
cascade_NotS <- HCV$cascade_name[-6]
InfectedPop_all <- popResults_MidYear(HCV, bestResults,Population = HCV$popNames,
                                      Disease_prog = NULL, 
                                      Cascade = cascade_NotS, 
                                      param = paramResults ,endYear = 50)%>%
  as.data.frame()%>%group_by(year, population)%>%
  summarise_at(.vars = c(colnames(popS)[-1]),sum) 



InfectedPop_all_range <-popResults_range(HCV, InfectedPop_all, 
                                         Population = HCV$popNames,
                                         Disease_prog = NULL, 
                                         Cascade = NULL, end_Y = 50)

InfectedPop_allPlot <- indicatorPlot(InfectedPop_all_range, 
                                     ylabel = "Number of HCV infections",
                                     xlimits = c(HCV$startYear, 
                                                 HCV$startYear+21, 1),
                                     calibration_Y = 2004,
                                     rangeun = "y", 
                                     groupPlot = NULL, 
                                     facetPlot = population, 
                                     simulateYear = HCV$simY) +
  ggtitle("Number of HCV infections by MSM population") 


InfectedPop_allPlot <- InfectedPop_allPlot + 
  facet_custom (~population,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(InfectedPop_all_range, HCV$startYear +21)[1],
                                                     limit_facet(InfectedPop_all_range, HCV$startYear +21)[5]))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(InfectedPop_all_range, HCV$startYear +21)[2],
                                                     limit_facet(InfectedPop_all_range, HCV$startYear +21)[6]))),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(InfectedPop_all_range, HCV$startYear +21)[3],
                                                     limit_facet(InfectedPop_all_range, HCV$startYear +21)[7]))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(InfectedPop_all_range, HCV$startYear +21)[4],
                                                     limit_facet(InfectedPop_all_range, HCV$startYear +21)[8]
                                                   )))))

##### Number of cured #####
###### total ######
Cured_all <- popResults_MidYear(HCV, bestResults,Population = NULL,
                                Disease_prog = NULL, 
                                Cascade = "cured", 
                                param = paramResults ,endYear = 50)%>%
  as.data.frame()%>%group_by(year)%>%
  summarise_at(.vars = c(colnames(popS)[-1]),sum) 


Cured_all_range <-popResults_range(HCV, Cured_all, Population = NULL,
                                   Disease_prog = NULL , 
                                   Cascade = NULL, end_Y = 50)

Cured_allPlot <- indicatorPlot(Cured_all_range, 
                               ylabel = "Number of HCV cured",
                               xlimits = c(HCV$startYear, 
                                           HCV$startYear+21, 1),
                               calibration_Y = 2004,
                               rangeun = "y", 
                               groupPlot = NULL, 
                               facetPlot = NULL, 
                               simulateYear = HCV$simY) +
  ggtitle("Number of HCV cured") + expandy(Cured_all_range, HCV$startYear + 21)

###### sub-population ######
CuredPop <- popResults_MidYear(HCV, bestResults,Population = HCV$popNames,
                                Disease_prog = NULL, 
                                Cascade = "cured", 
                                param = paramResults ,endYear = 50)%>%
  as.data.frame()%>%group_by(year, population)%>%
  summarise_at(.vars = c(colnames(popS)[-1]),sum) 


CuredPop_range <-popResults_range(HCV, CuredPop, Population = HCV$popNames,
                                   Disease_prog = NULL , 
                                   Cascade = NULL, end_Y = 50)

CuredPopPlot <- indicatorPlot(CuredPop_range , 
                               ylabel = "Number of HCV cured",
                               xlimits = c(HCV$startYear, 
                                           HCV$startYear+21, 1),
                               calibration_Y = 2004,
                               rangeun = "y", 
                               groupPlot = NULL, 
                               facetPlot = population, 
                               simulateYear = HCV$simY) +
  ggtitle("Number of HCV cured") + 
  facet_custom (~population,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(CuredPop_range, HCV$startYear +21)[1],
                                                     limit_facet(CuredPop_range, HCV$startYear +21)[5]))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(CuredPop_range, HCV$startYear +21)[2],
                                                     limit_facet(CuredPop_range, HCV$startYear +21)[6]))),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(CuredPop_range, HCV$startYear +21)[3],
                                                     limit_facet(CuredPop_range, HCV$startYear +21)[7]))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(limit_facet(CuredPop_range, HCV$startYear +21)[4],
                                                     limit_facet(CuredPop_range, HCV$startYear +21)[8]
                                                   )))))


##### Number of susceptible ##### 
###### total ######
Sus_all <- popResults_MidYear(HCV, bestResults,Population = NULL,
                              Disease_prog = "s", 
                              Cascade = NULL, 
                              param = paramResults ,endYear = 50)%>%
  as.data.frame()%>%group_by(year)%>%
  summarise_at(.vars = c(colnames(popS)[-1]),sum)%>%ungroup()


Sus_all_range <-popResults_range(HCV, Sus_all , Population = NULL,
                                 Disease_prog = NULL , 
                                 Cascade = NULL, end_Y = 50)

Sus_allPlot <- indicatorPlot(Sus_all_range, 
                             ylabel = "Number of susceptible",
                             xlimits = c(HCV$startYear, 
                                         HCV$startYear+21, 1),
                             calibration_Y = 2004,
                             rangeun = "y", 
                             groupPlot = NULL, 
                             facetPlot = NULL, 
                             simulateYear = HCV$simY) +
  ggtitle("Number of susceptible") + expandy(Sus_all_range , HCV$startYear + 21)

###### sub-population ######
SusPop <- popResults_MidYear(HCV, bestResults,Population = HCV$popNames,
                              Disease_prog = "s", 
                              Cascade = NULL, 
                              param = paramResults ,endYear = 50)%>%
  as.data.frame()%>%group_by(year, population)%>%
  summarise_at(.vars = c(colnames(popS)[-1]),sum)%>%ungroup()
gc()

SusPop_range <-popResults_range(HCV, SusPop, Population = HCV$popNames,
                                 Disease_prog = NULL, 
                                 Cascade = NULL, end_Y = 50)

SusPopPlot <- indicatorPlot(SusPop_range, 
                             ylabel = "Number of susceptible",
                             xlimits = c(HCV$startYear, 
                                         HCV$startYear+21, 1),
                             calibration_Y = 2004,
                             rangeun = "y", 
                             groupPlot = NULL, 
                             facetPlot = population, 
                             simulateYear = HCV$simY) +
  ggtitle("Number of susceptible")

##### by all components ##### 
allstate <- popResults_MidYear(HCV, bestResults, Population = NULL,
                              Disease_prog = HCV$progress_name, 
                              Cascade = HCV$cascade_name, 
                              param = paramResults, endYear = 50)%>%
  as.data.frame()

allstate_range <- popResults_range(HCV, allstate, Population = NULL, 
                                   Disease_prog = HCV$progress_name, 
                                   Cascade = HCV$cascade_name, end_Y = 50)

# reorder
allstate_range$state <- factor(allstate_range$state,
                                     levels = HCV$component_nameWOS)

# plot 
allstatePlot <- indicatorPlot(allstate_range, ylabel = "number of each disease stage", 
                              facetPlot = NULL, 
                              calibration_Y = 2004,
                              rangeun = "t", 
                              xlimits = c(HCV$startYear, 
                                          HCV$startYear+21, 1),
                              groupPlot = NULL, 
                              observationData = NULL, 
                              simulateYear = HCV$simY) + 
  ggtitle("Numbers of people in different HCV stage")

allstatePlot <- allstatePlot +facet_wrap(state~., scale = "free", nrow = 6) +
  scale_x_continuous(limits = c(1, 22), breaks = seq(1, 22,5), 
                     labels = seq(2004, 2025, 5))

# sub-population 

popstate <- popResults_MidYear(HCV, bestResults, Population = HCV$popNames,
                               Disease_prog =HCV$progress_name, 
                               Cascade = HCV$cascade_name, 
                               param = paramResults, endYear = 50)%>%
  as.data.frame()


popstate_range <- popResults_range(HCV, popstate, Population = HCV$popNames, 
                                   Disease_prog = HCV$progress_name, 
                                   Cascade = HCV$cascade_name, end_Y = 50)%>%
  filter(state != "s") 

popstate_range$state <- factor(popstate_range$state,
                               levels = HCV$component_nameWOS)

popstate_rangeSp <- list() 

# setting length of list 
Lnlst <- lapply(seq(1,length(HCV$component_nameWOS), by =HCV$ncascade ), 
                function(i) HCV$component_nameWOS[i:(i+ HCV$ncascade-1)])
Lnlst[[4]]
for (i in 1:length(Lnlst)){ 
  popstate_rangeSp[[i]] <- popstate_range%>%filter(state%in%Lnlst[[i]]) 
  }

popstate_rangeSp[[1]]
# plot 
popstatePlot <- list()

popstatePlot <- lapply(popstate_rangeSp, function(x) { 
  a <- indicatorPlot(x, 
                     ylabel = "Numbers at each disease stage by population", 
                     facetPlot = NULL, 
                     calibration_Y = 2004,
                     rangeun = "t", 
                     xlimits = c(HCV$startYear, 
                                 HCV$startYear+21, 1),
                     groupPlot = NULL, 
                     observationData = NULL, 
                     simulateYear = HCV$simY) + 
    ggtitle("Numbers of people in each HCV stage")
  
  a <- a + facet_wrap(.~state + population, scale ="free", ncol = 4)
  
})



##### by cascade ##### 
allcascade <- popResults_MidYear(HCV, bestResults, Population = NULL,
                               Disease_prog = NULL, 
                               Cascade = HCV$cascade_name, 
                               param = paramResults, endYear = 50)%>%
  as.data.frame()

cascade_range <- popResults_range(HCV, allcascade, Population = NULL, 
                                   Disease_prog = NULL, 
                                   Cascade = HCV$cascade_name, end_Y = 50)%>%
  group_by(year, cascade)

# reorder
cascade_range$cascade <- factor(cascade_range$cascade,
                             levels = HCV$cascade_name)


# plot 
allcascade_rangePlot <- indicatorPlot(cascade_range, ylabel = "", 
                              facetPlot = NULL, 
                              calibration_Y = 2004,
                              rangeun = "t", 
                              xlimits = c(HCV$startYear, 
                                          HCV$startYear+21, 1),
                              groupPlot = NULL, 
                              observationData = NULL, 
                              simulateYear = HCV$simY) + 
  ggtitle("Numbers of people in different HCV stage")

allcascade_rangePlot <- allcascade_rangePlot + facet_wrap (~cascade,
                                                             scales = "free", ncol = 2)

# sub group 

popCascade <- popResults_MidYear(HCV, bestResults, Population = HCV$popNames,
                   Disease_prog = NULL, 
                   Cascade = HCV$cascade_name, 
                   param = paramResults, endYear = 50)%>%
  as.data.frame()

popCascade_range <- popResults_range(HCV, popCascade, Population = HCV$popNames, 
                                  Disease_prog = NULL, 
                                  Cascade = HCV$cascade_name, end_Y = 50)%>%
  group_by(year, population, cascade)

# reorder
popCascade_range$cascade <- factor(popCascade_range$cascade,
                                levels = HCV$cascade_name)

popCascade_range[is.na(popCascade_range)] <-0

str(popCascade_range)
popCascade_rangePlot <- indicatorPlot(popCascade_range, ylabel = "", 
                                      facetPlot = NULL, 
                                      calibration_Y = 2004,
                                      rangeun = "t", 
                                      xlimits = c(HCV$startYear, 
                                                  HCV$startYear+21, 1),
                                      groupPlot = NULL, 
                                      observationData = NULL, 
                                      simulateYear = HCV$simY) + 
  ggtitle("Numbers of people in different HCV stage by population")

popCascade_rangePlot <- popCascade_rangePlot + 
  facet_wrap(.~ cascade + population ,ncol = 4, scales = "free")





##### by disease progress #####
# all population 

DisProg <- popResults_MidYear(HCV, bestResults, Population = NULL,
                   Disease_prog = c("s",HCV$progress_name), 
                   Cascade = NULL, param = paramResults, endYear = 50)%>%
  as.data.frame()

DisProg_range <- popResults_range(HCV, DisProg, Population = NULL, 
                               Disease_prog = c("s",HCV$progress_name),
                               Cascade = NULL, end_Y = 50)%>%
  group_by(year, disease_prog)

# reorder
DisProg_range$disease_prog <- factor(DisProg_range$disease_prog,
                                     levels = c("s", HCV$progress_name)) 




DisProgPlot <- indicatorPlot(DisProg_range, ylabel = "number of each disease stage", 
                              facetPlot = NULL, 
                              calibration_Y = 2004,
                              rangeun = "t", 
                              xlimits = c(HCV$startYear, 
                                          HCV$startYear+21, 1),
                              groupPlot = NULL, 
                              observationData = NULL, 
                              simulateYear = HCV$simY) + 
  ggtitle("Numbers of people in different HCV stage")

DisProgPlot <- DisProgPlot + facet_wrap(.~disease_prog, scale ="free_y")


                                               


# sub-population
popDisProg <- popResults_MidYear(HCV, bestResults, Population = HCV$popNames,
                   Disease_prog = HCV$progress_name, 
                   Cascade = NULL, param = paramResults, endYear = 50)%>%
  as.data.frame()%>%group_by(year, population, disease_prog)


popDisProg_range <- popResults_range(HCV, popDisProg, Population = HCV$popNames, 
                                  Disease_prog = HCV$progress_name,
                                  Cascade = NULL, end_Y = 50)

# reorder
popDisProg_range$disease_prog <- factor(popDisProg_range$disease_prog,
                                     levels = HCV$progress_name)
popDisProg_rangeSp <- list()
popDisProg_rangeSp[[1]] <- popDisProg_range%>%
  filter(disease_prog%in%HCV$progress_name[1:3])
popDisProg_rangeSp[[2]] <- popDisProg_range%>%
  filter(disease_prog%in%HCV$progress_name[4:6])
popDisProg_rangeSp[[3]] <- popDisProg_range%>%
  filter(disease_prog%in%HCV$progress_name[7:10])

# plot 
popDisProgPlot <- list()

popDisProgPlot <- lapply(popDisProg_rangeSp, function(x) { 
  a <- indicatorPlot(x, 
                     ylabel = "number of each disease stage", 
                     facetPlot = NULL, 
                     calibration_Y = 2004,
                     rangeun = "t", 
                     xlimits = c(HCV$startYear, 
                                 HCV$startYear+21, 1),
                     groupPlot = NULL, 
                     observationData = NULL, 
                     simulateYear = HCV$simY) + 
    ggtitle("Numbers of people in different HCV stage")
  
  a <- a + facet_wrap(.~disease_prog + population, ncol = 4, scales = "free")
  
})



##### flow  #####

# all population
# remove allpop, poptransition, inflow, and outflow 
Incprog <- names(bestResults)
Dlst <- Incprog[!Incprog%in%c("allPops","newS", "newpop_tran", 
                              "inflow", "outflow","death_hcv", "HCVdeathState", 
                              "newDeathState")]

cascadeflow <- lapply(Dlst, function(x) 
  indicatorResults(HCV, bestResults, x, 
                   pop="all",
                   paramR = paramResults, range = "y",
                   endY = 50) )

cascadeflowPlot <- lapply(seq_along(cascadeflow), 
                          function(i){ 
                            a <- indicatorPlot(cascadeflow[[i]],
                                               ylabel = Dlst[i], 
                                               facetPlot = NULL,
                                               calibration_Y = 2004,
                                               rangeun = "y", 
                                               xlimits = c(HCV$startYear, 
                                                           HCV$startYear+21, 1), 
                                               groupPlot = NULL, 
                                               observationData = NULL, 
                                               simulateYear = HCV$simY) + 
                              ggtitle(paste0("Numbers of", " ", Dlst[i]))
                            
                             a <- a + expandy(cascadeflow[[i]], HCV$startYear+21) 
                            })

# sub-pop # 

cascadePopflow <- lapply(Dlst, function(x) 
  indicatorResults(HCV, bestResults, x, 
                   pop=HCV$popNames,
                   paramR = paramResults, range = "y",
                   endY = 50) )

# y-scale limint

limFacet <- lapply(cascadePopflow, function(x) limit_facet(x, HCV$startYear+21))


# list of plot
cascadePopflowPlot <- lapply(
  seq_along(cascadePopflow), 
  function(i){
    a <- indicatorPlot(cascadePopflow[[i]],
                       ylabel = Dlst[i],
                       facetPlot = NULL,
                       calibration_Y = 2004,
                       rangeun = "y", 
                       xlimits = c(HCV$startYear, 
                                   HCV$startYear+21, 1), 
                       groupPlot = NULL, 
                       observationData = NULL,
                       simulateYear = HCV$simY) +
      ggtitle(paste0("Numbers of", " ", Dlst[i], 
                     " ", "by population"))
    a <- a + 
      facet_custom (~population,
                    scales = "free", ncol = 2,
                    scale_overrides = 
                      list(
                        scale_new(1,
                                  scale_y_continuous(limits = 
                                                       c(limFacet[[i]][1],
                                                         limFacet[[i]][5]))),
                        scale_new(2,
                                  scale_y_continuous(limits = 
                                                       c(limFacet[[i]][2],
                                                         limFacet[[i]][6]))),
                        scale_new(3,
                                  scale_y_continuous(limits = 
                                                       c(limFacet[[i]][3],
                                                         limFacet[[i]][7]))),
                        scale_new(4,
                                  scale_y_continuous(limits = 
                                                       c(limFacet[[i]][4],
                                                         limFacet[[i]][8])))
                        ) )
  })




  
cascadeflow_lay <- do.call("grid.arrange", c(cascadeflowPlot, ncol=4))  

cascadeflow_lay1 <- do.call("grid.arrange", c(cascadePopflowPlot[c(1:4)], ncol=4))
cascadeflow_lay2 <- do.call("grid.arrange", c(cascadePopflowPlot[c(5:8)], ncol=4))  
cascadeflow_lay3 <- do.call("grid.arrange", c(cascadePopflowPlot[c(9:10)], ncol=3)) 


#### Section 4. output setting ####
##### grid.arrange for calibration plots #####
lay <- rbind(c(1,1,1,1,2,2,3,3),
             c(1,1,1,1,2,2,3,3),
             c(4,4,5,5,6,6,7,7),
             c(4,4,5,5,6,6,7,7))
Plot_CaliPop <- grid.arrange(MSMPopPlot, HIVPrevPlot,HIVIncPlot, HIVDiagPlot, 
                             HIVDiagnumPlot,PrEPHIVPlot, PrEPHIVIncPlot, 
                             layout_matrix = lay) 


lay2 <- rbind(c(1,1,2,2,3,3,4,4),
             c(1,1,2,2,3,3, 4,4),
             c(5,5,6,6,7,7,8,8),
             c(5,5,6,6,7,7,8,8))
             
Plot_CaliHCV <- grid.arrange(totalPrevPlot, popPrevPlot, HCVIncPlot, 
                             HCVPopIncPlot, totalCumTreatmentPlot, 
                             PopCumTreatmentPlot, totalCumCurePlot, 
                             PopCumCurePlot, 
                             layout_matrix = lay2) 

lay3 <- rbind(c(1,1,2,2,3,3),
              c(4,4,5,5,6,6),
              c(4,4,5,5,6,6))

Plot_SIC <-grid.arrange(Sus_allPlot, Infected_allPlot, Cured_allPlot, 
                        SusPopPlot,InfectedPop_allPlot, CuredPopPlot,
                        layout_matrix = lay3)

allstateNum <- grid.arrange(allstatePlot)
allCascadeNum <- grid.arrange(allcascade_rangePlot)
allDisprogNum <- grid.arrange(DisProgPlot) 

#### HCV death ####
HCdeath <- indicatorResults(HCV, bestResults, "newHCVdeaths", 
                            pop="all",
                            paramR = paramResults, range = "y",
                            endY = 50)
HCdeath <- HCdeath%>%select(!c(min, max, Med, Mu, q5, q25, q75, q95))
tempHCVmor_all <- cbind(year = seq(HCV$startYear , 50-1 ,1),
                      as.data.frame((HCdeath[, -1]/(tempTotal_all[, -1] - 
                                       tempNOTInfected_all[, -1])*100)))%>%
  tibble::as_tibble() 

tempHCVmor_all_range <- popResults_range(HCV, tempHCVmor_all, Population = NULL,
                                       Disease_prog = NULL , 
                                       Cascade = NULL, end_Y = 50)

HCdeathPop <- indicatorResults(HCV, bestResults, "newHCVdeaths", 
                               pop=HCV$popNames,
                               paramR = paramResults, range = "y",
                               endY = 50)

HCdeathPop <- HCdeathPop%>%select(!c(min, max, Med, Mu, q5, q25, q75, q95))

tempHCVmor <- cbind(year = rep(seq(HCV$startYear , 50-1 ,1), each = HCV$npops),
                  population = HCV$popNames,
                  
                  as.data.frame(100*HCdeathPop[, -c(1,2)]/(tempTotal[, -c(1,2)] - 
                                       tempNOTInfected[ ,-c(1,2)])))%>%
  tibble::as_tibble() 


# summarise q5/q95 and median for HCV prevalence in MSM subgroups

tempHCVmor_range <- popResults_range(HCV, tempHCVmor, Population = HCV$popNames,
                                   Disease_prog = NULL , 
                                   Cascade = NULL, end_Y = 50)



#### 
#### treatment coverage: cumulative(newtreat)/[people living with HCV_y + cure (except for acute)] ####

statepop <- popResults_MidYear(HCV, bestResults, Population = NULL,
                               Disease_prog = HCV$progress_name, 
                               Cascade = HCV$cascade_name, 
                               param = paramResults, endYear = 50, 
                               YearCut = "End")%>%
  as.data.frame()

statepop <- statepop%>%filter(!state%in% c("s", "a_cured"))%>%group_by(year)%>%
  summarise_at(.vars = c(colnames(popS)[-1]),sum) 



livingHCVCure <- cbind(year = seq(HCV$startYear , 50-1 ,1),
                       as.data.frame(TreatmentInit_cumulative[,-1]/(statepop[, -1] )*100))


statepopPop <- popResults_MidYear(HCV, bestResults, Population = HCV$popNames,
                               Disease_prog = HCV$progress_name, 
                               Cascade = HCV$cascade_name, 
                               param = paramResults, endYear = 50, 
                               YearCut = "End")%>%
  as.data.frame()
statepopPop <- statepopPop%>%filter(!state%in% c("s", "a_cured"))%>%group_by(year, population)%>%
  summarise_at(.vars = c(colnames(popS)[-1]),sum) 


livingHCVCure_pop <- cbind(year = rep(seq(HCV$startYear , 50-1 ,1), each = HCV$npops),
      population = HCVInfectPop$population,
      as.data.frame(TreatmentInitPop_cumulative[,-c(1,2)]/(statepopPop[, -c(1,2)] )*100))%>%
  tibble::as_tibble() 
# grid object put in the list #
# list of plot put in c()
#### output ####
g1 <- marrangeGrob(grobs=c(list(Plot_CaliPop,Plot_CaliHCV,Plot_SIC),
                           list(allCascadeNum,allDisprogNum),
                           list(cascadeflow_lay, cascadeflow_lay1,
                                cascadeflow_lay2, cascadeflow_lay3, allstateNum),
                           c(popstatePlot,popDisProgPlot)), 
                   nrow=1, ncol=1)
ggsave(file="calibrationPlot20220901.pdf", g1, width = 15 , height = 10.2) 



 save(popS_range, male_pop, MSMPop_range, MSMHIVnum, HIV_range, HIVPrev, HIVdiag,
     HIVDiagRate_range, PrEPDropInc_range, HIVPrEPdrop, PrEPtoHIVP_range, 
     PrEPHIVInc_range, MSMHIVInc, HIVInc_range, 
     tempNOTInfected_all, tempTotal_all, tempPrev_all_range,
     tempNOTInfected, tempTotal, tempPrev_rang, MSMHCVPrev,
     HCVInfect_all, HCVInc_range, HCVPopInc_range, HCVInfectPop, 
     HCVPopInc_range, MSMHCVInc, 
     TreatmentInit,  TreatmentInit_range, TreatmentInitPop_range, 
     TreatmentInitPop_range, TreatmentInit_cumulativeRange, HCVTreat, 
     TreatmentInitPop_cumulative, TreatmentInitPop_cumulativeRange, CureInit, 
     CureInitPop, CureInit_cumulativeRange, HCVsvr, CureInitPop_cumulativeRange, 
     tempHCVmor_all_range, tempHCVmor_range,HCdeath, HCdeathPop,livingHCVCure,
     livingHCVCure_pop, allstate,
     file = file.path(projectFolder,
                      paste0(project_name, "SumResults", ".rda")))
