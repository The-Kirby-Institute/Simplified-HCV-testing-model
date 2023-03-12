# expected to be a function to update the model parameters, and rerun the HCVMSM 
# model for the best estimate parameters. Then output the model plot. 
library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)

basePath <- getwd()
dataPath <- file.path(basePath, "01. Data/model input") 
Rcode <- file.path(basePath, "03. Code")
urrTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
runSamples <- TRUE
saveAsBase <- TRUE  # if TRUE doesn't append time to results and overwrites
# a base file. Useful for storing main results or
# testing

source(file.path(Rcode, "/Functions/HCV_model.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R"))

# extend to same length 

best_estimates <- as.data.frame(matrix(0, ncol = ncol(constantsDf),
                                       nrow = HCV$npts))
colnames(best_estimates) <- colnames(constantsDf)

for (var in colnames(constantsDf)) {
  
  tempValues <- constantsDf[, var]
  
  for (year in 1:(HCV$nyears-1)) {
    
    indices <- ((year-1)* 1/HCV$timestep + 1): (year * 1/HCV$timestep + 1)
    
    yearValues <- seq(tempValues[year], tempValues[year+1], 
                      length = (1 + 1/HCV$timestep))
    
    best_estimates[indices, var] <- yearValues
  }  
}

#### manually change the entry number  
MSMentry <- read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"),
                               "MSM_entry.csv"),header = TRUE)%>%
  as.data.frame()

MSMentry[c(18:2000),1] <- seq(2021, 4003,1)

MSMentry2021 <- MSMentry[c(1:17),2]

MSMafrer2021 <- rep(MSMentry[17,2], 1983)

MSM_Entry <- MSMentry%>%
  mutate(number=c(MSMentry2021, MSMafrer2021))%>%
  mutate(number = number*0.052)

Ent <- as.data.frame(matrix(0, ncol = 1,
                            nrow = HCV$npts)) 


for (year in 1:(HCV$nyears-1)){
  
  indices <- ((year-1)* 1/HCV$timestep + 1): (year * 1/HCV$timestep + 1)
  
  yearvalue <- seq(MSM_Entry$number[year], MSM_Entry$number[year+1], 
                   length = (1 + 1/HCV$timestep))
  
  Ent[indices, 1] <-yearvalue 
  
}


best_estimates <- best_estimates%>%mutate(entry = Ent$V1) 

#### finding steady stage ####

tic <- proc.time()

steady<- HCVMSM(HCV, HCV$timestep ,best_estimates, best_initial_pop,
                disease_progress, fib, pop_array, dfList, 
                modelrun="steady")

toc <- proc.time() - tic 


#### extract infected proportion as the infected population allocation #### 

df_list <- lapply(steady, as.data.frame.table)

allpop <- df_list$allPops%>%mutate(time = rep(seq(HCV$startYear,
                                                  (HCV$endYear - HCV$timestep),
                                                  0.1), 
                                              each=HCV$ncomponent*HCV$npops),
                                   Frequency=round(Freq, digits = 3))

popPro_extract <- allpop%>%
  mutate(time = rep(seq(HCV$startYear, (HCV$endYear - HCV$timestep), 0.1), 
                    each=HCV$ncomponent*HCV$npops),
         Frequency=round(Freq, digits = 3))%>%
  filter(time==1999.9)%>%
  mutate(cascade_status = sub("^[^_]*_", "", Var2), 
         dis_prog = sub("\\_.*", "", Var2),
         SI = ifelse(cascade_status%in%c("s", "cured"), "S","I"),
         parameter =Var2)%>%group_by(Var1 ,SI)%>%
  mutate(total = sum(Frequency),
         value = ifelse(Frequency==0,0, round(Frequency/total, digits = 4)))%>%
  ungroup()%>%group_by(Var1)%>%mutate(pop_prop = ifelse(
    Frequency==0,0, round(Frequency/sum(Frequency), digits = 4)))%>%
  ungroup()%>%select(Var1,parameter, value, SI)

basePath <- getwd()

write.csv(popPro_extract, 
          file.path(basePath,"01. Data/model input/Estimate_initial_pop.csv")) 

#### number of MSM in each population ####

estPops<- read.csv(file.path(DataFolder, "Estimate_initial_pop.csv"), 
                   header = TRUE)%>%select(-"X")

init_pop <- filter(initialPops, parameter == "init_pop")$value

pop_prop <- initialPops%>%filter(parameter%in% c("pop_prop1", "pop_prop2", 
                                                 "pop_prop3", "pop_prop4"))%>%
  select(value)%>%unlist()%>%as.vector()

popProp <- init_pop*pop_prop 

init_prop_s <- initial_pop%>%filter(parameter%in% c("s1", "s2", "s3", "s4"))%>%
  select(value)%>%unlist()%>%as.vector()

init_prop_I <-c(1-init_prop_s)

estPops <- estPops%>%mutate(
  pop_group = rep(c(popProp),dim(estPops)[1]/HCV$npops),
  SIprop = ifelse(estPops$SI=="S", 
                  rep(init_prop_s, HCV$diseaseprogress_n*HCV$npops),
                  rep(init_prop_I, HCV$ncomponent*HCV$npops - 
                        HCV$diseaseprogress_n*HCV$npops)),
  est_pop = value*pop_group*SIprop)


best_est_pop <- as.matrix(as.data.frame(matrix(estPops$est_pop, 
                                               ncol = HCV$ncomponent,  
                                               nrow = HCV$npops)))

colnames(best_est_pop) <- c(HCV$component_name)

#### calibration ####

tic <- proc.time()


calibrateInit<- HCVMSM(HCV, HCV$timestep ,best_estimates, best_est_pop,
                       disease_progress, fib, pop_array, dfList, modelrun="UN")

toc <- proc.time() - tic 

#### create calibration plot ---------------------------------------------------
#### total population ####
# total 

popSize <- popResults_MidYear(HCV, calibrateInit, Population = NULL,
                              Disease_prog = NULL , 
                              Cascade = NULL)%>%as.data.frame()
# real world data

male_pop <-read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"), 
                              "Malepop.csv"), header = TRUE)%>%
  as.data.frame()%>%mutate(realPop = Total.male.population*0.052,
                           time = Year-2003,
                           low = Total.male.population*0.035,
                           up = Total.male.population*0.069)%>%
  dplyr::select(realPop, time, low, up)%>%as.data.frame()

totalPopPlot <- indicatorPlot(popSize, ylabel = "Population size", 
                              facetPlot = NULL, 
                              calibration_Y = 2004,
                              range = FALSE, 
                              xlimits = c(HCV$startYear, 
                                          HCV$startYear+21, 1),
                              groupPlot = NULL, 
                              observationData = male_pop) + 
  ggtitle("Overall population")  

# adjust y-axis limits 
totalPopPlot <- totalPopPlot + scale_y_continuous(limits = c(100000, 2000000))

# population by group

# notification data

MSMHIVnum <-read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"), 
                               "HIV_num.csv"), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year-2003, realPop = number,
                           low = lower,
                           up =upper)

MSMPop <- popResults_MidYear(HCV, calibrateInit,Population = HCV$popNames,
                             Disease_prog = NULL , 
                             Cascade = NULL)%>%as.data.frame() 

MSMPopPlot <- indicatorPlot(MSMPop, ylabel = "Population size by MSM", 
                            facetPlot = population,
                            calibration_Y = 2004,
                            range = FALSE, 
                            xlimits = c(HCV$startYear, 
                                        HCV$startYear+21, 1),
                            groupPlot = NULL, 
                            observationData = MSMHIVnum) + 
  geom_segment(data = MSMHIVnum, aes ( y = low, yend = up, x = time, xend = time))+
  geom_vline(xintercept=2021-2003, linetype = "dashed", size = 1)

# adjust y-axis limits in each facet 
MSMPopPlot <- MSMPopPlot + 
  facet_custom (~population, scales = "free", ncol = 2, 
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(
                                limits = c(limit_facet(MSMPop)[1],
                                           limit_facet(MSMPop)[1+4]))),
                    scale_new(2,
                              scale_y_continuous(
                                limits = c(limit_facet(MSMPop)[2],
                                           limit_facet(MSMPop)[6]))),
                    scale_new(3,
                              scale_y_continuous(
                                limits = c(limit_facet(MSMPop)[3], 
                                           limit_facet(MSMPop)[3+4]))),
                    scale_new(4,
                              scale_y_continuous(
                                limits = 
                                  c(limit_facet(MSMPop)[4],
                                    limit_facet(MSMPop)[4+4])))))

#### prevalence ####

tempNOTInfected_all <- popResults_MidYear(HCV, calibrateInit,Population = NULL,
                                          Disease_prog = c("s", "cured"), 
                                          Cascade = NULL)%>%group_by(year)

tempTotal_all <- popResults_MidYear(HCV, calibrateInit,Population = NULL,
                                    Disease_prog = NULL, 
                                    Cascade = NULL)

tempPrev_all <- data.frame(year = seq(HCV$startYear,HCV$endYear-1,1), 
                           best = (100 * (tempTotal_all$best-
                                            tempNOTInfected_all$best) / 
                                     tempTotal_all$best))%>%tbl_df()

totalPrevPlot <- indicatorPlot(tempPrev_all, 
                               ylabel = "Overall HCV prevalence (%)",
                               xlimits = c(HCV$startYear, 
                                           HCV$startYear+21, 1),
                               calibration_Y = 2004,
                               range = FALSE, 
                               groupPlot = NULL, 
                               facetPlot = NULL) +
  ggtitle("Overall HCV prevalence")

totalPrevPlot <- totalPrevPlot + expandy(totalPrevPlot)

# HCV prevalence by group 

tempNOTInfected <- popResults_MidYear(HCV, calibrateInit,
                                      Population = HCV$popNames,
                                      Disease_prog = c("s", "cured"), 
                                      Cascade = NULL)%>%
  group_by(year, population)%>%summarise(best = best)

tempTotal <- popResults_MidYear(HCV, calibrateInit,
                                Population = HCV$popNames,
                                Disease_prog = NULL, 
                                Cascade = NULL)

tempPrev <- data.frame(year = rep(seq(HCV$startYear,
                                      HCV$endYear-1,1),each=HCV$npops), 
                       population = HCV$popNames,
                       best = (100 *(tempTotal$best - tempNOTInfected$best) / 
                                 tempTotal$best))%>%tbl_df()

# population convert to factors 
# be careful about the NaN value 

tempPrev <- tempPrev%>%FactorPop(., HCV$popNames)

tempPrev[is.na(tempPrev)] <- 0
# real world data 

MSMHCVPrev <-read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"), 
                                "HCVPrev.csv"), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year-2003, 
                           realPop = HCV.prevalence*100,
                           up = upper*100,
                           low = lower*100)

popPrevPlot <- indicatorPlot(tempPrev, 
                             ylabel = "Population HCV prevalence (%)",
                             xlimits = c(HCV$startYear, 
                                         HCV$startYear+21, 1),
                             calibration_Y = 2004,
                             range = FALSE, 
                             groupPlot = NULL, 
                             facetPlot = population,
                             observationData = MSMHCVPrev) +
  ggtitle("HCV prevalence by population") 

popPrevPlot <- popPrevPlot + 
  facet_custom (~population,scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(
                                limits =c(0,15))),
                    scale_new(2,
                              scale_y_continuous(
                                limits = c(-1,1))),
                    scale_new(3, 
                              scale_y_continuous(
                                limits = c(0,20))),
                    scale_new(4,
                              scale_y_continuous(
                                limits = c(0,40))))) 


# New Infections 

NewI <- indicatorResults(HCV, calibrateInit, "newInfections", population="all")

NewIPlot <- indicatorPlot(NewI, ylabel = "Number of HCV new Infections", 
                          facetPlot = NULL,
                          calibration_Y = 2004,
                          range = FALSE, 
                          xlimits = c(HCV$startYear, 
                                      HCV$startYear+21, 1),
                          groupPlot = NULL, 
                          observationData = NULL) + 
  ggtitle("Overall HCV new infections")

NewIPlot <- NewIPlot + expandy(NewIPlot)

# by group 

NewIgroup <- indicatorResults(HCV, calibrateInit, "newInfections")

NewIgroupPlot <- indicatorPlot(NewIgroup, ylabel = "Number of HCV new Infections", 
                               facetPlot = population,
                               calibration_Y = 2004,
                               range = FALSE, 
                               xlimits = c(HCV$startYear, 
                                           HCV$startYear+21, 1),
                               groupPlot = NULL, 
                               observationData = NULL) + 
  ggtitle("HCV new infections by MSM group")

NewIgroupPlot <- NewIgroupPlot + 
  facet_custom (~population,scales = "free", ncol = 2, 
                scale_overrides = list(
                  scale_new(1,
                            scale_y_continuous(
                              limits = c(limit_facet(NewIgroup)[1],
                                         limit_facet(NewIgroup)[5]))),
                  scale_new(2,
                            scale_y_continuous(
                              limits = c(limit_facet(NewIgroup)[2], 
                                         limit_facet(NewIgroup)[6]))),
                  scale_new(3, scale_y_continuous(
                    limits = c(limit_facet(NewIgroup)[3], 
                               limit_facet(NewIgroup)[7]))),
                  scale_new(4,scale_y_continuous(
                    limits = c(limit_facet(NewIgroup)[4],
                               limit_facet(NewIgroup)[8])))))

# cumulative
NewIcum <- NewI%>%mutate(best = cumsum(best))

totalCumInfectionsPlot <- indicatorPlot(NewIcum, 
                                        ylabel = "Number of HCV new Infections", 
                                        facetPlot = NULL,
                                        calibration_Y = 2004,
                                        range = FALSE, 
                                        xlimits = c(HCV$startYear, 
                                                    HCV$startYear+21, 1),
                                        groupPlot = NULL, 
                                        observationData = NULL) + 
  ggtitle("Overall cumulative infections") 

totalCumInfectionsPlot <- totalCumInfectionsPlot + 
  expandy(totalCumInfectionsPlot)


# cumulative by group 
NewIcumGroup <- NewIgroup%>%group_by(population)%>%mutate(best = cumsum(best))

groupCumInfectionsPlot <- indicatorPlot(NewIcumGroup, 
                                        ylabel = "Number of HCV new Infections", 
                                        facetPlot = population,
                                        calibration_Y = 2004,
                                        range = FALSE, 
                                        xlimits = c(HCV$startYear, 
                                                    HCV$startYear+21, 1),
                                        groupPlot = NULL, 
                                        observationData = NULL) + 
  ggtitle("Cumulative infections by group") 

# exact limits for y-axis 

groupCumInfectionsPlot <- groupCumInfectionsPlot + 
  facet_custom (~population,scales = "free", ncol = 2, 
                scale_overrides = list(
                  scale_new(1,
                            scale_y_continuous(
                              limits = c(limit_facet(limit_y(NewIcumGroup, 21))[1],
                                         limit_facet(limit_y(NewIcumGroup, 21))[5]))),
                  scale_new(2,
                            scale_y_continuous(
                              limits = c(limit_facet(limit_y(NewIcumGroup, 21))[2], 
                                         limit_facet(limit_y(NewIcumGroup, 21))[6]))),
                  scale_new(3, scale_y_continuous(
                    limits = c(limit_facet(limit_y(NewIcumGroup, 21))[3], 
                               limit_facet(limit_y(NewIcumGroup, 21))[7]))),
                  scale_new(4,scale_y_continuous(
                    limits = c(limit_facet(limit_y(NewIcumGroup, 21))[4],
                               limit_facet(limit_y(NewIcumGroup, 21))[8])))))

# incidence  
NewIAll_incidence <- dplyr::select(NewI, -year)

NewIGroup_incidence <- NewIgroup%>%ungroup()%>%
  dplyr::select(., -year)%>%arrange(population)

popSize_incidence <- popSize%>%dplyr::select(., -year)

MSMPop_incidence <- MSMPop%>%dplyr::select(., -year)%>%
  arrange(population)%>%FactorPop(HCV$popNames)



TotalInc <- data.frame(year = head(HCV$years, -1),
                       best = (1e3 * NewIAll_incidence$best/
                                 popSize$best)) %>%
  tbl_df()

totalIncidencePlot <- indicatorPlot(TotalInc,
                                    ylabel = "Incidence per 1,000",
                                    xlimits = c(HCV$startYear,
                                                HCV$startYear+21, 1),
                                    calibration_Y = 2004,
                                    range = FALSE,
                                    groupPlot = NULL,
                                    facetPlot = NULL) +
  ggtitle("Overall HCV incidence") 


totalIncidencePlot <- totalIncidencePlot + expandy(totalIncidencePlot)

# by group 
## data 
HCVinc <- read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"),
                             "HCVInc.csv"),
                   header = TRUE)%>%as.data.frame()%>%
  mutate(realPop = HCV.incidence*1000,time = year-2003,
         up = upper*1000,
         low = lower*1000)%>%
  dplyr::select(realPop, time, population, up, low)%>%as.data.frame()

tempPopInc <- NewIgroup%>%ungroup()%>%mutate(Pop = MSMPop$best,
                                             incidentCase = best,
                                             best = incidentCase/Pop*1000)%>%
  tbl_df()


#data.frame(year = rep( head(HCV$years, -1), HCV$npops),
#                       population = savePops,
#                       best = (1e3 * NewIGroup_incidence$best /
#                                 MSMPop_incidence$best)) %>%
#tbl_df()


popIncidencePlot <- indicatorPlot(tempPopInc,
                                  ylabel = "Incidence per 1,000",
                                  xlimits = c(HCV$startYear,
                                              HCV$startYear+21, 1),
                                  calibration_Y = 2004,
                                  range = FALSE,
                                  groupPlot = NULL,
                                  facetPlot = population,
                                  observationData = HCVinc) +
  ggtitle("HCV incidence by MSM group") 


popIncidencePlot <- popIncidencePlot + 
  facet_custom (~population,scales = "free", ncol = 2,
                scale_overrides = list(
                  scale_new(1,
                            scale_y_continuous(
                              limits = c(0,10))),
                  scale_new(2,
                            scale_y_continuous(
                              limits = c(-1,1))),
                  scale_new(3, scale_y_continuous(
                    limits = c(0,20))),
                  scale_new(4,scale_y_continuous(
                    limits = c(0,40))))) 

NewSAll <- indicatorResults(HCV, calibrateInit, "newS", population="all") 

 NewSallPlot <- indicatorPlot(NewSAll, 
                             ylabel = "New entry",
                             xlimits = c(HCV$startYear, 
                                         HCV$startYear+21, 1),
                             calibration_Y = 2004,
                             range = FALSE, 
                             groupPlot = NULL, 
                             facetPlot = NULL,
                             observationData = NULL)  + 
  ggtitle("Number of new entry")

 NewSallPlot <- NewSallPlot + expandy(NewSallPlot)

# by group 
 
 MSMHIVnum <-read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"), 
                                "HIV_num.csv"), header = TRUE)%>%
   as.data.frame()%>%mutate(time = year-2003, realPop = number,
                            low = lower,
                            up = upper)
 
 MSMPop <- popResults_MidYear(HCV, calibrateInit,Population = HCV$popNames,
                              Disease_prog = NULL , 
                              Cascade = NULL)%>%as.data.frame() 
 
 MSMPopPlot <- indicatorPlot(MSMPop, ylabel = "Population size by MSM", 
                             facetPlot = population,
                             calibration_Y = 2004,
                             range = FALSE, 
                             xlimits = c(HCV$startYear, 
                                         HCV$startYear+21, 1),
                             groupPlot = NULL, 
                             observationData = MSMHIVnum) +
   ggtitle("number of HIV infected MSM")
 
 # adjust y-axis limits in each facet 
 MSMPopPlot <- MSMPopPlot + 
   facet_custom (~population, scales = "free", ncol = 2, 
                 scale_overrides = 
                   list(
                     scale_new(1,
                               scale_y_continuous(
                                 limits = c(limit_facet(MSMPop)[1],
                                            limit_facet(MSMPop)[1+4]))),
                     scale_new(2,
                               scale_y_continuous(
                                 limits = c(limit_facet(MSMPop)[2],
                                            limit_facet(MSMPop)[6]))),
                     scale_new(3,
                               scale_y_continuous(
                                 limits = c(limit_facet(MSMPop)[3], 
                                            limit_facet(MSMPop)[3+4]))),
                     scale_new(4,
                               scale_y_continuous(
                                 limits = 
                                   c(limit_facet(MSMPop)[4],
                                     limit_facet(MSMPop)[4+4])))))

 
 HIVInfected <- popResults_MidYear(HCV, calibrateInit,
                                   Population = c("HIV+", "HIV+d"),
                                   Disease_prog = NULL, 
                                   Cascade = NULL) 
 HIVInf <- HIVInfected%>%ungroup()
 
 HIVInf <-HIVInf%>% dplyr::group_by(year)%>%dplyr::summarise(HIVInf = sum(best))
 
 
 HIVdiag <- read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"), 
                               "HIVdiag.csv"), header = TRUE)%>%
   as.data.frame()%>%
   mutate(time = year - 2003, 
          realPop = HIV.diagnosis.rate*100,
          up = upper*100,
          low = lower*100)
 
 HIVInfected <- popResults_MidYear(HCV, calibrateInit,Population = c("HIV+d"),
                                   Disease_prog = NULL, 
                                   Cascade = NULL)
 
 HIVDiagRate <- data.frame(year = seq(HCV$startYear,
                                      HCV$endYear-1,1), 
                           best = (100 * (HIVInfected$best/HIVInf$HIVInf)))%>%
   tbl_df() 
 
 
 HIVDiagPlot <- indicatorPlot(HIVDiagRate, 
                              ylabel = "HIV Diagnosis among HIV positive MSM (%)",
                              xlimits = c(HCV$startYear, 
                                          HCV$startYear+21, 1),
                              calibration_Y = 2004,
                              range = FALSE, 
                              groupPlot = NULL, 
                              facetPlot = NULL,
                              observationData = HIVdiag) + 
   ggtitle("HIV Diagnosis (%)") 
 
 HIVDiagPlot <- HIVDiagPlot + expandy(HIVDiagPlot)
 
 
 
 
 
 #### HIV prevalence #### 
 MSMPop_HIV <- popResults_MidYear(HCV, calibrateInit,Population = HCV$popNames,
                              Disease_prog = NULL , 
                              Cascade = NULL)%>%as.data.frame()%>%
   filter(population%in%c("HIV+", "HIV+d"))%>%group_by(year)%>%
   summarise(total_inf = sum(best))
 
 
 tempPrev_HIV <- data.frame(year = seq(HCV$startYear,
                                       HCV$endYear-1,1), 
                            best = (100 *
                                      (MSMPop_HIV$total_inf / popSize$best))
                            )%>%tbl_df()
 
 HIVprev <- read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"), 
                               "HIVPrev.csv"), header = TRUE)%>%
   as.data.frame()%>%
   mutate(time = year - 2003, 
          realPop = HIV.prevalence*100,
          up = upper*100,
          low = lower*100)
 
 
 HIVprevPlot <- indicatorPlot(tempPrev_HIV, 
               ylabel = "HIV prevalence (%)",
               xlimits = c(HCV$startYear, 
                           HCV$startYear+21, 1),
               calibration_Y = 2004,
               range = FALSE, 
               groupPlot = NULL, 
               facetPlot = NULL,
               observationData = HIVprev) + 
   ggtitle("HIV Prevalence (%)") + 
   scale_y_continuous(limits = c(0,15))
 
 # using as caption 
 capBeta <- paste0("beta value:", paste0(best_estimates$beta1[1],", ",
                                     best_estimates$beta2[1],", ",
                                     best_estimates$beta3[1],", ",
                                     best_estimates$beta4[1], sep = ""), sep="")
 
 
 InitPro <-initialPops%>%filter(parameter%in% c("s1", "s2", "s3", "s4"))%>%
   select(value)%>%unlist()
 numIPro <- as.numeric(1-InitPro)%>%c()
 numIPro[2] <-0
 paste0(numIPro, sep =",")
 capInit <- paste0("Proportion of infected in Initial population:", 
                   paste0(numIPro[1],", ", numIPro[2],", ", numIPro[3], ", ",
                          numIPro[4]))
 
 ## using as table 
 table <- data.frame(beta = c(best_estimates$beta1[1], best_estimates$beta2[1],
                              best_estimates$beta3[1], best_estimates$beta4[1]),
                     Initial_pop = c(numIPro[1], numIPro[2], numIPro[3],
                                     numIPro[4])) 
 row.names(table) <- HCV$popNames
 table1_grob <- gridExtra::tableGrob(table)
 
 grid.arrange(popPrevPlot, popIncidencePlot, MSMPopPlot, HIVDiagPlot, 
              HIVprevPlot,table1_grob,
              ncol = 3)
 
 
 
 
 grid.arrange(popPrevPlot, popIncidencePlot,  
              nrow = 1,
              bottom = textGrob(capBeta,
                gp = gpar(fontface = 3, fontsize = 9),
                hjust = 1,
                x = 1),
              top = textGrob(capInit ,
                             gp = gpar(fontface = 3, fontsize = 9),
                             hjust = 1,
                             x = 1  
              ))


 
 
 

 

