# This script is preparing the data for figures and tables generation   
rm(list = ls()) 
gc()
library(here)
here()

# Load useful libraries
library("dplyr")
library("tidyr")
library("ggthemes")
library("stringr")
# Setup directories after setting working directory to source file 
# directory 

#file path of "TWHCV-model" project
codepath <- file.path(here() %>% dirname(), 'TWHCV-model/03. Code/Functions')

DataFolder <- file.path(here(), "01. DATA/model input")

# Rda file path 
# load the .rda file of base estimate 
rdapath <- file.path(here()%>%dirname(), "Taiwan-MSM-HCV-model")

# output file path 
# dir.create("02. Output") # create subdircetory 
# dir.create("02. Output/RDA") 
# dir.create("02. Output/Results")
# dir.create("02. Output/Figs")
outputdt <- here("02. Output/RDA")

outputfig <- here("02. Output/Figs")

# source 


source(file.path(codepath, "plotOptions.R"))

source(file.path(codepath, "plotFunctions.R"))

source(file.path(codepath, "plotManuscript.R"))

source(file.path(here(),"AggregateRes.R"))

# load rda files 
files <- list.files(path = paste0(outputdt, "/", sep = ""),
                    pattern = '^Outcome_.*\\.rda')
for( f in files){ 
  
  load(file.path(paste0(outputdt, "/" ,f, sep ="")))
  
  }

projectFile <- file.path(rdapath , paste0("HCVModel",".rda"))

projectVars <- load(projectFile)

discountR <- 0.03

Currency_converter <- 14.07

#### desired epi outcomes shown in figs and tables #### 
# (1) HCV incidence of testing strategies compared to base scenario overall MSM
# We have 4 scenarios in the main analysis means ouput 4 plots by scenarios  
# (2) HCV incidence of testing strategies compared to base scenario by targeted MSM
# We have 4 scenarios in the main analysis means ouput 4 plots, each plots contains number of subgroups targeted at 
# (3) number of HCV infections averted overall MSM (cumulative from year of testing scale up to target year) 
# (4) number of HCV infections averted by targeted MSM 

# names of each data set 
# combine results of different scenarios in a list 
Outcome_Base_epi <- list("PrEPsame" = Outcome_Base_PrEPsame$epi,
                         "PrEP" = Outcome_Base_PrEP$epi,
                         "HIVD" = Outcome_Base_HIVD$epi,
                         "PrEPnHIVD" = Outcome_Base_PrEPnHIVD$epi)

Outcome_Scen_epi <- list("PrEPsame" = Outcome_Scen_PrEPsame$epi,
                         "PrEP" = Outcome_Scen_PrEP$epi,
                         "HIVD" = Outcome_Scen_HIVD$epi,
                         "PrEPnHIVD" = Outcome_Scen_PrEPnHIVD$epi)

Outcome_Base_cost <- list("PrEPsame" = Outcome_Base_PrEPsame$costqaly,
                         "PrEP" = Outcome_Base_PrEP$costqaly,
                         "HIVD" = Outcome_Base_HIVD$costqaly,
                         "PrEPnHIVD" = Outcome_Base_PrEPnHIVD$costqaly)

Outcome_Scen_cost <- list("PrEPsame" = Outcome_Scen_PrEPsame$costqaly,
                         "PrEP" = Outcome_Scen_PrEP$costqaly,
                         "HIVD" = Outcome_Scen_HIVD$costqaly,
                         "PrEPnHIVD" = Outcome_Scen_PrEPnHIVD$costqaly)




dt_name <- names(Outcome_Base_epi)

Btt <- list()
Bttpop <- list()
Bttnumpop <- list()
Bttnumpop_reinf <- list()

for(n in dt_name){ 
  
  Btt[[n]] <- Outcome_Base_epi[[n]]$Inci
  
  Bttpop[[n]] <- Outcome_Base_epi[[n]]$Inci_pop
  
  Bttnumpop[[n]] <- Outcome_Base_epi[[n]]$infnum_pop
  
  Bttnumpop_reinf[[n]] <- Outcome_Base_epi[[n]]$reinfnum_pop
  
  
  }


Stt <- list()
Sttpop <- list()
Sttnumpop <- list()
Sttnumpop_reinf <- list()

for(n in dt_name){
  
  Stt[[n]] <- Outcome_Scen_epi[[n]]$Inci
  
  Sttpop[[n]] <- Outcome_Scen_epi[[n]]$Inci_pop
  
  Sttnumpop[[n]] <- Outcome_Scen_epi[[n]]$infnum_pop
  
  Sttnumpop_reinf[[n]] <- Outcome_Scen_epi[[n]]$reinfnum_pop

  }

Inci <- list()
Inci_pop <- list()

for(n in dt_name){ 
  Stt[[n]] <- popRange(Stt[[n]], pop = "y", test = NULL)%>%
    mutate(population = "All")
  Btt[[n]] <- popRange(Btt[[n]], pop = "y", test = NULL)%>%
    mutate(population = "All")
  
  Sttpop[[n]] <- popRange(Sttpop[[n]], pop = NULL, test = NULL, Casca = NULL)
  
  Bttpop[[n]] <- popRange(Bttpop[[n]], pop = NULL, test = NULL, Casca = NULL)
  
}

for(n in dt_name){ 
  Inci[[n]] <- rbind(Btt[[n]], Stt[[n]])%>%
    mutate(year = (year - HCV$cabY + 1),
           testing = factor(testing, levels = c("Status Quo", "POC_antibody",
                                                "Reflex_RNA", "POC_RNA")))%>%
    arrange(year, testing)%>%
    select(year, population , testing, best:ncol(.))
  
  
  Inci_pop[[n]] <- rbind(Bttpop[[n]], Sttpop[[n]])%>%
    mutate(year = (year - HCV$cabY + 1),
           testing = factor(testing, levels = c("Status Quo", "POC_antibody",
                                                "Reflex_RNA", "POC_RNA")))%>%
    arrange(year, population, testing)%>%
    select(year, population, testing, best:ncol(.))
}

# WHO elimination threshold 

thres <- list()
threspop <- list()

for(n in dt_name){ 
  
  threspop[[n]] <- Bttpop[[n]]%>%
    filter(year == 2015)%>%select(population, best)%>%mutate(thres = best*0.2)
  
  thres[[n]] <- Btt[[n]]%>%
    filter(year == 2015)%>%select(best)%>%mutate(thres = best*0.2)
  
  }

threspop_All <- list()

threspop_PrEP <- list()

threspop_HIVD <- list()

for(i in dt_name){ 
  
  
  threspop_All[[i]] <- rbind(thres[[i]], thres[[i]], 
                             thres[[i]], thres[[i]])%>%
    mutate(testing = unique(Inci[[i]]$testing))
  
  threspop_PrEP[[i]] <- rbind(threspop[[i]][threspop[[i]]$population =="HIV-" ,],
                              threspop[[i]][threspop[[i]]$population =="HIV-" ,], 
                              threspop[[i]][threspop[[i]]$population =="HIV-" ,], 
                              threspop[[i]][threspop[[i]]$population =="HIV-" ,])%>%
    mutate(testing = unique(Inci[[i]]$testing))
  
  threspop_HIVD[[i]] <- rbind(threspop[[i]][threspop[[i]]$population =="HIV+d" ,],
                              threspop[[i]][threspop[[i]]$population =="HIV+d" ,], 
                              threspop[[i]][threspop[[i]]$population =="HIV+d" ,], 
                              threspop[[i]][threspop[[i]]$population =="HIV+d" ,])%>%
    mutate(testing = unique(Inci[[i]]$testing))
  
  }




# (2) HCV incidence of testing strategies compared to base scenario by targeted MSM

labnam <- list(c("Testing scenarios"), c("Status Quo","Point-of-care antibody testing",
                                         "Clinic-based reflex RNA testing", 
                                         "Single visit point-of-care RNA testing"))

Fig_Inci_pop <- list()

gplotlist <- list()

for(i in dt_name){ 
  Fig_Inci_pop[[i]][["HIV-PrEP"]] <- ScePlot(dtt = Inci_pop[[n]]%>%filter(population == "HIV-PrEP"), 
                                             ylabel = "Incidence (per 1,000 PY)",
                                             ribbonarea = NULL,
                                             xlimits = c(2022, 2030, 1), 
                                             facetPlot = NULL,
                                             labelN = labnam) + 
   geom_hline(data = threspop_PrEP[[i]], 
                  aes(yintercept = thres), linetype = "dashed") + 
    geom_text(data = threspop_PrEP[[i]],
              aes(label = paste0("WHO elimination threshold"), 
                  y = threspop_PrEP[[i]]$thres, x = 2027), vjust = -1) + 
    theme_Publication() + ggtitle("MSM on PrEP") + 
    scale_y_continuous(limits = c(0, 5)) + 
    theme(legend.position = "")
    
    
  
  Fig_Inci_pop[[i]][["HIV+d"]] <- ScePlot(dtt = Inci_pop[[n]]%>%filter(population == "HIV+d"), 
                                             ylabel = "Incidence (per 1,000 PY)",
                                             ribbonarea = NULL,
                                             xlimits = c(2022, 2030, 1), 
                                             facetPlot = NULL,
                                             labelN = labnam) + 
    geom_hline(data = threspop_HIVD[[i]], 
               aes(yintercept = thres), linetype = "dashed") + 
    geom_text(data = threspop_HIVD[[i]],
              aes(label = paste0("WHO elimination threshold"), 
                  y = threspop_HIVD[[i]]$thres, x = 2027), vjust = -1) + 
    theme_Publication() + ggtitle("MSM diagnosed with HIV") + 
    scale_y_continuous(limits = c(0, 50)) + 
    theme(legend.position = "")
  # to put legend in the last plot 
  if(isTRUE(i == "PrEPnHIVD")){
    Fig_Inci_pop[[i]][["All"]] <- ScePlot(dtt = Inci[[i]], 
                                          ylabel = "Incidence (per 1,000 PY)",
                                          ribbonarea = NULL,
                                          xlimits = c(2022, 2030, 1), 
                                          facetPlot = NULL,
                                          labelN = labnam) + 
      geom_hline(data = thres[[i]], 
                 aes(yintercept = thres), linetype = "dashed") + 
      geom_text(data = threspop_HIVD[[i]],aes(label = paste0("WHO elimination threshold"), 
                                              y = 1, x = 2027), vjust = -1) + 
      theme_Publication() + ggtitle("Overall MSM") + 
      scale_y_continuous(limits = c(0, 20))
    
  }else{ 
    Fig_Inci_pop[[i]][["All"]] <- ScePlot(dtt = Inci[[i]], 
                                          ylabel = "Incidence (per 1,000 PY)",
                                          ribbonarea = NULL,
                                          xlimits = c(2022, 2030, 1), 
                                          facetPlot = NULL,
                                          labelN = labnam) + 
      geom_hline(data = thres[[i]], 
                 aes(yintercept = thres), linetype = "dashed") + 
      geom_text(data = threspop_HIVD[[i]],aes(label = paste0("WHO elimination threshold"), 
                                              y = 1, x = 2027), vjust = -1) + 
      theme_Publication() + ggtitle("Overall MSM") + 
      scale_y_continuous(limits = c(0, 20)) + 
      theme(legend.position = "")  
    
    }
  
    
}


Fig_Inci <- list()
for(i in names(Fig_Inci_pop)){ 
  
  Fig_Inci[[i]] <- list(Fig_Inci_pop[[i]]$All, 
                        Fig_Inci_pop[[i]]$`HIV+d`,
                        Fig_Inci_pop[[i]]$`HIV-PrEP`)
  }

taglab <- list("A", "B", "C", "D")
names(taglab) <- names(Sttpop)
for(i in names(Sttpop)){
  if(isTRUE(i == "PrEPnHIVD")){
    gplotlist[[i]] <- ggarrange(plotlist = Fig_Inci[[i]], 
                                labels = c(taglab[[i]]), nrow = 1, 
                                common.legend = TRUE, legend="bottom")
  }else{
    gplotlist[[i]] <- ggarrange(plotlist = Fig_Inci[[i]], 
                                labels = c(taglab[[i]]), nrow = 1, 
                                common.legend = FALSE)  
    
    }
  
}

# substr legend 

Flist <- ggarrange(gplotlist$PrEPsame, gplotlist$PrEP, gplotlist$HIVD, 
                   gplotlist$PrEPnHIVD, 
          common.legend = FALSE, nrow = length(names(taglab)),legend="bottom")

ggsave(path = outputfig, file="Inci_pop.pdf", Flist, height = 25, width = 16)

#### speperate primary infection and reinfection ####
Sttnumpop_range <- list()
Bttnumpop_range <- list()
Sttnumpop_reinf_range <- list()
Bttnumpop_reinf_range <- list()
NumInf_pop <- list()
NumInf_pop_range <- list()
NumreInf_pop <- list()
NumreInf_pop_range <- list()

for(n in dt_name){ 
  Sttnumpop_range[[n]] <- popRange(Sttnumpop[[n]], pop = NULL, test = NULL)
  Bttnumpop_range[[n]] <- popRange(Bttnumpop[[n]], pop = NULL, test = NULL)
  
  Sttnumpop_reinf_range[[n]] <- popRange(Sttnumpop_reinf[[n]], pop = NULL, test = NULL)
  Bttnumpop_reinf_range[[n]] <- popRange(Bttnumpop_reinf[[n]], pop = NULL, test = NULL)
  
  # number of infections
  ##combined status quo and scenarios 
  
  NumInf_pop[[n]] <- rbind(Sttnumpop[[n]], Bttnumpop[[n]])%>%
    mutate(testing = factor(testing, levels = c("Status Quo", "POC_antibody",
                                                "Reflex_RNA", "POC_RNA")))%>%
    arrange(year, population, testing)%>%
    select(year, population , testing, best:ncol(.))
  
  NumInf_pop_range[[n]] <- rbind(Sttnumpop_range[[n]], Bttnumpop_range[[n]])%>%
    mutate(testing = factor(testing, levels = c("Status Quo", "POC_antibody",
                                                "Reflex_RNA", "POC_RNA")))%>%
    arrange(year, population, testing)%>%
    select(year, population , testing, best:ncol(.))
  
  # number of reinfections
  ##combined status quo and scenarios 
  NumreInf_pop[[n]] <- rbind(Sttnumpop_reinf[[n]], Bttnumpop_reinf[[n]])%>%
    mutate(testing = factor(testing, levels = c("Status Quo", "POC_antibody",
                                                "Reflex_RNA", "POC_RNA")))%>%
    arrange(year, population, testing)%>%
    select(year, population , testing, best:ncol(.))
  
  NumreInf_pop_range[[n]] <- rbind(Sttnumpop_reinf_range[[n]], Bttnumpop_reinf_range[[n]])%>%
    mutate(testing = factor(testing, levels = c("Status Quo", "POC_antibody",
                                                "Reflex_RNA", "POC_RNA")))%>%
    arrange(year, population, testing)%>%
    select(year, population , testing, best:ncol(.))  
  
  

}

# get number of primary infection: total num of infection - num of reinfection 

NumpriInf_pop <- list()

NumpriInf_pop_range <- list()

for(n in dt_name){ 
  
  NumInf_pop[[n]] <- NumInf_pop[[n]]%>%arrange(year, population, testing)
  
  NumreInf_pop[[n]] <- NumreInf_pop[[n]]%>%arrange(year, population, testing)
  
  NumpriInf_pop[[n]] <- cbind(year = NumInf_pop[[n]]$year,
                              population = NumInf_pop[[n]]$population,
                              testing = NumInf_pop[[n]]$testing,
                              as.data.frame(NumInf_pop[[n]][ ,-c(1:3)] - 
                                              NumreInf_pop[[n]][, -c(1:3)])
                              )
  
  NumpriInf_pop_range[[n]] <- popRange(NumpriInf_pop[[n]], pop = NULL, 
                                       test = NULL)
}

# number of infection, reinfection and primary infection in overall population 
NumInf <- list()
NumreInf <- list()
NumpriInf <- list()

NumInf_range <- list()
NumreInf_range <- list()
NumpriInf_range <- list()

for(n in dt_name){ 
  
  NumInf[[n]] <- NumInf_pop[[n]]%>%group_by(year, testing)%>%
    summarize(across(c(best, c(paste0("set", seq(1, HCV$numberSamples,1)))), ~sum(.x)))%>%
    mutate(population = "All")%>%select(year, population, testing, best: ncol(.))
  
  NumreInf[[n]] <- NumreInf_pop[[n]]%>%group_by(year, testing)%>%
    summarize(across(c(best, c(paste0("set", seq(1, HCV$numberSamples,1)))), ~sum(.x)))%>%
    mutate(population = "All")%>%select(year, population, testing, best: ncol(.))
  
  NumpriInf[[n]] <- NumpriInf_pop[[n]]%>%group_by(year, testing)%>%
    summarize(across(c(best, c(paste0("set", seq(1, HCV$numberSamples,1)))), ~sum(.x)))%>%
    mutate(population = "All")%>%select(year, population, testing, best: ncol(.))
  
  
  NumInf_range[[n]] <- NumInf[[n]]%>%popRange(., pop = NULL, test = NULL)
  
  NumreInf_range[[n]] <- NumreInf[[n]]%>%popRange(., pop = NULL, test = NULL)
  
  NumpriInf_range[[n]] <- NumpriInf[[n]]%>%popRange(., pop = NULL, test = NULL)


}

#### plot numbers of infections, reinfection and primary infection 
# number of infections 

Fig_NumInf_pop <- list()

gplotlist <- list()

for(i in dt_name){ 
  Fig_NumInf_pop[[i]][["HIV-PrEP"]] <- ScePlot(dtt = NumInf_pop_range[[n]]%>%filter(population == "HIV-PrEP"), 
                                             ylabel = "Number of HCV infection",
                                             ribbonarea = NULL,
                                             xlimits = c(2022, 2030, 1), 
                                             facetPlot = NULL,
                                             labelN = labnam) + 
    theme_Publication() + ggtitle("MSM on PrEP") + 
    scale_y_continuous(limits = c(0, 300)) + 
    theme(legend.position = "")
  
  
  
  Fig_NumInf_pop[[i]][["HIV+d"]] <- ScePlot(dtt = NumInf_pop_range[[n]]%>%filter(population == "HIV+d"), 
                                          ylabel = "Number of HCV infection",
                                          ribbonarea = NULL,
                                          xlimits = c(2022, 2030, 1), 
                                          facetPlot = NULL,
                                          labelN = labnam) + 
    theme_Publication() + ggtitle("MSM diagnosed with HIV") + 
    scale_y_continuous(limits = c(0, 2000)) + 
    theme(legend.position = "")
  # to put legend in the last plot 
  if(isTRUE(i == "PrEPnHIVD")){
    Fig_NumInf_pop[[i]][["All"]] <- ScePlot(dtt = NumInf_range[[i]], 
                                          ylabel = "Number of HCV infection",
                                          ribbonarea = NULL,
                                          xlimits = c(2022, 2030, 1), 
                                          facetPlot = NULL,
                                          labelN = labnam) + 
      theme_Publication() + ggtitle("Overall MSM") + 
      scale_y_continuous(limits = c(0, 3000))
    
  }else{ 
    Fig_NumInf_pop[[i]][["All"]] <- ScePlot(dtt = NumInf_range[[i]], 
                                          ylabel = "Number of HCV infection",
                                          ribbonarea = NULL,
                                          xlimits = c(2022, 2030, 1), 
                                          facetPlot = NULL,
                                          labelN = labnam) + 
      theme_Publication() + ggtitle("Overall MSM") + 
      scale_y_continuous(limits = c(0, 3000)) + 
      theme(legend.position = "")  
    
  }
  
  
}


for(i in names(Fig_NumInf_pop)){ 
  
  Fig_NumInf[[i]] <- list(Fig_NumInf_pop[[i]]$All, 
                        Fig_NumInf_pop[[i]]$`HIV+d`,
                        Fig_NumInf_pop[[i]]$`HIV-PrEP`)
}

taglab <- list("A", "B", "C", "D")
names(taglab) <- names(Fig_NumInf_pop)
for(i in names(Fig_NumInf_pop)){
  if(isTRUE(i == "PrEPnHIVD")){
    gplotlist[[i]] <- ggarrange(plotlist = Fig_NumInf[[i]], 
                                labels = c(taglab[[i]]), nrow = 1, 
                                common.legend = TRUE, legend="bottom")
  }else{
    gplotlist[[i]] <- ggarrange(plotlist = Fig_NumInf[[i]], 
                                labels = c(taglab[[i]]), nrow = 1, 
                                common.legend = FALSE)  
    
  }
  
}

# substr legend 

Flist <- ggarrange(gplotlist$PrEPsame, gplotlist$PrEP, gplotlist$HIVD, 
                   gplotlist$PrEPnHIVD, 
                   common.legend = FALSE, nrow = length(names(taglab)),legend="bottom")

ggsave(path = outputfig, file="NumInf_pop.pdf", Flist, height = 25, width = 16)

##### number of reinfection ####

#### plot numbers of infections, reinfection and primary infection 
Fig_NumreInf_pop <- list()

gplotlist <- list()

for(i in dt_name){ 
  Fig_NumreInf_pop[[i]][["HIV-PrEP"]] <- ScePlot(dtt = NumreInf_pop_range[[n]]%>%filter(population == "HIV-PrEP"), 
                                               ylabel = "Number of HCV reinfection",
                                               ribbonarea = NULL,
                                               xlimits = c(2022, 2030, 1), 
                                               facetPlot = NULL,
                                               labelN = labnam) + 
    theme_Publication() + ggtitle("MSM on PrEP") + 
    scale_y_continuous(limits = c(0, 100)) + 
    theme(legend.position = "")
  
  
  
  Fig_NumreInf_pop[[i]][["HIV+d"]] <- ScePlot(dtt = NumreInf_pop_range[[n]]%>%filter(population == "HIV+d"), 
                                            ylabel = "Number of HCV reinfection",
                                            ribbonarea = NULL,
                                            xlimits = c(2022, 2030, 1), 
                                            facetPlot = NULL,
                                            labelN = labnam) + 
    theme_Publication() + ggtitle("MSM diagnosed with HIV") + 
    scale_y_continuous(limits = c(0, 1000)) + 
    theme(legend.position = "")
  # to put legend in the last plot 
  if(isTRUE(i == "PrEPnHIVD")){
    Fig_NumreInf_pop[[i]][["All"]] <- ScePlot(dtt = NumreInf_range[[i]], 
                                            ylabel = "Number of HCV reinfection",
                                            ribbonarea = NULL,
                                            xlimits = c(2022, 2030, 1), 
                                            facetPlot = NULL,
                                            labelN = labnam) + 
      theme_Publication() + ggtitle("Overall MSM") + 
      scale_y_continuous(limits = c(0, 1500))
    
  }else{ 
    Fig_NumreInf_pop[[i]][["All"]] <- ScePlot(dtt = NumreInf_range[[i]], 
                                            ylabel = "Number of HCV reinfection",
                                            ribbonarea = NULL,
                                            xlimits = c(2022, 2030, 1), 
                                            facetPlot = NULL,
                                            labelN = labnam) + 
      theme_Publication() + ggtitle("Overall MSM") + 
      scale_y_continuous(limits = c(0, 1500)) + 
      theme(legend.position = "")  
    
  }
  
  
}


Fig_NumreInf <- list()

for(i in names(Fig_NumreInf_pop)){ 
  
  Fig_NumreInf[[i]] <- list(Fig_NumreInf_pop[[i]]$All, 
                          Fig_NumreInf_pop[[i]]$`HIV+d`,
                          Fig_NumreInf_pop[[i]]$`HIV-PrEP`)
}

taglab <- list("A", "B", "C", "D")
names(taglab) <- names(Fig_NumreInf_pop)
for(i in names(Fig_NumreInf_pop)){
  if(isTRUE(i == "PrEPnHIVD")){
    gplotlist[[i]] <- ggarrange(plotlist = Fig_NumreInf[[i]], 
                                labels = c(taglab[[i]]), nrow = 1, 
                                common.legend = TRUE, legend="bottom")
  }else{
    gplotlist[[i]] <- ggarrange(plotlist = Fig_NumreInf[[i]], 
                                labels = c(taglab[[i]]), nrow = 1, 
                                common.legend = FALSE)  
    
  }
  
}

# substr legend 

Flist <- ggarrange(gplotlist$PrEPsame, gplotlist$PrEP, gplotlist$HIVD, 
                   gplotlist$PrEPnHIVD, 
                   common.legend = FALSE, nrow = length(names(taglab)),legend="bottom")

ggsave(path = outputfig, file="NumreInf_pop.pdf", Flist, height = 25, width = 16)



##### number of primary infection  ####

#### plot numbers of infections, reinfection and primary infection 
Fig_NumpriInf_pop <- list()

gplotlist <- list()

for(i in dt_name){ 
  Fig_NumpriInf_pop[[i]][["HIV-PrEP"]] <- ScePlot(dtt = NumpriInf_pop_range[[n]]%>%filter(population == "HIV-PrEP"), 
                                                 ylabel = "Number of HCV primaryinfection",
                                                 ribbonarea = NULL,
                                                 xlimits = c(2022, 2030, 1), 
                                                 facetPlot = NULL,
                                                 labelN = labnam) + 
    theme_Publication() + ggtitle("MSM on PrEP") + 
    scale_y_continuous(limits = c(0, 100)) + 
    theme(legend.position = "")
  
  
  
  Fig_NumpriInf_pop[[i]][["HIV+d"]] <- ScePlot(dtt = NumpriInf_pop_range[[n]]%>%filter(population == "HIV+d"), 
                                              ylabel = "Number of HCV primary infection",
                                              ribbonarea = NULL,
                                              xlimits = c(2022, 2030, 1), 
                                              facetPlot = NULL,
                                              labelN = labnam) + 
    theme_Publication() + ggtitle("MSM diagnosed with HIV") + 
    scale_y_continuous(limits = c(0, 1000)) + 
    theme(legend.position = "")
  # to put legend in the last plot 
  if(isTRUE(i == "PrEPnHIVD")){
    Fig_NumpriInf_pop[[i]][["All"]] <- ScePlot(dtt = NumpriInf_range[[i]], 
                                              ylabel = "Number of HCV primary infection",
                                              ribbonarea = NULL,
                                              xlimits = c(2022, 2030, 1), 
                                              facetPlot = NULL,
                                              labelN = labnam) + 
      theme_Publication() + ggtitle("Overall MSM") + 
      scale_y_continuous(limits = c(0, 1500))
    
  }else{ 
    Fig_NumpriInf_pop[[i]][["All"]] <- ScePlot(dtt = NumpriInf_range[[i]], 
                                              ylabel = "Number of HCV primaryinfection",
                                              ribbonarea = NULL,
                                              xlimits = c(2022, 2030, 1), 
                                              facetPlot = NULL,
                                              labelN = labnam) + 
      theme_Publication() + ggtitle("Overall MSM") + 
      scale_y_continuous(limits = c(0, 1500)) + 
      theme(legend.position = "")  
    
  }
  
  
}


Fig_NumpriInf <- list()

for(i in names(Fig_NumpriInf_pop)){ 
  
  Fig_NumpriInf[[i]] <- list(Fig_NumpriInf_pop[[i]]$All, 
                            Fig_NumpriInf_pop[[i]]$`HIV+d`,
                            Fig_NumpriInf_pop[[i]]$`HIV-PrEP`)
}

taglab <- list("A", "B", "C", "D")
names(taglab) <- names(Fig_NumpriInf_pop)
for(i in names(Fig_NumpriInf_pop)){
  if(isTRUE(i == "PrEPnHIVD")){
    gplotlist[[i]] <- ggarrange(plotlist = Fig_NumpriInf[[i]], 
                                labels = c(taglab[[i]]), nrow = 1, 
                                common.legend = TRUE, legend="bottom")
  }else{
    gplotlist[[i]] <- ggarrange(plotlist = Fig_NumpriInf[[i]], 
                                labels = c(taglab[[i]]), nrow = 1, 
                                common.legend = FALSE)  
    
  }
  
}

# substr legend 

Flist <- ggarrange(gplotlist$PrEPsame, gplotlist$PrEP, gplotlist$HIVD, 
                   gplotlist$PrEPnHIVD, 
                   common.legend = FALSE, nrow = length(names(taglab)),legend="bottom")

ggsave(path = outputfig, file="NumpriInf_pop.pdf", Flist, height = 25, width = 16)


#### desired cost outcomes shown in figs and tables ####
# (1) cost along the HCV care: 
  # (1) cost of HCV diagnosis: consisting of cost_ab + cost_RNA
  # (2) cost of HCV treatment: cosisting of cost_treat + cost_Pops[["Treat"]]
  # (3) cost of HCV management and others co-mordibility: cosisting of cost_Pops[[!= c("Treatment", "s")]]
# structure of Scen_cost_flow list: Scen_cost_flow[["PrEP or HIVD scenario"]][["cost along cascade"]][["testing scenarios"]]
# Structure of Scen_cost_box list: Scen_cost_box[["PrEP or HIVD scenario"]][["cost/QALY of components"]][["testing scenarios"]]
Base_cost_flow <- list()
Base_cost_box <- list()
Scen_cost_flow <- list()
Scen_cost_box <- list()

for(i in dt_name){ 
  
  Base_cost_flow[[i]] <- Outcome_Base_cost[[i]]$cost_flow
  
  Base_cost_box[[i]] <- Outcome_Base_cost[[i]]$cost_box
  
  
  for(n in names(Outcome_Scen_cost$PrEPsame)){ 
   
    Scen_cost_box[[i]][[n]] <- Outcome_Scen_cost[[i]][[n]]$cost_box
    
    Scen_cost_flow[[i]][[n]] <- Outcome_Scen_cost[[i]][[n]]$cost_flow
    
  }
  
}

# colnames 
colset <- Scen_cost_flow[[1]][[1]][[1]]%>%select(., contains("set"))%>%
  select(., str_sort(paste0("set", seq(1, HCV$numberSamples,1)),numeric = FALSE))

PrEPS <- names(Scen_cost_flow)

testingS <- names(Scen_cost_flow[[1]])

cost_HCVtreatPops_Base <- list()

cost_HCV_Base <- list()

cost_Treat_noagr <- list()

cost_Retreat_noagr <- list()

cost_HCVtreatPops_Scen <- list()

cost_HCV_Scen <- list()

Scen_cost_Treat_noagr <- list()

Scen_cost_Retreat_noagr <- list()

# stage cost of HCV test and treat in Status quo

for(i in PrEPS){ 
  
  cost_HCVtreatPops_Base[[i]] <- Base_cost_box[[i]][["costPops"]]%>%
    filter(cascade == "treat")%>%group_by(year, population)%>%
    summarize(across(c(best, colnames(colset)), ~sum(.x)))%>%
    select(year,population, best, colnames(colset))%>%as.data.frame()

  cost_Treat_noagr[[i]] <- Base_cost_flow[[i]][["costTreatment"]]%>%
    select(year, population, best, colnames(colset))%>%as.data.frame()
  
  cost_Retreat_noagr[[i]] <- Base_cost_flow[[i]][["costRetreat"]]%>%
    select(year, population, best, colnames(colset))%>%as.data.frame()
  
  cost_HCV_Base[[i]][["diag"]] <- 
    cbind(year = Base_cost_flow[[i]][["costTestingAb"]]$year, 
          population = Base_cost_flow[[i]][["costTestingAb"]]$population,
          as.data.frame(Base_cost_flow[[i]][["costTestingAb"]][, -c(1:2)] +  
                          Base_cost_flow[[i]][["costTestingAg"]][, -c(1:2)] + 
                          Base_cost_flow[[i]][["costnewTestingPOCT"]][, -c(1:2)] 
                                  ))%>%
    select(year, population, best, colnames(colset))
  
  
  cost_HCV_Base[[i]][["treat"]] <- 
    cbind(year = Base_cost_flow[[i]][["costTreatment"]]$year, 
          population = Base_cost_flow[[i]][["costTreatment"]]$population, 
          as.data.frame(cost_Treat_noagr[[i]][, -c(1:2)] + 
                          cost_HCVtreatPops_Base[[i]][ , -c(1:2)] + 
                          cost_Retreat_noagr[[i]][, -c(1:2)]
          ))%>%
    select(year, population, best, colnames(colset))
  
  
  cost_HCV_Base[[i]][["management"]] <- Base_cost_box[[i]][["costPops"]]%>%
    filter(!cascade %in% c("treat", "s"))%>%
    filter(!state %in% c("a_cured", "f0_cured", "f1_cured", "f2_cured"))%>%
    group_by(year, population)%>%
    summarize(across(c(best, colnames(colset)), ~sum(.x)))%>%
    select(year, population, best, colnames(colset))
}
  


# stage cost of HCV test and treat in scenarios 
# saving like xxx[[PrEP scenario]][[HCV testing scenarios]][[stage of diag/treatment/managemenmt]]

for(i in PrEPS){
  
  for(t in testingS){ 
    
    cost_HCVtreatPops_Scen[[i]][[t]] <- Scen_cost_box[[i]][[t]][["costPops"]]%>%
      filter(cascade == "treat")%>%group_by(year, population)%>%
      summarize(across(c(best, colnames(colset)), ~sum(.x)))%>%
      select(year,population, best, colnames(colset))%>%as.data.frame()
    
    Scen_cost_Treat_noagr[[i]][[t]] <- Scen_cost_flow[[i]][[t]][["costTreatment"]]%>%
      select(year, population, best, colnames(colset))%>%as.data.frame()
    
    Scen_cost_Retreat_noagr[[i]][[t]] <- Scen_cost_flow[[i]][[t]][["costRetreat"]]%>%
      select(year, population, best, colnames(colset))%>%as.data.frame()
    
    cost_HCV_Scen[[i]][[t]][["diag"]] <- 
      cbind(year = Scen_cost_flow[[i]][[t]][["costTestingAb"]]$year, 
            population = Scen_cost_flow[[i]][[t]][["costTestingAb"]]$population,
            as.data.frame(Scen_cost_flow[[i]][[t]][["costTestingAb"]][, -c(1:2)] + 
                            Scen_cost_flow[[i]][[t]][["costTestingAg"]][, -c(1:2)] + 
                            Scen_cost_flow[[i]][[t]][["costnewTestingPOCT"]][, -c(1:2)] 
            ))%>%
      select(year, population, best, colnames(colset))
    
    cost_HCV_Scen[[i]][[t]][["treat"]] <- 
      cbind(year = Scen_cost_flow[[i]][[t]][["costTreatment"]]$year, 
            population = Scen_cost_flow[[i]][[t]][["costTreatment"]]$population,
            as.data.frame(Scen_cost_Treat_noagr[[i]][[t]][, -c(1:2)] + 
                            cost_HCVtreatPops_Scen[[i]][[t]][ , -c(1:2)] + 
                            Scen_cost_Retreat_noagr[[i]][[t]][, -c(1:2)]  
            ))%>%
      select(year, population, best, colnames(colset))
    
    cost_HCV_Scen[[i]][[t]][["management"]] <- 
      Scen_cost_box[[i]][[t]][["costPops"]]%>%
      filter(!cascade %in% c("treat", "s"))%>%
      filter(!state %in% c("a_cured", "f0_cured", "f1_cured", "f2_cured"))%>%
      group_by(year, population)%>%
      summarize(across(c(best, colnames(colset)), ~sum(.x)))%>%
      select(year, population, best, colnames(colset))
    
    }
  }

# combind status quo cost into list 
Carecas <- names(cost_HCV_Scen[[1]][[1]])

for(i in dt_name){ 
  
  cost_HCVtreatPops_Scen[[i]][["Status Quo"]] <- cost_HCVtreatPops_Base[[i]]
  
  for(m in Carecas){ 
    
    cost_HCV_Scen[[i]][["Status Quo"]][[m]] <- cost_HCV_Base[[i]][[m]]
    
    }
  }

# create dt for the total cost plot 
# undiscounted, in NTD
cost_stage <- list()
for(i in dt_name){ 
  for(t in c(testingS, "Status Quo")){ 
    for(m in Carecas){ 
      
      cost_stage[[i]][[t]][[m]] <- cost_HCV_Scen[[i]][[t]][[m]]%>%
        mutate(testing = t, 
               Cascade = m)
    }
    cost_stage[[i]][[t]] <- do.call("rbind", cost_stage[[i]][[t]])
    
  }
  cost_stage[[i]] <- do.call("rbind", cost_stage[[i]])%>%
    select(year, population, testing, Cascade, best, 
           paste0("set", seq(1,HCV$numberSamples,1)))%>%
    popRange(.)
}


for(i in dt_name){ 
  cost_stage[[i]] <- cost_stage[[i]]%>%group_by(population, year,testing)%>%
    summarize(across(c(best, colnames(colset)), ~sum(.x)))%>%
    arrange(year, population, testing)
}

# total cost (accumulating from 2022 to 2090) 
HCVcost <- list()

HCVcost_dis <- list()

HCVcostCum <- list()

cost_HCV_Scent <- list()

cost_HCV_Scent_dis <- list() 

cost_HCV_Scent_dis_cum <- list()

HCVcostCum_dis <- list()

# cost still calculated in NTD 
for(i in dt_name){ 
  for(t in names(cost_HCV_Scen[[1]])){ 
   for(m in Carecas){ 
      cost_HCV_Scent[[i]][[t]][[m]] <- cost_HCV_Scen[[i]][[t]][[m]]%>%
        arrange(year, population)%>%as.data.frame()%>%
        mutate(Casca = m)%>%
        mutate(year = year + HCV$cabY - 1)%>%
        filter(year >= 2022)%>%
        mutate(id = year - 2022)%>%
        mutate(discount = (1 + discountR)^id)
      
      cost_HCV_Scent_dis[[i]][[t]][[m]] <-  cost_HCV_Scent[[i]][[t]][[m]]%>%
        mutate(across(c(best, colnames(colset)), function(x) x/discount))%>%
        arrange(year, population, Casca)
      
    }
    # combine stage cost 
    HCVcost[[i]][[t]] <- do.call("rbind", cost_HCV_Scent[[i]][[t]])%>%
      as_tibble()%>%mutate(testing = t)%>%
      select(year, population, testing, Casca, best,colnames(colset))
    
    # cumulative stage cost 
    HCVcostCum[[i]][[t]] <- HCVcost[[i]][[t]]%>%
      arrange(year, population, testing, Casca)%>%
      group_by(population, testing, Casca)%>%
      mutate(across(c(best, colnames(colset)), ~cumsum(.x))) 
    
    # discount cost 
    HCVcost_dis[[i]][[t]] <- do.call("rbind", cost_HCV_Scent_dis[[i]][[t]])%>%
      as_tibble()%>%mutate(testing = t)%>%
      select(year, population, testing, Casca, best,colnames(colset))
    # cumulative discount cost 
    HCVcostCum_dis[[i]][[t]] <- HCVcost_dis[[i]][[t]]%>%
      arrange(year, population, testing, Casca)%>%
      group_by(population, testing, Casca)%>%
      mutate(across(c(best, colnames(colset)), ~cumsum(.x)))%>%
      select(year, population, testing, Casca, best,colnames(colset))
  }
}

# plot the annual discounted cost 
HCVcost_dis_all <- list()
HCVcost_dis_year <- list()

for(i in PrEPS){ 
  
  HCVcost_dis[[i]] <- do.call("rbind", HCVcost_dis[[i]])%>%
    select(year, population, testing, Casca, best, 
           paste0("set", seq(1, HCV$numberSamples, 1)))
}
 
for(i in PrEPS){  
  HCVcost_dis_all[[i]] <-  HCVcost_dis[[i]]%>%
    group_by(year, testing, Casca)%>%
    summarize(across(c(best, colnames(colset)), ~sum(.x)))%>%
    mutate(population = "All")%>%
    select(year, population, testing, Casca, 
           best, paste0("set", seq(1, HCV$numberSamples, 1)))
  
}

for(i in PrEPS){
  HCVcost_dis_year[[i]] <- rbind(HCVcost_dis[[i]], HCVcost_dis_all[[i]])
    
  }






# combind testing into one dataframe   
HCVcost_ScenCum <- list()

HCVcost_ScenCum_dis <- list()

for(i in names(HCVcost)){ 
  
  HCVcost_ScenCum[[i]] <- do.call("rbind", HCVcostCum[[i]])
  
  HCVcost_ScenCum[[i]]$testing <- factor(HCVcost_ScenCum[[i]]$testing, 
                                         levels = c("Status Quo", "POC_antibody", 
                                                    "Reflex_RNA", "POC_RNA"))
  
  HCVcost_ScenCum_dis[[i]] <- do.call("rbind", HCVcostCum_dis[[i]])
  
  HCVcost_ScenCum_dis[[i]]$testing <- factor(HCVcost_ScenCum_dis[[i]]$testing, 
                                         levels = c("Status Quo", "POC_antibody", 
                                                    "Reflex_RNA", "POC_RNA"))
  
  
}

HCVcostCum_dis_all <- list()

HCVcostCum_dis_all_range <- list()

for(i in dt_name){ 
  HCVcostCum_dis_all[[i]] <-  HCVcost_ScenCum_dis[[i]]%>%
    group_by(year, testing, Casca)%>%
    summarize(across(c(best, colnames(colset)), ~sum(.x)))
  }
for(i in dt_name){ 
  HCVcostCum_dis_all_range[[i]] <-popRange(HCVcostCum_dis_all[[i]], 
                                           target_year = NULL, pop = "all", 
                                           test = NULL, Casca = "y")
}



#### QALY #### 
# sum up QALY in each year by pops 
QALYpop <- list()

HCVQALYCum <- list()

for(i in dt_name){
  
  Scen_cost_box[[i]][["Status Quo"]][["QALYPops"]] <- Base_cost_box[[i]][["QALYPops"]]
  
  for(t in names(cost_HCV_Scen[[1]])){
    
    QALYpop[[i]][[t]] <- Scen_cost_box[[i]][[t]][["QALYPops"]]%>%
      group_by(year, population)%>%
      summarize(across(c(best, colnames(colset)), ~sum(.x)))%>%
      mutate(year = HCV$cabY + year - 1)%>%
      filter(year >= 2022)%>%
      mutate(id = year - 2022)%>%
      mutate(discount = (1 + discountR)^id,
             testing = t)
    
  }
}
  for(i in dt_name){    
    # discounted
    HCVQALYCum[[i]] <-  do.call("rbind", QALYpop[[i]])%>%
      as_tibble()%>%arrange(year, population, testing)%>%
      mutate(across(c(best, colnames(colset)), function(x) x/discount))%>%
      arrange(year, population)%>%
      group_by(population, testing)%>%
      mutate(across(c(best, colnames(colset)), ~cumsum(.x)))%>%
      select(year, population, testing, best,colnames(colset))%>%
      arrange(year, population, testing)
    }
   
#### QALY gained plot #### 
QALY_year <- list()
for(i in dt_name){ 
  
  QALY_year[[i]] <- do.call("rbind", QALYpop[[i]])%>%
    as_tibble()%>%arrange(year, population, testing)%>%
    mutate(across(c(best, colnames(colset)), function(x) x/discount))%>%
    arrange(year, population)%>%
    group_by(population, testing)%>%
    select(year, population, testing, best,colnames(colset))%>%
    arrange(year, population, testing)
  
  }


save(Inci, Inci_pop, NumInf_pop, NumInf_pop_range,
     NumreInf_pop, NumreInf_pop_range, NumpriInf_pop, NumpriInf_pop_range,
     NumInf, NumInf_range,
     NumreInf, NumreInf_range, NumpriInf, NumpriInf_range,
     HCVcost, HCVcostCum, cost_HCV_Scent, 
     cost_HCV_Scent_dis, HCVcostCum_dis, 
     HCVcost_ScenCum_dis, 
     QALYpop, HCVQALYCum, 
     file = file.path(outputdt, paste0("measureOutcome.rda")))

