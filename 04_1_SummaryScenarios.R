# this script used for wrangling the dataset of scenarios results 
# to out put as prevalence, incidence and the results reported in the manuscript 
rm(list = ls())
library(tidyverse)
library(data.table)
library(formattable)
memory.limit(size = 5000000)
basePath <- getwd()
Rcode <- file.path(basePath, "03. Code")
DataFolder <- file.path(basePath, "01. DATA/model input" )
OutputFolder <- file.path(basePath, "04. Output" )
projectFolder <- file.path(basePath)
FigFolder <- file.path(basePath, "04. Output/Result_figure" )
ResultsFolder <- file.path(basePath, "04. Output/Results")

source(file.path(Rcode, "/Functions/HCV_model.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R"))
source(file.path(Rcode, "/Functions/plotManuscript.R"))


project_name <- "HCVModel"
load(file.path(ResultsFolder, paste0("results_", "2022-10-28_23-00-37"), 
               paste0("results_", "2022-10-28_23-00-37", ".rda")))
projectFile <- file.path(basePath,paste0(project_name, ".rda"))
projectVars <- load(projectFile)
load(file.path(projectFolder, paste0(project_name, "ScenResults",".rda")))


scenarioResults$main <- scenarioResults$Main
scenarioResults <- within(scenarioResults, rm(Main))
sceSenName <- names(scenarioResults) 
HCV$sParam <- 1000

scenarioResultsPar <- list()
scenarioResultsPar[["main"]] <- 
  get(load(file.path(projectFolder, 
                     paste0(project_name,"ScenResultsPar_main",".rda"))))
# outcome plot (prev, incidence, annunal incidence, people living with CHCV)
# status quo
bR <- list()
bPar <- list()
bR[["Status quo"]] <- list(bestResults)
bPar[["Status quo"]] <- list(paramResults)

# Prev
bttP <- scenario_Prevalence(bR[["Status quo"]],bPar[["Status quo"]],
                            pop = NULL, statusQ = "y")

ttP <- scenario_Prevalence(scenarioResults[["main"]], 
                           scenarioResultsPar[["main"]], 
                           pop = NULL)

dtp <- rbind(bttP, ttP)%>%mutate(Med = round(Med, 2))


# inci
tt <- scenario_Incidence(scenarioResults[["main"]], 
                         scenarioResultsPar[["main"]], 
                         pop = NULL)

btt <- scenario_Incidence(bR[["Status quo"]],bPar[["Status quo"]],
                            pop = NULL, statusQ = "y")

thres <- btt%>%mutate(year = 2003+ year)%>%filter(year ==2015)%>%
  select(best)%>%mutate(thres = best*0.2)

dtt <- rbind(btt, tt)


ttpop <- scenario_Incidence(scenarioResults[["main"]], 
                         scenarioResultsPar[["main"]], 
                         pop = "pop")

bttpop <- scenario_Incidence(bR[["Status quo"]],bPar[["Status quo"]],
                          pop = "pop", statusQ = "y")
threspop <- bttpop%>%mutate(year = 2003+ year)%>%filter(year ==2015)%>%
  select(best)%>%mutate(thres = best*0.2)
threspop[2,] <- threspop[1,]

dttpop <- rbind(bttpop, ttpop)

#number of people living with chronic HCV 
schc <- list()
schc <- scenario_CHCVnum(scenarioResults[["main"]], 
                         scenarioResultsPar[["main"]], 
                         pop = NULL, statusQ = NULL)

bchc <- scenario_CHCVnum(bR[["Status quo"]],bPar[["Status quo"]],
                           pop = NULL, statusQ = "y")

dchc <- rbind(bchc, schc)

# number of HCV new infections 
ttnum <- scenario_Incidence(scenarioResults[["main"]], 
                            scenarioResultsPar[["main"]], 
                            pop = NULL, indicator = "newInfections")

sttnum <- scenario_Incidence(bR[["Status quo"]],bPar[["Status quo"]],
                           pop = NULL,indicator = "newInfections", statusQ = "y")


dtnum <- rbind(sttnum , ttnum)

Prev_main <-dtp
Inci_main <-dtt
CHCV_main <- dchc
NumInf <- dtnum
save(Prev_main,Inci_main,CHCV_main,NumInf,
     file = file.path(projectFolder,
                      paste0(project_name, "ScePlot_main", ".rda")))
scenarioResultsPar <- within(scenarioResultsPar, rm(scenarioResultsPar[["main"]]))


a <- ScePlot(dtt = Prev_main, ylabel = "Prevalence",
             ribbonarea = NULL,
             xlimits = c(2015, 2030, 5), facetPlot = NULL,
             labelN = NULL) + scale_y_continuous(expand = c(0,0), 
                                                 limits = c(0,6), 
                                                 breaks = seq(0, 6,1)) + theme_Publication()
a
b <- ScePlot(dtt = Inci_main, ylabel = "Incidence",
             ribbonarea = NULL,
             xlimits = c(2015, 2030, 5), facetPlot = NULL,
             labelN = NULL) + scale_y_continuous(expand = c(0,0), limits = c(0,20), 
                                                 breaks = seq(0, 20,5)) + theme_Publication()


c <- ScePlot(dtt = CHCV_main, ylabel = "Living with chronic HCV",
             ribbonarea = NULL,
             xlimits = c(2015, 2030, 5), facetPlot = NULL,
             labelN = NULL) + scale_y_continuous(expand = c(0,0), 
                                                 limits = c(0,10000), 
                                                 breaks = seq(0, 10000,1000)) + theme_Publication()

c
d <- ScePlot(dtt = NumInf, ylabel = "Annual Incidence",
             ribbonarea = NULL,
             xlimits = c(2015, 2030, 5), facetPlot = NULL,
             labelN = NULL) + scale_y_continuous(expand = c(0,0), 
                                                 limits = c(0,3000), 
                                                 breaks = seq(0, 3000,500)) + theme_Publication()

a <- a + theme(legend.position = "none") + labs(tag = "A")
b <- b + theme(legend.position = "none") + labs(tag = "B")
c <- c + theme(legend.position = "none") + labs(tag = "C")
d <- d + theme(legend.position = "none") + labs(tag = "D")
fig2 <- ggarrange(a, b, c,d, common.legend = TRUE, legend="top") + 
  theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm")) 

ggsave(file="mainFig2.png", fig2, width = 16, height = 9)


tic <- proc.time()
#### outcome table and datset for cascade plot ####
SceRestbl <- list()

for(i in sceSenName){ 
  scenarioResultsPar <- list()
  
  scenarioResultsPar[[i]] <- 
    get(load(file.path(projectFolder, 
                       paste0(project_name,"ScenResultsPar_", 
                              i,".rda"))))
  
  SceRestbl[["tbl"]][[i]] <- Restbl(scenarioResults[[i]], 
                                                scenarioResultsPar[[i]], 
                                                pop = NULL)%>%
    mutate(SensiSce = i)
  
  rm(list=paste("scenarioResultsPar_",i,sep  = ""))
} 


for(i in sceSenName){ 
  scenarioResultsPar <- list()
  scenarioResultsPar[[i]] <- get(load(file.path(projectFolder, 
                                                paste0(project_name,
                                                       "ScenResultsPar_", 
                                                       i,".rda"))))
  SceRestbl[["cas"]][[i]] <- 
    lapply(seq_along(scenarioResults[[i]]), function(x){ 
      a <- HCVCas(HCV,scenarioResults[[i]][[x]], 
                  scenarioResultsPar[[i]][[x]], 
                  pop = NULL, scenario ="y", initY = 2017)})
  
  rm(list=paste("scenarioResultsPar_",i,sep  = ""))
} 






# function to update object in rda file   
resave <- function(..., list = character(), file) {
  previous  <- load(file)
  var.names <- c(list, as.character(substitute(list(...)))[-1L])
  for (var in var.names) assign(var, get(var, envir = parent.frame()))
  save(list = unique(c(previous, var.names)), file = file)
} 

save(SceRestbl,
     file = file.path(projectFolder,
                      paste0(project_name, "SceRestbl", ".rda")))

SqRestbl <- list()


SqRestbl[["tbl"]] <- Restbl(bR[[1]],
                            bPar[[1]], pop = NULL, statusQ = "y")

SqRestbl[["cas"]] <- HCVCas(HCV,bestResults, 
                            paramResults, 
                            pop = NULL, scenario ="y", initY = 2017)

save(SqRestbl,
     file = file.path(projectFolder,
                      paste0(project_name, "SqRestbl", ".rda")))




