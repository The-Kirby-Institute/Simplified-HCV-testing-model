#### sensitivity analysis #### 

# snapshot of sensitivity analysis 


# before 2010 (capture beta changing) 

# 2017 Epi outcomes 

# 2022 Epi outcomes 

# 2030 Epi outcomes 

# Lifetime QALY and cost 


#=============================== Example =======================================
# x1 <- rnorm(n = 10, mean = 120, sd = 130)
# x2 <- rnorm(n = 10, mean = 80, sd = 5)
# x3 <- rnorm(n = 10, mean = 40, sd = 20)
# y <- 2 + (0.5 * x1) - (1.7 * x2) + (0.2 * x3)

# dat.df01 <- data.frame(x1 = x1, x2 = x2, x3 = x3, y = y)
# epi.prcc(dat.df01, sided.test = 2, conf.level = 0.95)
#===============================================================================


library(here)
here()

# load libraries
library(dplyr)
library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)
library(doParallel)
library(sensitivity)
library(purrr)
library(ggpubr)
# we specify the number of cores/workers we want to use
registerDoParallel(cores = detectCores() - 1)

codep <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/TWHCV-model"

com_codeP <- "/Users/jjwu/Documents/Simplified-HCV-testing-model/03. Code"

epidatapath <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Taiwan-MSM-HCV-model"

project_codep <- "/Users/jjwu/Documents/Simplified-HCV-testing-model/Projects/02.CEA_TWHCVMSM"

#file path of "TWHCV-model" project
Rcode <- file.path(codep, '03. Code')

scenariopath <- file.path(epidatapath, '01. DATA/model input/Scenarios')

# save and cost data folder 
dtp <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/02.CEA_TWHCVMSM"

DataFolder <- file.path(dtp, "01. DATA/model input")

outputdt <- file.path(dtp, "02. Output/RDA")

outputfig <- file.path(dtp, "02. Output/Figs")

# source 

source(file.path(com_codeP, "Functions/plotFunctions.R"))

source(file.path(com_codeP, "Functions/plotOptions.R"))

source(file.path(com_codeP, "Functions/Scenarios_wrangling.R"))

source(file.path(com_codeP, "Functions/plotManuscript.R"))

source(file.path(project_codep, "AggregateRes.R"))

projectFile <- file.path(epidatapath , paste0("HCVModel",".rda"))

load(file.path(outputdt , "cost.rda"))

load(file.path(outputdt , "toft.rda"))

load(file.path(outputdt , "Param_dfExt.rda"))

load(file.path(outputdt , "mainS_PrEPHIV.rda"))

load(file.path(outputdt , "pop_rate.rda"))

load(file.path(outputdt, "measureOutcome.rda"))

load(file.path(rdapath , paste0("HCVModel", "cali",".rda")))

load(file.path(rdapath , paste0("HCVModel", "param",".rda")))

source(file.path(here(), "AggregateRes.R"))
source(file.path(common_codepath, "plotOptions.R"))
# load HCV incidence dt 

files <- list.files(path = paste0(outputdt, "/", sep = ""),
                    pattern = '^Outcome_Base_.*\\.rda')
for( f in files){ 
  
  load(file.path(paste0(outputdt, "/" ,f, sep ="")))
  
}


Outcome_Base <- list() 

Outcome_Base[["Incidence"]] <- list(PrEPsame = Outcome_Base_PrEPsame$epi$Inci, 
                                   PrEP = Outcome_Base_PrEP$epi$Inci,
                                   HIVD = Outcome_Base_HIVD$epi$Inci,
                                   PrEPnHIVD = Outcome_Base_PrEPnHIVD$epi$Inci)


cost_y <- list()
qaly_y <- list()


for(i in names(HCVcost_ScenCum_dis)){
  
  cost_y[[i]] <- HCVcost_ScenCum_dis[[i]]%>%ungroup()%>%
    filter(testing == "Status Quo" & year == 2090)%>%
    summarize(across(c(best, c(paste0("set", seq(1, 1000,1)))), ~sum(.x)))
              
              qaly_y[[i]] <- HCVQALYCum[[i]]%>%
                filter( testing == "Status Quo" & year == 2090)%>%ungroup()%>%
                summarize(across(c(best, c(paste0("set", seq(1, 1000,1)))), ~sum(.x)))
  
  
}


for(i in names(HCVcost_ScenCum_dis)){
  Outcome_Base[["Lifecost"]][[i]] <- cost_y[[i]]
  
  Outcome_Base[["LifeQALY"]][[i]] <- qaly_y[[i]]
}



# timepoint 

# year 2015
# year 2022
# year 2030
# year 2090 

calibrateY <- 2004
Y15 <- 2015 
Y22 <- 2022
Y30 <- 2030
Y90 <- 2090
Y15_pt <- length(seq(calibrateY, Y15, by = HCV$timestep))
Y22_pt <- length(seq(calibrateY, Y22, by = HCV$timestep))
Y30_pt <- length(seq(calibrateY, Y30, by = HCV$timestep))
Y90_pt <- length(seq(calibrateY, Y90, by = HCV$timestep))
# status quo in each PrEP & HIV treatment scenario

#### outcomes ####
Y <- list()

Y_cost <- list()

Y_QALY <- list()


for(n in names(Outcome_Base[["Incidence"]])){ 
  Y[[n]] <- Outcome_Base[["Incidence"]][[n]]%>%
    filter(year == Y22 - HCV$cabY + 1)%>%
    select(paste0("set", seq(1, 1000,1)))%>%t()
  
  
  Y_cost[[n]] <- Outcome_Base$Lifecost[[n]]%>%
    select(paste0("set", seq(1, 1000,1)))%>%t()
  
  Y_QALY[[n]] <- Outcome_Base$LifeQALY[[n]]%>%
    select(paste0("set", seq(1, 1000,1)))%>%t()
  
}


#### parameters ####
# Param_estimatesOff
# Param_Pops >| changing to % of subpops 
# Param_dfListExtend
# Param_costList[[set]][[state]], QALY, flow, POCab, ReflexRNA, POCRNA


# subpop% 
popP <- list()
coname <- c("HIVD", "HIVunD", "HIVN", "HIVPrEP")

subpop_percent_range$`HIV-`$PrEPsame <- subpop_percent_range$`HIV-`$Base
for(i in names(subpop_percent_range)){
  subpop_percent_range[[i]][["PrEPsame"]] <- subpop_percent_range[[i]][["Base"]]
  
  for(m in names(Outcome_Base[["Incidence"]])){ 
    
    popP[[i]][[m]] <- subpop_percent_range[[i]][[m]]%>%
      filter(year == (Y22 - HCV$cabY + 1))%>%
      select(paste0("set", seq(1,HCV$numberSamples,1)))%>%t()
    
    colnames(popP[[i]][[m]]) <-i
    
  }
  
  
}


popP <- popP%>%list_transpose() # transpose list 

pop_testing <- list()

for(m in names(popP)){ 
  
  pop_testing[[m]] <- do.call("cbind", popP[[m]])
  
  
  }

for(m in names(popP)){
  
  colnames(pop_testing[[m]]) <- coname
}


# create dataframe 
aa <- list()
sen_estimate <- list()
for(set in 1:HCV$numberSamples){ 
  
  for( i in names(Param_dfListExtend[[1]])){ 
    aa[[i]][[set]] <- unique(as.vector(Param_dfListExtend[[set]][[i]][, , Y22_pt])) 
    
    names(aa[[i]][[set]]) <- c(paste0(i, ".",1:length(aa[[i]][[set]])))
  }
  
  sen_estimate[[set]] <- Param_estimatesOff[[set]][Y22_pt, ]  
}


sen_estimate_bind <- do.call("rbind", sen_estimate)

tt <- list()

  for( i in names(Param_dfListExtend[[1]])){
     tt[[i]] <- do.call("rbind", aa[[i]])
      }

tt <- do.call("cbind", tt)

tn <- combinParam(Param_fib, HCV)
tx <- combinParam(Param_disease_progress, HCV)

#### cost and QALY parameters #### 
costqaly <- list()
for(set in 1:HCV$numberSamples){ 
  
  costqaly[["state"]][[set]] <- matrix(Param_costList[[set]]$state[, , Y90_pt], nrow = 1)
  
  names(costqaly[["state"]][[set]]) <- c(paste0(rep(dimnames(Param_costList[[1]]$state)[[2]], 
                                                 each = HCV$npops), "_costS_",1:HCV$npops))
  
  costqaly[["state"]][[set]] <- costqaly[["state"]][[set]][!duplicated(costqaly[["state"]][[set]])]
  
  # QALY
  costqaly[["QALY"]][[set]] <- as.vector(Param_costList[[set]]$QALY[, , Y90_pt])
  
  names(costqaly[["QALY"]][[set]]) <- c(paste0(rep(dimnames(Param_costList[[set]]$QALY)[[2]], 
                                             each = HCV$npops), "_qaly_",1:HCV$npops))
  
  costqaly[["QALY"]][[set]] <- costqaly[["QALY"]][[set]][!duplicated(costqaly[["QALY"]][[set]])]
  
  
  # costflow
  costqaly[["costflow"]][[set]] <- as.vector(Param_costList[[set]]$flow[, , Y90_pt])
  
  names(costqaly[["costflow"]][[set]]) <- c(paste0(rep(dimnames(Param_costList[[set]]$flow)[[2]], 
                                                 each = HCV$npops), "_costS_",1:HCV$npops)) 
  
  
  costqaly[["costflow"]][[set]] <- costqaly[["costflow"]][[set]][!duplicated(costqaly[["costflow"]][[set]])]
  
}

for( n in names(costqaly)){ 
  
  costqaly[[n]] <- do.call("rbind", costqaly[[n]])
  
  }

test_t <- list()
ttnm <- list()

for (i in names(Outcome_Base[["Incidence"]])){ 
  test_t[["Inci"]][[i]] <- cbind(pop_testing[[i]], tt, sen_estimate_bind, tx)%>%
    as.data.frame()
  
  test_t[["cost"]][[i]] <- cbind(pop_testing[[i]], tt, sen_estimate_bind, tx, 
                                 costqaly$state,costqaly$costflow)%>%
    as.data.frame()
  
  test_t[["qaly"]][[i]] <- cbind(pop_testing[[i]], tt, sen_estimate_bind, tx,
                                 costqaly$QALY)%>%as.data.frame()
  
  
  ttnm[["Inci"]][[i]] <- test_t[["Inci"]][[i]][,colSums(test_t[["Inci"]][[i]], 
                                                        na.rm = TRUE) != 0]%>%
    select(-c("HCVP1", "HCVP3", "HCVP4", "HIV_undiag_diseaseprog",
              "HIV_diag_diseaseprog"))
  
  ttnm[["cost"]][[i]] <- test_t[["cost"]][[i]][,colSums(test_t[["cost"]][[i]], 
                                                        na.rm = TRUE) != 0]%>%
    select(-c("HCVP1", "HCVP3", "HCVP4", "HIV_undiag_diseaseprog",
              "HIV_diag_diseaseprog"))
  
  ttnm[["qaly"]][[i]] <- test_t[["qaly"]][[i]][,colSums(test_t[["qaly"]][[i]], 
                                                        na.rm = TRUE) != 0]%>%
    select(-c("HCVP1", "HCVP3", "HCVP4", "HIV_undiag_diseaseprog",
              "HIV_diag_diseaseprog"))
  
  
  }

# HCV testing 
# remove columns with all 0


# Global sensitivity analysis: Partial rank coefficient correlation 
tic <- proc.time()
tx <- as.data.frame(tx)
test_tx <- tx[, c(2,5)]
colnames(test_tx) <- c("af01", "f_f")

str(ttnm[["cost"]][[1]])


tic <- proc.time()

PRCC_inci <- pcc(ttnm[["Inci"]][[1]], Y[[1]], semi = FALSE, 
                 rank = TRUE, nboot = 1000, conf = 0.95)

save(PRCC_inci, 
     file = file.path(outputdt, "PRCC_inci.rda"))
    
PRCC_cost <- pcc(ttnm[["cost"]][[1]], Y_cost[[1]], semi = FALSE, 
                 rank = TRUE, nboot = 1000, conf = 0.95)

save(PRCC_cost, 
     file = file.path(outputdt, "PRCC_cost.rda"))
    
PRCC_qaly <- pcc(ttnm[["qaly"]][[1]], Y_QALY[[1]], semi = FALSE, 
                 rank = TRUE, nboot = 1000, conf = 0.95)

save(PRCC_qaly, 
     file = file.path(outputdt, "PRCC_qaly.rda"))
  
toc <- proc.time() - tic 

toc

load(file.path(outputdt , "PRCC_inci.rda"))

load(file.path(outputdt , "PRCC_cost.rda"))

load(file.path(outputdt , "PRCC_qaly.rda"))

PRCC_inci$PRCC

#### 
sen_ana <- list()

sen_ana[["Inci"]] <- PRCC_inci$PRCC%>%
  mutate(parameter = rownames(.))

sen_ana[["cost"]] <- PRCC_cost$PRCC%>%
  mutate(parameter = rownames(.))

sen_ana[["qaly"]] <- PRCC_qaly$PRCC%>%
  mutate(parameter = rownames(.)) 


# confidence interval calculation 




# reorder by abstract value of original 
# extract first 10 
sen_ana_rank <- list()
for (i in names(sen_ana)){ 
  
  sen_ana_rank[[i]] <- sen_ana[[i]]%>%arrange(., desc(abs(original)))%>%
    head(., 11)
  }
sen_ana_rank$Inci <- sen_ana_rank$Inci%>%filter(parameter!= "a_f04")%>%
  mutate(parameter = factor(parameter,
                            levels = c("a_f01", "f4_hcc1", "spc1", "HIVD", 
                                       "f3_hcc1", "f4_dc1", "HIVPrEP",
                                       "Cure_morhcc_Reduction", "spc3", "f3_hcc4"),
                            labels = c("disease progression: Acute → F0",
                                       "disease progression: F4→ HCC (HIV-)",
                                       "Spontaneous clearance (HIV-)",
                                       "proportion of MSM diagnosed with HIV and on treatment",
                                       "disease progression: F3→ HCC (HIV-)",
                                       "disease progression: F4→ DC (HIV-)",
                                       "proportion of MSM on PrEP",
                                       "HCC-related death (post-cured)",
                                       "Spontaneous clearance (HIV+)",
                                       "disease progression: F3→ HCC (diagnosed with HIV and on treatment)"
                                       )))

sen_ana_rank$cost <- sen_ana_rank$cost[c(1:10), ]%>%
  mutate(parameter = 
           factor(parameter, 
                  levels = c("cured.2", "HIVD", "a_f01", "spc4", "mordc", 
                             "morhcc", "morlt", "HIVd_mor", "dc_hcc4", "spc2"),
                  labels = c("SVR rate", 
                             "proportion of MSM diagnosed with HIV and on treatment",
                             "disease progression: Acute → F0",
                             "Spontaneous clearance (HIV+)",
                             "DC related death",
                             "HCC related death",
                             "liver transplant related death",
                             "Mortality of MSM diagnosed with HIV and on treatment",
                             "disease progression: DC → HCC (HIV diagnosed with HIV and on treatment)",
                             "Spontaneous clearance (HIV-)")))


sen_ana_rank$qaly <- sen_ana_rank$qaly[c(1:10), ]%>%
  mutate(parameter = 
           factor(parameter, 
                  levels = c("HIVD", "cured.2", "hcc_lt1", "spc2", "f4_hcc4", 
                             "f3_f42", "f4_hcc3", "hcc_lt2", "HIVN", "dc_lt3"),
                  labels = c("proportion of MSM diagnosed with HIV and on treatment",
                             "SVR rate", 
                             "disease progression: HCC → liver transplant (HIV- not on PrEP)",
                             "Spontaneous clearance (HIV-)",
                             "disease progression: F4 → HCC (HIV diagnosed on treatment)",
                             "disease progression: F3 → F4 (HIV-)",
                             "disease progression: F4 → HCC (HIV undiagnosed)",
                             "disease progression: HCC → liver transplant (HIV- on PrEP)",
                             "proportion of MSM living with HIV and not on PrEP",
                             "disease progression: DC → liver transplant (HIV undiagnosed)")))

tornado_p[[1]] <- sen_ana_rank[[1]]%>%arrange(original)

tornado_p <- list()
### unsolved: rename parameters 
title_name <- c("HCV incidence", "Lifetime cost", "Lifetime QALYs")
for(i in 1:length(names(sen_ana_rank))){ 
  
  tornado_p[[i]] <- sen_ana_rank[[i]]%>%arrange(abs(original))%>%
     ggplot(data = ., aes(x = parameter, 
                              y = original))+ 
    geom_segment( aes(xend=parameter, y=`min. c.i.`, yend = `max. c.i.`)) +
    geom_point(aes(y = original), size=1, color="orange") +
    geom_hline( aes(yintercept = 0), linetype = "dashed") +
    coord_flip() +
    scale_y_continuous(limits = c(-1, 1)) +
    theme_bw() +
    xlab("Parameters") +
    ylab("Partial rank correlation coefficients") + 
    
    ggtitle(title_name[i]) + 
    theme_Publication(base_size = 10) + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5,hjust = 0.5,
                               face = "bold",size = rel(1)))
    
  
  }
tornado_p[[1]]
tornado_p[[2]]
tornado_p[[3]]

# need to edit the label of the parameters afterward 

ggsave(path = outputfig, file="PRCC_coefficient_inci.png", tornado_p[[1]], height = 4, 
       width = 8, dpi = 800)
ggsave(path = outputfig, file="PRCC_coefficient_cost.png", tornado_p[[2]], height = 4, 
       width = 8, dpi = 800)

ggsave(path = outputfig, file="PRCC_coefficient_qaly.png", tornado_p[[2]], height = 4, 
       width = 8, dpi = 800)
