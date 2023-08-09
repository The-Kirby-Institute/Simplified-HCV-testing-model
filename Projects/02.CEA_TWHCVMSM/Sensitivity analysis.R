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
# we specify the number of cores/workers we want to use
registerDoParallel(cores = detectCores() - 1)

#file path of "TWHCV-model" project
codepath <- file.path(here() %>% dirname(), 'TWHCV-model/03. Code/Functions')

DataFolder <- file.path(here(), "01. DATA/model input")

# Rda file path 
# load the .rda file of base estimate 
rdapath <- file.path(here()%>%dirname(), "Taiwan-MSM-HCV-model")

outputdt <- here("02. Output/RDA")

projectFile <- file.path(rdapath , paste0("HCVModel",".rda"))

projectVars <- load(projectFile)

load(file.path(rdapath , paste0("HCVModel", "cali",".rda")))

load(file.path(outputdt , "cost.rda"))

load(file.path(outputdt , "toft.rda"))

load(file.path(outputdt , "Param_dfExt.rda"))

load(file.path(outputdt , "mainS_PrEPHIV.rda"))

load(file.path(outputdt , "pop_rate.rda"))

load(file.path(outputdt, "measureOutcome.rda"))

load(file.path(rdapath , paste0("HCVModel", "cali",".rda")))

load(file.path(rdapath , paste0("HCVModel", "param",".rda")))

source(file.path(here(),"AggregateRes.R"))

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
#### 
sen_ana <- list()

sen_ana[["Inci"]] <- PRCC_inci$PRCC%>%
  mutate(parameter = rownames(.))

sen_ana[["cost"]] <- PRCC_cost$PRCC%>%
  mutate(parameter = rownames(.))

sen_ana[["qaly"]] <- PRCC_qaly$PRCC%>%
  mutate(parameter = rownames(.))

# reorder by abstract value of original 
# extract first 10 
sen_ana_rank <- list()
for (i in names(sen_ana)){ 
  
  sen_ana_rank[[i]] <- sen_ana[[i]]%>%arrange(., desc(abs(original)))%>%
    head(., 10)
  }
sen_ana_rank$qaly
tornado_p <- list()
### unsolved: rename parameters 

for(i in names(sen_ana_rank)){ 
  
  tornado_p[[i]] <- sen_ana_rank[[i]]%>%arrange(abs(original))%>%
    mutate(param = factor(parameter, levels = parameter))%>%
    ggplot(data = ., aes(x = param, 
                              y = original))+ 
    geom_segment( aes(xend=param, yend=0)) +
    geom_point( size=2.5, color="orange") +
    coord_flip() +
    scale_y_continuous(limits = c(-1, 1)) +
    theme_bw() +
    xlab("") +
    
    ggtitle(i)
  
  }

tornado_p$cost
