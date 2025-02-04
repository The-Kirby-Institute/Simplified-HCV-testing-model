# this script imports the datasets tidy up in Resultssummary.R 
# {Res_dt.rda} and {epiRes_timestep.rda} 
# we aggregate the number to annual numbers and apply cost discount in this script. 
# three datasets we work on in this script: {Num_box}, {Resflow_dt}, {Rescost_dt}
gc()
rm(list = ls())
tic <- proc.time()

project_name <- "POC_AU"

codefun_path <- paste("/Users/jjwu/Documents/Simplified-HCV-testing-model")

data_path <- paste("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/", 
                   project_name, sep = "")
# Load useful libraries

library("readr")
library("dplyr")
library("tidyr")
library("purrr")
library("parallel")
library("pacman")
library("doMC")
library("openxlsx")

Rcode <- file.path(codefun_path, "03. Code")

DataFolder <- file.path(data_path, "01. DATA/model input" )
OutputFolder <- file.path(data_path, "02. Output")
# project specific code path 
Proj_code <- file.path(codefun_path, paste0("projects/", project_name))

load(file.path(OutputFolder, paste0(project_name, ".rda")))


# load rda into list of list 
rda2list <- function(file) {
  e <- new.env()
  load(file, envir = e)
  as.list(e)
}

files <- list.files(OutputFolder, pattern = paste0(project_name,"Res_dt_"))

Res_dt <- Map(rda2list, file.path(OutputFolder, files))

name_file <- sub("POC_AURes_dt_", "", files)

names(Res_dt) <- tools::file_path_sans_ext(name_file)

source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Proj_code, "/model_timestep.R")) 

AUdiscount <- 0.05
unitC_DAA <- 35956.37
unitC_secline_DAA <- 44613.66
cap <- 200000000
endY <- 100
# summarize the number of flow for each pop in each year 
par_col <- c("best", paste0("set", seq(1, POC_AU$numberSamples,1)))
names(Res_dt) <- c(names(Res_dt)[-7], "Status quo")
Res_numbox <- list()

for(i in names(Res_dt)){ 
  Res_numbox[[i]] <- Res_dt[[i]]$Num_box[[i]]
}

save(Res_numbox,
     file = file.path(OutputFolder,
                      paste0(project_name,"Res_numbox" ,".rda"))) 



Resflow_year_pop <- list()

for(i in names(Res_dt)){ 
  Resflow_year_pop[[i]] <- Res_dt[[i]]$Resflow_dt[[i]]
}

for(i in names(Res_dt)){
  for(indic in names(Resflow_year_pop[[1]])){
    Resflow_year_pop[[i]][[indic]] <- Resflow_year_pop[[i]][[indic]]%>%
      as_tibble()%>%ungroup()%>%arrange(population, year)%>%
      group_by(year, population)%>%
      summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%
      arrange(year)
  }
}

Resflow_year_all <- list()
for(i in names(Resflow_year_pop)){
  for(indic in names(Resflow_year_pop[[1]])){
    
    Resflow_year_all[[i]][[indic]] <- Resflow_year_pop[[i]][[indic]]%>%
      ungroup()%>%group_by(year)%>%
      summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%
      arrange(year)
    
  }
}

# scenarios 
Resflow_sc_year_pop <- list()

for(i in names(Res_dt)){ 
  Resflow_sc_year_pop[[i]] <- Res_dt[[i]]$Resflow_sc_dt[[i]]
}

for(i in names(Res_dt)){
  for(indic in names(Resflow_sc_year_pop[[1]])){
    Resflow_sc_year_pop[[i]][[indic]] <- Resflow_sc_year_pop[[i]][[indic]]%>%
      as_tibble()%>%ungroup()%>%arrange(population, year)%>%
      group_by(year, population)%>%
      summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%
      arrange(year)
    
    
  }
}

Resflow_sc_year_all <- list()
for(i in names(Resflow_sc_year_pop)){
  for(indic in names(Resflow_sc_year_pop[[1]])){
    
    Resflow_sc_year_all[[i]][[indic]] <- Resflow_sc_year_pop[[i]][[indic]]%>%
      ungroup()%>%group_by(year)%>%
      summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%
      arrange(year)
    
  }
}




# separate cost of DAA and other costs related to treatment initiation  
# treatment and retreat 
# data sum in yearly
RescostDAA <- list()
for(i in names(Res_dt)){ 
  RescostDAA[[i]] <- Res_dt[[i]]$Rescost_dt[[i]]
}

for(i in names(Res_dt)){
  for(indic in names(RescostDAA[[i]])){
    RescostDAA[[i]][[indic]] <- RescostDAA[[i]][[indic]]%>%
      as_tibble()%>%ungroup()%>%arrange(timestep,population)%>%
      mutate(year = c(rep(seq(1, endY - 2, 1), each = POC_AU$npops*1/POC_AU$timestep),
                      c(rep(endY - 1, POC_AU$npops*(1/POC_AU$timestep)))))%>%
      group_by(year, population)%>%
      summarise(across(c(par_col),~ sum(.x, na.rm = TRUE)))%>%
      arrange(year)
    }
}

cost_TreatDAA <- list()
cost_RetreatDAA <- list()
cost_TreatDAA_sc <- list()
cost_TreatOther <- list()
cost_RetreatOther <- list()
cost_TreatOther_sc <- list()
cost_totalDAA <- list()
index_col <- list()

for(i in names(RescostDAA)){
  
  index_col[[i]] <- cbind.data.frame(year = RescostDAA[[i]]$cost_Treatment$year,
                          population = RescostDAA[[i]]$cost_Treatment$population)
  cost_TreatDAA[[i]] <-  cbind(as.data.frame(index_col[[i]]), 
                               Resflow_year_pop[[i]]$Treatment[ ,c(par_col)]*unitC_DAA)
    
  cost_RetreatDAA[[i]] <- cbind(as.data.frame(index_col[[i]]), 
                                Resflow_year_pop[[i]]$Retreat[, c(par_col)]*unitC_secline_DAA)
  
  cost_TreatDAA_sc[[i]] <- cbind(as.data.frame(index_col[[i]]), 
                                 Resflow_sc_year_pop[[i]]$Treatment_sc[, c(par_col)]*unitC_DAA)
}

index_col <- list()
for(i in names(RescostDAA)){
  
  index_col[[i]] <- cbind.data.frame(year = RescostDAA[[i]]$cost_Treatment$year, 
                          population = RescostDAA[[i]]$cost_Treatment$population)

  cost_TreatOther[[i]] <- cbind(as.data.frame(index_col[[i]]), 
                                (RescostDAA[[i]][["cost_Treatment"]][, c(par_col)] - 
                                   cost_TreatDAA[[i]][, c(par_col)] - 
                                  RescostDAA[[i]][["cost_Treatment_sc"]][, c(par_col) ]))
  
  cost_RetreatOther[[i]] <- cbind(as.data.frame(index_col[[i]]), 
                                  (RescostDAA[[i]][["cost_Retreat"]][, c(par_col)] - 
                                    cost_RetreatDAA[[i]][, c(par_col)]))
  
  cost_TreatOther_sc[[i]] <- cbind(as.data.frame(index_col[[i]]), 
                                   (RescostDAA[[i]][["cost_Treatment_sc"]][, c(par_col)] - 
                                      cost_TreatDAA_sc[[i]][, c(par_col)]))
  
  cost_totalDAA[[i]] <- cbind(as.data.frame(index_col[[i]]), 
                              (cost_TreatDAA[[i]][, c(par_col)] + 
                                 cost_RetreatDAA[[i]][, c(par_col)] + 
                                 cost_TreatDAA_sc[[i]][, c(par_col)]))
}  


for(i in names(RescostDAA)){
 RescostDAA[[i]] <- append(RescostDAA[[i]], 
                   list("cost_TreatDAA" = cost_TreatDAA[[i]], 
                        "cost_RetreatDAA" = cost_RetreatDAA[[i]],
                        "cost_TreatDAA_sc" = cost_TreatDAA_sc[[i]],
                        "cost_TreatOther" = cost_TreatOther[[i]],
                        "cost_RetreatOther" = cost_RetreatOther[[i]],
                        "cost_TreatOther_sc" = cost_TreatOther_sc[[i]],
                        "cost_totalDAA" = cost_totalDAA[[i]])
                   )
}

RescostDAA_totalpop <- list()

for(i in names(RescostDAA)){
  for(indic in names(RescostDAA[[1]])){
    RescostDAA_totalpop[[i]][[indic]] <- RescostDAA[[i]][[indic]]%>%as_tibble()%>%
      ungroup()%>%arrange(year)%>%group_by(year)%>%
      summarise(across(c(par_col),~ sum(.x, na.rm = TRUE)))%>%ungroup()
  }
  
}

RescostDAA_sc <- list()
for(i in names(RescostDAA)){
  
  RescostDAA_sc[[i]] <- list("cost_TreatDAA" = cost_TreatDAA_sc[[i]],
                             "cost_TreatOther" = cost_TreatOther_sc[[i]],
                             "cost_totalDAA" = cost_TreatDAA_sc[[i]])  
}

RescostDAA_sc_totalpop <- list()

for(i in names(RescostDAA_sc)){
  for(indic in names(RescostDAA_sc[[1]])){ 
    RescostDAA_sc_totalpop[[i]][[indic]] <- RescostDAA_sc[[i]][[indic]]%>%ungroup()%>%
      arrange(year)%>%group_by(year)%>%
      summarise(across(c(par_col),~ sum(.x, na.rm = TRUE)))%>%ungroup()
    }
  }


# remove the original cost_Treatment and cost_Retreat since we breakdown the cost of treatment and retreat 
# into DAA,and nonDAA cost.

Rescost_DAA_dt <- RescostDAA

for(i in names(RescostDAA)){ 
  Rescost_DAA_dt[[i]]$cost_Treatment <- NULL
  Rescost_DAA_dt[[i]]$cost_Retreat <- NULL

}


Rescost_DAA_sc_dt <- list()

for(i in names(RescostDAA)){ 
  Rescost_DAA_sc_dt[[i]] <- list("cost_ab_sc" = RescostDAA[[i]]$cost_ab,
                                 "cost_RNA_sc" = RescostDAA[[i]]$cost_RNA_sc,
                                 "cost_POCT_sc" = RescostDAA[[i]]$cost_POCT_sc,
                                 "cost_TreatDAA" = cost_TreatDAA_sc[[i]],
                                 "cost_TreatOther" = cost_TreatOther_sc[[i]],
                                 "cost_totalDAA" = cost_TreatDAA_sc[[i]])  
 
 
}

# cost_total DAA, exculdeDAA replaced cost_Treatment & cost_Retreat 
# valide with cost_total, 
# calculating total_cost_cap  
# extract the columns to replace in cost_treatment& cost_retreat 
# "cost_totalDAA_y": "cost_totalDAAcap_discum"
Rescost_year <- list()

for(i in names(Rescost_DAA_dt)){ 
  for(indic in names(Rescost_DAA_dt[[1]])){ 
    Rescost_year[[i]][[indic]] <- Rescost_DAA_dt[[i]][[indic]]%>%
      as_tibble()%>%mutate(year = year + POC_AU$cabY - 1,
                           id = year - POC_AU$simY, 
                           discount = ifelse(id>=0, (1 + AUdiscount)^id, NA))
     
      
    }
  
  }


Rescost_sc_year <- list() 

for(i in names(Rescost_DAA_sc_dt)){ 
  
  for(indic in names(Rescost_DAA_sc_dt[[1]])){ 
    Rescost_sc_year[[i]][[indic]] <- Rescost_DAA_sc_dt[[i]][[indic]]%>%
      as_tibble()%>%mutate(year = year + POC_AU$cabY - 1,
                           id = year - POC_AU$simY, 
                           discount = ifelse(id>=0, (1 + AUdiscount)^id, NA))
    }
  }


# cumulative yearly cost by population  & discount value 
Rescost_yearcum_pop <- list() 

for(i in names(Rescost_year)){ 
  
  for(indic in names(Rescost_year[[1]])){ 
    Rescost_yearcum_pop[[i]][[indic]] <- Rescost_year[[i]][[indic]]%>%
      as_tibble()%>%arrange(population)%>%group_by(population)%>%
      mutate(across(c(par_col), cumsum, .names = "{col}"))
    
  }
  
}

Rescost_yearcum_sc_pop <- list()

for(i in names(Rescost_sc_year)){ 
  for(indic in names(Rescost_sc_year[[1]])){
    Rescost_yearcum_sc_pop[[i]] <- Rescost_sc_year[[i]][[indic]]%>%
      as_tibble()%>%arrange(population)%>%group_by(population)%>%
      mutate(across(c(par_col), cumsum, .names = "{col}"))
      
  }
}

Rescost_year_all <- list()

for(i in names(Rescost_year)){ 
  
  for(indic in names(Rescost_year[[1]])){ 
    Rescost_year_all[[i]][[indic]] <- Rescost_year[[i]][[indic]]%>%
      as.data.frame()%>%
      ungroup()%>%
      group_by(year)%>%
      summarise(across(c(par_col),~ sum(.x, na.rm = TRUE)))%>%ungroup()
  }
}


for(i in names(Rescost_year)){ 
  Rescost_year_all[[i]][["cost_totalDAA_Cap"]] <- 
    Rescost_year_all[[i]][["cost_totalDAA"]]%>%
    mutate(across(c(par_col), ~ ifelse(.>=cap, cap, . ),.names = "{col}"))
  
  Rescost_year_all[[i]][["cost_total_Cap"]] <- 
    cbind(year =Rescost_year_all[[i]]$cost_totalDAA$year,
          as.data.frame(Rescost_year_all[[i]][["cost_compartment"]][, c(par_col)] + 
                          Rescost_year_all[[i]][["cost_ab"]][, c(par_col)] + 
                          Rescost_year_all[[i]][["cost_RNA"]][, c(par_col)] + 
                          Rescost_year_all[[i]][["cost_POCT"]][, c(par_col)] +
                          Rescost_year_all[[i]][["cost_totalDAA_Cap"]][, c(par_col)] + 
                          Rescost_year_all[[i]][["cost_TreatOther"]][, c(par_col)] + 
                          Rescost_year_all[[i]][["cost_RetreatOther"]][, c(par_col)] + 
                          Rescost_year_all[[i]][["cost_Cured"]][, c(par_col)]))
}

Rescost_year_sc_all <- list()

for(i in names(Rescost_sc_year)){ 
  for(indic in names(Rescost_sc_year[[1]])){ 
    
    Rescost_year_sc_all[[i]][[indic]] <- Rescost_sc_year[[i]][[indic]]%>%
      as.data.frame()%>%
      ungroup()%>%
      group_by(year)%>%
      summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%ungroup()
    
  }
}

Rescost_disyear_all <- list()

for(i in names(Rescost_year_all)){ 
  for(indic in names(Rescost_year_all[[1]])){ 
    
    Rescost_disyear_all[[i]][[indic]] <- Rescost_year_all[[i]][[indic]]%>%
      as.data.frame()%>%
      ungroup()%>%
      mutate(id = year - POC_AU$simY, 
             discount = ifelse(id>=0, (1 + AUdiscount)^id, NA))%>%
      mutate(across(c(par_col), ~./discount,
                    .names = "{col}"))
    
  }
}

Rescost_disyear_sc_all <- list()

for(i in names(Rescost_sc_year)){ 
  for(indic in names(Rescost_sc_year[[1]])){ 
    
    Rescost_disyear_sc_all[[i]][[indic]] <- Rescost_year_sc_all[[i]][[indic]]%>%
      as.data.frame()%>%
      ungroup()%>%
      mutate(id = year - POC_AU$simY, 
             discount = ifelse(id>=0, (1 + AUdiscount)^id, NA))%>%
      mutate(across(c(par_col), ~./discount,
                    .names = "{col}"))
    
  }
}


# cumulative cost each year for overall pops
Rescost_yearcum_all <- list()

for(i in names(Rescost_year_all)){ 
  for(indic in names(Rescost_year_all[[1]])){
    Rescost_yearcum_all[[i]][[indic]] <- Rescost_year_all[[i]][[indic]]%>%
      ungroup()%>%
      mutate(across(c(par_col), cumsum,
                    .names = "{col}"))%>%
      mutate(id = year - POC_AU$simY, 
             discount = ifelse(id>=0, (1 + AUdiscount)^id, NA))
    }
  }

# cumulative cost each year for overall pops
Rescost_yearcum_sc_all <- list()

for(i in names(Rescost_year_sc_all)){ 
  for(indic in names(Rescost_year_sc_all[[1]])){
    Rescost_yearcum_sc_all[[i]][[indic]] <- Rescost_year_sc_all[[i]][[indic]]%>%
      ungroup()%>%
      mutate(across(c(par_col), cumsum,
                    .names = "{col}"))%>%
      mutate(id = year - POC_AU$simY, 
             discount = ifelse(id>=0, (1 + AUdiscount)^id, NA))
    
  }
}

# discount cumulative cost each year for overall pops

Rescost_discum_all <- list()

for(i in names(Rescost_year_all)){ 
  for(indic in names(Rescost_year_all[[1]])){
    Rescost_discum_all[[i]][[indic]] <- Rescost_yearcum_all[[i]][[indic]]%>%
      ungroup()%>%
      mutate(across(c(par_col), ~./discount,
                    .names = "{col}"))
  }
}

# discount cumulative cost each year for overall pops
Rescost_discum_sc_all <- list()

for(i in names(Rescost_year_sc_all)){ 
  for(indic in names(Rescost_year_sc_all[[1]])){
    Rescost_discum_sc_all[[i]][[indic]] <- Rescost_yearcum_sc_all[[i]][[indic]]%>%
      ungroup()%>%
      mutate(across(c(par_col), ~ ./discount,
                    .names = "{col}"))
    
  }
}





# save rda files 
save(Resflow_year_pop, Resflow_year_all, 
     Resflow_sc_year_pop, Resflow_sc_year_all,
     RescostDAA, RescostDAA_totalpop, 
     Rescost_year, Rescost_sc_year,
     Rescost_yearcum_pop, Rescost_yearcum_sc_pop, 
     Rescost_year_all, Rescost_year_sc_all, 
     Rescost_disyear_all, Rescost_disyear_sc_all,
     Rescost_yearcum_all, Rescost_yearcum_sc_all,
     Rescost_discum_all, Rescost_discum_sc_all,
     AUdiscount,
     unitC_DAA,
     unitC_secline_DAA,
     cap,
     file = file.path(OutputFolder,
                      paste0(project_name,"Res_flowcost" ,".rda"))) 




