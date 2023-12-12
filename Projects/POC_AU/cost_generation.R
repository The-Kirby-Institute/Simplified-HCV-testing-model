# apply function {cost_model} to calculate the cost of each scenario 
rm(list = ls())

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

Rcode <- file.path(codefun_path, "03. Code")

DataFolder <- file.path(data_path, "01. DATA/model input" )
OutputFolder <- file.path(data_path, "02. Output")
# create the Subfolder to save figs 
# dir.create(file.path(data_path, "02. Output/Figs"))
OutputFig <- file.path(data_path, "02. Output/Figs")
load(file.path(OutputFolder, paste0(project_name, ".rda")))
load(file.path(OutputFolder, paste0(project_name, "cali.rda")))
load(file.path(OutputFolder, paste0(project_name, "cali_timev.rda")))
load(file.path(OutputFolder, paste0(project_name, "Sce_np.rda")))
load(file.path(OutputFolder, paste0(project_name, "S_NP_test.rda")))
load(file.path(OutputFolder, paste0(project_name, "S_NPscale_test.rda")))
load(file.path(OutputFolder, paste0(project_name, "Sce_npscale.rda")))
load(file.path(OutputFolder, paste0(project_name, "coverage_sq.rda")))
load(file.path(OutputFolder, paste0(project_name, "coverage_np.rda")))
load(file.path(OutputFolder, paste0(project_name, "coverage_npscale.rda")))
source(file.path(Rcode, "/Functions/HCV_model.R"))

source(file.path(Rcode, "/Functions/plotManuscript.R"))

source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Rcode, "/Functions/check_steady.R"))
source(file.path(Rcode, "/Functions/check_steady.R"))

source(file.path(codefun_path,paste0("Projects/", project_name, "/cost_model.R")))
# import cost 
#### cost #### 
files <- list.files(path = paste0(DataFolder, 
                                  "/cost/", sep =  ""), pattern = '*.csv')


costdfList <- lapply(files, function(f) {
  
  df <- read.csv(file.path(paste0(DataFolder, "/cost/", f, sep = "")), header = TRUE)
  
  df <- df[, -1]
  
  df <- df%>%as_tibble()
  
  df <- as.matrix(df, nrow = npops, ncol = length(.) + 1)
})

names(costdfList) <- c(gsub("^|.csv", "", files)) # ^: from beginning, \ end before .csv

c_sq <- cost_model(POC_AU, costdfList, coverage = coverage_sq, 
                 dfList, dfList_NP = dfList, 
                 calibrateInit,
                 S_Yint = 2022, S_Yfir = NULL, S_Ymid = NULL, S_Yend = 2051, 
                 censor_Y = 2030)

c_np <- cost_model(POC_AU, costdfList, coverage = coverage_np, 
                    dfList, dfList_NP = dfList_NP, 
                    Sce_np,
                    S_Yint = 2022, S_Yfir = NULL, S_Ymid = NULL, S_Yend = 2051,
                   censor_Y = 
                     2030)

c_npscale <- cost_model(POC_AU, costdfList, coverage = coverage_npscale, 
                    dfList, dfList_NP = dfList_NPscale, 
                    Sce_npscale,
                    S_Yint = 2022, S_Yfir = 2024, S_Ymid = 2027, S_Yend = 2051,
                    censor_Y = 2030)



indicators <- c("costTestingAb", "costTestingRNA", "costTestingPOCT", 
                "costTreatment", "costCured", "costRetreat") 

Sce_cost <- list("status quo" = c_sq, 
                 "National program" = c_np, 
                 "National program scale up" = c_npscale)
# trimp timepoints 
for(n in names(Sce_cost)){ 
  for(i in indicators){ 
    Sce_cost[[n]][[i]] <- Sce_cost[[n]][[i]][, c(1:dim(Sce_cost[[n]][["allPops"]])[3])]
  }
  Sce_cost[[n]][["costPops"]] <- Sce_cost[[n]][["costPops"]][, , c(1:dim(Sce_cost[[n]][["allPops"]])[3])]
}


x_np <- list()
for(n in names(Sce_cost)){
  for(i in indicators){ 
    x_np[[n]][[i]] <- indicatorResults(POC_AU, Sce_cost[[n]],
                                  i, 
                                  pop = NULL, range = NULL, 
                                  endY = 100, scenario = NULL)%>%
      mutate(population = rep(POC_AU$popNames, length(population)/POC_AU$npops))
  }
  
  x_np[[n]][["cost_state"]] <- popResults_MidYear(POC_AU, Sce_cost[[n]], 
                                                  Population = POC_AU$popNames,
                                                  Disease_prog = POC_AU$progress_name, 
                                                  Cascade = c("s",POC_AU$cascade_name), 
                                                  param = NULL,
                                                  endYear = 100, 
                                                  allp = "costPops")%>%
    as.data.frame()
  x_np[[n]][["costPops"]] <- x_np[[n]][["cost_state"]]%>%
    group_by(year, population)%>%summarize(best = sum(best))
  
  
  
}

discountR <- 0.03
cost_np <- list()
cost_stage <- list()
for(n in names(Sce_cost)){
  for(s in names(x_np[[1]])){ 
    
    cost_stage[[n]][[s]] <- x_np[[n]][[s]]%>%
      mutate(year = year + POC_AU$cabY - 1)%>%
      mutate(id = year - POC_AU$simY)%>%
      mutate(discount = (1 + discountR)^id)%>%filter(id>= 0)
    
    cost_stage[[n]][[s]] <- cost_stage[[n]][[s]]%>%mutate(discountValue = best/discount)%>%
      group_by(year)%>%summarise(discountValue = sum(discountValue))%>%ungroup()%>%
      mutate(discountValue_cum = cumsum(discountValue))
    
    }
  cost_np[[n]] <- cbind(year = (x_np[[n]]$costTestingAb$year + POC_AU$cabY - 1), 
                        population = x_np[[n]]$costTestingAb$population,
                        best = (x_np[[n]]$costTestingAb$best+ 
                                  x_np[[n]]$costTestingRNA$best +
                                  x_np[[n]]$costTestingPOCT$best +
                                  x_np[[n]]$costTreatment$best +
                                  x_np[[n]]$costCured$best +
                                  x_np[[n]]$costRetreat$best +
                                  x_np[[n]]$costPops$best))%>%as_tibble()
  
  cost_np[[n]]$year <- as.numeric(cost_np[[n]]$year)
  
  cost_np[[n]]$best <- as.numeric(cost_np[[n]]$best)
  
  cost_np[[n]] <- cost_np[[n]]%>%
    mutate(id = year - POC_AU$simY)%>%
    mutate(discount = (1 + discountR)^id)%>%filter(id>= 0)
  
  cost_np[[n]] <- cost_np[[n]]%>%mutate(discountValue = best/discount)%>%
    group_by(year)%>%summarise(discountValue = sum(discountValue))%>%
    ungroup()%>%
    mutate(discountValue_cum = cumsum(discountValue))
  
}

# cost of diagnosis 
c_stage <- list()
for(n in names(cost_stage)){ 
  c_stage[["diagnosis"]][[n]] <- 
    cbind(year = (cost_stage[[n]]$costTestingAb$year + POC_AU$cabY - 1), 
          discountValue = (cost_stage[[n]]$costTestingAb$discountValue + 
                             cost_stage[[n]]$costTestingRNA$discountValue +
                             cost_stage[[n]]$costTestingPOCT$discountValue),
          discountValue_cum = (cost_stage[[n]]$costTestingAb$discountValue_cum+ 
                                 cost_stage[[n]]$costTestingRNA$discountValue_cum +
                                 cost_stage[[n]]$costTestingPOCT$discountValue_cum ))%>%
    as_tibble()
  
  c_stage[["Treatment"]][[n]] <- 
    cbind(year = (cost_stage[[n]]$costTestingAb$year + POC_AU$cabY - 1), 
          discountValue = (cost_stage[[n]]$costTreatment$discountValue + 
                             cost_stage[[n]]$costRetreat$discountValue +
                             cost_stage[[n]]$costCured$discountValue),
          discountValue_cum = (cost_stage[[n]]$costTreatment$discountValue_cum+ 
                                 cost_stage[[n]]$costRetreat$discountValue_cum +
                                 cost_stage[[n]]$costCured$discountValue_cum ))%>%
    as_tibble()
  
  c_stage[["Management"]][[n]] <- 
    cbind(year = (cost_stage[[n]]$costTestingAb$year + POC_AU$cabY - 1), 
          discountValue = (cost_stage[[n]]$costPops$discountValue),
          discountValue_cum = (cost_stage[[n]]$costPops$discountValue_cum))%>%
    as_tibble()
  
}
  
c_stage$diagnosis <-
  dplyr::bind_rows(c_stage$diagnosis, .id = "scenario")
c_stage$Treatment<-
  dplyr::bind_rows(c_stage$Treatment, .id = "scenario")
c_stage$Management <-
  dplyr::bind_rows(c_stage$Management, .id = "scenario")

c_stage <- dplyr::bind_rows(c_stage, .id = "cost indicator")%>%
  mutate(year = year - POC_AU$cabY + 1)

ggplot(data = c_stage) + 
  geom_line(aes(x = year, y = discountValue_cum, colour = scenario)) + 
  facet_wrap(.~`cost indicator`, scale ="free") + 
  scale_x_continuous(limits = c(2022, 2050), expand = c(0, 0))

xt <- list()
discountR <- 0.03
for(n in names(x_np)){ 
  
  xt[[n]] <- x_np[[n]][["cost_state"]]%>%
    mutate(year = year + POC_AU$cabY - 1)%>%
    mutate(id = year - POC_AU$simY)%>%
    mutate(discount = (1 + discountR)^id)%>%filter(id>= 0)%>%
    group_by(year, state)%>%mutate(discountValue = best/discount)%>%
    ungroup()%>%group_by(state)%>%
    mutate(discountValue_cum = cumsum(discountValue))
    
  }

xxt <- dplyr::bind_rows(xt, .id= "scenario")

c_dab <- lapply(cost_stage, function(x) x$costTestingAb)%>%
  dplyr::bind_rows(, .id= "scenario")%>%as.data.frame()

c_drna <- lapply(cost_stage, function(x) x$costTestingRNA)%>%
  dplyr::bind_rows(, .id= "scenario")%>%as.data.frame()

c_drnaonly <- lapply(cost_stage, function(x) x$costTestingPOCT)%>%
  dplyr::bind_rows(, .id= "scenario")%>%as.data.frame()

c_t <- lapply(cost_stage, function(x) x$costRetreat)%>%
  dplyr::bind_rows(, .id= "scenario")%>%as.data.frame()

c_tp <- ggplot(data = c_t) + 
  geom_line(aes( x = year, y = discountValue_cum, colour = scenario))+
  scale_x_continuous(limit= c(2022, 2050)) + 
  ggtitle("Retreat")

cost_stage$`status quo`$costRetreat
diag <- list("Ab test"= c_dab, "RNA test" = c_drna, "RNAonly/POCT" = c_drnaonly)

diag_p <- list()
for(i in 1: length(diag)){ 
  diag_p[[i]] <- ggplot(data = diag[[i]]) + 
    geom_line(aes( x = year, y = discountValue_cum, colour = scenario))+
    scale_x_continuous(limit= c(2022, 2050)) + 
    ggtitle(names(diag)[i])
    
}

diag_p[[1]]

ggplot(data = xxt) + 
  geom_line(aes( x = year, y = discountValue_cum, colour = scenario)) + 
  facet_wrap(.~ state, scale = "free") + 
  scale_x_continuous(limit= c(2022, 2050))
save(cost_np, x_np,cost_stage, c_stage,
     file = file.path(OutputFolder ,
                      paste0(project_name,"cost_np" ,".rda")))

cost_np$`status quo` <- cost_np$`status quo`%>%
  mutate(scenario = "Status quo")
##### cost-saving 
cost_np$`status quo` <- cost_np$`status quo`%>%
  mutate(sq_cumval = cost_np$`status quo`$discountValue_cum)%>%
  mutate(c_saving = sq_cumval - discountValue_cum, 
         scenario = "Status quo")

cost_np$`National program` <- cost_np$`National program`%>%
  mutate(sq_cumval = cost_np$`status quo`$discountValue_cum)%>%
  mutate(c_saving = sq_cumval - discountValue_cum, 
         scenario = "National program")


cost_np$`National program scale up` <- cost_np$`National program scale up`%>%
  mutate(sq_cumval = cost_np$`status quo`$discountValue_cum)%>%
  mutate(c_saving = sq_cumval - discountValue_cum, 
         scenario = "National program scale up")

cost_saving <- rbind(cost_np$`National program`, 
                     cost_np$`National program scale up`)



cost_saving_p <- cost_saving%>%
  filter(year%in% c(2030, 2040, 2050))%>%
           ggplot(., aes(x = as.character(year), y = c_saving)) +
  geom_histogram(stat='identity', aes(fill = scenario), 
           position="dodge") + 
  scale_y_continuous(limits = c(-15000000, 15000000), 
                     breaks = seq(-15000000, 15000000, 1000000),
                     labels = paste0(seq(-15,15,1), "M"))+ 
  theme_linedraw() +
  
  scale_fill_brewer(palette="Paired") +
  labs(x = "Year", y = "Cost saving (AUD)")

ggsave(path = OutputFig, file="cost_saving.png", cost_saving_p, 
       height = 6, width = 6, dpi = 800)


# total cost 
cost_total <- rbind(cost_np$`status quo`, 
                    cost_np$`National program`, 
                     cost_np$`National program scale up`)
cost_total_p <- cost_total%>%
  filter(year%in% c(2030, 2040, 2050))%>%
  ggplot(., aes(x = as.character(year), y = discountValue_cum)) +
  geom_histogram(stat='identity', aes(fill = scenario), 
                 position="dodge") + 
  scale_y_continuous(limits = c(0, 700000000), 
                     breaks = seq(0, 700000000, 50000000),
                     labels = paste0(seq(0,700,50), "M"))+ 
  theme_linedraw() +

  scale_fill_brewer(palette="Paired") +
  labs(x = "Year", y = "Cumulative cost (AUD)")

ggsave(path = OutputFig, file="cost_saving.png", cost_saving_p, 
       height = 6, width = 6, dpi = 800)

# cost of management/ treatment/ management 
# three scenario 
# stacked 
# also extract 2030 2040 2050 

c_stage_p <- c_stage%>%filter(year %in% c(2030,2040,2050))%>%
  mutate(`cost indicator` = factor(`cost indicator`, 
                                   levels = c("diagnosis", "Management", "Treatment"),
                                   labels = c("Diagnosis", "Managment", "Treatment")))%>%
 ggplot(data = ., aes(x = year, y = discountValue_cum)) + 
  geom_histogram(stat='identity', aes(fill = scenario), 
                 position="dodge") + 
  facet_wrap(.~`cost indicator`, scale ="free") +
  theme_linedraw() +
  scale_fill_brewer(palette="Paired") +
  labs(x = "Year", y = "Cumulative cost (AUD)") 

c_stage_p <- c_stage_p + 
  facet_custom (~`cost indicator`,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 120000000), 
                                                 breaks = seq(0, 120000000, 
                                                              10000000),
                                                 labels = paste0(seq(0, 120, 10), "M"))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 200000000), 
                                                 breaks = seq(0, 200000000, 
                                                              10000000),
                                                 labels = paste0(seq(0, 200, 10), "M"))),
                    
                    scale_new(3,
                              scale_y_continuous(limits =  c(0, 400000000), 
                                                 breaks = seq(0, 400000000, 
                                                              50000000),
                                                 labels = paste0(seq(0, 400, 50), "M")))))

ggsave(path = OutputFig, file="cost_stage.png", c_stage_p, 
       height = 6, width = 9, dpi = 800)



# advanced liver disease cost 
ad_cost <- xxt%>%
  group_by(year, scenario, population, disease_prog)%>%
  summarize(discountValue = sum(discountValue),
            discountValue_cum = sum(discountValue_cum))%>%ungroup()%>%
  filter(disease_prog %in% c("dc", "hcc", "lt", "plt"))

ad_cost_p <- ad_cost%>%filter(year %in% c(2030,2040,2050))%>%
  mutate(disease_prog = factor(disease_prog, 
                               levels = c("dc", "hcc", "lt", "plt"),
                               labels = c("Decompensated cirrhosis", 
                                          "Hepatocellular carcinoma", 
                                          "Liver transplant", 
                                          "Post-liver transplant")))%>%
  ggplot(data = ., aes(x = year, y = discountValue_cum)) + 
  geom_histogram(stat='identity', aes(fill = scenario), 
                 position="dodge") + 
  facet_wrap(.~disease_prog, scale ="free") +
  theme_linedraw() +
  scale_fill_brewer(palette="Paired") +
  labs(x = "Year", y = "Cumulative cost (AUD)") 

ggsave(path = OutputFig, file="ad_costsaving.png", ad_cost_p, 
       height = 6, width = 8, dpi = 800)
# cost saving in the advanced liver diseases 
