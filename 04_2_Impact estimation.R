# Cascade% calculation 

rm(list = ls())

library(data.table)
library(formattable)
library("readr")
library("dplyr")
library("tidyr")
library("purrr")
library(stringr) # for str_squish()
library(flextable)
library(officer)
library(forcats)
library(expss)
library(gt)
library(tidyverse)



basePath <- getwd()
projectFolder <- file.path(basePath)
project_name <- "HCVModel"
DataFolder <- file.path(basePath, "01. DATA/model input" )
ResultsFolder <- file.path(basePath, "04. Output/Results")
Rcode <- file.path(basePath, "03. Code")
  # if TRUE doesn't append time to results and overwrites
# a base file. Useful for storing main results or
# testing
projectFile <- file.path(basePath,paste0(project_name, ".rda"))
projectVars <- load(projectFile)
source(file.path(Rcode, "/Functions/HCV_model.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R"))
source(file.path(Rcode, "/Functions/plotManuscript.R"))
load(file.path(projectFolder, paste0(project_name,"SceRestbl",".rda")))
load(file.path(projectFolder, 
               paste0(project_name,"SqRestbl",".rda")))

# combine in the same data frame 
# rearrange columns 
coName <- colnames(SceRestbl$tbl$tarpop)
coNameORD <- c(coName[1], coName[13], coName[2], coName[12])
coName <- coName[!coName%in% coNameORD]
coName <- c(coNameORD, coName)
SenscenImp <-do.call("rbind", SceRestbl$tbl)%>%
  mutate(year = year + HCV$cabY - 1)%>%
  select(coName)

# aligned col with scenarios 
SqName <- colnames(SqRestbl$tbl)
SqnoN <- coName <- coName[!coName%in%SqName]
names(SqnoN) <- SqnoN
ta <- SqRestbl$tbl%>%mutate(!!sym(SqnoN) := "Status Quo", 
                            year = year  + HCV$cabY - 1) # to pass string in dplyr !!sym(string):= 

ttbl <- rbind(ta, SenscenImp)

SenImpact <- ttbl%>%
  mutate(SensiSce = factor(SensiSce, levels = c("Status Quo", "main","timeL", "timeU", 
                                                "diagL", "diagU", "reinfL", 
                                                "reinfU", "tarpop"), 
                           labels = c("Status Quo", "Base case", 
                                      "Pessimistic scale-up (5 years)",
                                      "Optimized scale-up (1 year)",
                                      "Pessimistic impact on HCV diagnosis",
                                      "Optimistic impact on HCV diagnosis",
                                      "No reinfection", 
                                      "Higher reinfection rate (1.5*primary infection)",
                                      "Prioritized MSM subgroups who regularly engaged with HIV prevention and care services")),
         testing = factor(testing, 
                          levels = c("Status Quo", "POC_antibody",
                                     "DBS", "Reflex_RNA", 
                                     "POC_RNA"),
                          labels = c("Status Quo", "POC antibody test",
                                     "Dried blood spot", "Reflex RNA", 
                                     "POC RNA")))%>%
  group_by(SensiSce, testing)%>%
  mutate(best = format(round(best, 1), nsmall = 1, trim = TRUE),
         q5 = format(round(q5,1), nsmall = 1, trim = TRUE),
         q95 = format(round(q95,1),nsmall = 1, trim = TRUE))%>%
  mutate(val = paste(best,"(", q5 ,"-", q95, ")", 
                        sep = ""))%>%
  select(year, SensiSce, testing, indicator, val)%>%spread(indicator, val)

SenImpact <- SenImpact%>%arrange(SensiSce, testing)%>%filter(year == 2030)%>%
  select(SensiSce, testing, Prevalence, Incidence, CHCVNum, InciNum, HCVdeath)


# combine with cascade coverage 
Scename <- unique(SceRestbl$tbl$tarpop$testing)

testcoverage <- SqRestbl$cas

scecoverage <- SceRestbl$cas

cascade <- c(names(testcoverage))

SqCover <- list()
SqCover <- lapply(seq_along(testcoverage), function(x){ 
  a <- testcoverage[[x]]%>%as.data.frame()%>%
    mutate(scenario = "Status Quo",cascade = cascade[x])
})
SqCover <-do.call("rbind", SqCover)
SqCoverplot <- SqCover%>%as.data.frame()%>%
  select(year, scenario, cascade, 
         c(cum_best,min, max, Med,Mu, q5, q25, q75, q95))%>%
  filter(year == 2030)
SenSceCover <- list()
for(i in names(scecoverage)){
  
  SenSceCover[[i]] <- lapply(seq_along(scecoverage[[i]]), function(x){
    m <- lapply(seq_along(scecoverage[[i]][[1]]), function(y){
      a <- scecoverage[[i]][[x]][[y]]%>%as.data.frame()%>%
        mutate(scenario = Scename[x],cascade = cascade[y])
    })
    aa <- do.call("rbind", m)
  })  
  
  SenSceCover[[i]] <- do.call("rbind", SenSceCover[[i]])
  
  
}
casplot <- list()
casplot <- lapply(seq_along(SenSceCover), function(x){ 
  
  a <- SenSceCover[[x]]%>%as.data.frame()%>%
    select(year, scenario, cascade, 
           c(cum_best,min, max, Med,Mu, q5, q25, q75, q95))%>%
    filter(year == 2030)
  
  a <- rbind(SqCoverplot, a)
})


### reorder dataset 
casF <- unique(casplot[[1]]$cascade)
sceF <- unique(casplot[[1]]$scenario)
names(SenSceCover)
casplot <- lapply(seq_along(SenSceCover), function(x){
  
  a <- casplot[[x]]%>%
    mutate(cascade = factor(cascade, levels =c(casF),
                            labels = c("Infected", "Diagnosed",
                                       "Treatment initiated","Cured")),
           scenario = factor(scenario, levels = c(sceF),
                             labels = c("Status Quo", "POC antibody test",
                                        "Dried blood spot", "Reflex RNA",
                                        "POC RNA")),
           SensiSce = names(SenSceCover)[x]
                             
           )
    
})

ScenaCover <- do.call("rbind", casplot)%>%filter(cascade!="Infected")%>%
  select(year, scenario, cascade, cum_best, SensiSce)%>%
  mutate(cum_best =format(round(cum_best, 1), nsmall = 1, trim = TRUE))%>%
  spread(key = cascade, value =cum_best)%>%
  mutate(Cascade = paste0("(", Diagnosed,"%" ,", ", `Treatment initiated`, "%",
                          ", ", Cured, "%", ")"))%>%
  mutate(SensiSce = factor(SensiSce, levels = c("Status Quo", "main","timeL", "timeU", 
                                                "diagL", "diagU", "reinfL", 
                                                "reinfU", "tarpop"), 
                           labels = c("Status Quo", "Base case", 
                                      "Pessimistic scale-up (5 years)",
                                      "Optimized scale-up (1 year)",
                                      "Pessimistic impact on HCV diagnosis",
                                      "Optimistic impact on HCV diagnosis",
                                      "No reinfection", 
                                      "Higher reinfection rate (1.5*primary infection)",
                                      "Prioritized MSM subgroups who regularly engaged with HIV prevention and care services")),
         testing = scenario)%>%
  select(!c( Diagnosed, `Treatment initiated`, Cured, scenario, year))

# clean status quo 

ScenaCoversq <- ScenaCover%>%filter(testing =="Status Quo")
ScenaCoversq$SensiSce <- ScenaCoversq$testing
ScenaCoversen <- ScenaCover%>%filter(testing !="Status Quo") 

ScenaCover <- rbind(ScenaCoversq[1,] , ScenaCoversen)

SenImpactX <- merge(SenImpact, ScenaCover, by = c("testing", "SensiSce"))%>%
  arrange(SensiSce, testing)

tib <- as_grouped_data(x = SenImpactX, groups = c("SensiSce"), columns = NULL)
tib%>%flextable::as_flextable()

a <- tib%>% 
  flextable::as_flextable( ) %>% 
  flextable::compose(
    i = ~ !is.na(SensiSce), # when var_group not NA
    j = "testing", # on column "var"
    # create a paragraph containing a chunk containing value of `var_group`
    value = as_paragraph(as_chunk(SensiSce))) %>% 
  hline(i = ~ !is.na(SensiSce), border = officer::fp_border() ) %>% 
  autofit()%>%
  set_header_labels(., Prevalence = "Prevalence, %",
                    Incidence = "Incidence, per 1000PY*",
                    CHCVNum = "Number of people living with chronic HCV",
                    InciNum = "Number of annual new infections",
                    HCVdeath = "Number of annual HCV-related deaths",
                    Cascade = "Cascade \n(Diagnosed, Treatment initiated, Cured)",
                    testing = "POCT scenario")%>%
  bold( bold = TRUE,part=c("all"))
str(a)
# arrange scenario orders                              


#### table #### 

str(SceRestbl$tbl$tarpop)
# fig 3 impact of poct in different scenarios 
SenCuInf <- list()
for( i in names(SceRestbl$tbl)){  
  SenCuInf[[i]] <- SceRestbl$tbl[[i]]%>%filter(indicator =="InciNum")%>%
    mutate(year = year + HCV$cabY - 1 )%>%filter(year>=2022)%>%
    group_by(testing)%>%
    mutate(across(c(best, Mu, Med, max, min, q5, q25, q75, q95), cumsum))%>%
    filter(year == 2030)
}
tt <- SqRestbl$tbl%>%filter(indicator =="InciNum")%>%
  mutate(year = year + HCV$cabY - 1)%>%filter(year>=2022)%>%
  mutate_at(3:ncol(.) ,cumsum)%>%filter(year == 2030)
cumInf <- paste0(format(round(tt$best,0),nsmall = 0, trim = TRUE), 
                 " (", format(round(tt$q5,0),nsmall = 0,trim = TRUE), 
                 " - ", 
                 format(round(tt$q95,0),nsmall = 0, trim = TRUE),
                 ")")
senR <- list()
senRR <- list()
senRRlow <- list()
senRRmed <- list()
senRRup <- list()
senMed <- list()
senIQR <- list()
senQ1 <- list()
senQ3 <- list()
testing = c(unique(SceRestbl$tbl$main$testing))
for(m in names(SceRestbl$tbl)){
  for(i in testing){ 
    
    senRR[[m]][[i]] <- (SenCuInf[[m]][SenCuInf[[m]]$testing ==i, "best"] - tt[1, "best"])/tt[1, "best"]*100
    senRRlow[[m]][[i]] <- (SenCuInf[[m]][SenCuInf[[m]]$testing ==i, "q25"] - tt[1, "q25"])/tt[1, "q25"]*100
    senRRup[[m]][[i]] <- (SenCuInf[[m]][SenCuInf[[m]]$testing ==i, "q75"] - tt[1, "q75"])/tt[1, "q75"]*100
    senRRmed[[m]][[i]] <- (SenCuInf[[m]][SenCuInf[[m]]$testing ==i, "Med"] - tt[1, "Med"])/tt[1, "Med"]*100
    senR[[m]][[i]] <- paste0(format(abs(senRR[[m]][[i]]),nsmall = 1, trim = TRUE), "%", sep ="")
    
    senMed[[m]][[i]] <- SenCuInf[[m]][SenCuInf[[m]]$testing ==i, "Med"]
    
    senQ1[[m]][[i]] <- SenCuInf[[m]][SenCuInf[[m]]$testing ==i, "q25"]
    
    senQ3[[m]][[i]] <- SenCuInf[[m]][SenCuInf[[m]]$testing ==i, "q75"]
    
    senIQR[[m]][[i]] <- paste0(SenCuInf[[m]][SenCuInf[[m]]$testing ==i, "q25"], 
                               " - ", 
                               SenCuInf[[m]][SenCuInf[[m]]$testing ==i, "q75"] )
    
    
  } }

ReduRang <- list()
minmaxred <- list()
for(i in testing){
  
  ReduRang[[i]] <- cbind(Sen = c(names(senRR)),
                         red = rbind(senRR[[1]][[i]],
                                     senRR[[2]][[i]],
                                     senRR[[3]][[i]],
                                     senRR[[4]][[i]],
                                     senRR[[5]][[i]],
                                     senRR[[6]][[i]],
                                     senRR[[7]][[i]],
                                     senRR[[8]][[i]]))
  
  minmaxred[[i]] <- c(min(ReduRang[[i]]$best), max(ReduRang[[i]]$best))
}

reSenName <- c("main","timeL", "timeU","diagL", "diagU", "reinfL", 
               "reinfU", "tarpop") 

reSenNamelv <- c("Base case", 
                 "Pessimistic scale-up (5 years)",
                 "Optimized scale-up (1 year)",
                 "Pessimistic impact on HCV diagnosis",
                 "Optimistic impact on HCV diagnosis",
                 "No reinfection", 
                 "Higher reinfection rate (1.5*primary infection)",
                 "Targeted MSM engaged in HIV services")
SenCUMIQR <- list()

for(i in reSenName){ 
  
  SenCUMIQR[[i]] <- c(paste0(
    format(round(senMed[[i]][[1]]$Med,0),nsmall = 0, trim = TRUE), " (",
    format(round(senQ1[[i]][[1]]$q25,0),nsmall = 0, trim = TRUE), " - ", 
    format(round(senQ3[[i]][[1]]$q75,0),nsmall = 0, trim = TRUE), ")"),
    paste0(
      format(round(senMed[[i]][[2]]$Med,0),nsmall = 0, trim = TRUE), " (",
      format(round(senQ1[[i]][[2]]$q25,0),nsmall = 0, trim = TRUE), " - ", 
      format(round(senQ3[[i]][[2]]$q75,0),nsmall = 0, trim = TRUE), ")"),
    paste0(
      format(round(senMed[[i]][[3]]$Med,0),nsmall = 0, trim = TRUE), " (",
      format(round(senQ1[[i]][[3]]$q25,0),nsmall = 0, trim = TRUE), " - ", 
      format(round(senQ3[[i]][[3]]$q75,0),nsmall = 0, trim = TRUE), ")"),
    paste0(
      format(round(senMed[[i]][[4]]$Med,0),nsmall = 0, trim = TRUE), " (",
      format(round(senQ1[[i]][[4]]$q25,0),nsmall = 0, trim = TRUE), " - ", 
      format(round(senQ3[[i]][[4]]$q75,0),nsmall = 0, trim = TRUE), ")")
    
  )
  
} 
testingName = c("Point-of-care antibody testing","Dried blood spot testing", 
                "Clinic-based reflex RNA testing", "Single visit Point-of-care RNA testing")
table <- list()
for(i in reSenName){
  table[[i]] <- cbind(Scenarios = c(i, testingName),
                      `Cumulative new infections (IQR)` = 
                        c("", SenCUMIQR[[i]]),
                      
                      `Reduction (%)` = c("", c(do.call("rbind",senRRmed[[i]])$Med)),
                      
                      `CI low` = c("", c(do.call("rbind",senRRlow[[i]])$q25)),
                      
                      `CI high` = c("", c(do.call("rbind",senRRup[[i]])$q75)))%>%
    as.data.frame()
  
}

for(i in seq_along(reSenName)){ 
  table[[i]]$Scenarios <- c(reSenNamelv[i], c(testingName))
  
  
}


table<-  do.call("rbind", table)


# indent the subgroup if there is a number in the placebo column
table <- table%>%
  mutate(
    `Reduction (%)` = round(as.numeric(`Reduction (%)`),1), 
    `CI low` = round(as.numeric(`CI low`),1),
    `CI high` = round(as.numeric(`CI high`),1)
  )
table$Scenarios <- ifelse(is.na(table$`Reduction (%)`), table$Scenarios,
                          paste0("   ", table$Scenarios))

# remove indent of the first row
table
# use forester to create the table with forest plot
a <- foresterR(left_side_data = table[,1:2],
               estimate = table$`Reduction (%)`,
               estimate_col_name = "Reduction (%, IQR)",
               ci_low = table$`CI low`,
               ci_high = table$`CI high`,
               display = TRUE,
               file_path = ("forester_plot.png"),
               font_family = "sans",
               null_line_at = 0,
               xlim = c(-70, 40),
               xbreaks = seq( -70, 40, 10),
               render_as = "png") + 
  scale_x_continuous(limits = c(-70, 40), breaks = seq(-70,40,10),
                     labels = c("-70%", "-60%","-50%","-40%","-30%","-20%","-10%",
                                "0%","10%","20%","30%","40%"))





...#### cascade plot #### 
cascade_dt <- ScenaCover%>%
  mutate(Infected = 100,
         Diagnosed = as.numeric(substr(Cascade, 2, 5)),
         Treatment_initiated = as.numeric(substr(Cascade, 8, 12)),
         Cured = as.numeric(substr(Cascade, 15, 19)))
cascade_dtmain <- cascade_dt%>%filter(SensiSce%in%c("Status Quo", "Base case"))%>%
  select(!c("SensiSce", "Cascade"))%>%
  gather(cascade, value, -(testing))%>%
  mutate(cascade = factor(cascade, levels = c("Infected", "Diagnosed",
                                              "Treatment_initiated","Cured"), 
                          labels = c("Infected", "Diagnosed",
                                     "Treatment initiated","Cured")),
         scenario = factor(testing, levels = c("Status Quo", "POC antibody test",
                                               "Dried blood spot", "Reflex RNA",
                                               "POC RNA"),
                           labels = c("Status Quo", "POC antibody test",
                                      "Dried blood spot", "Reflex RNA",
                                      "POC RNA"),
         )
  )
fig4 <-ggplot(cascade_dtmain , aes(fill = scenario, y = value, x = cascade)) + 
  geom_bar(position = "dodge", stat = "identity",aes(fill = scenario)) + 
  labs(y = "%") + 
  scale_fill_OkabeIto(name  ="Testing scenarios",
                      breaks=c(unique(casplot[[1]]$scenario)),
                      labels=c( 
                        "Status Quo",
                        "POC antibody test", 
                        "Dried blood spot", 
                        "Reflex RNA", 
                        "POC RNA")) + plotOpts + 
  theme(axis.text.x = element_text(face = "bold", 
                                   angle = 0,
                                   size=14,
                                   colour="black", hjust = 0.5, vjust = 0.5)) + 
  geom_text(aes(label= paste0(round(value, digit = 0), "%", sep = "")), 
            position=position_dodge(width=0.9), vjust=-0.25, 
            
            size = 5) + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 110), 
                     breaks = seq(0, 100,10),
                     labels = c(paste0(seq(0, 100,10), "%", sep = ""))) 
ggsave(file="mainFig4.png", fig4, width = 15 , height = 9)
#### dataset for boxplot of cumulative new infection 
# status quo 

HCVInfect_all <- indicatorResults(HCV, bestResults, "newInfections", pop="all",
                                  paramR = paramResults, range = "y",
                                  endY = 30)%>%mutate(year = year + HCV$cabY - 1)%>%
  filter(year >=2022)
sqcum <- HCVInfect_all%>%
  mutate(across(c(best, paste0("set", seq(1, HCV$numberSamples,1))),
                cumsum))%>%
  select(c(year, best, paste0("set", seq(1, HCV$numberSamples,1))))%>%filter(year == 2030)

sqcum <- rbind(sqcum, sqcum, sqcum, sqcum)


NumInfcum <- HCVInfect_all%>%filter(year>=2021)%>%
  mutate_if(is.numeric, funs(c(first(.), (. - first(.))[-1]))) %>%
  filter(year==10)

bestcolsq <- which(colnames(NumInfcum) == "cum_best")
popSizeQuansq <- NumInfcum%>%gather(., "simulation", "estimate", 
                              bestcol:ncol(NumInfcum))%>%
  filter(simulation!= "cum_best")

popQsq <- popSizeQuansq%>%group_by(year) %>%
  summarise(min = min(estimate, na.rm = TRUE),
            max = max(estimate, na.rm = TRUE),
            Med = median(estimate, na.rm = TRUE),
            Mu = mean(estimate, na.rm = TRUE),
            q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
            q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
            q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
            q95 = quantile(estimate, prob = 0.975, na.rm = TRUE),
            lw = quantile(estimate, prob = 0.25, na.rm = TRUE) - 
              1.5*((quantile(estimate, prob = 0.75, na.rm = TRUE) - 
                      quantile(estimate, prob = 0.25, na.rm = TRUE)) ),
            uw = quantile(estimate, prob = 0.75, na.rm = TRUE) + 
              1.5*((quantile(estimate, prob = 0.75, na.rm = TRUE) - 
                      quantile(estimate, prob = 0.25, na.rm = TRUE)) )
  )%>%
  ungroup()%>%arrange(year)%>% 
  select(everything())
NumInfcum[, c("min", "max", "Med", "Mu", 
          "q5", "q25", "q75", "q95", "lw", "uw")] <- 
  popQsq%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95, lw, uw)) 


NumInfcum$q25

ttt <- SceResAvert$main[[1]]
ttt <- ttt%>%arrange(year, testing)%>%mutate( year = year + HCV$cabY - 1)


tttaf <- ttt%>%group_by(testing)%>%filter(year>=2021)%>%
  mutate_if(is.numeric, funs(c(first(.), (. - first(.))[-1]))) %>%
  filter(year == 10)


bindAvert <- cbind(year =tttaf$year, testing = tttaf$testing, as.data.frame(
  tttaf[ ,-c(1,2)] - sqcum[, -1]))

sq <- NumInf%>%filter(testing =="Status Quo")
sq <- rbind(sq, sq,sq,sq,sq)
scenn <- cbind(year = NumInf$year, testing = NumInf$testing,
               as.data.frame(NumInf[ ,-c(1,2)] - sq[ ,-c(1,2)]))
              


tttff <- list()
for(i in unique(tttaf$testing)){ 
  tttff[[i]] <- tttaf%>%filter(testing == i) 
  
  
  tttff[[i]] <- cbind(year = tttff[[i]]$year, testing = tttff[[i]]$testing, 
                      as.data.frame(tttff[[i]][ , -c(1,2)] - NumInfcum[ , -c(1)], na.rm = TRUE))

  }


tttal <- do.call("rbind", tttff)%>%arrange(year, testing)
colnames(tttal)


bestcol <- which(colnames(tttaf) == "cum_best")
popSizeQuan <- tttaf%>%gather(., "simulation", "estimate", 
                             bestcol:ncol(tttaf))%>%
    filter(simulation!= "cum_best")

    popQ <- popSizeQuan%>%group_by(year, testing) %>%
      summarise(min = min(estimate, na.rm = TRUE),
                max = max(estimate, na.rm = TRUE),
                Med = median(estimate, na.rm = TRUE),
                Mu = mean(estimate, na.rm = TRUE),
                q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
                q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
                q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
                q95 = quantile(estimate, prob = 0.975, na.rm = TRUE),
                lw = quantile(estimate, prob = 0.25, na.rm = TRUE) - 
                  1.5*((quantile(estimate, prob = 0.75, na.rm = TRUE) - 
                          quantile(estimate, prob = 0.25, na.rm = TRUE)) ),
                uw = quantile(estimate, prob = 0.75, na.rm = TRUE) + 
                  1.5*((quantile(estimate, prob = 0.75, na.rm = TRUE) - 
                          quantile(estimate, prob = 0.25, na.rm = TRUE)) )
                )%>%
      ungroup()%>%arrange(year, testing)%>% 
      select(everything())
    tttaf[, c("min", "max", "Med", "Mu", 
           "q5", "q25", "q75", "q95", "lw", "uw")] <- 
      popQ%>%select(., c(min, max, Med, Mu, q5, q25, q75, q95, lw, uw)) 
str(tttal)
tttaf%>%filter(year == 10)%>%ggplot(data = ., aes(x = testing,
        ymin = min,
        lower = q25,
        middle = Mu,
        upper = q75,
        ymax = max)) +
  geom_boxplot(stat = "identity") + xlab("Year") + ylab("Cumulative new infections") + 
  scale_y_continuous(limits = c(0,60000), breaks = seq(0, 60000, 5000)) +
  geom_hline(yintercept = NumInfcum$q25, linetype ="dashed") + 
  geom_hline(yintercept = NumInfcum$q75, linetype ="dashed")


cumboxplot <- list()
SenSceTag <-c("Status Quo", "Base case",
              "Pessimistic scale-up (5 years)",            
              "Optimized scale-up (1 year)",
              "Pessimistic impact on HCV diagnosis",
              "Optimistic impact on HCV diagnosis",
              "No reinfection", 
              "Higher reinfection rate (1.5*primary infection)",
              "Prioritized MSM PrEP users")

cumboxplot <- lapply(seq_along(SenSceCum), function(x){  
  labelN = list(c("Testing scenarios"), c("Status Quo", 
                                          "POC antibody test", 
                                          "Dried blood spot", 
                                          "Reflex RNA", 
                                          "POC RNA"))
  a <-SenSceCum[[x]]%>%filter(year == 2030)%>%ggplot(.,                              # Draw ggplot2 boxplot
       aes(x = testing,
           ymin = lower,
           lower = q25,
           middle = Med,
           upper = q75,
           ymax = upper)) +
  geom_boxplot(stat = "identity") + xlab("Year") + ylab("Cumulative new infections") + 
    plotOpts + 
    scale_color_OkabeIto(name  =labelN[[1]],
                         breaks=c(unique(SenSceCum[[x]]$testing)),
                         labels=labelN[[2]]) +
    scale_fill_OkabeIto(name  =labelN[[1]],
                        breaks=c(unique(SenSceCum[[x]]$testing)),
                        labels=labelN[[2]]) + 
    
    ggtitle(SenSceTag[x])
  })






try <- SceRestbl$tbl$tarpop%>%filter(indicator =="InciNum")%>%group_by(testing)%>%
  mutate(across(where(is.numeric),
                cumsum, .names = "cum_{col}"))%>%select(!cum_year)%>%
  mutate(year = year + HCV$cabY - 1)
trysq <- SqRestbl$tbl%>%filter(indicator =="InciNum")%>%group_by(testing)%>%
  mutate(across(where(is.numeric),
                cumsum, .names = "cum_{col}"))%>%select(!cum_year)%>%
  mutate(year = year + HCV$cabY - 1)


tcum <- rbind(trysq, try)%>%as.data.frame()%>%
  mutate(upper = cum_q75 + 1.5*(cum_q75-cum_q25), 
         lower = cum_q25 + 1.5*(cum_q75-cum_q25))



scenarioResults$Main