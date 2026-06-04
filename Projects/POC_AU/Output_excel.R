# this script is to generate the model outcomes to excel sheet to check the detail 
# 



library(openxlsx)
load(file.path(OutputFolder, paste0(project_name, "ResSum.rda")))
load(file.path(OutputFolder, paste0(project_name, "ResSum_cost.rda")))
load(file.path(OutputFolder, paste0(project_name, "Sce_sq.rda")))
load(file.path(OutputFolder, paste0(project_name, "Sce_np.rda")))

# components 
Sce <- list("status quo" = Sce_sq, 
            "National program" = Sce_np)

# compartment number and cost in each year 

Num_component_dt <- list()
Cost_component_dt <- list()
component <- list()
for(n in names(Sce)){ 
  Num_component_dt[[n]] <- modres.t(POC_AU, Sce[[n]], endYear = 100)%>%
    as.data.frame()%>%
    select(year, timestep, population, state, cascade, disease_prog, best) 
  
  Cost_component_dt[[n]] <- modres.t(POC_AU, Sce[[n]], endYear = 100, allp = "costPops")%>%
    as.data.frame()%>%
    select(year, timestep, population, state, cascade, disease_prog, best) 
  
  component[[n]] <- cbind(as.data.frame(Num_component_dt[[n]]), cost = Cost_component_dt[[n]]$best)%>%
    mutate(year = year + POC_AU$cabY - 1)
  
}

component_dt <- dplyr::bind_rows(component, .id = "scenario")

 

# cascade 
cascade_fun <- function(pj,dlst, cost = NULL){ 
  cascade_dt <- lapply(names(dlst), function(x) dlst[[x]]%>%
                         group_by(timestep, population, cascade)%>%
                         summarize(best = round(sum(best)), digits = 2)%>%
                         spread(cascade, best)%>%
                         mutate(year = timestep%/%1 - 1)%>%
                         select(year,timestep, s, undiag, diag_ab, 
                                diag_RNA, treat, 
                                treat_f, cured)%>%
                         mutate(year = pj$cabY + year))
  
  if(!is.null(cost)){ 
    
  cascade_dt <- lapply(cascade_dt, setNames,c("population", "year","timestep",  "cost_s", 
                                              "cost_undiag", "cost_diag_ab", 
                                              "cost_diag_RNA", "cost_treat", 
                                              "cost_treat_f", "cost_cured"))
    
    
      
  }
  names(cascade_dt) <- names(dlst)
  cascade_dt <-   dplyr::bind_rows(cascade_dt, .id = "scenario")
  
  return(cascade_dt)
}

Num_cascade_dt <- cascade_fun(POC_AU, Num_component_dt)
Cost_cascade_dt <- cascade_fun(POC_AU, Cost_component_dt, cost= "y")


cascade_dt <- cbind(as.data.frame(Num_cascade_dt), as.data.frame(Cost_cascade_dt[, 5:11]))
component_dt$disease_prog
# disease_prog
# cascade 
dis_prog_fun <- function(pj,dlst, cost = NULL){ 
  dis_prog_dt <- lapply(names(dlst), function(x) dlst[[x]]%>%
                         group_by(timestep, population, disease_prog)%>%
                         summarize(best = round(sum(best)), digits = 2)%>%
                         spread(disease_prog, best)%>%
                         mutate(year = timestep%/%1 - 1)%>%
                         select(year,timestep, unique(dlst[[x]]$disease_prog))%>%
                         mutate(year = pj$cabY + year))
  
  if(!is.null(cost)){ 
    
    dis_prog_dt <- lapply(dis_prog_dt, setNames,c("population", "year","timestep",  
                                                 c(paste0("cost_", c(unique(dlst[[1]]$disease_prog))))))
    
    
    
  }
  names(dis_prog_dt) <- names(dlst)
  dis_prog_dt <-   dplyr::bind_rows(dis_prog_dt, .id = "scenario")
  
  return(dis_prog_dt)
}

Num_disprog_dt <- dis_prog_fun(POC_AU, Num_component_dt)
Cost_disprog_dt <- dis_prog_fun(POC_AU, Cost_component_dt, cost= "y")


disprog_dt <- cbind(as.data.frame(Num_disprog_dt), as.data.frame(Cost_disprog_dt[, 5:15]))

# flows
casflow.index <- c("newTestingAb", "newTestingAg", "newTestingPOCT",
                   "newTreatment", "newRetreat",  "newCured",
                   "costTestingAb", "costTestingAg","costTestingPOCT",  
                   "costTreatment", "costRetreat",  "costCured")

cas_flow_dt <- list()
for(n in names(Sce)){ 
  for(i in casflow.index){
    cas_flow_dt[[n]][[i]] <- 
      modres.flow.t(POC_AU, Sce[[n]], allp = i, endYear = 100)%>%
      as.data.frame()
    }
  cas_flow_dt[[n]] <- dplyr::bind_rows(cas_flow_dt[[n]], .id = "indicator")
}

cas_flow_dt <- dplyr::bind_rows(cas_flow_dt, .id = "scenario")%>%
  mutate(indicator = factor(indicator, levels = unique(cas_flow_dt$`status quo`$indicator)))

# spread indicators to columns 
cas_flow_dt <- cas_flow_dt%>%spread(indicator, best)

# newinf, newhcvdeath 
epiflow.index <- c("newInfections", "newHCVdeaths", "newreinfection")
epi_flow_dt <- list()
for(n in names(Sce)){ 
  for(i in epiflow.index){
    epi_flow_dt[[n]][[i]] <- 
      modres.flow.t(POC_AU, Sce[[n]], allp = i, endYear = 100)%>%
      as.data.frame()
  }
  epi_flow_dt[[n]] <- dplyr::bind_rows(epi_flow_dt[[n]], .id = "indicator")
}
unique(epi_flow_dt$`status quo`$indicator)
epi_flow_dt <- dplyr::bind_rows(epi_flow_dt, .id = "scenario")%>%
  mutate(indicator = factor(indicator, levels = unique(epi_flow_dt$`status quo`$indicator)))

# spread indicators to columns 
epi_flow_dt <- epi_flow_dt%>%spread(indicator, best)



modres_lst <- list("compartments" = component_dt, "cascade" = cascade_dt,
                   "disprog" = disprog_dt, 
                   "cascade flows" = cas_flow_dt,
                   "epi flows" = epi_flow_dt
                   )

write.xlsx(modres_lst, file = file.path(OutputFolder, paste0("num_component.xlsx")))


# yearly number 

#epi: prevalence, incidence 
#epi: number of advanced liver diseases 
tempPrevRNA_subpop <- dplyr::bind_rows(tempPrevRNA_subpop, .id = "scenario")
HCVInc_subpop <- dplyr::bind_rows(HCVInc_subpop, .id = "scenario")

epi <- rbind(tempPrevRNA_subpop, HCVInc_subpop)%>%
  mutate( year = year + POC_AU$cabY - 1)%>%spread(indicator, best)

ind <- lapply(Sce_flow_sub, function(x){ 
  a <- x[!names(x)%in%
           c("newS", "newEntry", "newDeath", "newLeave",
             "inflow", "outflow", "costTestingAg")]
  
  a <- dplyr::bind_rows(a, .id = 'indicators')%>%
    select(scenario, year, population, indicators, best)
  })

ind$`status quo`%>%filter(year %in% c(7:10))

epi.flow.sub <- dplyr::bind_rows(ind, .id = 'scenario')%>%
  select(scenario, year, population, indicators, best)%>%
  mutate(year = year + POC_AU$cabY - 1)%>%
  spread(indicators, best)
  

write.xlsx(preinc, file = file.path(OutputFolder, paste0("preinc.xlsx")), 
           append=TRUE)



# cost: 
# total cost 
# components cost + cost_testingab, RNA, poct, treatment, retreatment, cured 

box_cost <- component_dt%>%group_by(scenario, year, population)%>%
  summarise(cost = sum(cost))

flow_cost <- cas_flow_dt%>%mutate(year = lag(year,POC_AU$npops))%>%
  select(scenario, year, population, costTestingAb,costTestingAg,costTestingPOCT,
           costTreatment,costRetreat,costCured)%>%
  ungroup()%>%group_by(scenario, year, population)%>%
  mutate(flow_cost = costTestingAb + costTestingAg + costTestingPOCT + 
           costTreatment + costRetreat + costCured)%>%na.omit()%>%
  summarise_at(.var = c("costTestingAb","costTestingAg","costTestingPOCT", 
                        "costTreatment","costRetreat","costCured", "flow_cost"),
               .funs = c(sum = sum
                         ))%>%
  mutate(year = year + POC_AU$cabY)

cost_all <- cbind(as.data.frame(flow_cost), costCompartment_sum = box_cost$cost)%>%
  mutate(cost_total = costCompartment_sum + flow_cost_sum)%>%
  mutate(id = year - POC_AU$simY)%>%
  mutate(discount = (1 + 0.05)^id)%>%filter(id>= 0)%>%
  mutate(discountValue = cost_total/discount)

for(n in names(Sce_cost)){
  for(s in names(x_np[[1]])){ 
    
    x_np[[n]][[s]] <- x_np[[n]][[s]]%>%
      mutate(year = year + POC_AU$cabY - 1)%>%
      mutate(id = year - POC_AU$simY)%>%
      mutate(discount = (1 + 0.05)^id)%>%filter(id>= 0)%>%
      mutate(discountValue = best/discount)%>%
      group_by(year)%>%summarise(discountValue = sum(discountValue),
                                 Value = sum(best))%>%ungroup()%>%
      mutate(discountValue_cum = cumsum(discountValue), 
             Value_cum = cumsum(Value))
  }
  }

x_np <- lapply(x_np, function(x) dplyr::bind_rows(x, .id = "categories")%>%
                 filter(categories != "cost_state"))


x_np_total <- lapply(x_np, function(x) x%>%filter(categories != "cost_state")%>%
                       group_by(year)%>%
                       summarise(discountValue = sum(discountValue),
                                 Value = sum(Value),
                                 Value_cum = sum(Value_cum), 
                                 discountValue_cum = sum(discountValue_cum))%>%
                       mutate(categories = "Total")%>%
                       select(categories, year, discountValue, Value, 
                              discountValue_cum, Value_cum))

cost_detail <- list()
for(i in 1:2){ 
  
  cost_detail[[i]] <- rbind(x_np_total[[i]], x_np[[i]])%>%
    pivot_wider(id_cols = year, 
                names_from = categories, 
                values_from = c("discountValue", "Value", "discountValue_cum", 
                                "Value_cum"))
  
  }


write.xlsx(cost_detail, file = file.path(OutputFolder, paste0("tota_cost.xlsx")), 
           append=TRUE) 



ROI <- ROI%>%mutate(Year = year, 
                    Scenario = factor(scenario, levels = c(unique(scenario)),
                                      labels = c("National Program", "scale up national program")),
                    `Cost (discounted, AUD)` = discountValue,
                    `Cumulative Cost (discounted, AUD)` = discountValue_cum,
                    `Return of Investment` = round(ROI, digits = 6))%>%
  select(-c(year, scenario, discountValue, discountValue_cum, ROI))

library(ggplot2)
ROI_p <- ggplot(data = ROI%>%filter(Year == 2100), aes(x = Scenario, y = `Return of Investment`, group = Scenario)) + 
  geom_line(aes(colour = Scenario)) + 
  scale_y_continuous(limits = c(-1, 1)) + 
  scale_color_manual(values = c("#A6CEE3", "#1F78B4")) +
  theme_bw()+
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") 


x_np$`status quo`$costTreatment%>%filter(year == 7)

indicatorResults(POC_AU, calibrateInit, "newTreatment", 
                 pop=POC_AU$popNames,
                 paramR = NULL, range = NULL,
                 endY = 100)%>%filter(year == 7
                                      )


