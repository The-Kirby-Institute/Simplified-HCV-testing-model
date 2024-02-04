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






modres_lst <- list("compartments" = component_dt, "cascade" = cascade_dt,
                   "disprog" = disprog_dt)

write.xlsx(modres_lst, file = file.path(OutputFolder, paste0("num_component.xlsx")))

write.xlsx(list_of_datasets, file = "writeXLSX2.xlsx")


# ad liver disease 

Num_dc_all <- lapply(Num_dc_all, function(x) x%>%
                       mutate(best = round(best, digits = 2))%>%
                       mutate(disease_prog = "DC")
                     )

Num_hcc_all <- lapply(Num_hcc_all, function(x) x%>%
                        mutate(best = round(best, digits = 2))%>%
                        mutate(disease_prog = "HCC")
)

Num_lt_all <- lapply(Num_lt_all, function(x) x%>%
                        mutate(best = round(best, digits = 2))%>%
                        mutate(disease_prog = "liver transplant")
)

Num_plt_all <- lapply(names(Num_component), function(x)
  Num_component[[x]]%>%
  filter(disease_prog%in%c("plt"))%>%
  group_by(year)%>%summarize(best = sum(best))%>%
  mutate(scenario = x, disease_prog = "post-liver transplant"))
names(Num_plt_all) <- names(Num_component)

num_ad <- list()
for(i in 1:3){ 
  
  num_ad[[i]] <- rbind(Num_dc_all[[i]], Num_hcc_all[[i]], Num_lt_all[[i]], 
                       Num_plt_all[[i]])%>%spread(disease_prog, best)%>%
    mutate(year = year + POC_AU$cabY - 1 )%>%
    select(!scenario)
  
  
  
  
  }


write.xlsx(num_ad, file = file.path(OutputFolder, paste0("num_ad.xlsx")), append=TRUE)  



tempPrevRNA_subpop <- lapply(tempPrevRNA_subpop, function(x) x%>%
                               mutate(indicator = "RNA prevlence"))

HCVInc_subpop <- lapply(HCVInc_subpop, function(x) x%>%
                               mutate(indicator = "HCV incidence"))

preinc <- list()
for(i in 1:3){ 
  
  preinc[[i]] <- rbind(tempPrevRNA_subpop[[i]], HCVInc_subpop[[i]])%>%
                          spread(indicator, best)%>%
    mutate(year = year + POC_AU$cabY - 1 )
  
  
  
  
}

write.xlsx(preinc, file = file.path(OutputFolder, paste0("preinc.xlsx")), 
           append=TRUE)

# number of flow 
ind_v <- c("newInfections", "newHCVdeaths",  "newTestingAb", "newTestingRNA", 
           "newTestingPOCT", "newTreatment", "newRetreat")

indflow <- lapply(Sce_flow_all, function(x){ 
  
  a <- x[names(x)%in%ind_v]
  for(n in ind_v){ 
    
    a[[n]] <- a[[n]]%>%mutate(indicator = n)
    
  }
  
  mm <- do.call("rbind", a)
  
  mm <- mm%>%
    spread(indicator, best)%>%
    mutate(year = year + POC_AU$cabY - 1 )
  return(mm)
}
  
  )
write.xlsx(indflow, file = file.path(OutputFolder, paste0("flow_num.xlsx")), 
           append=TRUE)


# cost: 
# total cost 

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
for(i in 1:3){ 
  
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


