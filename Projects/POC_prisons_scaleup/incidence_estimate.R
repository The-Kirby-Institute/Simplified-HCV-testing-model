#### incidence ##### 

name_parset <- c("best")
all_scenarios <- c(list("Status quo" = Sce_sq), Sce_np)
subpop_N <- list()
pop_N <- list()
total_N <- list()
commu_N <- list()
prison_N <- list()
prisonPWID_N <- list()
HCVInfect_subpop <- list()
HCVInfect_subpop_P <- list()
HCVInfect_subpop_C <- list()
pop_N_P <- list()
HCVInc_subpop_P <- list()
HCVInc_subpop_C <- list()
HCVInc_subpop <- list()
HCVInfect_setting <- list()
HCVInc_setting <- list()
for(i in names(all_scenarios)){ 
  subpop_N[[i]] <- lapply(POC_AU$popNames, function(x){ 
    
    a <- popResults_MidYear(POC_AU, all_scenarios[[i]],
                            Population = x,
                            Disease_prog = NULL, 
                            Cascade = NULL, param = NULL, 
                            endYear = endY)%>%ungroup() 
  })
  
  names(subpop_N[[i]]) <- POC_AU$popNames
  
  pop_N[[i]] <- dplyr::bind_rows(subpop_N[[i]], .id = 'population')
  
  total_N[[i]] <- pop_N[[i]]%>%group_by(year)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))
  
  #### N: community ####
  commu_N[[i]] <- pop_N[[i]]%>%filter(population %in% c("C_PWID", "C_fPWID"))%>%
    dplyr::group_by(year)%>%summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))
  
  #### N: prison ####
  prison_N[[i]] <- pop_N[[i]]%>%filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID"))%>%
    dplyr::group_by(year)%>%summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))
  
  #### N: prisonPWID ####
  prisonPWID_N[[i]] <- pop_N[[i]]%>%filter(population %in% c("P_PWID", "P_fPWID"))%>%
    dplyr::group_by(year)%>%summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))
  
  pop_N[[i]] <- pop_N[[i]]%>%arrange(year, population)
  prison_N[[i]] <- prison_N[[i]]%>%arrange(year)
  prisonPWID_N[[i]] <- prisonPWID_N[[i]]%>%arrange(year)
  commu_N[[i]] <- commu_N[[i]]%>%arrange(year)
  
  
  HCVInfect_subpop[[i]] <- indicatorResults(POC_AU, all_scenarios[[i]], "newInfections", 
                                       pop=POC_AU$popNames,
                                       paramR = NULL, range = NULL,
                                       endY = endY)
  
  HCVInfect_subpop_P[[i]] <- HCVInfect_subpop[[i]]%>%
    filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID"))
  
  
  HCVInfect_subpop_C[[i]] <- HCVInfect_subpop[[i]]%>%
    filter(population %in% c("C_PWID", "C_fPWID"))
  
  
  
  pop_N_P[[i]] <- pop_N[[i]]%>%
    filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID"))
  
  
  pop_N_C[[i]] <- pop_N[[i]]%>%
    filter(population %in% c("C_PWID", "C_fPWID"))
  
  
  HCVInc_subpop_C[[i]] <- cbind(year = HCVInfect_subpop_C[[i]]$year,
                           population = HCVInfect_subpop_C[[i]]$population,
                           as.data.frame(100*HCVInfect_subpop_C[[i]][, name_parset]/
                                           pop_N_C[[i]][ ,name_parset]))
  
  HCVInc_subpop_P[[i]] <- cbind(year = HCVInfect_subpop_P[[i]]$year,
                           population = HCVInfect_subpop_P[[i]]$population,
                           as.data.frame(100*HCVInfect_subpop_P[[i]][, name_parset]/
                                           (2*pop_N_P[[i]][ ,name_parset])))
  
  HCVInc_subpop[[i]] <- rbind(HCVInc_subpop_C[[i]], HCVInc_subpop_P[[i]])%>%arrange(year, population)%>%
    tibble::as_tibble()
  
  
  HCVInfect_setting[[i]] <- list()
  
  HCVInfect_setting[[i]][["commu"]] <- HCVInfect_subpop_C[[i]]%>%group_by(year)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)
  
  HCVInfect_setting[[i]][["prisons"]] <- HCVInfect_subpop_P[[i]]%>%group_by(year)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)
  
  HCVInfect_setting[[i]][["prisonsPWID"]] <- HCVInfect_subpop[[i]]%>%
    filter(population %in% c("P_PWID", "P_fPWID"))%>%group_by(year)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)
  
  
  HCVInc_setting[[i]] <- list()
  HCVInc_setting[[i]][["commu"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                     as.data.frame(1000*(HCVInfect_setting[[i]][["commu"]][ , name_parset]/ 
                                                           commu_N[[i]][ ,name_parset])))%>%
    tibble::as_tibble()
  HCVInc_setting[[i]][["prisons"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                       as.data.frame(1000*(HCVInfect_setting[[i]][["prisons"]][ , name_parset]/ 
                                                             (2*prison_N[[i]][ ,name_parset]))))%>%
    tibble::as_tibble()
  
  HCVInc_setting[[i]][["prisonsPWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                           as.data.frame(1000*(HCVInfect_setting[[i]][["prisonsPWID"]][ , name_parset]/ 
                                                                 (2*prisonPWID_N[[i]][ ,name_parset]))))%>%
    tibble::as_tibble()
  
  
  
  
  
  
  }

lab_name <- c("Status Quo", "Prison_testing_I",
              "Prison_testing_II", 
              "Prison_testing_III",
              "Program sustained",
              "Program scale-up")
col_pal <- c(paletteer_d("nationalparkcolors::Acadia"))
View(HCVInc_setting$Prison_testing_II$commu)
HCVInc_setting <- lapply(HCVInc_setting, function(x) bind_rows(x,
                                                               .id = 'setting'))

HCVInc_setting <- bind_rows(HCVInc_setting, .id = "Scenarios")

HCVInc_setting <- HCVInc_setting%>%
  mutate(setting = factor(setting, levels = c("commu", "prisons", "prisonsPWID") , 
                          labels = c("Community", "Prisons_total", "prisons_excluded_nonPWID")),
         
         Scenarios = factor(Scenarios, levels = unique(HCVInc_setting$Scenarios), 
                            labels = lab_name )
         )
inc_setting <- ggplot(HCVInc_setting, aes(x = year, y = best, colour = Scenarios)) +
  geom_line() +
   geom_hline(yintercept = unlist(HCVInc_setting[1, "best"]*0.1), linetype = "dashed") + 
  scale_x_continuous(limits = c(1, 51), breaks = seq(1, 51, 5), 
                     labels = seq(1, 51, 5) + POC_AU$cabY - 1) + 
  
  scale_color_manual(name = "Scenarios", values = col_pal )
  
inc_setting <- inc_setting  +  facet_custom (~setting,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 15))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 100))),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 150)))
                  )) + 
  theme_Publication_facet() + 
  labs(y = "HCV incidence 1000-PY")
  
  
HCVInc_subpop <- HCVInc_subpop%>%bind_rows(., .id = "Scenarios")
HCVInc_subpop <- HCVInc_subpop%>%
  mutate(population = factor(population, levels = POC_AU$popNames,
                             labels = c("PWID in community", "former PWID in community",
                                        "PWID in prisons", "former PWID in prisons", 
                                        "nonPWID in prisons")),
         Scenarios = factor(Scenarios, levels = unique(HCVInc_subpop$Scenarios),
                            labels = lab_name))

hline_df <- data.frame(
  population = c("PWID in community", "PWID in prisons"),  # only these will get hline
  yintercept = 0.2)
)

Inc_subpop <- ggplot(HCVInc_subpop, aes(x = year, y = best, colour = Scenarios)) +
  geom_line() +
  geom_hline(data = hline_df, 
             aes(yintercept = yintercept), 
             linetype = "dashed") + 
  
  scale_x_continuous(limits = c(1, 51), breaks = seq(1, 51, 5), 
                     labels = seq(1, 51, 5) + POC_AU$cabY - 1) + 
  
  scale_color_manual(name = "Scenarios", values = col_pal )

Inc_subpop <-  Inc_subpop + facet_custom (~population,
                               scales = "free", ncol = 3,
                               scale_overrides = 
                                 list(
                                   scale_new(1,
                                             scale_y_continuous(limits = 
                                                                  c(0, 1))),
                                   scale_new(2,
                                             scale_y_continuous(limits = 
                                                                  c(0, 4))),
                                   scale_new(3,
                                             scale_y_continuous(limits = 
                                                                  c(0, 1))),
                                   scale_new(4,
                                             scale_y_continuous(limits = 
                                                                  c(0, 5))),
                                   scale_new(5,
                                             scale_y_continuous(limits = 
                                                                  c(0, 30)))
                                 )) + 
  theme_Publication_facet() + 
  labs(y = "HCV incidence 1000-PY")


ggsave(file = file.path(OutputFig, paste0("Inc_subpop", ".png")), 
       Inc_subpop, 
       width = 12, height = 8, bg = "white", dpi = 300) 

ggsave(file = file.path(OutputFig, paste0("Inc_setting", ".png")), 
       inc_setting , 
       width = 12, height = 8, bg = "white", dpi = 300) 
