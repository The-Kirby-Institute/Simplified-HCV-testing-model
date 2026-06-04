# this script is organizing the model output at each time step and export to spread sheet 
# for checking bugs and details 

# new infections
# from 2022-2042
inf_cum <- epi.flow.sub%>%ungroup()%>%group_by(scenario, year)%>%
  summarise(newinf = sum(newInfections))%>%ungroup()%>%
  mutate(cum_inf = cumsum(newinf))%>%select(!cum_inf)%>%
  spread(scenario, newinf)%>%
  mutate(`Averted cases` = `status quo` - `National program`)%>%
  gather(scenario, newinf, -year)%>%
  mutate(scenario = factor(scenario, levels = unique(scenario)))

p_inf_cum  <- ggplot(inf_cum, aes(x=year, y=newinf)) +
  geom_segment( aes(x=year, xend=year, y=0, yend=newinf), color="grey") +
  geom_point( color="orange", size=2) + 
  facet_wrap(.~scenario, scales = "free") + 
  scale_x_continuous(limits = c(2022, 2040))+ 
  theme_linedraw()  + 
  ylab("Cases per year ") + 
  ggtitle("New infections")

ggsave(file=file.path(OutputFig, "newinf.png"), p_inf_cum, 
       width = 8 , height = 6, bg = "white")

# HCV related death 

HCVdeath_cum <- epi.flow.sub%>%ungroup()%>%group_by(scenario, year)%>%
  summarise(newHCVdeaths = sum(newHCVdeaths))%>%ungroup()%>%
  mutate(cum_hcvd = cumsum(newHCVdeaths))%>%select(!cum_hcvd)%>%
  spread(scenario, newHCVdeaths)%>%
  mutate(`Averted cases` = `status quo` - `National program`)%>%
  gather(scenario, newHCVdeaths, -year)%>%
  mutate(scenario = factor(scenario, levels = unique(scenario)))

p_HCVdeath  <- ggplot(HCVdeath_cum, aes(x=year, y=newHCVdeaths)) +
  geom_segment( aes(x=year, xend=year, y=0, yend=newHCVdeaths), color="grey") +
  geom_point( color="orange", size=2) + 
  facet_wrap(.~scenario, scales = "free") + 
  scale_x_continuous(limits = c(2022, 2040))+ 
  theme_linedraw()  + 
  ylab("Cases per year ") + 
  ggtitle("HCV releated death") 

p_HCVdeath <- p_HCVdeath + facet_custom (~scenario,
                                         scales = "free", ncol = 3,
                                         scale_overrides = 
                                           list(
                                             scale_new(1,
                                                       scale_y_continuous(limits = 
                                                                            c(0, 600),
                                                                          breaks = seq(0, 600, 100),
                                                                          labels = seq(0, 600, 100))),
                                             scale_new(2,
                                                       scale_y_continuous(limits = 
                                                                            c(0, 600),
                                                                          breaks = seq(0, 600, 100),
                                                                          labels = seq(0, 600, 100))),
                                             scale_new(3,
                                                       scale_y_continuous(limits = 
                                                                            c(0, 5),
                                                                          breaks = seq(0, 5, 1),
                                                                          labels = seq(0, 5, 1)))))

ggsave(file=file.path(OutputFig, "hcvd_case.png"), p_HCVdeath , width = 8 , height = 6, bg = "white")




ggsave(file=file.path(OutputFig, "newinf.png"), p_inf_cum, 
       width = 8 , height = 6, bg = "white")

# advanced liver disease 
ad_case <- Num_disprog_dt%>%group_by(scenario, year)%>%
  summarise(s = sum(s),
            a = sum(a), 
            f0 = sum(f0), 
            f1 = sum(f1), 
            f2 = sum(f2), 
            f3 = sum(f3),
            f4 = sum(f4),
            dc = sum(dc), 
            hcc =sum(hcc), 
            lt = sum(lt), 
            plt = sum(plt))%>%
  mutate(best= dc + hcc + lt + plt)%>%
  select(scenario, year, best)%>%
  group_by(scenario)%>%
  spread(scenario, best)%>%
  mutate(`Averted cases` = `status quo` - `National program`)%>%
  gather(scenario, best, -year)
ad_case <- ad_case%>%mutate(scenario = factor(scenario, levels = unique(scenario)))

p_ad_case <- ggplot(ad_case, aes(x=year, y=best)) +
  geom_segment( aes(x=year, xend=year, y=0, yend=best), color="grey") +
  geom_point( color="orange", size=2) + 
  facet_wrap(.~scenario, scales = "free") + 
  scale_x_continuous(limits = c(2022, 2040))+ 
  theme_linedraw()  + 
  ylab("Cases per year ") + 
  ggtitle("Advanced liver disease cases")

p_ad_case <- p_ad_case +
  facet_custom (~scenario,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 80000),
                                                 breaks = seq(0, 80000, 10000),
                                                 labels = seq(0, 80000, 10000))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 80000),
                                                 breaks = seq(0, 80000, 10000),
                                                 labels = seq(0, 80000, 10000))),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 100),
                                                 breaks = seq(0, 100, 10),
                                                 labels = seq(0, 100, 10)))))

ggsave(file=file.path(OutputFig, "ad_case.png"), p_ad_case, width = 8 , height = 6, bg = "white")


c_all <- cost_all%>%group_by(scenario, year)%>%
  summarise(best = sum(discountValue))%>%
  ungroup()%>%
  group_by(scenario)%>%
  mutate(best = cumsum(best))%>%
  spread(scenario, best)%>%
  mutate(`addtional cost` = `National program` - `status quo`)%>%
  gather(scenario, best, -year)
c_all <- c_all%>%mutate(scenario = factor(scenario, levels = unique(scenario)))

total_cost <- ggplot(c_all, aes(x=year, y=best)) +
  geom_segment( aes(x=year, xend=year, y=0, yend=best), color="grey") +
  geom_point( color="orange", size=2) + 
  facet_wrap(.~scenario, scales = "free") + 
  scale_x_continuous(limits = c(2022, 2040))+ 
  theme_linedraw()  + 
  ylab("Cumulative Cost (Millions)")


total_cost <- total_cost + 
  ggtitle("Total cost") +
  facet_custom (~scenario,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 4500000000),
                                                 breaks = seq(0, 4500000000, 500000000),
                                                 labels = seq(0, 4500000000, 500000000)/1000000)),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 4500000000),
                                                 breaks = seq(0, 4500000000, 500000000),
                                                 labels = seq(0, 4500000000, 500000000)/1000000)),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(-10000000, 50000000),
                                                 breaks = seq(-10000000, 50000000, 5000000),
                                                 labels = seq(-10000000, 50000000, 5000000)/1000000))))

ggsave(file=file.path(OutputFig, "total_cost.png"), total_cost, width = 8 , height = 6, bg = "white")


# cost stage: diagnosis, treatment, management 


c_stage <- cost_all%>%
  mutate(diag_sum = costTestingAb_sum + costTestingAg_sum + costTestingPOCT_sum, 
         treatment_sum = costTreatment_sum + costRetreat_sum + costCured_sum, 
         management_sum = costCompartment_sum)%>%
  select(scenario, year, population, diag_sum , treatment_sum, 
         management_sum, id, discount)%>%
  group_by(scenario, year)%>%
  summarise(diag_sum = sum(diag_sum),
            treatment_sum = sum(treatment_sum), 
            management_sum = sum(management_sum))%>%
  mutate(id = year - POC_AU$simY)%>%
  mutate(discount = (1 + 0.05)^id)%>%filter(id>= 0)%>%
  mutate(dis_diag= diag_sum/discount, 
         dis_treat= treatment_sum/discount,
         dis_manage= management_sum/discount)%>%
  ungroup()%>%
  group_by(scenario)%>%
  mutate(dis_diag = cumsum(dis_diag), 
         dis_treat = cumsum(dis_treat),
         dis_manage = cumsum(dis_manage))%>%
  ungroup()

c_sc <- c_stage%>%filter(scenario =="status quo")
c_np <- c_stage%>%filter(scenario =="National program")  
c_add <-  cbind(scenario =  "addtional cost", 
                year = c_sc$year, 
                id = c_sc$id,
                discount = c_sc$discount, 
                as.data.frame(c_np[, c(3,4,5, 8,9,10)] - 
                                c_sc[, c(3,4,5, 8,9,10)]))
                
c_stage <- rbind(c_stage, c_add)               
  
c_stage_long <- gather(c_stage, stage, best, 
                       -c(scenario, year, id, discount, diag_sum, treatment_sum, management_sum)) 
 
c_stage_long <- c_stage_long%>%mutate(scenario = factor(scenario, levels = unique(scenario)))

p_stage_cost <- list()
for(i in unique(c_stage_long$stage)){ 
  
  p_stage_cost[[i]] <- 
    ggplot(c_stage_long%>%filter(stage == i), aes(x=year, y=best)) +
    geom_segment( aes(x=year, xend=year, y=0, yend=best), color="grey") +
    geom_point( color="orange", size=2) + 
    facet_wrap(.~scenario, scales = "free") + 
    scale_x_continuous(limits = c(2022, 2040))+ 
    theme_linedraw()  + 
    ylab("Cumulative Cost (Millions)") 
  }
p_stage_cost[[1]] <- p_stage_cost[[1]] + 
  ggtitle("Cost: diagnosis") + 
  facet_custom (~scenario,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 400000000),
                                                 breaks = seq(0, 400000000, 50000000),
                                                 labels = seq(0, 400000000, 50000000)/1000000)),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 400000000),
                                                 breaks = seq(0, 400000000, 50000000),
                                                 labels = seq(0, 400000000, 50000000)/1000000)),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 5000000),
                                                 breaks = seq(0, 5000000, 1000000),
                                                 labels = seq(0, 5000000, 1000000)/1000000))))



p_stage_cost[[2]] <- p_stage_cost[[2]] + 
  ggtitle("Cost: treatment") + 
  facet_custom (~scenario,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 3000000000),
                                                 breaks = seq(0, 3000000000, 500000000),
                                                 labels = seq(0, 3000000000, 500000000)/1000000)),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 3000000000),
                                                 breaks = seq(0, 3000000000, 500000000),
                                                 labels = seq(0, 3000000000, 500000000)/1000000)),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(-10000000, 35000000),
                                                 breaks = seq(-10000000, 35000000, 5000000),
                                                 labels = seq(-10000000, 35000000, 5000000)/1000000))))

p_stage_cost[[3]] <- p_stage_cost[[3]] + 
  ggtitle("Cost: management") + 
  facet_custom (~scenario,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 2000000000),
                                                 breaks = seq(0, 2000000000, 500000000),
                                                 labels = seq(0, 2000000000, 500000000)/1000000)),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 2000000000),
                                                 breaks = seq(0, 2000000000, 500000000),
                                                 labels = seq(0, 2000000000, 500000000)/1000000)),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(-10000000, 0),
                                                 breaks = seq(-10000000, 0, 5000000),
                                                 labels = seq(-10000000, 0, 5000000)/1000000))))

names(p_stage_cost) <- c("diagnosis", "treatment", "management")
for(i in names(p_stage_cost)){ 
  ggsave(file=file.path(OutputFig, paste0("stage_cost_",i ,".png")), 
         p_stage_cost[[i]], width = 8 , height = 6, bg = "white")  
  }

  


epi.flow.sub
