 # figure generation for the model test  

total_population <- bestPopSizes%>%group_by(time_step)%>%
  summarize(total = sum(best))%>%as.data.frame()

ggplot(data = total_population, aes(x = time_step, y = total)) +
  geom_line() + 
  xlab("Year") + ylab("Total population") +
  plotOpts 

ggsave("total_population.png", width = 6, height = 4, 
       path = paste0(basePath, "/04. Output"))


## population by MSM group overtime (midYear value)
subpopulation <- bestPopSizes%>%group_by(population, time_step)%>%
  summarize(total = sum(best))%>%as.data.frame

ggplot(data = subpopulation, aes(x = time_step, y = total)) +
  geom_line() + 
  xlab("Year") + ylab("Population size") +
  scale_x_continuous(expand = c(0,0), limits = c(0,HCV$npts), 
                     breaks = seq(0,HCV$npts,100), labels = seq(1,HCV$npts,100)) +
  facet_grid(.~population) +
  plotOpts 
ggsave("population_size.png", width = 8, height = 4, 
       path = paste0(basePath, "/04. Output"))

## prevalence overtime (total infected number )
prev <- bestPopSizes%>%filter(!state %in% c("s", "a_cured", "f0_cured", "f1_cured",
                                    "f2_cured", "f3_cured", "f4_cured", 
                                    "dc_cured", "hcc_cured", "lt_cured", 
                                    "plt_cured"))%>%group_by(time_step)%>%
  summarize(total = sum(best))%>%as.data.frame()
prev1 <- cbind(prev, total_population$total)%>%
  mutate(prev = (total/total_population$total)*100)

ggplot(data = prev1, aes(x = time_step, y = prev)) +
  geom_line() + 
  xlab("Year") + ylab("Overall prevalence (%)") +
  scale_x_continuous(expand = c(0,0), limits = c(0,HCV$npts), 
                     breaks = seq(0,HCV$npts,100), labels = seq(0,HCV$npts,100)) +
  
  plotOpts 
ggsave("overallPrev.png", width = 6, height = 4, 
       path = paste0(basePath, "/04. Output"))


## prevalence by MSM groups

prev_group <- bestPopSizes%>%filter(!state %in% c("s", "a_cured", "f0_cured", "f1_cured",
                                         "f2_cured", "f3_cured", "f4_cured", 
                                         "dc_cured", "hcc_cured", "lt_cured", 
                                         "plt_cured"))%>%
  group_by(population, time_step)%>%
  summarize(total = sum(best))%>%as.data.frame()
prev_group1 <- cbind(prev_group, subpopulation$total)%>%
  mutate(prev = (total/subpopulation$total)*100)

ggplot(data = prev_group1, aes(x = time_step, y = prev)) +
  geom_line() + 
  xlab("Year") + ylab("Prevalence by MSM group") +
  scale_x_continuous(expand = c(0,0), limits = c(0,396), 
                     breaks = seq(0,HCV$npts,100), labels = seq(0,HCV$npts,100)) +
  facet_grid(.~population) +
  plotOpts 
ggsave("groupPrev.png", width = 8, height = 4, 
       path = paste0(basePath, "/04. Output"))
stage_flow%>%filter(dis_prog=="f0")%>%filter(best!=0)

## each stage 
stage_flow <-bestPopSizes%>%
  mutate(dis_prog =ifelse(.$state=="s","s", gsub("_.*","",.$state)))%>%
  group_by(population, dis_prog, time_step)%>%
  summarize(total = sum(best))%>%as.data.frame()

ggplot(data = stage_flow, aes(x = time_step, y = total)) +
  geom_line(aes(color = dis_prog), position = position_dodge(width = 0.1)) + 
  xlab("Year") + ylab("Population size") +
  scale_x_continuous(expand = c(0,0), limits = c(0,396), 
                     breaks = seq(3,395,40), labels = seq(1,99,10)) +
  facet_grid(.~population) +
  plotOpts 


str(stage_flow$dis_prog)
                                            