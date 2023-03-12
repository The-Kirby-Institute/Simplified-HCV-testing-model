#### check steady state #### 
## use plot and find the steady state by visual 
## might need to find the steady state time point automatically
## My understand of steady state: y(t)-y(t-1)=0????? 

check_steady <-function(model_result){
  source(file.path(Rcode, "/Functions/plotOptions.R"))
  
  df_list <- lapply(model_result, as.data.frame.table)
  allpop <- df_list$allPops%>%mutate(time = rep(seq(1,1999.9,0.1), each=244),
                                     Frequency=round(Freq, digits = 3))
  allpop_cascade <- allpop%>%mutate(cascade_status = sub("^[^_]*_", "", Var2), 
                                    dis_prog = sub("\\_.*", "", Var2))%>%
    group_by(Var1, dis_prog, cascade_status, time)%>%
    summarise(statenum = sum(Frequency))
  
  state_change <- allpop_cascade%>%
    mutate(state = ifelse(cascade_status=="cured", "cured",dis_prog))%>%
    group_by(Var1, state, time)%>%summarise(state_num = sum(statenum))%>%
    filter(time<=1500)%>%
    ggplot(data = ., aes(x = time, y = state_num )) + 
    geom_line(aes(colour = state)) + 
    xlab("Year") + ylab("number") +labs(tag = "") + 
    scale_x_continuous(expand = c(0,0), limits = c(0, 1500), 
                       breaks = seq(0, 1500, 100))  + 
    facet_wrap(.~Var1, ncol=2, scales = "free") + 
    plotOpts 
  
  return(state_change)
}


check_steady(steady)
ta <- round(steady$allPops, digits = 3)%>%as.data.frame()%>%t()
taint <- round(steady$newInfections, digits = 5)%>%as.data.frame()%>%t()
