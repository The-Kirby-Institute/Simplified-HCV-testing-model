#### check steady state #### 
## use plot and find the steady state by visual 
## might need to find the steady state time point automatically
## My understand of steady state: y(t)-y(t-1)=0?????  

check_steady <-function(model_result, endY,timestep, Ncomp, Tequilibrium){
  # endY: endY
  # timestep
  #Ncomp = number of components*number of pops 
  source(file.path(Rcode, "/Functions/plotOptions.R"))
  
  df_list <- lapply(model_result, as.data.frame.table)
  allpop <- df_list$allPops%>%
    mutate(time = rep(seq(1,(endY - timestep),timestep), each=Ncomp),
                                     Frequency=round(Freq, digits = 3))
  allpop_cascade <- allpop%>%mutate(cascade_status = sub("^[^_]*_", "", Var2), 
                                    dis_prog = sub("\\_.*", "", Var2))%>%
    group_by(Var1, dis_prog, cascade_status, time)%>%
    summarise(statenum = sum(Frequency))
  
  state_change <- allpop_cascade%>%
    mutate(state = ifelse(cascade_status=="cured", "cured",dis_prog))%>%
    group_by(Var1, state, time)%>%summarise(state_num = sum(statenum))%>%
    filter(time<=Tequilibrium)%>%
    ggplot(data = ., aes(x = time, y = state_num )) + 
    geom_line(aes(colour = state)) + 
    xlab("Year") + ylab("number") +labs(tag = "") + 
    scale_x_continuous(expand = c(0,0), limits = c(0, Tequilibrium), 
                       breaks = seq(0, Tequilibrium, 250))  + 
    facet_wrap(.~Var1, ncol=2, scales = "free") + 
    plotOpts 
  
  return(state_change)
}



