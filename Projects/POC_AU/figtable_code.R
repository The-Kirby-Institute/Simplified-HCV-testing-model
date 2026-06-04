#### this script contains functions for figures and tables 

plot_pocau <- function(pj,dt,indicator, type, year_obs, UI = NULL){ 
  # type: new, cum, and avert  
  if(length(unique(dt[[indicator]]$scenario)) == 2){ 
    col_pal <- c("#000000", "#E69F00")
    
  } 
  else{col_pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#264653")
  }
  
  if(isTRUE(type == "new")){ 
    
    a <- dt[[indicator]]%>%ggplot(data = ., 
                                  aes(x = year, colour = scenario)) + 
      geom_line(aes(y = best, colour = scenario, linetype = scenario), size = 1) + 
      scale_x_continuous(expand = c(0,0), limits = c(pj$simY - 1,2050), breaks = c(seq(pj$simY - 1,2050, 5), 2050))+
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt[[indicator]]$scenario)) - 1))) +
      theme_Publication()
    
  }
  else if(isTRUE(type =="cum")){ 
    a <- dt[[indicator]]%>%filter(year %in% year_obs)%>%
      ggplot(data =. , aes(x = as.character(year),  colour = scenario)) + 
      geom_errorbar(aes(ymin = q5, ymax = q95, colour = scenario), width = 0.5,
                    size = 1, 
                    position = "dodge") + 
      geom_point(aes(y = Mu, colour = scenario),position = position_dodge(width = 0.5),size = 1.5) +
      labs(x = "", y = indicator) + 
      scale_x_discrete(labels = c(paste0((year_obs - pj$simY + 1), "-Year", sep = ""))) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_fill_manual(name = "Scenarios", values = col_pal ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt[[indicator]]$scenario)) - 1))) +
      theme_Publication()
    
    
  }
  
  else if(isTRUE(type =="avert")){ 
    a <- dt[[indicator]]%>%filter( year %in% year_obs)%>%
      ggplot(data =. , aes(x = as.character(year), colour = scenario)) + 
      geom_errorbar(aes(ymin = q5, ymax = q95, colour = scenario),
                    size = 1, position = "dodge" , width = 0.5) + 
      geom_point(aes(y = Mu, colour = scenario), position = position_dodge(width = 0.5) , size = 1.5) + 
      labs(x = "", y = indicator) + 
      scale_x_discrete(labels = c(paste0((year_obs - pj$simY + 1), "-Year", sep = ""))) + 
      scale_color_manual(name = "Scenarios", values = col_pal[-1] ) + 
      scale_fill_manual(name = "Scenarios", values = col_pal[-1] ) + 
      
      theme_Publication()
    
  } 
  a <- a + theme(legend.key.size = unit(1,"line"))
  
  return(a)
}
lim_ident <- function(dt, year_range){ 
  
  lim <- dt%>%as.data.frame()%>%ungroup()%>%
    filter(year %in% year_range)%>%
    summarise(x = max(q95, na.rm = TRUE))%>%
    mutate(lim = case_when( 
      x <1 ~5, 
      x>=1 & x <10 ~ 10, 
      x>=10 & x< 30 ~ 30, 
      x>=30 & x<60 ~ 60, 
      x>=60 & x<80 ~ 80, 
      x>80 & x<=100 ~100,
      x>100 & x<=1000 ~ (x%/%100 + 1)*100,
      x >1000 & x<=10000 ~ (x%/%1000 + 1)*1000,
      x >10000 & x <= 100000  ~ (x%/%10000 + 1)*10000,
      x >100000 & x <= 1000000  ~ (x%/%100000 + 1)*100000,
      x >1000000 & x <= 10000000 ~ (x%/%1000000 + 1)*1000000
    ))
  
  return(lim)  
}

# plot function for generating cascade numbers 
Cas_num_plot <- function(pj, dt, obdt =NULL, xlimits, UI = NULL, population = NULL){ 
  # @ UI the name of scenario for showing UI range
  if(length(unique(dt$scenario)) == 2){ 
    col_pal <- c("#000000", "#E69F00")
    
  } 
  else{col_pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
  }
  
  if(is.null(obdt) & is.null(UI) & !is.null(population)){
    traj_plot <- ggplot(dt, aes(x = year, y = best)) + 
      geom_line(aes(colour = scenario, linetype = scenario), size = 1) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt$scenario)) - 1))) + 
      coord_cartesian(xlim = xlimits[1:2]) +
      scale_x_continuous(expand = c(0.01, 0), limits =c(xlimits[1], xlimits[2]) ,
                         breaks = seq(xlimits[1], xlimits[2], 
                                      by = xlimits[3]),
                         labels = seq(pj$cabY + xlimits[1] - 1, 
                                      (pj$cabY + xlimits[2] - 1),
                                      xlimits[3])) + 
      theme(panel.spacing = unit(2, "lines")) + theme_Publication() + 
      theme(legend.key.size = unit(1,"line"))
  }
  
  else if(is.null(obdt) & is.null(UI)){ 
    traj_plot <- ggplot(dt, aes(x = year, y = best)) + 
      geom_line(aes(colour = scenario, linetype = scenario)) + 
      
      facet_wrap(~ population, scale ="free", ncol = 2 ) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt$scenario)) - 1))) + 
      coord_cartesian(xlim = xlimits[1:2]) +
      scale_x_continuous(expand = c(0.01, 0), limits =c(xlimits[1], xlimits[2]) ,
                         breaks = seq(xlimits[1], xlimits[2], 
                                      by = xlimits[3]),
                         labels = seq(pj$cabY + xlimits[1] - 1, 
                                      (pj$cabY + xlimits[2] - 1),
                                      xlimits[3])) +
      theme(panel.spacing = unit(2, "lines")) + theme_Publication_facet() + 
      theme(legend.key.size = unit(1,"line"))
  }
  else if(is.null(obdt) & !is.null(UI)){ 
    traj_plot <- ggplot(dt, aes(x = year, y = best)) + 
      geom_line(aes(colour = scenario, linetype = scenario)) + 
      geom_ribbon(data =dt%>%filter(scenario == UI) ,aes(ymin = q5, ymax = q95), fill = "#000000", alpha = 0.2) +
      facet_wrap(~ population, scale ="free", ncol = 2 ) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_fill_manual(name = "Scenarios", values = "#000000" ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt$scenario)) - 1))) + 
      coord_cartesian(xlim = xlimits[1:2]) +
      scale_x_continuous(expand = c(0.01, 0), limits =c(xlimits[1], xlimits[2]) ,
                         breaks = seq(xlimits[1], xlimits[2], 
                                      by = xlimits[3]),
                         labels = seq(pj$cabY + xlimits[1] - 1, 
                                      (pj$cabY + xlimits[2] - 1),
                                      xlimits[3])) + 
      theme(panel.spacing = unit(2, "lines")) + theme_Publication_facet() + 
      theme(legend.key.size = unit(1,"line"))
    
  }
  
  else if(!is.null(obdt) & is.null(UI)){
    traj_plot <- ggplot(dt, aes(x = year, y = best)) + 
      geom_line(aes(colour = scenario, linetype = scenario)) + 
      facet_wrap(~ population, scale ="free", ncol = 2 ) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt$scenario)) - 1))) + 
      coord_cartesian(xlim = xlimits[1:2]) +
      scale_x_continuous(expand = c(0.01, 0), limits =c(xlimits[1], xlimits[2]) ,
                         breaks = seq(xlimits[1], xlimits[2], 
                                      by = xlimits[3]),
                         labels = seq(pj$cabY + xlimits[1] - 1, 
                                      (pj$cabY + xlimits[2] - 1),
                                      xlimits[3])) + 
      geom_point(data=obdt, aes(y=realPop, x = year), 
                 colour = "black", size = 1) +
      geom_segment(data = obdt, 
                   aes ( y = low, yend = up, x = year, xend = year)) +
      theme(panel.spacing = unit(2, "lines")) + theme_Publication_facet() + 
      theme(legend.key.size = unit(1,"line"))
    
    
  }
  else if(!is.null(obdt) & !is.null(UI) & is.null(population)){
    traj_plot <- ggplot(dt, aes(x = year, y = best)) + 
      geom_line(aes(colour = scenario, linetype = scenario)) + 
      geom_ribbon(data =dt%>%filter(scenario == UI), aes(ymin = q5, ymax = q95), fill = "#000000", alpha = 0.2) +
      facet_wrap(~ population, scale ="free", ncol = 2 ) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_fill_manual(name = "Scenarios", values = "#000000" ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt$scenario)) - 1))) + 
      coord_cartesian(xlim = xlimits[1:2]) +
      scale_x_continuous(expand = c(0.01, 0), limits =c(xlimits[1], xlimits[2]) ,
                         breaks = seq(xlimits[1], xlimits[2], 
                                      by = xlimits[3]),
                         labels = seq(pj$cabY + xlimits[1] - 1, 
                                      (pj$cabY + xlimits[2] - 1),
                                      xlimits[3])) + 
      geom_point(data=obdt, aes(y=realPop, x = year), 
                 colour = "black", size = 1) +
      geom_segment(data = obdt, 
                   aes ( y = low, yend = up, x = year, xend = year)) +
      theme(panel.spacing = unit(2, "lines")) + theme_Publication_facet() + 
      theme(legend.key.size = unit(1,"line"))
  }else if(!is.null(obdt) & !is.null(UI) & !is.null(population)){
    traj_plot <- ggplot(dt, aes(x = year, y = best)) + 
      geom_line(aes(colour = scenario, linetype = scenario)) + 
      geom_ribbon(data =dt%>%filter(scenario == UI), aes(ymin = q5, ymax = q95), fill = "#000000", alpha = 0.2) +
      
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_fill_manual(name = "Scenarios", values = "#000000" ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt$scenario)) - 1))) + 
      coord_cartesian(xlim = xlimits[1:2]) + 
      geom_point(data=obdt, aes(y=realPop, x = year), 
                 colour = "black", size = 1)  +
      theme(panel.spacing = unit(2, "lines")) + theme_Publication() + 
      theme(legend.key.size = unit(1,"line"))
  }
  
  return(traj_plot)
}
