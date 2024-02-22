# each timesteps 
# compartments 
modres.t <- function(pg, Best, endYear, allp = NULL) {
  df_list <- lapply(Best, as.data.frame.table)
  if(!is.null(allp)){ 
    allpop <- df_list[[allp]]
  }
  else{ 
    allpop <- df_list$allPops
    
  }
  allpop <- allpop%>%
    mutate(time = rep(seq(1.0,(endYear - pg$timestep), pg$timestep), 
                      each=pg$ncomponent*pg$npops),
           cascade = sub("^[^_]*_", "", Var2), 
           disease_prog = sub("\\_.*", "", Var2))%>%
    dplyr::select(-Var3)%>%ungroup()
  
  allpop <- allpop%>%
    mutate(time = c(rep(seq(pg$startYear, endYear- 1*pg$timestep,
                            pg$timestep),each = pg$npops*pg$ncomponent)))
  
  names(allpop) <- c("population", "state", "best","timestep", "cascade",
                     "disease_prog")
  
  ## MidyearIndex
  timelong <- seq(pg$startYear, endYear, pg$timestep) 
  
  allpop <- allpop%>%filter(timestep%in% timelong)%>%
    mutate(year = timestep%/%1)
}

# for flows 
modres.flow.t <- function(pg, Best, endYear, allp = NULL) {
  flowpop <- as.data.frame.table(Best[[allp]])
  
  flowpop  <- flowpop %>%
    mutate(dt = rep(seq(1.0,(endYear - pg$timestep), pg$timestep), 
                      each=pg$npops))%>%ungroup()%>%
    mutate(year = dt%/%1 - 1, 
           population = Var1, 
           best = Freq)%>%
    mutate(timestep = dt)%>%
    select(year, timestep, population, best)
  return(flowpop)

  }
  

# extract the population-specific resutls 
N_pop_sum <- function(dt, pop = NULL){ 
  
  if(is.null(pop)){ 
    a <- dt%>%group_by(timestep, population)%>%
      summarise_at("best", sum)%>%arrange(timestep)
  }
  else{ 
    a <- dt%>%filter(population %in% pop)%>%
      group_by(timestep)%>%
      summarise_at("best", sum)%>%arrange(timestep)
    
    }
  
  
  return(a)
  }

# extract the population-specific resutls 
N_pop_casdisprog <- function(dt, pop = NULL, cas = NULL, disprog = NULL){ 
  if(!is.null(disprog)){ 
    x <- dt%>%filter(disease_prog %in% disprog)
    
  }
  else{ 
    
    x <- dt
    
    }
  
  
  if(is.null(pop) & is.null(cas)){ 
    a <- x%>%group_by(timestep, population)%>%
      summarise_at("best", sum)%>%arrange(timestep, population)
  }
  else if(!is.null(pop) & is.null(cas)){ 
    a <- x%>%filter(population %in% pop)%>%
      group_by(timestep, cascade)%>%
      summarise_at("best", sum)%>%arrange(timestep)
    
  }
  else if(is.null(pop) & !is.null(cas)){
    a <- x%>%filter(cascade %in% cas)%>%
      group_by(timestep, population)%>%
      summarise_at("best", sum)%>%arrange(timestep, population)
    
  }else if(!is.null(pop) & !is.null(cas)){
    a <- x%>%filter(cascade %in% cas & population %in% pop)%>%
      group_by(timestep)%>%
      summarise_at("best", sum)%>%arrange(timestep, population)
    
  }
  
  
  return(a)
}
