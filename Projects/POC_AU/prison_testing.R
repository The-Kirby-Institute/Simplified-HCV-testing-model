# prisons testing scenarios 
fm <- list()
fm[["2022"]] <- c(1.1,1.1,  2.4,  2.4 , 0.01)
fm[["2023"]] <- c(1.1, 1.1, 18.5, 18.5, 1)
fm[["2024"]] <- c(1.1, 1.1, 18.5, 18.5, 1)
# fm[["2025"]] <- c(0.1, 0.1, 15, 15, 1)
fm[["2025"]] <- c(1.1, 1.1, 15, 15, 1)
fm[["2026"]] <- c(1.1, 1.1, 15, 15, 1)
#fm[["2026"]] <- c(0.1, 0.1, 15, 15, 1)
# fm[["2027"]] <- c(0.1, 0.1, 13, 13, 1)
fm[["2027"]] <- c(0.8, 0.8, 13, 13, 1)
fm[["2028"]] <- c(0.8, 0.8, 15, 15, 1)
fm[["2029"]] <- c(0.8, 0.8, 13, 13, 1)
fm[["2030"]] <- c(0.8, 0.8, 15, 15, 1)


#### prison_testing_I: program sustained by reducing testing in community ####
# 2027

odd_num_test <- c(0.8, 1)

Ccal[[2027]] <- list("C" = Ccal[[2026]]$C*odd_num_test[1],
                     "P" = Ccal[[2026]]$P*odd_num_test[2])

frac_test[[2027]] <- frac_test[[2026]]

frac_ab[["2027"]] <- c(unlist(as.numeric(frac_test[[2027]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2027]]$P$reflex)))

dfList_NP_2027 <- dfList_NPPhaseII
for(i in param_var){  
  dfList_NP_2027[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2027, index = i, 
                                   frac_testing = frac_test[[2027]],
                                   S_Yint = 2027, S_Yend = 2028, r_Yend = 2028, NPlst = NPlst, 
                                   fp = c(Ccal[[2027]]$C*fm[["2027"]][1], 
                                          Ccal[[2027]]$C*fm[["2027"]][2], 
                                          Ccal[[2027]]$P*fm[["2027"]][3], 
                                          Ccal[[2027]]$P*fm[["2027"]][4], 
                                          Ccal[[2027]]$P*fm[["2027"]][5]))
  
}

for(i in param_var){
  # begining of 2027
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPPhaseII[[i]])[3]
  dfList_NP_2027[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2027"]] <- c(Ccal[[2027]]$C, Ccal[[2027]]$P)
n_ab_np[["2027"]] <- n_ab_np[["2026"]]*odd_num_test
xfs[["2027"]] <- fs_estimate(num_ab = n_ab_np[["2027"]], 
                             cov_np = coverage_np[["2027"]], 
                             frac_ab = frac_ab[["2027"]], 
                             fp = fm[["2027"]], year = 2027, endY = 100,
                             modsim = Sce_sq)

View(xfs[["2027"]][[1]])
ini_dt <- (2027 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2027 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2027"]] <- fs[["2026"]]
fs[["2027"]][1, ini_dt: end_dt] <-   xfs[["2027"]][[1]][1, ini_dt: end_dt]/fm[["2027"]][1]
fs[["2027"]][2, ini_dt: end_dt] <-   xfs[["2027"]][[1]][2, ini_dt: end_dt]/fm[["2027"]][2]
fs[["2027"]][3, ini_dt: end_dt] <-   xfs[["2027"]][[1]][3, ini_dt: end_dt]/fm[["2027"]][3]
fs[["2027"]][4, ini_dt: end_dt] <-   xfs[["2027"]][[1]][4, ini_dt: end_dt]/fm[["2027"]][4]
fs[["2027"]][5, ini_dt: end_dt] <-   xfs[["2027"]][[1]][5, ini_dt: end_dt]/xfs[["2027"]][[1]][5, ini_dt: end_dt]




# 2028 

odd_num_test <- c(0.9, 1.05)

Ccal[[2028]] <- list("C" = Ccal[[2027]]$C*odd_num_test[1],
                     "P" = Ccal[[2027]]$P*odd_num_test[2])

frac_test[[2028]] <- frac_test[[2026]]

frac_ab[["2028"]] <- c(unlist(as.numeric(frac_test[[2028]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2028]]$P$reflex)))

dfList_NP_2028 <- dfList_NP_2027
for(i in param_var){  
  dfList_NP_2028[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2028, index = i, 
                                   frac_testing = frac_test[[2028]],
                                   S_Yint = 2028, S_Yend = 2029, r_Yend = 2029, NPlst = NPlst, 
                                   fp = c(Ccal[[2028]]$C*fm[["2028"]][1], 
                                          Ccal[[2028]]$C*fm[["2028"]][2], 
                                          Ccal[[2028]]$P*fm[["2028"]][3], 
                                          Ccal[[2028]]$P*fm[["2028"]][4], 
                                          Ccal[[2028]]$P*fm[["2028"]][5]))
  
}

for(i in param_var){
  # begining of 2026
  b_pt <- (2029 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2028[[i]])[3]
  dfList_NP_2028[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2028"]] <- c(Ccal[[2028]]$C, Ccal[[2028]]$P)
n_ab_np[["2028"]] <- n_ab_np[["2027"]]*odd_num_test
xfs[["2028"]] <- fs_estimate(num_ab = n_ab_np[["2028"]], 
                             cov_np = coverage_np[["2028"]], 
                             frac_ab = frac_ab[["2028"]], 
                             fp = fm[["2028"]], year = 2028, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2028 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2028"]] <- fs[["2027"]]
fs[["2028"]][1, ini_dt: end_dt] <-   xfs[["2028"]][[1]][1, ini_dt: end_dt]/fm[["2028"]][1]
fs[["2028"]][2, ini_dt: end_dt] <-   xfs[["2028"]][[1]][2, ini_dt: end_dt]/fm[["2028"]][2]
fs[["2028"]][3, ini_dt: end_dt] <-   xfs[["2028"]][[1]][3, ini_dt: end_dt]/fm[["2028"]][3]
fs[["2028"]][4, ini_dt: end_dt] <-   xfs[["2028"]][[1]][4, ini_dt: end_dt]/fm[["2028"]][4]
fs[["2028"]][5, ini_dt: end_dt] <-   xfs[["2028"]][[1]][5, ini_dt: end_dt]/xfs[["2028"]][[1]][5, ini_dt: end_dt]

# 2029

odd_num_test <- c(0.9, 1)

Ccal[[2029]] <- list("C" = Ccal[[2028]]$C*odd_num_test[1],
                     "P" = Ccal[[2028]]$P*odd_num_test[2])

frac_test[[2029]] <- frac_test[[2028]]

frac_ab[["2029"]] <- c(unlist(as.numeric(frac_test[[2029]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2029]]$P$reflex)))

dfList_NP_2029 <- dfList_NP_2028
for(i in param_var){  
  dfList_NP_2029[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2029, index = i, 
                                   frac_testing = frac_test[[2029]],
                                   S_Yint = 2029, S_Yend = 2030, r_Yend = 2030, NPlst = NPlst, 
                                   fp = c(Ccal[[2029]]$C*fm[["2029"]][1], 
                                          Ccal[[2029]]$C*fm[["2029"]][2], 
                                          Ccal[[2029]]$P*fm[["2029"]][3], 
                                          Ccal[[2029]]$P*fm[["2029"]][4], 
                                          Ccal[[2029]]$P*fm[["2029"]][5]))
  
}

for(i in param_var){
  # begining of 2029
  b_pt <- (2030 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2029[[i]])[3]
  dfList_NP_2029[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2029"]] <- c(Ccal[[2029]]$C, Ccal[[2029]]$P)
n_ab_np[["2029"]] <- n_ab_np[["2028"]]*odd_num_test
xfs[["2029"]] <- fs_estimate(num_ab = n_ab_np[["2029"]], 
                             cov_np = coverage_np[["2029"]], 
                             frac_ab = frac_ab[["2029"]], 
                             fp = fm[["2029"]], year = 2029, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2029 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2029 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2029"]] <- fs[["2028"]]
fs[["2029"]][1, ini_dt: end_dt] <-   xfs[["2029"]][[1]][1, ini_dt: end_dt]/fm[["2029"]][1]
fs[["2029"]][2, ini_dt: end_dt] <-   xfs[["2029"]][[1]][2, ini_dt: end_dt]/fm[["2029"]][2]
fs[["2029"]][3, ini_dt: end_dt] <-   xfs[["2029"]][[1]][3, ini_dt: end_dt]/fm[["2029"]][3]
fs[["2029"]][4, ini_dt: end_dt] <-   xfs[["2029"]][[1]][4, ini_dt: end_dt]/fm[["2029"]][4]

fs[["2029"]][5, ini_dt: end_dt] <-   xfs[["2029"]][[1]][5, ini_dt: end_dt]/xfs[["2029"]][[1]][5, ini_dt: end_dt]


# 2030

odd_num_test <- c(0.9, 1)

Ccal[[2030]] <- list("C" = Ccal[[2029]]$C*odd_num_test[1],
                     "P" = Ccal[[2029]]$P*odd_num_test[2])

frac_test[[2030]] <- frac_test[[2029]]

frac_ab[["2030"]] <- c(unlist(as.numeric(frac_test[[2030]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2030]]$P$reflex)))

Prison_testing_I <- dfList_NP_2029
for(i in param_var){  
  Prison_testing_I[[i]] <- Param_cal(pj = POC_AU, dlist = Prison_testing_I, index = i, 
                                        frac_testing = frac_test[[2030]],
                                        S_Yint = 2030, S_Yend = 2031, r_Yend = 2031, NPlst = NPlst, 
                                        fp = c(Ccal[[2030]]$C*fm[["2030"]][1], 
                                               Ccal[[2030]]$C*fm[["2030"]][2], 
                                               Ccal[[2030]]$P*fm[["2030"]][3], 
                                               Ccal[[2030]]$P*fm[["2030"]][4], 
                                               Ccal[[2030]]$P*fm[["2030"]][5]))
  
}

for(i in param_var){
  # begining of 2029
  b_pt <- (2031 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(Prison_testing_I[[i]])[3]
  Prison_testing_I[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2030"]] <- c(Ccal[[2030]]$C, Ccal[[2030]]$P)
n_ab_np[["2030"]] <- n_ab_np[["2029"]]*odd_num_test
xfs[["2030"]] <- fs_estimate(num_ab = n_ab_np[["2030"]], 
                             cov_np = coverage_np[["2030"]], 
                             frac_ab = frac_ab[["2030"]], 
                             fp = fm[["2030"]], year = 2030, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2030 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2030 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["Prison_testing_I"]] <- fs[["2029"]]
fs[["Prison_testing_I"]][1, ini_dt: end_dt] <-   xfs[["2030"]][[1]][1, ini_dt: end_dt]/fm[["2030"]][1]
fs[["Prison_testing_I"]][2, ini_dt: end_dt] <-   xfs[["2030"]][[1]][2, ini_dt: end_dt]/fm[["2030"]][2]
fs[["Prison_testing_I"]][3, ini_dt: end_dt] <-   xfs[["2030"]][[1]][3, ini_dt: end_dt]/fm[["2030"]][3]
fs[["Prison_testing_I"]][4, ini_dt: end_dt] <-   xfs[["2030"]][[1]][4, ini_dt: end_dt]/fm[["2030"]][4]
fs[["Prison_testing_I"]][5, ini_dt: end_dt] <-   xfs[["2030"]][[1]][5, ini_dt: end_dt]/xfs[["2030"]][[1]][5, ini_dt: end_dt]



#### prison_testing_II ####
# 2027

# 2027

odd_num_test <-  c(1.03, 1.114)

Ccal[[2027]] <- list("C" = Ccal[[2026]]$C*odd_num_test[1],
                     "P" = Ccal[[2026]]$P*odd_num_test[2])

frac_test[[2027]] <- frac_test[[2026]]

frac_ab[["2027"]] <- c(unlist(as.numeric(frac_test[[2027]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2027]]$P$reflex)))

dfList_NP_2027 <- dfList_NPPhaseII
for(i in param_var){  
  dfList_NP_2027[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2027, index = i, 
                                   frac_testing = frac_test[[2027]],
                                   S_Yint = 2027, S_Yend = 2028, r_Yend = 2028, NPlst = NPlst, 
                                   fp = c(Ccal[[2027]]$C*fm[["2027"]][1], 
                                          Ccal[[2027]]$C*fm[["2027"]][2], 
                                          Ccal[[2027]]$P*fm[["2027"]][3], 
                                          Ccal[[2027]]$P*fm[["2027"]][4], 
                                          Ccal[[2027]]$P*fm[["2027"]][5]))
  
}

for(i in param_var){
  # begining of 2027
  b_pt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NPPhaseII[[i]])[3]
  dfList_NP_2027[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2027"]] <- c(Ccal[[2027]]$C, Ccal[[2027]]$P)
n_ab_np[["2027"]] <- n_ab_np[["2026"]]*odd_num_test
xfs[["2027"]] <- fs_estimate(num_ab = n_ab_np[["2027"]], 
                             cov_np = coverage_np[["2027"]], 
                             frac_ab = frac_ab[["2027"]], 
                             fp = fm[["2027"]], year = 2027, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2027 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2027 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2027"]] <- fs[["2026"]]
fs[["2027"]][1, ini_dt: end_dt] <-   xfs[["2027"]][[1]][1, ini_dt: end_dt]/fm[["2027"]][1]
fs[["2027"]][2, ini_dt: end_dt] <-   xfs[["2027"]][[1]][2, ini_dt: end_dt]/fm[["2027"]][2]
fs[["2027"]][3, ini_dt: end_dt] <-   xfs[["2027"]][[1]][3, ini_dt: end_dt]/fm[["2027"]][3]
fs[["2027"]][4, ini_dt: end_dt] <-   xfs[["2027"]][[1]][4, ini_dt: end_dt]/fm[["2027"]][4]
fs[["2027"]][5, ini_dt: end_dt] <-   xfs[["2027"]][[1]][5, ini_dt: end_dt]/xfs[["2027"]][[1]][5, ini_dt: end_dt]





# 2028 
odd_num_test <-  c(1, 0.6)
Ccal[[2028]] <- list("C" = Ccal[[2027]]$C*odd_num_test[1],
                     "P" = Ccal[[2027]]$P*odd_num_test[2])
fm[["2028"]] <- c(0.8, 0.8, 12, 12, 1)

Ccal[[2028]] <- lapply(Ccal[[2027]],function(x) x*odd_num_test)

frac_test[[2028]] <- frac_test[[2026]]

frac_ab[["2028"]] <- c(unlist(as.numeric(frac_test[[2028]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2028]]$P$reflex)))

dfList_NP_2028 <- dfList_NP_2027
for(i in param_var){  
  dfList_NP_2028[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2028, index = i, 
                                   frac_testing = frac_test[[2028]],
                                   S_Yint = 2028, S_Yend = 2029, r_Yend = 2029, NPlst = NPlst, 
                                   fp = c(Ccal[[2028]]$C*fm[["2028"]][1], 
                                          Ccal[[2028]]$C*fm[["2028"]][2], 
                                          Ccal[[2028]]$P*fm[["2028"]][3], 
                                          Ccal[[2028]]$P*fm[["2028"]][4], 
                                          Ccal[[2028]]$P*fm[["2028"]][5]))
  
}

for(i in param_var){
  # begining of 2026
  b_pt <- (2029 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2028[[i]])[3]
  dfList_NP_2028[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2028"]] <- c(Ccal[[2028]]$C, Ccal[[2028]]$P)
n_ab_np[["2028"]] <- n_ab_np[["2027"]]*odd_num_test
xfs[["2028"]] <- fs_estimate(num_ab = n_ab_np[["2028"]], 
                             cov_np = coverage_np[["2028"]], 
                             frac_ab = frac_ab[["2028"]], 
                             fp = fm[["2028"]], year = 2028, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2028 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2028 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2028"]] <- fs[["2027"]]
fs[["2028"]][1, ini_dt: end_dt] <-   xfs[["2028"]][[1]][1, ini_dt: end_dt]/fm[["2028"]][1]
fs[["2028"]][2, ini_dt: end_dt] <-   xfs[["2028"]][[1]][2, ini_dt: end_dt]/fm[["2028"]][2]
fs[["2028"]][3, ini_dt: end_dt] <-   xfs[["2028"]][[1]][3, ini_dt: end_dt]/fm[["2028"]][3]
fs[["2028"]][4, ini_dt: end_dt] <-   xfs[["2028"]][[1]][4, ini_dt: end_dt]/fm[["2028"]][4]
fs[["2028"]][5, ini_dt: end_dt] <-   xfs[["2028"]][[1]][5, ini_dt: end_dt]/xfs[["2028"]][[1]][5, ini_dt: end_dt]

# 2029
# fm[["2029"]] <- c(0.8, 0.8, 10, 10, 1)
# fm[["2030"]] <- c(0.8, 0.8, 10, 10, 1)

odd_num_test <- c(1,1.137)
Ccal[[2029]] <- list("C" = Ccal[[2028]]$C*odd_num_test[1],
                     "P" = Ccal[[2028]]$P*odd_num_test[2])



frac_test[[2029]] <- frac_test[[2028]]

frac_ab[["2029"]] <- c(unlist(as.numeric(frac_test[[2029]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2029]]$P$reflex)))

dfList_NP_2029 <- dfList_NP_2028
for(i in param_var){  
  dfList_NP_2029[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP_2029, index = i, 
                                   frac_testing = frac_test[[2029]],
                                   S_Yint = 2029, S_Yend = 2030, r_Yend = 2030, NPlst = NPlst, 
                                   fp = c(Ccal[[2029]]$C*fm[["2029"]][1], 
                                          Ccal[[2029]]$C*fm[["2029"]][2], 
                                          Ccal[[2029]]$P*fm[["2029"]][3], 
                                          Ccal[[2029]]$P*fm[["2029"]][4], 
                                          Ccal[[2029]]$P*fm[["2029"]][5]))
  
}


for(i in param_var){
  # begining of 2029
  b_pt <- (2030 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(dfList_NP_2029[[i]])[3]
  dfList_NP_2029[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2029"]] <- c(Ccal[[2029]]$C, Ccal[[2029]]$P)
n_ab_np[["2029"]] <- n_ab_np[["2028"]]*odd_num_test
xfs[["2029"]] <- fs_estimate(num_ab = n_ab_np[["2029"]], 
                             cov_np = coverage_np[["2029"]], 
                             frac_ab = frac_ab[["2029"]], 
                             fp = fm[["2029"]], year = 2029, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2029 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2029 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["2029"]] <- fs[["2028"]]
fs[["2029"]][1, ini_dt: end_dt] <-   xfs[["2029"]][[1]][1, ini_dt: end_dt]/fm[["2029"]][1]
fs[["2029"]][2, ini_dt: end_dt] <-   xfs[["2029"]][[1]][2, ini_dt: end_dt]/fm[["2029"]][2]
fs[["2029"]][3, ini_dt: end_dt] <-   xfs[["2029"]][[1]][3, ini_dt: end_dt]/fm[["2029"]][3]
fs[["2029"]][4, ini_dt: end_dt] <-   xfs[["2029"]][[1]][4, ini_dt: end_dt]/fm[["2029"]][4]
fs[["2029"]][5, ini_dt: end_dt] <-   xfs[["2029"]][[1]][5, ini_dt: end_dt]/xfs[["2029"]][[1]][5, ini_dt: end_dt]


# 2030

odd_num_test <- c(1, 1.065)
Ccal[[2030]] <- list("C" = Ccal[[2029]]$C*odd_num_test[1],
                     "P" = Ccal[[2028]]$P*odd_num_test[2])


frac_test[[2030]] <- frac_test[[2029]]

frac_ab[["2030"]] <- c(unlist(as.numeric(frac_test[[2030]]$C$reflex)),
                       unlist(as.numeric(frac_test[[2030]]$P$reflex)))

Prison_testing_II <- dfList_NP_2029
for(i in param_var){  
  Prison_testing_II[[i]] <- Param_cal(pj = POC_AU, dlist = Prison_testing_II, index = i, 
                                        frac_testing = frac_test[[2030]],
                                        S_Yint = 2030, S_Yend = 2031, r_Yend = 2031, NPlst = NPlst, 
                                        fp = c(Ccal[[2030]]$C*fm[["2030"]][1], 
                                               Ccal[[2030]]$C*fm[["2030"]][2], 
                                               Ccal[[2030]]$P*fm[["2030"]][3], 
                                               Ccal[[2030]]$P*fm[["2030"]][4], 
                                               Ccal[[2030]]$P*fm[["2030"]][5]))
  
}

for(i in param_var){
  # begining of 2029
  b_pt <- (2031 - POC_AU$cabY)/POC_AU$timestep + 1 
  
  # length of the time points
  dim_length <- dim(Prison_testing_II[[i]])[3]
  Prison_testing_II[[i]][, , b_pt: dim_length] <- dfList_NP[[i]][, , b_pt: dim_length]
} 

coverage_np[["2030"]] <- c(Ccal[[2030]]$C, Ccal[[2030]]$P)
n_ab_np[["2030"]] <- n_ab_np[["2029"]]*odd_num_test
xfs[["2030"]] <- fs_estimate(num_ab = n_ab_np[["2030"]], 
                             cov_np = coverage_np[["2030"]], 
                             frac_ab = frac_ab[["2030"]], 
                             fp = fm[["2030"]], year = 2030, endY = 100,
                             modsim = Sce_sq)


ini_dt <- (2030 - POC_AU$cabY)/POC_AU$timestep + 1 
end_dt <- ((2030 + 1 ) - POC_AU$cabY)/POC_AU$timestep
# fs[["2024"]][1, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][1, ini_dt:end_dt ]/fm[["2024"]][1]
# fs[["2024"]][2, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][2, ini_dt:end_dt ]/fm[["2024"]][2]
# fs[["2024"]][3, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][3, ini_dt:end_dt ]/fm[["2024"]][3]
# fs[["2024"]][4, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][4, ini_dt:end_dt ]/fm[["2024"]][4]
# fs[["2024"]][5, ini_dt:end_dt ] <-   xfs[["2024"]][[1]][5, ini_dt:end_dt ]/xfs[["2024"]][[1]][5, ini_dt:end_dt ] 

fs[["Prison_testing_II"]] <- fs[["2029"]]
fs[["Prison_testing_II"]][1, ini_dt: end_dt] <-   xfs[["2030"]][[1]][1, ini_dt: end_dt]/fm[["2030"]][1]
fs[["Prison_testing_II"]][2, ini_dt: end_dt] <-   xfs[["2030"]][[1]][2, ini_dt: end_dt]/fm[["2030"]][2]
fs[["Prison_testing_II"]][3, ini_dt: end_dt] <-   xfs[["2030"]][[1]][3, ini_dt: end_dt]/fm[["2030"]][3]
fs[["Prison_testing_II"]][4, ini_dt: end_dt] <-   xfs[["2030"]][[1]][4, ini_dt: end_dt]/fm[["2030"]][4]
fs[["Prison_testing_II"]][5, ini_dt: end_dt] <-   xfs[["2030"]][[1]][5, ini_dt: end_dt]/xfs[["2030"]][[1]][5, ini_dt: end_dt]








#===============================================================================
#
#                               run simulations
#
#===============================================================================
endY <- 100
Sce_np <- list()
tic <- proc.time()
Sce_np <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                disease_progress, pop_array,
                dfList,  
                param_cascade_sc = Prison_testing_II, 
                fib = fib, 
                modelrun = "UN", proj = "POC_AU", end_Y = endY, 
                cost = NULL, costflow = NULL, 
                costflow_Neg = NULL, 
                fc = fs[["Prison_testing_II"]])


   
  
  toc <- proc.time() - tic

  print(toc)


  test <- list()
  cl_ext <- names(Sce_np)[c(10:22)]
  for(i in cl_ext){
    
    test[[i]] <- modres.flow.t(POC_AU, Sce_np, endYear = 100, 
                               allp = i)%>%
      ungroup()%>%
      group_by(year, population)%>%
      summarise(best = sum(best))
    
    
  }
  
  test_sq <- list()
  for(i in cl_ext){
    
    test_sq[[i]] <- modres.flow.t(POC_AU, Sce_sq, endYear = 100, 
                                  allp = i)%>%
      ungroup()%>%
      group_by(year, population)%>%
      summarise(best = sum(best))
    
    
  }
  
  
  test <- dplyr::bind_rows(test, .id = 'index')%>%group_by(year, population)%>%spread(index, best)
  test_sq <- dplyr::bind_rows(test_sq, .id = 'index')%>%group_by(year, population)%>%spread(index, best)
  
  
  test_fscal <- test%>%mutate(Ab = newTestingAb_sc + newTestingAb_sc_neg, 
                              RNA = (newTestingAg_sc+ newTestingAg_sc_neg + newTestingPOCT_sc + 
                                       newTestingPOCT_sc_neg ))%>%select(year, population, Ab, 
                                                                         RNA)%>%
    mutate(setting = ifelse(population %in% c("C_PWID", "C_fPWID"), "C", "P"))%>%
    ungroup()%>%
    select(-c(population))%>%
    gather(index, value, -c(year, setting))%>%
    group_by(year, setting, index)%>%summarise(value = sum(value))%>%
    mutate(scenario = "prisons_testing_I")
  dtp <- c(6667,6667, 11476,11476,20393, 20393, 12000, 13000, 11000, 14000, 10000, 15000,9000, 16000,8000, 17000,7000, 18000 )
  # View(test_fscal%>%filter(year%in% c(7:15)))
  View(test_fscal%>%filter(year%in% c(7:15))%>%group_by(year, setting)%>%
         summarise(tot_NP_test = sum(value))%>%
         mutate(year = year + 2015))
  test%>%mutate(setting = ifelse(population %in% POC_AU$popNames[3:5], "P", "C"))%>%
    group_by(setting, year)%>%summarise(best = sum(newTreatment_sc))%>%
    filter(year%in% c(7:12))


