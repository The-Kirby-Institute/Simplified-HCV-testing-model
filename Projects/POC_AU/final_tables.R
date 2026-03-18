#### manuscript numbers and tables #### 
library(openxlsx)
library(readxl)
dt_path <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/POC_AU/02. Output"
Prev_dt <- read_excel(file.path(paste0(dt_path, "/PrevInc_epi.xlsx")), sheet = "tempPrevRNA_setting")
View(Prev_dt)
Prev_dt_tab <- Prev_dt%>%filter(setting != "prisonsPWID")%>%
  mutate(year = (year + POC_AU$cabY - 1),
         scenario = factor(scenario, levels = sce_level,
                           labels = sce_label),
         setting = factor(setting, levels = c("commu", "prisons"),
                          labels = c("Community", "Prison"))
         )

Prev_dt_tab <- Prev_dt_tab%>%select(scenario, setting, year, best, q5, q95)%>%
  filter(year <= 2030)%>%
  mutate(best = round(best, digits = 1), 
         q5 = round(q5, digits = 1),
         q95 = round(q95, digits = 1))

Prev_dt_tab$scenario
scenario_order <- c("(1) No national program", "(2) Foundational implementation", 
                    "(3) Program succession", "(4) Program sustained", "(5) Program scale-up")

# Get reference values (No National Program)
ref_data <- Prev_dt_tab %>%
  filter(scenario == scenario_order[1]) %>%
  filter(year <= 2030)%>%
  select(setting, year, ref_best = best, ref_q5 = q5, ref_q95 = q95)

# Calculate prevalence and reduction
table_data <- Prev_dt_tab %>% filter(year <= 2030)%>%
  left_join(ref_data, by = c("setting", "year")) %>%
  mutate(
    # Prevalence with 95% CrI
    prevalence = sprintf("%.1f (%.1f, %.1f)", best, q5, q95),
    
    # Relative reduction (%)
    red_best = (ref_best- best) / ref_best * 100,

    
    # Reduction with 95% CrI (show "-" for reference)
    reduction = if_else(
      scenario == "No national program",
      "-",
      sprintf("%.1f", red_best)
    ))%>%
  select(setting, year, scenario, prevalence, reduction) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(prevalence, reduction),
    names_glue = "{scenario}_{.value}"
  ) %>%
  arrange(setting, year)






library(tidyverse)
library(gt)
table_data %>%
  gt(groupname_col = "setting") %>%
  tab_header(
    title = "HCV RNA prevalence"
  ) %>%
  
  # Spanners for each scenario
  tab_spanner(label = "(1) No National Program", 
              columns = starts_with(scenario_order[1])) %>%
  tab_spanner(label = "(2) Foundational Implementation", 
              columns = starts_with(scenario_order[2])) %>%
  tab_spanner(label = "(3) Program Succession", 
              columns = starts_with(scenario_order[3])) %>%
  tab_spanner(label = "(4) Program Sustained", 
              columns = starts_with(scenario_order[4])) %>%
  tab_spanner(label = "(5) Program scale-up", 
              columns = starts_with(scenario_order[5])) %>%
  
  # Rename columns
  cols_label(
    year = "Year",
    `(1) No national program_prevalence` = "Prevalence",
    `(1) No national program_reduction` = "Reduction (%)",
    `(2) Foundational implementation_prevalence` = "Prevalence",
    `(2) Foundational implementation_reduction` = "Reduction (%)",
    `(3) Program succession_prevalence` = "Prevalence",
    `(3) Program succession_reduction` = "Reduction (%)",
    `(4) Program sustained_prevalence` = "Prevalence",
    `(4) Program sustained_reduction` = "Reduction (%)",
    `(5) Program scale-up_prevalence` = "Prevalence",
    `(5) Program scale-up_reduction` = "Reduction (%)"
  ) %>%
  
  # Reorder columns
  cols_move(
    columns = c(
      `(1) No national program_prevalence`, `(1) No national program_reduction`,
      `(2) Foundational implementation_prevalence`, `(2) Foundational implementation_reduction`,
      `(3) Program succession_prevalence`, `(3) Program succession_reduction`,
      `(4) Program sustained_prevalence`, `(4) Program sustained_reduction`,
      `(5) Program scale-up_prevalence`, `(5) Program scale-up_reduction`
    ),
    after = year
  ) %>%
  
  # Styling
  cols_align(align = "center", columns = -year) %>%
  cols_align(align = "left", columns = year) %>%
  
  tab_footnote(
    footnote = "Values are presented as estimate (95% credible interval). Reduction is relative change compared to No National Program."
  ) %>%
  
  tab_options(
    table.width = pct(100),
    table.font.size = px(11),
    heading.title.font.size = px(14),
    heading.title.font.weight = "bold",
    heading.align = "left",
    column_labels.font.weight = "bold",
    column_labels.border.top.width = px(2),
    column_labels.border.top.color = "black",
    column_labels.border.bottom.width = px(1),
    column_labels.border.bottom.color = "black",
    row_group.font.weight = "bold",
    row_group.border.top.width = px(1),
    row_group.border.bottom.width = px(0),
    table_body.border.bottom.width = px(2),
    table_body.border.bottom.color = "black",
    table_body.hlines.width = px(0),
    footnotes.font.size = px(10)
  )%>%
  gtsave(., file = file.path(OutputFig, paste0("HCV_RNA_prevalence_scenarios.docx")))


#### new infection #### 
# files <- list.files(
#  path = file.path(dt_path,"Figs" ,urrTime),
#  pattern = "Resflow_year_all_range_.*\\.xlsx$",
#  full.names = TRUE  # returns full path
#)

#flow_dt <- lapply(files, function(f) {
#  read_excel(f, sheet = "newInfections")
# })
# names(flow_dt) <- c(sce_label[2:5], sce_label[1])

# flow_dt <- bind_rows(flow_dt, .id = "scenario")
View(Resflow_all_lst$Resflow_cum_avert$newInfections)

ref_data <- Resflow_all_lst$Resflow_year$newInfections%>%
  filter(scenario == scenario_order[1]) %>%
  filter(year <= 2030)%>%
  mutate(ref_best = best,
         ref_q5 = q5, 
         ref_q95 = q95)%>%
  arrange(year) %>%
  select(year, ref_best, ref_q5, ref_q95)

# Calculate prevalence and reduction
table_data <- Resflow_all_lst$Resflow_year$newInfections  %>% filter(year <= 2030) %>%
  
  arrange(scenario, year) %>%group_by(scenario)%>%
 
  ungroup() %>%
  left_join(ref_data, by = "year") %>%
  mutate(
    # Prevalence with 95% CrI
    newinfection = sprintf("%.0f (%.0f, %.0f)", best, q5, q95),
    
    # Relative reduction (%)
    red_best = (ref_best- best) / ref_best * 100,
    
    
    # Reduction with 95% CrI (show "-" for reference)
    reduction = if_else(
      scenario == scenario_order[1],
      "-",
      sprintf("%.0f", red_best)
    ))%>%
  select(year, scenario, newinfection, reduction) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(newinfection, reduction),
    names_glue = "{scenario}_{.value}"
  ) %>%
  arrange(year)

View(table_data)




library(tidyverse)
library(gt)
table_data %>%
  gt() %>%
  tab_header(
    title = "HCV new infections"
  ) %>%
  
  # Spanners for each scenario
  tab_spanner(label = "(1) No National Program", 
              columns = starts_with(scenario_order[1])) %>%
  tab_spanner(label = "(2) Foundational Implementation", 
              columns = starts_with(scenario_order[2])) %>%
  tab_spanner(label = "(3) Program Succession", 
              columns = starts_with(scenario_order[3])) %>%
  tab_spanner(label = "(4) Program Sustained", 
              columns = starts_with(scenario_order[4])) %>%
  tab_spanner(label = "(5) Program scale-up", 
              columns = starts_with(scenario_order[5])) %>%
  
  # Rename columns
  cols_label(
    year = "Year",
    `(1) No national program_newinfection` = "HCV new infections",
    `(1) No national program_reduction` = "Reduction (%)",
    `(2) Foundational implementation_newinfection` = "HCV new infections",
    `(2) Foundational implementation_reduction` = "Reduction (%)",
    `(3) Program succession_newinfection` = "HCV new infections",
    `(3) Program succession_reduction` = "Reduction (%)",
    `(4) Program sustained_newinfection` = "HCV new infections",
    `(4) Program sustained_reduction` = "Reduction (%)",
    `(5) Program scale-up_newinfection` = "HCV new infections",
    `(5) Program scale-up_reduction` = "Reduction (%)"
  ) %>%
  
  # Reorder columns
  cols_move(
    columns = c(
      `(1) No national program_newinfection`, `(1) No national program_reduction`,
      `(2) Foundational implementation_newinfection`, `(2) Foundational implementation_reduction`,
      `(3) Program succession_newinfection`, `(3) Program succession_reduction`,
      `(4) Program sustained_newinfection`, `(4) Program sustained_reduction`,
      `(5) Program scale-up_newinfection`, `(5) Program scale-up_reduction`
    ),
    after = year
  ) %>%
  
  # Styling
  cols_align(align = "center", columns = -year) %>%
  cols_align(align = "left", columns = year) %>%
  
  tab_footnote(
    footnote = "Values are presented as estimate (95% credible interval). Reduction is relative change compared to No National Program."
  ) %>%
  
  tab_options(
    table.width = pct(100),
    table.font.size = px(11),
    heading.title.font.size = px(14),
    heading.title.font.weight = "bold",
    heading.align = "left",
    column_labels.font.weight = "bold",
    column_labels.border.top.width = px(2),
    column_labels.border.top.color = "black",
    column_labels.border.bottom.width = px(1),
    column_labels.border.bottom.color = "black",
    row_group.font.weight = "bold",
    row_group.border.top.width = px(1),
    row_group.border.bottom.width = px(0),
    table_body.border.bottom.width = px(2),
    table_body.border.bottom.color = "black",
    table_body.hlines.width = px(0),
    footnotes.font.size = px(10)
  )%>%
  gtsave(., file = file.path(OutputFig, paste0("HCV_newinfections_scenarios.docx"))) 

#### screened #### 
col_est <- c("best", paste0("set", seq(1,1000,1)))
num_screened <-data.frame()
num_screened <- xt_screened%>%group_by(year,scenario)%>%
  dplyr::mutate(across(all_of(c("best",set_cols)),sum, na.rm = TRUE))%>%
  mutate(scenario = factor(scenario, level = unique(xt_screened$scenario), label = scenario_order))%>%
  slice(1)%>%
  ungroup()%>%select(!c(population, NP))%>%mutate(across(all_of(set_cols), ~na_if(., 0)))


ref_data <- num_screened%>%
  arrange(scenario, year)%>%
  filter(scenario == scenario_order[1]) %>%
  filter(year <= 2030)%>%
  rowwise() %>%
  mutate(
    min    = min(c_across(starts_with("set")),na.rm = TRUE),
    max    = max(c_across(starts_with("set")),na.rm = TRUE),
    median = median(c_across(starts_with("set")),na.rm = TRUE),
    mean   = mean(c_across(starts_with("set")),na.rm = TRUE),
    q2.5   = quantile(c_across(starts_with("set")), 0.025,na.rm = TRUE),
    q97.5  = quantile(c_across(starts_with("set")), 0.975,na.rm = TRUE)
  )%>% 
  select(year, ref_best = best, ref_q5 = q2.5, ref_q95 = q97.5)

table_data <- num_screened%>% filter(year <= 2030) %>%
  
  arrange(scenario, year) %>%group_by(scenario)%>%
  rowwise() %>%
  mutate(
    min    = min(c_across(starts_with("set")),na.rm = TRUE),
    max    = max(c_across(starts_with("set")),na.rm = TRUE),
    median = median(c_across(starts_with("set")),na.rm = TRUE),
    mean   = mean(c_across(starts_with("set")),na.rm = TRUE),
    q2.5   = quantile(c_across(starts_with("set")), 0.025,na.rm = TRUE),
    q97.5  = quantile(c_across(starts_with("set")), 0.975,na.rm = TRUE)
  )%>% 
  
  ungroup() %>%
  left_join(ref_data, by = "year") %>%
  mutate(
    
    newinfection = sprintf("%.0f (%.0f, %.0f)", best, q2.5, q97.5),
    
    # Relative increase (%)
    red_best = (best - ref_best) / ref_best * 100,
    
    
    #increa with 95% CrI (show "-" for reference)
    increase = if_else(
      scenario == scenario_order[1],
      "-",
      sprintf("%.0f", red_best)
    ))%>%
  select(year, scenario, newinfection, increase ) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(newinfection, increase),
    names_glue = "{scenario}_{.value}"
  ) %>%
  arrange(year)

View(table_data)

table_data %>%
  gt() %>%
  tab_header(
    title = "HCV screened"
  ) %>%
  
  # Spanners for each scenario
  tab_spanner(label = "(1) No National Program", 
              columns = starts_with(scenario_order[1])) %>%
  tab_spanner(label = "(2) Foundational Implementation", 
              columns = starts_with(scenario_order[2])) %>%
  tab_spanner(label = "(3) Program Succession", 
              columns = starts_with(scenario_order[3])) %>%
  tab_spanner(label = "(4) Program Sustained", 
              columns = starts_with(scenario_order[4])) %>%
  tab_spanner(label = "(5) Program scale-up", 
              columns = starts_with(scenario_order[5])) %>%
  
  # Rename columns
  cols_label(
    year = "Year",
    `(1) No National Program_newinfection` = "HCV screened",
    `(1) No National Program_increase` = "Increase (%)",
    `(2) Foundational implementation_newinfection` = "HCV screened",
    `(2) Foundational implementation_increase` = "Increase (%)",
    `(3) Program succession_newinfection` = "HCV screened",
    `(3) Program succession_increase` = "Increase (%)",
    `(4) Program sustained_newinfection` = "HCV screened",
    `(4) Program sustained_increase` = "Increase (%)",
    `(5) Program scale-up_newinfection` = "HCV screened",
    `(5) Program scale-up_increase` = "Increase (%)"
  ) %>%
  
  # Reorder columns
  cols_move(
    columns = c(
      `(1) No National Program_newinfection`, `(1) No National Program_increase`,
      `(2) Foundational implementation_newinfection`, `(2) Foundational implementation_increase`,
      `(3) Program succession_newinfection`, `(3) Program succession_increase`,
      `(4) Program sustained_newinfection`, `(4) Program sustained_increase`,
      `(5) Program scale-up_newinfection`, `(5) Program scale-up_increase`
    ),
    after = year
  ) %>%
  
  # Styling
  cols_align(align = "center", columns = -year) %>%
  cols_align(align = "left", columns = year) %>%
  
  tab_footnote(
    footnote = "Values are presented as estimate (95% credible interval). Reduction is relative change compared to No National Program."
  ) %>%
  
  tab_options(
    table.width = pct(100),
    table.font.size = px(11),
    heading.title.font.size = px(14),
    heading.title.font.weight = "bold",
    heading.align = "left",
    column_labels.font.weight = "bold",
    column_labels.border.top.width = px(2),
    column_labels.border.top.color = "black",
    column_labels.border.bottom.width = px(1),
    column_labels.border.bottom.color = "black",
    row_group.font.weight = "bold",
    row_group.border.top.width = px(1),
    row_group.border.bottom.width = px(0),
    table_body.border.bottom.width = px(2),
    table_body.border.bottom.color = "black",
    table_body.hlines.width = px(0),
    footnotes.font.size = px(10)
  )%>%
  gtsave(., file = file.path(OutputFig, paste0("HCV_screened_scenarios.docx"))) 

#### screen cumulative ####

col_est <- c("best", paste0("set", seq(1,1000,1)))

num_screened_cum <- num_screened%>%filter(year >= 2022)%>%
  group_by(scenario)%>%
  dplyr::mutate(across(all_of(col_est), ~ cumsum(ifelse(is.na(.x), 0, .x))))%>%
  ungroup()%>%mutate(across(all_of(col_est), ~na_if(., 0)))


ref_data <- num_screened_cum%>%
  arrange(scenario, year)%>%
  filter(scenario == scenario_order[1]) %>%
  filter(year <= 2030)%>%
  rowwise() %>%
  mutate(
    min    = min(c_across(starts_with("set")),na.rm = TRUE),
    max    = max(c_across(starts_with("set")),na.rm = TRUE),
    median = median(c_across(starts_with("set")),na.rm = TRUE),
    mean   = mean(c_across(starts_with("set")),na.rm = TRUE),
    q2.5   = quantile(c_across(starts_with("set")), 0.025,na.rm = TRUE),
    q97.5  = quantile(c_across(starts_with("set")), 0.975,na.rm = TRUE)
  )%>% 
  select(year, ref_best = best, ref_q5 = q2.5, ref_q95 = q97.5)

table_data <- num_screened_cum%>% filter(year <= 2030) %>%
  
  arrange(scenario, year) %>%group_by(scenario)%>%
  rowwise() %>%
  mutate(
    min    = min(c_across(starts_with("set")),na.rm = TRUE),
    max    = max(c_across(starts_with("set")),na.rm = TRUE),
    median = median(c_across(starts_with("set")),na.rm = TRUE),
    mean   = mean(c_across(starts_with("set")),na.rm = TRUE),
    q2.5   = quantile(c_across(starts_with("set")), 0.025,na.rm = TRUE),
    q97.5  = quantile(c_across(starts_with("set")), 0.975,na.rm = TRUE)
  )%>% 
  
  ungroup() %>%
  left_join(ref_data, by = "year") %>%
  mutate(
    
    newinfection = sprintf("%.0f (%.0f, %.0f)", best, q2.5, q97.5),
    
    # Relative increase (%)
    red_best = (best - ref_best) / ref_best * 100,
    
    
    #increa with 95% CrI (show "-" for reference)
    increase = if_else(
      scenario == scenario_order[1],
      "-",
      sprintf("%.0f", red_best)
    ))%>%
  select(year, scenario, newinfection, increase ) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(newinfection, increase),
    names_glue = "{scenario}_{.value}"
  ) %>%
  arrange(year)

View(table_data)

table_data %>%
  gt() %>%
  tab_header(
    title = "HCV screened (cumulative)"
  ) %>%
  
  # Spanners for each scenario
  tab_spanner(label = "(1) No National Program", 
              columns = starts_with(scenario_order[1])) %>%
  tab_spanner(label = "(2) Foundational Implementation", 
              columns = starts_with(scenario_order[2])) %>%
  tab_spanner(label = "(3) Program Succession", 
              columns = starts_with(scenario_order[3])) %>%
  tab_spanner(label = "(4) Program Sustained", 
              columns = starts_with(scenario_order[4])) %>%
  tab_spanner(label = "(5) Program scale-up", 
              columns = starts_with(scenario_order[5])) %>%
  
  # Rename columns
  cols_label(
    year = "Year",
    `(1) No National Program_newinfection` = "HCV screened (cumulative)",
    `(1) No National Program_increase` = "Increase (%)",
    `(2) Foundational implementation_newinfection` = "HCV screened (cumulative)",
    `(2) Foundational implementation_increase` = "Increase (%)",
    `(3) Program succession_newinfection` = "HCV screened (cumulative)",
    `(3) Program succession_increase` = "Increase (%)",
    `(4) Program sustained_newinfection` = "HCV screened (cumulative)",
    `(4) Program sustained_increase` = "Increase (%)",
    `(5) Program scale-up_newinfection` = "HCV screened (cumulative)",
    `(5) Program scale-up_increase` = "Increase (%)"
  ) %>%
  
  # Reorder columns
  cols_move(
    columns = c(
      `(1) No National Program_newinfection`, `(1) No National Program_increase`,
      `(2) Foundational implementation_newinfection`, `(2) Foundational implementation_increase`,
      `(3) Program succession_newinfection`, `(3) Program succession_increase`,
      `(4) Program sustained_newinfection`, `(4) Program sustained_increase`,
      `(5) Program scale-up_newinfection`, `(5) Program scale-up_increase`
    ),
    after = year
  ) %>%
  
  # Styling
  cols_align(align = "center", columns = -year) %>%
  cols_align(align = "left", columns = year) %>%
  
  tab_footnote(
    footnote = "Values are presented as estimate (95% credible interval). Reduction is relative change compared to No National Program."
  ) %>%
  
  tab_options(
    table.width = pct(100),
    table.font.size = px(11),
    heading.title.font.size = px(14),
    heading.title.font.weight = "bold",
    heading.align = "left",
    column_labels.font.weight = "bold",
    column_labels.border.top.width = px(2),
    column_labels.border.top.color = "black",
    column_labels.border.bottom.width = px(1),
    column_labels.border.bottom.color = "black",
    row_group.font.weight = "bold",
    row_group.border.top.width = px(1),
    row_group.border.bottom.width = px(0),
    table_body.border.bottom.width = px(2),
    table_body.border.bottom.color = "black",
    table_body.hlines.width = px(0),
    footnotes.font.size = px(10)
  )%>%
  gtsave(., file = file.path(OutputFig, paste0("HCV_screened_cumulative_scenarios.docx"))) 


#### treatment #### 

ref_data <- Resflow_all_lst$Resflow_year$Tot_Treatment %>%
  filter(scenario == scenario_order[1]) %>%
  filter(year <= 2030)%>%
  select(year, ref_best = best, ref_q5 = q5, ref_q95 = q95)

ref_data_all <- Resflow_all_lst$Resflow_year$Tot_Treatment %>%
  filter(year <= 2030, scenario == scenario_order[1]) %>%
  select(year, best, starts_with("set"))
alt_scn <- unique(Resflow_all_lst$Resflow_year$Tot_Treatment$scenario)[4]
ref_scn <- unique(Resflow_all_lst$Resflow_year$Tot_Treatment$scenario)[1]
draw_cols <- c(paste0("set", seq(1,1000,1)))  # set1..set1000
cols_to_diff <- c("best", draw_cols)

treat_data_incre <- Resflow_all_lst$Resflow_year$Tot_Treatment%>%filter(scenario == alt_scn & year <= 2030)
inc_data <- cbind(
  year = treat_data_incre$year,
  scenario = alt_scn, 
  as.data.frame(treat_data_incre[,cols_to_diff] - ref_data_all[, cols_to_diff]))%>%
  rowwise() %>%
  mutate(
    inc_q2.5  = quantile(c_across(all_of(draw_cols)), 0.025, na.rm = TRUE),
    inc_q97.5 = quantile(c_across(all_of(draw_cols)), 0.975, na.rm = TRUE)
  ) %>%
  ungroup() 

inc_data <- inc_data%>% mutate(incre_treat =  sprintf("%.0f (%.0f, %.0f)", best, inc_q2.5, inc_q97.5))
View(inc_data%>%select(incre_treat))
# Calculate prevalence and reduction
table_data <- Resflow_all_lst$Resflow_year$Tot_Treatment %>% filter(year <= 2030)%>%
  arrange(scenario, year) %>%group_by(scenario)%>%

  left_join(ref_data, by = c( "year")) %>%
  mutate(
    # Prevalence with 95% CrI
    newinfection = sprintf("%.0f (%.0f, %.0f)", best, q5, q95),
    
    # Relative reduction (%)
    red_best = ( best - ref_best) / ref_best * 100,
    
    
    # Reduction with 95% CrI (show "-" for reference)
    reduction = if_else(
      scenario == "No national program",
      "-",
      sprintf("%.0f", red_best)
    ))%>%
  select(year, scenario, newinfection, reduction) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(newinfection, reduction),
    names_glue = "{scenario}_{.value}"
  ) %>%
  arrange(year)






library(tidyverse)
library(gt)
table_data %>%
  gt() %>%
  tab_header(
    title = "Treatment initiation"
  ) %>%
  
  # Spanners for each scenario
  tab_spanner(label = "(1) No National Program", 
              columns = starts_with(scenario_order[1])) %>%
  tab_spanner(label = "(2) Foundational Implementation", 
              columns = starts_with(scenario_order[2])) %>%
  tab_spanner(label = "(3) Program Succession", 
              columns = starts_with(scenario_order[3])) %>%
  tab_spanner(label = "(4) Program Sustained", 
              columns = starts_with(scenario_order[4])) %>%
  tab_spanner(label = "(5) Program scale-up", 
              columns = starts_with(scenario_order[5])) %>%
  
  # Rename columns
  cols_label(
    year = "Year",
    `(1) No national program_newinfection` = "Treatment initiation",
    `(1) No national program_reduction` = "Increase (%)",
    `(2) Foundational implementation_newinfection` = "Treatment initiation",
    `(2) Foundational implementation_reduction` = "Increase (%)",
    `(3) Program succession_newinfection` = "Treatment initiation",
    `(3) Program succession_reduction` = "Increase (%)",
    `(4) Program sustained_newinfection` = "Treatment initiation",
    `(4) Program sustained_reduction` = "Increase (%)",
    `(5) Program scale-up_newinfection` = "Treatment initiation",
    `(5) Program scale-up_reduction` = "Increase (%)"
  ) %>%
  
  # Reorder columns
  cols_move(
    columns = c(
      `(1) No national program_newinfection`, `(1) No national program_reduction`,
      `(2) Foundational implementation_newinfection`, `(2) Foundational implementation_reduction`,
      `(3) Program succession_newinfection`, `(3) Program succession_reduction`,
      `(4) Program sustained_newinfection`, `(4) Program sustained_reduction`,
      `(5) Program scale-up_newinfection`, `(5) Program scale-up_reduction`
    ),
    after = year
  ) %>%
  
  # Styling
  cols_align(align = "center", columns = -year) %>%
  cols_align(align = "left", columns = year) %>%
  
  tab_footnote(
    footnote = "Values are presented as estimate (95% credible interval). Reduction is relative change compared to No National Program."
  ) %>%
  
  tab_options(
    table.width = pct(100),
    table.font.size = px(11),
    heading.title.font.size = px(14),
    heading.title.font.weight = "bold",
    heading.align = "left",
    column_labels.font.weight = "bold",
    column_labels.border.top.width = px(2),
    column_labels.border.top.color = "black",
    column_labels.border.bottom.width = px(1),
    column_labels.border.bottom.color = "black",
    row_group.font.weight = "bold",
    row_group.border.top.width = px(1),
    row_group.border.bottom.width = px(0),
    table_body.border.bottom.width = px(2),
    table_body.border.bottom.color = "black",
    table_body.hlines.width = px(0),
    footnotes.font.size = px(10)
  )%>%
  gtsave(., file = file.path(OutputFig, paste0("HCV_treatment_scenarios.docx"))) 




#### HCV death ####
ref_data <- Resflow_all_lst$Resflow_year$HCVdeath%>%
  filter(scenario == scenario_order[1]) %>%
  filter(year <= 2030)%>%
  mutate(ref_best = best,
         ref_q5 = q5, 
         ref_q95 = q95)%>%
  arrange(year) %>%
  select(year, ref_best, ref_q5, ref_q95)

# Calculate prevalence and reduction
table_data <- Resflow_all_lst$Resflow_year$HCVdeath  %>% filter(year <= 2030) %>%
  
  arrange(scenario, year) %>%group_by(scenario)%>%
  
  ungroup() %>%
  left_join(ref_data, by = "year") %>%
  mutate(
    # Prevalence with 95% CrI
    newinfection = sprintf("%.0f (%.0f, %.0f)", best, q5, q95),
    
    # Relative reduction (%)
    red_best = (ref_best- best) / ref_best * 100,
    
    
    # Reduction with 95% CrI (show "-" for reference)
    reduction = if_else(
      scenario == scenario_order[1],
      "-",
      sprintf("%.0f", red_best)
    ))%>%
  select(year, scenario, newinfection, reduction) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(newinfection, reduction),
    names_glue = "{scenario}_{.value}"
  ) %>%
  arrange(year)

View(table_data)




library(tidyverse)
library(gt)
table_data %>%
  gt() %>%
  tab_header(
    title = "HCV deaths"
  ) %>%
  
  # Spanners for each scenario
  tab_spanner(label = "(1) No National Program", 
              columns = starts_with(scenario_order[1])) %>%
  tab_spanner(label = "(2) Foundational Implementation", 
              columns = starts_with(scenario_order[2])) %>%
  tab_spanner(label = "(3) Program Succession", 
              columns = starts_with(scenario_order[3])) %>%
  tab_spanner(label = "(4) Program Sustained", 
              columns = starts_with(scenario_order[4])) %>%
  tab_spanner(label = "(5) Program scale-up", 
              columns = starts_with(scenario_order[5])) %>%
  
  # Rename columns
  cols_label(
    year = "Year",
    `(1) No national program_newinfection` = "HCV deaths",
    `(1) No national program_reduction` = "Reduction (%)",
    `(2) Foundational implementation_newinfection` = "HCV deaths",
    `(2) Foundational implementation_reduction` = "Reduction (%)",
    `(3) Program succession_newinfection` = "HCV deaths",
    `(3) Program succession_reduction` = "Reduction (%)",
    `(4) Program sustained_newinfection` = "HCV deaths",
    `(4) Program sustained_reduction` = "Reduction (%)",
    `(5) Program scale-up_newinfection` = "HCV deaths",
    `(5) Program scale-up_reduction` = "Reduction (%)"
  ) %>%
  
  # Reorder columns
  cols_move(
    columns = c(
      `(1) No national program_newinfection`, `(1) No national program_reduction`,
      `(2) Foundational implementation_newinfection`, `(2) Foundational implementation_reduction`,
      `(3) Program succession_newinfection`, `(3) Program succession_reduction`,
      `(4) Program sustained_newinfection`, `(4) Program sustained_reduction`,
      `(5) Program scale-up_newinfection`, `(5) Program scale-up_reduction`
    ),
    after = year
  ) %>%
  
  # Styling
  cols_align(align = "center", columns = -year) %>%
  cols_align(align = "left", columns = year) %>%
  
  tab_footnote(
    footnote = "Values are presented as estimate (95% credible interval). Reduction is relative change compared to No National Program."
  ) %>%
  
  tab_options(
    table.width = pct(100),
    table.font.size = px(11),
    heading.title.font.size = px(14),
    heading.title.font.weight = "bold",
    heading.align = "left",
    column_labels.font.weight = "bold",
    column_labels.border.top.width = px(2),
    column_labels.border.top.color = "black",
    column_labels.border.bottom.width = px(1),
    column_labels.border.bottom.color = "black",
    row_group.font.weight = "bold",
    row_group.border.top.width = px(1),
    row_group.border.bottom.width = px(0),
    table_body.border.bottom.width = px(2),
    table_body.border.bottom.color = "black",
    table_body.hlines.width = px(0),
    footnotes.font.size = px(10)
  )%>%
  gtsave(., file = file.path(OutputFig, paste0("HCV_HCVdeath_scenarios.docx"))) 










#### CEA tables####
set_cols <- paste0("set", 1:1000)   # or: grep("^set\\d+$", names(df), value = TRUE)

x_catcost_CEA <- lapply(list(cost_disydaacap_categories_bind,cost_disydaanocap_categories_bind), function(x) x%>%
                          filter(year>= 2022)%>%
                          group_by(scenario, sensitivity, Categories)%>%
                          mutate(across(c(par_col, "min", "max", "Med", "Mu", "q5", 
                                          "q25", "q75", "q95"), cumsum, .names = "{col}"))%>%ungroup()%>%
                          arrange(scenario, sensitivity))

names(x_catcost_CEA) <- c("discount_cap", "discount_nocap")
View(x_catcost_CEA$discount_cap)
CEA_cost <- x_catcost_CEA$discount_cap%>%group_by(scenario, sensitivity, year)%>%
  filter(year == 2041)%>%
  mutate(
    across(c("best", all_of(set_cols)), ~ sum(.x))
  ) %>%
  mutate(
    min    = min(c_across(all_of(set_cols)),na.rm = TRUE),
    max    = max(c_across(all_of(set_cols)),na.rm = TRUE),
    median = median(c_across(all_of(set_cols)),na.rm = TRUE),
    mean   = mean(c_across(all_of(set_cols)),na.rm = TRUE),
    q2.5   = quantile(c_across(all_of(set_cols)), 0.025,na.rm = TRUE),
    q97.5  = quantile(c_across(all_of(set_cols)), 0.975,na.rm = TRUE)
  )

CEA_cost <- CEA_cost%>%select(!Categories)%>%slice(1)%>%
  select(sensitivity, scenario, year, best, all_of(set_cols), min, max, median, mean, q2.5, q97.5)

bench_scn <- "(1) No national program"

CEA_incre <- CEA_cost  %>%
  group_by(sensitivity, year) %>%
  left_join(
    CEA_cost  %>%
      filter(scenario == bench_scn) %>%
      select(sensitivity, year, best, all_of(set_cols)) %>%
      rename(best_bench = best) %>%
      rename_with(~ paste0(.x, "_bench"), all_of(set_cols)),
    by = c("sensitivity", "year")
  ) %>%
  mutate(
    inc_best = best - best_bench,
    across(
      all_of(set_cols),
      ~ .x - get(paste0(cur_column(), "_bench")),
      .names = "inc_{col}"
    )
  ) %>%
  rowwise() %>%
  mutate(
    inc_min    = min(c_across(starts_with("inc_set")),na.rm = TRUE),
    inc_max    = max(c_across(starts_with("inc_set")),na.rm = TRUE),
    inc_median = median(c_across(starts_with("inc_set")),na.rm = TRUE),
    inc_mean   = mean(c_across(starts_with("inc_set")),na.rm = TRUE),
    inc_q2.5   = quantile(c_across(starts_with("inc_set")), 0.025,na.rm = TRUE),
    inc_q97.5  = quantile(c_across(starts_with("inc_set")), 0.975,na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(-all_of(set_cols), -ends_with("_bench"))

View(CEA_incre)


CEA_incre <- CEA_incre%>%
  select(sensitivity, scenario.x, best, q2.5, q97.5, inc_best, paste0("inc_set", seq(1,1000,1)), 
         inc_min, inc_max, inc_median, inc_mean, inc_q2.5, inc_q97.5)%>%
  rename(scenario = scenario.x)

#### QAYL ####
QALY_flat <- map_dfr(names(CEAanalysis), function(sens) {
  map_dfr(names(CEAanalysis[[sens]][["20y"]][["QALY"]]), function(scn) {
    CEAanalysis[[sens]][["20y"]][["QALY"]][[scn]] %>%
      mutate(sensitivity = sens, scenario = scn)
  })
})
QALY_incre <- QALY_flat %>%
  left_join(
    QALY_flat %>%
      filter(scenario == bench_scn) %>%
      select(sensitivity, year, best, all_of(set_cols)) %>%
      rename_with(~ paste0(.x, "_bench"), c("best", all_of(set_cols))),
    by = c("sensitivity", "year")
  ) %>%
  mutate(
    inc_best = best - best_bench,
    across(
      all_of(set_cols),
      ~ .x - get(paste0(cur_column(), "_bench")),
      .names = "inc_{col}"
    )
  ) %>%
  rowwise() %>%
  mutate(
    inc_min    = min(c_across(starts_with("inc_set")),na.rm = TRUE),
    inc_max    = max(c_across(starts_with("inc_set")),na.rm = TRUE),
    inc_median = median(c_across(starts_with("inc_set")),na.rm = TRUE),
    inc_mean   = mean(c_across(starts_with("inc_set")),na.rm = TRUE),
    inc_q2.5   = quantile(c_across(starts_with("inc_set")), 0.025,na.rm = TRUE),
    inc_q97.5  = quantile(c_across(starts_with("inc_set")), 0.975,na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(-all_of(set_cols), -ends_with("_bench"))






####

ICER <- CEA_incre%>%select(sensitivity, scenario, inc_best, all_of(paste0("inc_set", 1:1000))) %>%
  rename_with(~ paste0(.x, "_cost"), c("inc_best", paste0("inc_set", 1:1000))) %>%
  left_join(
    QALY_incre %>%
      select(sensitivity, scenario, inc_best, all_of(paste0("inc_set", 1:1000))) %>%
      rename_with(~ paste0(.x, "_qaly"), c("inc_best", paste0("inc_set", 1:1000))),
    by = c("sensitivity", "scenario")
  )%>%
  mutate({
    cost_mat <- as.matrix(pick(all_of(paste0("inc_set", 1:1000, "_cost"))))
    qaly_mat <- as.matrix(pick(all_of(paste0("inc_set", 1:1000, "_qaly"))))
    icer_mat <- cost_mat / qaly_mat
    colnames(icer_mat) <- paste0("icer_set", 1:1000)
    as_tibble(icer_mat)
  }) %>%
  mutate(
    icer_best   = inc_best_cost / inc_best_qaly,
    icer_mat    = as.matrix(pick(all_of(paste0("icer_set", 1:1000)))),
    icer_min    = apply(icer_mat, 1, min,      na.rm = TRUE),
    icer_max    = apply(icer_mat, 1, max,      na.rm = TRUE),
    icer_median = apply(icer_mat, 1, median,   na.rm = TRUE),
    icer_mean   = rowMeans(icer_mat,            na.rm = TRUE),
    icer_q2.5   = apply(icer_mat, 1, quantile, 0.025, na.rm = TRUE),
    icer_q97.5  = apply(icer_mat, 1, quantile, 0.975, na.rm = TRUE)
  ) %>%
  select(-icer_mat, -starts_with("icer_set"),
         -ends_with("_cost"), -ends_with("_qaly"))


#### table ####
rformat_icer <- function(x) {
  if_else(x < 0, "dominant", format(round(x), big.mark = ",", scientific = FALSE))
}

format_currency <- function(val, lo, hi, unit = "A$", scale = 1e6) {
  sprintf("%s%s\n(%s%s - %s%s)",
          unit, format(round(val / scale), big.mark = ","),
          unit, format(round(lo  / scale), big.mark = ","),
          unit, format(round(hi  / scale), big.mark = ","))
}

cost_str <- CEA_incre %>%
  arrange(sensitivity, scenario) %>%
  mutate(cost_str = format_currency(best, q2.5, q97.5, unit = "A$", scale = 1e6)) %>%
  select(sensitivity, scenario, cost_str)

qaly_str <- QALY_incre %>%
  arrange(sensitivity, scenario) %>%
  mutate(qaly_str = sprintf("%s\n(%s - %s)",
                            format(round(best / 1e3), big.mark = ","),
                            format(round(q5   / 1e3), big.mark = ","),
                            format(round(q95  / 1e3), big.mark = ","))) %>%
  select(sensitivity, scenario, qaly_str)

inc_cost_str <- CEA_incre %>%
  arrange(sensitivity, scenario) %>%
  mutate(inc_cost_str = if_else(scenario == bench_scn, "-",
                                sprintf("A$%s\n(A$%s - A$%s)",
                                        format(round(inc_best  / 1e3), big.mark = ","),
                                        format(round(inc_q2.5  / 1e3), big.mark = ","),
                                        format(round(inc_q97.5 / 1e3), big.mark = ",")))) %>%
  select(sensitivity, scenario,inc_cost_str)

inc_qaly_str <- QALY_incre %>%
  arrange(sensitivity, scenario) %>%
  mutate(inc_qaly_str = if_else(scenario == bench_scn, "-",
                                sprintf("%s\n(%s - %s)",
                                        format(round(inc_best),  big.mark = ","),
                                        format(round(inc_q2.5),  big.mark = ","),
                                        format(round(inc_q97.5), big.mark = ",")))) %>%
  select(sensitivity, scenario, inc_qaly_str)

icer_str <- ICER %>%
  arrange(sensitivity, scenario) %>%
  mutate(icer_str = if_else(scenario == bench_scn, "-",
                            sprintf("%s\n(%s - %s)",
                                    rformat_icer(icer_best),
                                    rformat_icer(icer_q2.5),
                                    rformat_icer(icer_q97.5)))) %>%
  select(sensitivity, scenario, icer_str)

cea_gt_df <- bind_cols(cost_str, qalys = qaly_str$qaly_str, inc_cost = inc_cost_str$inc_cost_str, 
                       inc_qaly = inc_qaly_str$inc_qaly_str, ICER = icer_str$icer_str)


cea_gt_df %>%
  gt(groupname_col = "sensitivity") %>%
  tab_header(title = "Cost-Effectiveness Analysis") %>%
  cols_label(
    scenario = "Scenarios",
    cost_str = html("Costs, million\n(discounted)"),
    qalys    = html("QALYs, thousand\n(discounted)"),
    inc_cost = html("Incremental costs, thousand\n(discounted)"),
    inc_qaly = html("QALYs gained\n(discounted)"),
    ICER     = "ICER"
  ) %>%
  cols_align(align = "center", columns = -scenario) %>%
  cols_align(align = "left",   columns = scenario) %>%
  # NO fmt() or fmt_markdown - keep raw \n strings
  tab_footnote(
    footnote = "Values are presented as estimate (95% credible interval). Incremental values relative to No National Program. Dominant = intervention both saves costs and gains QALYs."
  ) %>%
  tab_options(
    table.width                          = pct(100),
    table.font.size                      = px(11),
    heading.title.font.size              = px(14),
    heading.title.font.weight            = "bold",
    heading.align                        = "left",
    column_labels.font.weight            = "bold",
    column_labels.border.top.width       = px(2),
    column_labels.border.top.color       = "black",
    column_labels.border.bottom.width    = px(1),
    column_labels.border.bottom.color    = "black",
    row_group.font.weight                = "bold",
    row_group.border.top.width           = px(1),
    row_group.border.bottom.width        = px(0),
    table_body.border.bottom.width       = px(2),
    table_body.border.bottom.color       = "black",
    table_body.hlines.width              = px(0),
    footnotes.font.size                  = px(10)
  ) %>%
  gtsave(file = file.path(OutputFig, "HCV_ICER_scenarios.docx"))



#### ROI #### 

#### ROI #### 
na0 <- function(x) { x[is.na(x)] <- 0; x }
program_cost <- list()
for(i in names(Rescost_year_all)){ 
  for(m in names(Rescost_year_all[[1]])){ 
    
    program_cost[[i]][[m]]<- 
      cbind(year = Rescost_year_all[[i]][[m]]$cost_ab_sc$year, 
            as.data.frame(na0(Rescost_year_all[[i]][[m ]]$cost_ab_sc[, par_col])+ 
                            na0(Rescost_year_all[[i]][[m ]]$cost_RNA_sc[, par_col]) +
                            na0(Rescost_year_all[[i]][[m ]]$cost_POCT_sc[, par_col]) + 
                            na0(Rescost_year_all[[i]][[m ]]$cost_TreatOther_sc[, par_col])))
    
    
    program_cost[[i]][[m]][program_cost[[i]][[m]] == 0] <- NA
    
  }
  
}

program_cost_disy <- list()

for(i in names(program_cost)){ 
  for(m in names(program_cost[[1]])){
    program_cost_disy[[i]][[m]] <- program_cost[[i]][[m]]%>%
      mutate(id = year - POC_AU$simY, 
             discount = ifelse(id>=0, (1 + AUdiscount)^id, NA))%>%
      mutate(across(c(par_col), ~./discount,
                    .names = "{col}"))
    
  }
}

program_cost_discumy <- list()

for(i in names(program_cost_disy)){ 
  for(m in names(program_cost_disy[[1]])){
    program_cost_discumy[[i]][[m]] <- program_cost_disy[[i]][[m]]%>%
      mutate(across(
        all_of(par_col),                          # select columns by name
        ~cumsum(replace_na(., 0)),                # replace NA then cumsum
        .names = "{.col}"                     # note the dot in {.col}
      ))%>%mutate(across(
        all_of(par_col),
        ~na_if(., 0)                  # turn 0 back to NA
      ))%>%
      
      popResults_range(POC_AU, ., end_Y = 100-1)%>%
      as_tibble()
  }
}

set_cols <- paste0("set", 1:1000)   # or: grep("^set\\d+$", names(df), value = TRUE)

x_catcost_CEA <- lapply(list(cost_disydaacap_categories_bind,cost_disydaanocap_categories_bind), function(x) x%>%
                          filter(year>= 2022)%>%
                          group_by(scenario, sensitivity, Categories)%>%
                          mutate(across(c(par_col, "min", "max", "Med", "Mu", "q5", 
                                          "q25", "q75", "q95"), cumsum, .names = "{col}"))%>%ungroup()%>%
                          arrange(scenario, sensitivity))

names(x_catcost_CEA) <- c("discount_cap", "discount_nocap")
View(x_catcost_CEA$discount_cap)
CEA_cost <- x_catcost_CEA$discount_cap%>%group_by(scenario, sensitivity, year)%>%
  filter(year == 2041)%>%
  mutate(
    across(c("best", all_of(set_cols)), ~ sum(.x))
  ) %>%
  mutate(
    min    = min(c_across(all_of(set_cols)),na.rm = TRUE),
    max    = max(c_across(all_of(set_cols)),na.rm = TRUE),
    median = median(c_across(all_of(set_cols)),na.rm = TRUE),
    mean   = mean(c_across(all_of(set_cols)),na.rm = TRUE),
    q2.5   = quantile(c_across(all_of(set_cols)), 0.025,na.rm = TRUE),
    q97.5  = quantile(c_across(all_of(set_cols)), 0.975,na.rm = TRUE)
  )

CEA_cost <- CEA_cost%>%select(!Categories)%>%slice(1)%>%
  select(sensitivity, scenario, year, best, all_of(set_cols), min, max, median, mean, q2.5, q97.5)

bench_scn <- "(1) No national program"

CEA_incre <- CEA_cost  %>%
  group_by(sensitivity, year) %>%
  left_join(
    CEA_cost  %>%
      filter(scenario == bench_scn) %>%
      select(sensitivity, year, best, all_of(set_cols)) %>%
      rename(best_bench = best) %>%
      rename_with(~ paste0(.x, "_bench"), all_of(set_cols)),
    by = c("sensitivity", "year")
  ) %>%
  mutate(
    inc_best = best - best_bench,
    across(
      all_of(set_cols),
      ~ .x - get(paste0(cur_column(), "_bench")),
      .names = "inc_{col}"
    )
  ) %>%
  rowwise() %>%
  mutate(
    inc_min    = min(c_across(starts_with("inc_set")),na.rm = TRUE),
    inc_max    = max(c_across(starts_with("inc_set")),na.rm = TRUE),
    inc_median = median(c_across(starts_with("inc_set")),na.rm = TRUE),
    inc_mean   = mean(c_across(starts_with("inc_set")),na.rm = TRUE),
    inc_q2.5   = quantile(c_across(starts_with("inc_set")), 0.025,na.rm = TRUE),
    inc_q97.5  = quantile(c_across(starts_with("inc_set")), 0.975,na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(-all_of(set_cols), -ends_with("_bench"))
CEA_incre <- CEA_incre %>%
  rename(scenario = scenario.x) %>%
  filter(scenario != bench_scn)%>%
  select(-scenario.y)

program_cost_2030 <- map_dfr(names(program_cost_discumy), function(i) {
  map_dfr(names(program_cost_discumy[[i]]), function(m) {
    program_cost_discumy[[i]][[m]] %>%
      filter(year == 2030) %>%
      rowwise() %>%
      
      mutate(
        pc_best = best,
        pc_q5   = quantile(c_across(all_of(set_cols)), 0.025, na.rm = TRUE),
        pc_q95  = quantile(c_across(all_of(set_cols)), 0.975, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      select(pc_best,all_of(set_cols),
             pc_q5, pc_q95) %>%
      mutate(sensitivity = i, scenario = m)
  })
})
program_cost_2030_join <- program_cost_2030 %>%
  rename_with(~ paste0("pc_", .x), all_of(set_cols)) %>%
  select(sensitivity, scenario, pc_best, all_of(paste0("pc_", set_cols)))%>%
  arrange(sensitivity, scenario)

pc_set_cols <- paste0("pc_", set_cols)

roi_table <- CEA_incre %>%
  mutate(scenario = as.character(scenario)) %>%   # match type
  filter(scenario != bench_scn) %>%
  select(-c(best,min, max, median, mean, q2.5, q97.5))%>%arrange(sensitivity, scenario)%>%
  left_join(program_cost_2030_join, by = c("sensitivity", "scenario")) %>%
  mutate(roi_best = -inc_best / pc_best) %>%
  mutate(across(
    all_of(paste0("inc_", set_cols)),
    ~ -(.x / get(paste0("pc_", sub("inc_", "", cur_column())))),  # strip inc_ prefix
    .names = "roi_{col}"
  )) %>%
  rowwise() %>%
  rename_with(~ sub("roi_inc_", "roi_", .x), starts_with("roi_inc_"))%>%
  mutate(
    roi_min   = min(c_across(starts_with("roi_set")),            na.rm = TRUE),
    roi_q5    = quantile(c_across(starts_with("roi_set")), 0.025, na.rm = TRUE),
    roi_q95   = quantile(c_across(starts_with("roi_set")), 0.975, na.rm = TRUE),
    roi_lower = ifelse(roi_q5 < roi_best, roi_min, roi_q5)
  ) %>%
  ungroup() %>%
  select(sensitivity, scenario, roi_best, roi_q5, roi_q95, -starts_with("roi_set"))


# Format helpers
format_currency_roi <- function(val, lo, hi, unit = "A$", scale = 1e6) {
  sprintf("%s%s\n(%s%s - %s%s)",
          unit, format(round(val / scale), big.mark = ","),
          unit, format(round(lo  / scale), big.mark = ","),
          unit, format(round(hi  / scale), big.mark = ","))
}

format_roi <- function(val, lo, hi) {
  sprintf("%.2f\n(%.2f - %.2f)", val, lo, hi)
}

# Incremental cost string
inc_cost_str <- CEA_incre %>%
  filter(year == 2041) %>%
  arrange(sensitivity, scenario) %>%
  mutate(inc_cost_str = format_currency_roi(-inc_best, -inc_q97.5,-inc_q2.5,  
                                            unit = "A$", scale = 1e6)) %>%
  select(sensitivity, scenario, inc_cost_str)

# Program cost string
pc_str <- program_cost_2030 %>%
  arrange(sensitivity, scenario) %>%
  mutate(pc_str = format_currency_roi(pc_best, pc_q5, pc_q95, 
                                      unit = "A$", scale = 1e6)) %>%
  select(sensitivity, scenario, pc_str)

# ROI string
roi_str <- roi_table %>%
  arrange(sensitivity, scenario) %>%
  mutate(roi_str = format_roi(roi_best, roi_q5, roi_q95)) %>%
  select(sensitivity, scenario, roi_str)

# Combine
roi_gt_df <- inc_cost_str %>%
  left_join(pc_str,  by = c("sensitivity", "scenario")) %>%
  left_join(roi_str, by = c("sensitivity", "scenario"))
library(gt)
# GT table
roi_gt_df %>%
  gt(groupname_col = "sensitivity") %>%
  tab_header(title = "Return on Investment Analysis") %>%
  cols_label(
    scenario     = "Scenarios",
    inc_cost_str = html("Net savings, million<br>(discounted)"),
    pc_str       = html("Program Cost, million<br>(discounted)"),
    roi_str      = html("ROI<br>(net savings / program cost)")
  ) %>%
  cols_align(align = "center", columns = -scenario) %>%
  cols_align(align = "left",   columns = scenario) %>%
  tab_footnote(
    footnote = "Values are presented as estimate (95% credible interval). ROI = net savings divided by program cost."
  ) %>%
  tab_options(
    table.width                          = pct(100),
    table.font.size                      = px(11),
    heading.title.font.size              = px(14),
    heading.title.font.weight            = "bold",
    heading.align                        = "left",
    column_labels.font.weight            = "bold",
    column_labels.border.top.width       = px(2),
    column_labels.border.top.color       = "black",
    column_labels.border.bottom.width    = px(1),
    column_labels.border.bottom.color    = "black",
    row_group.font.weight                = "bold",
    row_group.border.top.width           = px(1),
    row_group.border.bottom.width        = px(0),
    table_body.border.bottom.width       = px(2),
    table_body.border.bottom.color       = "black",
    table_body.hlines.width              = px(0),
    footnotes.font.size                  = px(10)
  ) %>%
  gtsave(file = file.path(OutputFig%>%dirname(), "HCV_ROI_scenarios.docx"))



