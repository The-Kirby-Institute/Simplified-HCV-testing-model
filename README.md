# Simplified HCV testing model

This project contains the code for investigating the impact of point-of-care testing and other simplified HCV diagnosis strategies for HCV prevention on HCV incidence and prevalence across the different key populations. 


[![DOI](https://zenodo.org/badge/612906740.svg)](https://zenodo.org/badge/latestdoi/612906740)

Model developer, coder and maintainer of this repository: [Joyce Huei-Jiuan Wu](https://github.com/ninowwss)

ORCiD ID: [0000-0002-5085-1976](https://orcid.org/my-orcid?orcid=0000-0002-5085-1976)


Affiliation: The Kirby Institute, UNSW Sydney, NSW, Australia


# Installation
You need the following software & associated packages to run this model:

* R, a free statistical program to run and analyze the model results using the provided scripts.
* (Optional) RStudio, a useful user interface for the R program.
* R packages associated with TWHCV model: 

    `gt_0.8.0`, `expss_0.11.4`, `maditr_0.8.3`, `officer_0.6.1`, `flextable_0.8.6`, `ggthemes_4.2.4`,
 
    `ggpubr_0.6.0`, `xtable_1.8-4`, `gridExtra_2.3`, `cowplot_1.1.2`, `colorblindr_0.1.0`, `colorspace_2.1-1`, 
 
    `reshape2_1.4.4`, `formattable_0.2.1`, `data.table_1.14.8`, `lubridate_1.9.2`, `forcats_1.0.0`, `stringr_1.5.0`,    
 
    `dplyr_1.1.0`, `purrr_1.0.1`, `readr_2.1.4`, `tidyr_1.3.0`, `tibble_3.2.0`, `ggplot2_3.4.1`,  
 
    `tidyverse_2.0.0`  
    
* Clone or download the code from this repository into a convenient location on your computer. 

# Project structure 
```bash

├── 01.DATA/            # model input data & data for calibraiton  
│   ├── model input/           # model input data
│       ├── cost
│       ├── parameter_varied_stages
│       ├── Scenarios           
│       └── Used 
├── 02.Documents/	# model structure diagram 
├── 03.Code/
│    └── Functions             # functions for data wrangling, model simulation and results generation(plots) 
├── 04.Output/		#  models outputs
│    ├── Result_data   
│    ├── Result_figure
│    └── Results               # the associated simulation results are stored as .rda files in this subfolder
├── 00_CreatProject.R
├── 01.SetupModel.R
├── 02_CalibrateHCVmodel.R
├── 02_2_Uncertainty.R
├── 02_3_scenarioSetup.R
├── 03_RunningHCVModel.R
├── 04. SummaryResults_update.R
├── 04_1_SummaryScenarios.R
├── 04_2_Impact estimation.R
```
Note: Files in `01.DATA`,`02.Documents` and `04.Output` folder are stored locally and are not available in this online repository.

# Set up and run the TWHCV model
All the project files are stored in the main directory and 4 main sub-directories. This main README file describes the overall project. 

The project code is written in `R` version 4.2.2 as R or R markdown scripts with `Rstudio` version 2022.12.0.353. 


All model inputs stored as either `.csv` or `.xlsx` and outputs stored as `.rda` files.



`00_CreatProject.R`: This script is used to setup a project based on the specifications entered by the user. 
The user has specify a project name, a start and end year for simulations, the simulation time step, and population names (which determines the number of populations in the model). 
The script creates a project directory in the "projects" directory and creates a project specifications `.rda` file and relevant sub-directories. 
Template input `.csv` files are also created for the user to input parameter values:
 * `initial_populations.csv`: For specifying initial population sizes.
 * `parameters_constants.csv`: For specifying parameters which constant over time.
 * `diseaseProgress.csv`: For specifying parameters regarding to disease progression which constant over time.
 * `curedProgress.csv`: For specifying parameters regarding to disease progression after cured from HCV which constant over time.
 * parameter_variedstage_set (sub-folder): For specifying the base value of parameters regarding to HCV testing and treatment which change over time.
    7 csv files included in this sub-folder. 
      *    `tau_ab.csv`
      *    `tau_ag.csv`
      *    `tau_poct.csv`
      *    `eta.csv`
      *    `rho.csv`
      *    `lota.csv`
      *    `cured.csv`
 * `population_transitions.csv`: Specifies the rate populations moves from one to another (e.g. acquiring HIV). 
    If population X moves into population Y at rate r then r is entered in the cell corresponding to column named XY.


`01.SetupModel.R`: This script is used to read in all the specified projects input `.csv` files and create the parameter sets that will be run in the model. 
  This script needs to be run for the actual model simulations to be performed. 
  The project parameter sets consist of a "best estimates" parameter set and an ensemble of parameter sets generated through sampling from the input parameter ranges. 
  All parameter sets are appended to the project's .rda file. 
  If the user changes or updates the .csv input files then this script need to be rerun for the correct parmaters to be used.
  
Then opend `FitInitCondtion.R` (in sub-folder 03.Code/Functions) for finding the intital condition of the model. 
  This script runs the model simulation, reallocates the number of population in each components and produce plots of the model output for visual inspection. 
 you need to rerun `01.SetupModel.R` to save the final parameters in the project.

Once you finished the initiation condtion fitting, run `02_CalibrateHCVmodel.R` to save the final parameters in the project. 

Then run `calibrate_timevarying.R`. This R file contains a function which is useful for calibrating the model. 
A user can change the value of time-varying parameters, run simulation and save to a `.rda` file. 


Then open `04. SummaryResults_update.R` to produce plots of the model output for visual inspection. 
If visual inspection is not satisfied, go back to `calibrate_timevarying.R` changed value of parameters then rerun model simulation and 
open `04. SummaryResults_update.R` to produce plots of the model output for visual inspection. A user may repeat these steps for several times to finish calibration.

Once a user finish calibration, can further produce parameter sets of uncertainty by using `02_2_Uncertainty.R`.  

`02_3_scenarioSetup.R`: This script is for producing the parameter sets for each scenario. 
 

`03_RunningHCVModel.R`: This script is for HCV model simulation and save the results in a `.rda` file. 

`04_1_SummaryScenarios.R`: This script is used to produce figures of summary main results for a specified project.

`04_2_Impact estimation.R`: This script is used to produce figures and table of results regarding to main scenarios and senstivity analysis for a specified project.

Note as these are scripts, not functions, care should be taken to ensure the project specifications are correct before running a script especially after updates have been pulled from the repository.






# Publication 
The following publication is associated with this project and used the code in this repository to generate all of the results and figures.

* 



# Disclaimer
This study is made possible as part of a research-funded PhD being undertaken by [Joyce Huei-Jiuan Wu](https://github.com/ninowwss) under the University of New South Wales (UNSW) Scientia scholarship and associated with the Rapid Point of Care Research Consortium for infectious disease in the Asia Pacific (RAPID), which is funded by an NHMRC Centre for Research Excellence.
The Kirby Institute is funded by the Australian Government Department of Health and is affiliated with the Faculty of Medicine, UNSW Sydney, Australia. 
The views expressed in this publication do not necessarily represent the position of the Australian Government.

