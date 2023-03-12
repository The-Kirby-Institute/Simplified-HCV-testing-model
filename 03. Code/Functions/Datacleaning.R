#### Taiwan CDC data-- living HIV cases 
library("readxl")
CDCData <- read_excel(paste0(projectFolder, "/01. DATA", 
                  "/HIV AIDS data recode_Taiwan CDCdata_living.xlsx", sep = ""), 
           sheet = 1)%>%filter(sex ==1 & `risk factor`%in%c("4","5","6","7"))

# number of HIV+d in each year 
# HIV diagnosed before June of that year and still alive/die after middle of that year. 
# remove female 
# upper estimated: risk factor%in% c(4,5,6,7)
# lower estimated: risk factor%in%c(5,6)
# estimated point: risk factor%in% c(5,6,7)
CDCData$death[is.na(CDCData$death)] <- "2021-12-31"
CDCY <- list()

crieteria <- c(paste(seq(2004, 2015,1), "-06-30", sep=""))

CDCData%>%filter(death>crieteria[1])

for (i in 1: length(crieteria)){ 
  CDCY[[i]] <- CDCData%>%filter(`HIV date` <= crieteria[i] & death>crieteria[i] )
  
}



upper <- lapply(CDCY, function(x){nrow(x)})%>%unlist()
lower <- lapply(CDCY, function(x){nrow(filter(x, `risk factor`%in%c("5","6"))
                                       )})%>%unlist()
point <- lapply(CDCY, function(x){nrow(filter(x, `risk factor`%in%c("5","6", "7")
                                              ))})%>%unlist()

CDCHIVD <-data.frame(year = seq(2004,2015,1), 
                     realPop = point,
                     low = lower,
                     up =upper)%>%tbl_df()

write.csv(CDCHIVD, 
          file.path(basePath,"01. Data/model input/ECDCHIVD.csv")) 



# death per year 


CDCD <- list()

crieteriaD1 <- c(paste(seq(2004, 2016,1), "-01-01", sep=""))




for (i in 1: length(crieteriaD1)){ 
  CDCD[[i]] <- CDCData%>%filter(death>=crieteriaD1[i] & death<crieteriaD1[i+1])
  
}

upperD <- lapply(CDCD, function(x){nrow(x)})%>%unlist()
lowerD <- lapply(CDCD, function(x){nrow(filter(x, `risk factor`%in%c("5","6"))
)})%>%unlist()
pointD <- lapply(CDCD, function(x){nrow(filter(x, `risk factor`%in%c("5","6", "7")
))})%>%unlist()

CDCHIVdeath <-data.frame(year = seq(2004,2016,1), 
                     realPop = pointD,
                     low = lowerD,
                     up =upperD)%>%tbl_df()



CDCHIVdeath <- CDCHIVdeath%>%filter(year!=2016)%>%
  mutate(time =year - (HCV$cabY -1 ) )

write.csv(CDCHIVdeath, 
          file.path(basePath,"01. Data/CDCHIVdeath.csv")) 
