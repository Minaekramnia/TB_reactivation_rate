#Logistic regression for the risk factors for two cycles of NHANES dataset

rm(list=ls())

getwd()
setwd("/Users/mie368/Documents/Harvard/Data/Riskfactors/NHANES/2011-12/")
if (!require(survey)) {install.packages("survey");require(survey)}
if (!require(foreign)){install.packages("foreign");require(foreign)}
if (!require(dplyr))   {install.packages("dplyr");require(dplyr)}
if (!require(tidyr)) {install.packages("tidyr");require(tidyr)}

if (!require(mgcv)) {install.packages("mgcv");require(mgcv)}
if (!require(MASS)) {install.packages("MASS");require(MASS)}

if (!require(lme4)) {install.packages("lme4");require(lme4)}
if (!require(merTools)) {install.packages("merTools");require(merTools)}

if (!require(brms)) {install.packages("brms");require(brms)}
library(rlang)
if (!require(tidybayes)) {install.packages("tidybayes");require(tidybayes)}
if (!require(blme)) {install.packages("blme");require(blme)}
library("plyr")  #for append
#sasexport
if(!require('janitor')) {
  install.packages('janitor')
  library('janitor')
}
library("SASxport")
library("Hmisc")
library(dplyr)
#Packages Needed for Cross-Tab
req <- substitute(require(x, character.only = TRUE))
libs<-c("sjPlot")
sapply(libs, function(x) eval(req) || {install.packages(x); eval(req)})
##

rawdata.path<-"/Users/mie368/Documents/Harvard/Data/Riskfactors/NHANES/2011-12/"
options(survey.lonely.psu = "adjust")

DownloadImport<-function(ftp.filepath){
  tf<-tempfile()
  download.file(ftp.filepath,tf,mode="wb")
  read.xport(tf) #By putting read.xport as the final line, this function will return on r data frame
}

#MINA code for risk factor analysis NHANES
require(SASxport)
## Survey Analysis
# Define a function to call svymean and unweighted count
getSummary <- function(varformula, byformula, design){
  # Get mean, stderr, and unweighted sample size
  c <- svyby(varformula, byformula, design, unwtd.count ) 
  p <- svyby(varformula, byformula, design, FUN = svyciprop, deff=FALSE, vartype=c("se","ci"), method = "logit")
  outSum <- left_join(select(c,-se), p) 
  outSum
}

###############  NHANES 2011-2012  ###########################
## DATA IN
tstdat11 <- sasxport.get("/Users/mie368/Documents/Harvard/Data/NHANES_2011-12_data/TBX_G_2011.XPT")
igradat11 <- sasxport.get("/Users/mie368/Documents/Harvard/Data/NHANES_2011-12_data/TB_G_2011.XPT")
demogdat11 <- read.xport("/Users/mie368/Documents/Harvard/Data/NHANES_2011-12_data/DEMO_G_2011.XPT")

sub.demo1112<-demogdat11[,c("SEQN","SDDSRVYR","RIAGENDR","RIDAGEYR","RIDRETH1","RIDRETH3","DMDBORN4","WTINT2YR","WTMEC2YR","SDMVPSU","SDMVSTRA")]
sub.demo1112[,"racenha"]<-NA
sub.demo1112[which(sub.demo1112$RIDRETH3==1|sub.demo1112$RIDRETH3==2),"racenha"]<-"Hispanic"
sub.demo1112[which(sub.demo1112$RIDRETH3==3),"racenha"]<-"NHW"
sub.demo1112[which(sub.demo1112$RIDRETH3==4),"racenha"]<-"NHB"
sub.demo1112[which(sub.demo1112$RIDRETH3==6),"racenha"]<-"NHA"
sub.demo1112[which(sub.demo1112$RIDRETH3==7),"racenha"]<-"others"
sub.demo1112[,"frbn"]<-NA
sub.demo1112[which(sub.demo1112$DMDBORN4==1),"frbn"]<-0
sub.demo1112[which(sub.demo1112$DMDBORN4==2),"frbn"]<-1
#sub.demo1112<-  sub.demo1112[,-c(which(colnames(sub.demo1112)=="RIDRETH3"),which(colnames(sub.demo1112)=="DMDBORN4"))]

#merge IGRA, TST and Demographic datasets for 2011
DatNH11 <- merge(x=sub.demo1112,y=tstdat11,by.x="SEQN", by.y='seqn',all=F)
DatNH11 <- merge(x=DatNH11,y=igradat11,by.x="SEQN", by.y="seqn",all=F)

# datdesign<-
#   svydesign(
#     id=~SDMVPSU,
#     strata=~SDMVSTRA,
#     nest=T,
#     weights=~WTMEC2YR,
#     data=DatNH11
#   )
# 
# datdesign <- subset(datdesign, RIDAGEYR>=1 & RIDAGEYR<=80 & RIAGENDR >=1 & RIAGENDR <=2 & DMDBORN4 >=1 & DMDBORN4 <=2 )

#Age groups 6-80
#building our age/sex/US born categories
#newdata11 <- DatNH11 %>% filter(RIDAGEYR>=1 & RIDAGEYR<=80 & RIAGENDR >=1 & RIAGENDR <=2 & DMDBORN4 >=1 & DMDBORN4 <=2 )
#we drop 5 people because of the DBDBORN4
newdata11 <- DatNH11
newdata11$count <- 1
newdata11$AGECAT <- 0
newdata11$AGECAT[newdata11$RIDAGEYR>=1 & newdata11$RIDAGEYR<5] <-  "age_01_04"
newdata11$AGECAT[newdata11$RIDAGEYR>=5 & newdata11$RIDAGEYR<15] <- "age_5_14"
newdata11$AGECAT[newdata11$RIDAGEYR>=15 & newdata11$RIDAGEYR<25] <- "age_15_24"
newdata11$AGECAT[newdata11$RIDAGEYR>=25 & newdata11$RIDAGEYR<35] <- "age_25_34"
newdata11$AGECAT[newdata11$RIDAGEYR>=35 & newdata11$RIDAGEYR<45] <- "age_35_44"
newdata11$AGECAT[newdata11$RIDAGEYR>=45 & newdata11$RIDAGEYR<55] <- "age_45_54"
newdata11$AGECAT[newdata11$RIDAGEYR>=55 & newdata11$RIDAGEYR<65] <- "age_55_64"
newdata11$AGECAT[newdata11$RIDAGEYR>=65 & newdata11$RIDAGEYR<75] <- "age_65_74"
newdata11$AGECAT[newdata11$RIDAGEYR>=75 & newdata11$RIDAGEYR<=80] <- "age_75_80"

#turning tbdruind in to binary 
newdata11$tbdruind_bin[newdata11$tbdruind>10 ] <- 1
newdata11$tbdruind_bin[newdata11$tbdruind<=10 ] <- 0

#for 2011, replacing 3 (indeterminate result) to NA. 2 to 0 
newdata11$lbxtbin<-replace(newdata11$lbxtbin, newdata11$lbxtbin==3, 0)
newdata11$lbxtbin<-replace(newdata11$lbxtbin, newdata11$lbxtbin==2, 0)

#check if there is duplicates
length(which(table(newdata11$SEQN)>1))

#INPUT RISK FACTORS
####Diabete : DIQ010
dat_diab_11<-read.xport("/Users/mie368/Documents/Harvard/Data/Riskfactors/NHANES/2011-12/DIQ_G.XPT")

#Replacing ALL no, borderline, refused, don't know diabs to 0
dat_diab_11$DIQ010<-replace(dat_diab_11$DIQ010, dat_diab_11$DIQ010==2, 0)
dat_diab_11$DIQ010<-replace(dat_diab_11$DIQ010, dat_diab_11$DIQ010==3, 0)
dat_diab_11$DIQ010<-replace(dat_diab_11$DIQ010, dat_diab_11$DIQ010==7, 0)
dat_diab_11$DIQ010<-replace(dat_diab_11$DIQ010, dat_diab_11$DIQ010==9, 0)
#dat_diab_11$DIQ010<-replace(dat_diab_11$DIQ010, dat_diab_11$DIQ010==., 0)
#DIQ175C - Age

#HIV antibody test result: LBDHI
hiv_1112_loc<-"https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/HIV_G.XPT"
hiv_1112<-DownloadImport(hiv_1112_loc)

hiv_1112$LBDHI<-replace(hiv_1112$LBDHI, hiv_1112$LBDHI==2, 0)

#TBQ040 - Ever told you had active TB
#TBQ060: Ever lived in the same household with someone while that person was sick with tuberculosis or TB?
TBQ_1112_loc<-"https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/TBQ_G.XPT"
TBQ_1112<-DownloadImport(TBQ_1112_loc)

TBQ_1112$TBQ040<-replace(TBQ_1112$TBQ040, TBQ_1112$TBQ040==2, 0)
TBQ_1112$TBQ040<-replace(TBQ_1112$TBQ040, TBQ_1112$TBQ040==7, 0)
TBQ_1112$TBQ040<-replace(TBQ_1112$TBQ040, TBQ_1112$TBQ040==9, 0)

TBQ_1112$TBQ060<-replace(TBQ_1112$TBQ060, TBQ_1112$TBQ060==2, 0)
TBQ_1112$TBQ060<-replace(TBQ_1112$TBQ060, TBQ_1112$TBQ060==7, 0)
TBQ_1112$TBQ060<-replace(TBQ_1112$TBQ060, TBQ_1112$TBQ060==9, 0)

#Merge IGRA, TST and Demographic with Diab datasets for 2011
#newdata11 <- merge(x=newdata11,y=dat_diab ,by.x="SEQN",by.y = "SEQN",all=T)
newdata11_diab <- merge(x=newdata11,y=dat_diab_11, by="SEQN", all.x=T) #keep all rows of x dataframe

#Merge IGRA, TST and Demographic with HIV datasets for 2011
newdata11_hiv <- merge(x=newdata11_diab, y= hiv_1112, all.x=T)

#Merge IGRA, TST and Demographic with TBQ datasets for 2011
newdata11_tbq <- merge(x=newdata11_hiv, y= TBQ_1112, all.x=T)

#Merge IGRA, TST and Demographic with RENAL datasets for 2011 - ESRD
#datdemoncr<-get(load(file='/Users/mie368/Documents/Harvard/Data/Riskfactors/NHANES/2011-12/datcr_renal.rda'))
#updated results of ESRD from Yunfei. cleaned/merged data (datdemoncr_ver2) that has SDDSRVYR and esrd every year. 
datdemoncr_v2<-get(load("/Users/mie368/Documents/Harvard/reactivation_rates/datdemoncr_v2.rda"))

newdata11_renal_v2 <- merge(x=newdata11_tbq,y= datdemoncr_v2,all.x=T)

#Immunosuppression
#Prescription Medications (RXQ_RX_G) for Prednisone usage
RXQ_RX_1112_loc<-"https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/RXQ_RX_G.xpt"
RXQ_RX_1112<-DownloadImport(RXQ_RX_1112_loc)

RXQ_RX_1112$immunsup<-0
RXQ_RX_1112$immunsup[RXQ_RX_1112$RXDDRUG=="PREDNISONE"]<-1
RXQ_RX_1112_immunesup <- RXQ_RX_1112 %>% filter(RXQ_RX_1112$immunsup==1)

newdata11_renal_v2$PRED<-0
newdata11_renal_v2$PRED[newdata11_renal_v2$SEQN%in%RXQ_RX_1112_immunesup$SEQN] <-1

dat.premedication<-get(load("/Users/mie368/Documents/Harvard/reactivation_rates/dat.premedication.rda"))
newdata11_immuntnf <- merge(x=newdata11_renal_v2,y= dat.premedication,all.x=T)

#Cross tabulation Diabetes, IGRA
# #Cross tabulation HIV, IGRA
# sjPlot::tab_xtab(var.row = newdata11$tbdruind_bin, var.col = newdata11$LBDHI, title = "Table Title", show.row.prc = TRUE)
# sjPlot::tab_xtab(var.row = newdata11$lbxtbin , var.col = newdata11$LBDHI, title = "Table Title", show.row.prc = TRUE)
# 
# tabyl(newdata11,tbdruind_bin,LBDHI)
# tabyl(newdata11,lbxtbin,LBDHI)
# #Cross tabulation TBQ: Ever told you had active TB, IGRA
# sjPlot::tab_xtab(var.row = newdata11$tbdruind_bin, var.col = newdata11$TBQ040, title = "Table Title", show.row.prc = TRUE)
# sjPlot::tab_xtab(var.row = newdata11$lbxtbin , var.col = newdata11$TBQ040, title = "Table Title", show.row.prc = TRUE)
# 
# tabyl(newdata11,tbdruind_bin,TBQ040)
# tabyl(newdata11,lbxtbin,TBQ040)
# 
# #Cross tabulation TBQ: Household TB , IGRA
# tabyl(newdata11,tbdruind_bin,TBQ060)
# tabyl(newdata11,lbxtbin,TBQ060)
# 
# #Cross tabulation Renal: Household TB , IGRA
# sjPlot::tab_xtab(var.row = newdata11$tbdruind_bin, var.col = newdata11$esrd[newdata11$SDDSRVYR.x==7], title = "Table Title", show.row.prc = TRUE)
# sjPlot::tab_xtab(var.row = newdata11$lbxtbin , var.col = newdata11$esrd[newdata11$SDDSRVYR.x==7], title = "Table Title", show.row.prc = TRUE)
# 
# tabyl(newdata11,tbdruind_bin,TBQ060)
# tabyl(newdata11,lbxtbin,TBQ060)
# #Cross tabulation TBQ: Household TB , IGRA
# sjPlot::tab_xtab(var.row = newdata11$tbdruind_bin, var.col = newdata11$immunsup, title = "Table Title", show.row.prc = TRUE)
# sjPlot::tab_xtab(var.row = newdata11$lbxtbin , var.col = newdata11$immunsup, title = "Table Title", show.row.prc = TRUE)
# tabyl(newdata11,tbdruind_bin,immunsup)
# tabyl(newdata11,lbxtbin,immunsup)

###############  1999-2000  ##################################
tstdat99 <- sasxport.get("/Users/mie368/Documents/Harvard/Data/NHANES_1999_data/TB_nhanes99-20_lab_exp.XPT")
demogdat99 <- read.xport("/Users/mie368/Documents/Harvard/Data/NHANES_1999_data/DEMO_99_20.XPT")

#merge datasets for 1999
DatNH99 <- merge(x=demogdat99,y=tstdat99,by.x="SEQN", by.y='seqn',all=F)

#Race/Ethnicity
# sub.demo9900<-demogdat99[,c("SEQN","SDDSRVYR","RIAGENDR","RIDAGEYR","RIDRETH1","DMDBORN","WTINT2YR","WTMEC2YR","SDMVPSU","SDMVSTRA")]
# sub.demo9900[,"racenha"]<-"NOTM"
# sub.demo9900[,"frbn"]<-NA
# sub.demo9900[which(sub.demo9900$DMDBORN==1),"frbn"]<-0
# sub.demo9900[which(sub.demo9900$DMDBORN==2|sub.demo9900$DMDBORN==3),"frbn"]<-1
# sub.demo9900<-sub.demo9900[,-which(colnames(sub.demo9900)=="DMDBORN")]

#Age groups 1-85
#building our age/sex/US born categoreis
#newdata99 <- DatNH99 %>% filter(RIDAGEYR>=1 & RIDAGEYR<=85 & RIAGENDR >=1 & RIAGENDR <=2 & DMDBORN >=1 & DMDBORN <=3 )
newdata99 <- DatNH99
newdata99$count <- 1
newdata99$AGECAT <- 0
newdata99$AGECAT[newdata99$RIDAGEYR>=01 & newdata99$RIDAGEYR<05] <- "age_01_04"
newdata99$AGECAT[newdata99$RIDAGEYR>=05 & newdata99$RIDAGEYR<15] <- "age_05_14"
newdata99$AGECAT[newdata99$RIDAGEYR>=15 & newdata99$RIDAGEYR<25] <- "age_15_24"
newdata99$AGECAT[newdata99$RIDAGEYR>=25 & newdata99$RIDAGEYR<35] <- "age_25_34"
newdata99$AGECAT[newdata99$RIDAGEYR>=35 & newdata99$RIDAGEYR<45] <- "age_35_44"
newdata99$AGECAT[newdata99$RIDAGEYR>=45 & newdata99$RIDAGEYR<55] <- "age_45_54"
newdata99$AGECAT[newdata99$RIDAGEYR>=55 & newdata99$RIDAGEYR<65] <- "age_55_64"
newdata99$AGECAT[newdata99$RIDAGEYR>=65 & newdata99$RIDAGEYR<75] <- "age_65_74"
newdata99$AGECAT[newdata99$RIDAGEYR>=75 & newdata99$RIDAGEYR<=85] <- "age_75_85"

#for 99 place of born, replacing 3 to 2 to have one category for all non US born
newdata99$DMDBORN<-replace(newdata99$DMDBORN, newdata99$DMDBORN==3, 2)

#TST proportion positive - 1999
newdata99$tbdppds_bin[newdata99$tbdppds>10 ] <- 1
newdata99$tbdppds_bin[newdata99$tbdppds<=10 ] <- 0

#INPUT RISK FACTORS
#look for duplicates in SEQN
length(which(table(newdata99$SEQN)>1))

####Diabete : DIQ010
dat_diab_9900_loc<-"https://wwwn.cdc.gov/Nchs/Nhanes/1999-2000/DIQ.XPT"
dat_diab_9900<-DownloadImport(dat_diab_9900_loc)

#Replacing ALL no, borderline, refused, don't know diabs to 0
dat_diab_9900$DIQ010<-replace(dat_diab_9900$DIQ010, dat_diab_9900$DIQ010==2, 0)
dat_diab_9900$DIQ010<-replace(dat_diab_9900$DIQ010, dat_diab_9900$DIQ010==3, 0)
dat_diab_9900$DIQ010<-replace(dat_diab_9900$DIQ010, dat_diab_9900$DIQ010==7, 0)
dat_diab_9900$DIQ010<-replace(dat_diab_9900$DIQ010, dat_diab_9900$DIQ010==9, 0)
#DIQ175C - Age

#TBQ040 - Ever told you had active TB
#TBQ060 - Lived in household TB sick person
TB_9900_loc<-"https://wwwn.cdc.gov/Nchs/Nhanes/1999-2000/TBQ.XPT"
TB_9900<-DownloadImport(TB_9900_loc)

#Replacing ALL no, borderline, refused, don't know diabs to 0
TB_9900$TBQ040<-replace(TB_9900$TBQ040, TB_9900$TBQ040==2, 0)
TB_9900$TBQ040<-replace(TB_9900$TBQ040, TB_9900$TBQ040==7, 0)
TB_9900$TBQ040<-replace(TB_9900$TBQ040, TB_9900$TBQ040==9, 0)

TB_9900$TBQ060<-replace(TB_9900$TBQ060, TB_9900$TBQ060==2, 0)
TB_9900$TBQ060<-replace(TB_9900$TBQ060, TB_9900$TBQ060==7, 0)
TB_9900$TBQ060<-replace(TB_9900$TBQ060, TB_9900$TBQ060==9, 0)

#HIV antibody test result: LBDHI
hiv_9900_loc<-"https://wwwn.cdc.gov/Nchs/Nhanes/1999-2000/LAB03.XPT"
hiv_9900<-DownloadImport(hiv_9900_loc)

hiv_9900$LBDHI<-replace(hiv_9900$LBDHI, hiv_9900$LBDHI==2, 0)
hiv_9900$LBDHI<-replace(hiv_9900$LBDHI, hiv_9900$LBDHI==3, 0)

#Merge TST and Demographic with Diab datasets for 1999-2000
newdata99_diab <- merge(x=newdata99,y=dat_diab_9900, all.x=T)

#Merge TST and Demographic with TB questioner for 1999-2000
newdata99_hiv <- merge(x=newdata99_diab,y=hiv_9900, all.x=T)

#Merge TST and Demographic with TB questionair for 1999-2000
newdata99_tbq <- merge(x=newdata99_hiv,y=TB_9900, all.x=T)

#datdemoncr99<-get(load(file='/Users/mie368/Documents/Harvard/Data/Riskfactors/NHANES/1999-2000/datcr_renal.rda'))
newdata99_renal_v2 <- merge(x=newdata99_tbq,y= datdemoncr_v2, all.x=T)

#Immunosuppresion
RXQ_RX_9920_loc<-"https://wwwn.cdc.gov/Nchs/Nhanes/1999-2000/RXQ_RX.xpt"
RXQ_RX_9920<-DownloadImport(RXQ_RX_9920_loc)

RXQ_RX_9920$immunsup<-0
RXQ_RX_9920$immunsup[RXQ_RX_9920$RXD240B=="PREDNISONE"]<-1
RXQ_RX_9920_immunesup <- RXQ_RX_9920 %>% filter(RXQ_RX_9920$immunsup==1)
newdata99_immunsup <- merge(x= newdata99_renal_v2, y = RXQ_RX_9920_immunesup, all.x=T) #left join ===> edit to create similar to 11

newdata99_immunsup$PRED<-0
newdata99_immunsup$PRED[newdata99_immunsup$SEQN%in%RXQ_RX_9920_immunesup$SEQN] <-1

newdata99_immuntnf<-merge(x=newdata99_immunsup,y= dat.premedication,all.x=T)

# table(dat.premedication$tnf,dat.premedication$SDDSRVYR==1)
# 
# FALSE  TRUE
# 0 91177  9942
# 1    88     0
#For 1999 we dont have any TNF. only 5 for 2011

#newdata99_immune<-merge(x= newdata99, y = RXQ_RX_9920_immunesup, by = intersect(names(newdata99), names(RXQ_RX_9920_immunesup)), all = FALSE, all.x = T, all.y = T)

# #TNF Alpha
# dat.demodrugn4_99 <- dat.demodrugn4 %>% filter(dat.demodrugn4$SDDSRVYR==1)
# newdata99<- merge(x=newdata99,y=dat.demodrugn4_99,all.x=T)

####Subset dataframes 2011-1999
newdata11_sub<-newdata11_immuntnf[, c('SEQN','SDMVPSU','SDDSRVYR','WTMEC2YR','SDMVSTRA','AGECAT','count','racenha','RIDAGEYR','RIAGENDR','DMDBORN4','LBDHI','DIQ010','TBQ060', 'TBQ040','esrd','tnf','immu','PRED','tbdruind_bin','lbxtbin')]
newdata99_sub<-newdata99_immuntnf[, c('SEQN','SDMVPSU','SDDSRVYR','WTMEC2YR','SDMVSTRA','AGECAT','count','RIDAGEYR','RIAGENDR','DMDBORN','LBDHI','DIQ010','TBQ060', 'TBQ040','esrd','tnf','immu','PRED','tbdppds_bin')]

#creating similar columns in both dataframes
newdata99_sub$lbxtbin<-NA #IGRA
colnames(newdata11_sub)[which(colnames(newdata11_sub) == "DMDBORN4")] = "DMDBORN"
colnames(newdata99_sub)[which(colnames(newdata99_sub) == "tbdppds_bin")] = "tbdruind_bin"

newdata99_sub$racenha<-NA 
#see difference in colnames
setdiff(colnames(newdata11_sub), colnames(newdata99_sub))
setdiff(colnames(newdata99_sub), colnames(newdata11_sub))

#the same column names to stack 
#stack the 1 and 2 df, if they have same column names. 
newdata99_11<-rbind(newdata11_sub,newdata99_sub[,colnames(newdata11_sub)])###orders

#renaming colnames to more meaningful names
names(newdata99_11)[names(newdata99_11) == "RIAGENDR"] <- "Sex"
names(newdata99_11)[names(newdata99_11) == "DMDBORN"] <- "USNUSB"
names(newdata99_11)[names(newdata99_11) == "LBDHI"] <- "HIV"
names(newdata99_11)[names(newdata99_11) == "DIQ010"] <- "Diab"
names(newdata99_11)[names(newdata99_11) == "TBQ060"] <- "HHTB"
names(newdata99_11)[names(newdata99_11) == "TBQ040"] <- "EverActiveTB"
#names(newdata99_11)[names(newdata99_11) == "tbdruind_bin"] <- "TST"

datdesign<-
  svydesign(
    id=~SDMVPSU,
    strata=~SDMVSTRA,
    nest=T,
    weights=~WTMEC2YR,
    data=newdata99_11
  )
#subset -> omitting 99-77
datdesign <- subset(datdesign, RIDAGEYR>=1 & RIDAGEYR<=80 & Sex >=1 & Sex <=2 & USNUSB >=1 & USNUSB <=2 )

#renaming colnames to more meaningful names for 2011
names(newdata11_sub)[names(newdata11_sub) == "RIAGENDR"] <- "Sex"
names(newdata11_sub)[names(newdata11_sub) == "DMDBORN"] <- "USNUSB"
names(newdata11_sub)[names(newdata11_sub) == "LBDHI"] <- "HIV"
names(newdata11_sub)[names(newdata11_sub) == "DIQ010"] <- "Diab"
names(newdata11_sub)[names(newdata11_sub) == "TBQ060"] <- "HHTB"
names(newdata11_sub)[names(newdata11_sub) == "TBQ040"] <- "EverActiveTB"

#For 2011 only. 
datdesign_11<-
  svydesign(
    id=~SDMVPSU,
    strata=~SDMVSTRA,
    nest=T,
    weights=~WTMEC2YR,
    data=newdata11_sub
  )
#subset
datdesign_11<-subset(datdesign_11, RIDAGEYR>=1 & RIDAGEYR<=80 & Sex >=1 & Sex <=2 & USNUSB >=1 & USNUSB <=2 )

#the number of people that have both Igra+ and risk factor
x_diab_igrapos<-subset(datdesign_11, Diab==1 & lbxtbin==1)
svytotal(x = ~ Diab, design = x_diab_igrapos)
svytable(~Diab+lbxtbin, design = datdesign_11)

x_HHTB_igrapos<-subset(datdesign_11,HHTB==1 & lbxtbin==1)
svytotal(x = ~ HHTB, design = x_HHTB_igrapos)
svytable(~HHTB+lbxtbin, design = datdesign_11)

x_everTB_igrapos<-subset(datdesign_11,EverActiveTB==1 & lbxtbin==1)
svytotal(x = ~ EverActiveTB, design = x_everTB_igrapos)
svytable(~EverActiveTB+lbxtbin, design = datdesign_11)

x_tnf_igrapos<-subset(datdesign_11,tnf==1 & lbxtbin==1)
svytotal(x = ~ tnf, design = x_tnf_igrapos)
svytable(~tnf+lbxtbin, design = datdesign_11)

#renaming colnames to more meaningful names for 1999
names(newdata99_sub)[names(newdata99_sub) == "RIAGENDR"] <- "Sex"
names(newdata99_sub)[names(newdata99_sub) == "DMDBORN"] <- "USNUSB"
names(newdata99_sub)[names(newdata99_sub) == "LBDHI"] <- "HIV"
names(newdata99_sub)[names(newdata99_sub) == "DIQ010"] <- "Diab"
names(newdata99_sub)[names(newdata99_sub) == "TBQ060"] <- "HHTB"
names(newdata99_sub)[names(newdata99_sub) == "TBQ040"] <- "EverActiveTB"

#For 1999 only. 
datdesign_99<-
  svydesign(
    id=~SDMVPSU,
    strata=~SDMVSTRA,
    nest=T,
    weights=~WTMEC2YR,
    data=newdata99_sub
  )
#subset
datdesign_99<-subset(datdesign_99, RIDAGEYR>=1 & RIDAGEYR<=80 & Sex >=1 & Sex <=2 & USNUSB >=1 & USNUSB <=2 )

#Logistic Regression
#race
model1_race<-svyglm(tbdruind_bin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(racenha)+factor(HIV),family=quasibinomial,design=datdesign,na.action=na.omit)
mod1_race<-summary(model1_race)
round(mod1_race$coefficients,3)

out.coef.mod1r<-orfun(mod=model1_11_race)
z9<-round(out.coef.mod1r,3)
write.table(z9,file="zexport9.csv",sep=",")
#
model1_11_race<-svyglm(lbxtbin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(racenha)+factor(EverActiveTB),family=quasibinomial,design=datdesign_11,na.action=na.omit)
mod1_race<-summary(model1_11_race)
round(mod1_race$coefficients,3)

out.coef.mod1r<-orfun(mod=model1_11_race)
z9<-round(out.coef.mod1r,3)
write.table(z9,file="zexport9.csv",sep=",")

model1_11_race<-svyglm(lbxtbin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(racenha):Diab,family=quasibinomial,design=datdesign_11,na.action=na.omit)
mod1_race<-summary(model1_11_race)
round(mod1_race$coefficients,3)

out.coef.mod1r<-orfun(mod=model1_11_race)
z9<-round(out.coef.mod1r,3)
write.table(z9,file="zexport9.csv",sep=",")

#TST
model1<-svyglm(tbdruind_bin~factor(SDDSRVYR)+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(HIV)+factor(Diab)+factor(HHTB)+factor(EverActiveTB)+factor(esrd)+factor(immu)+factor(PRED)+factor(tnf),family=quasibinomial,design=datdesign,na.action=na.omit)
mod1<-summary(model1)
round(mod1$coefficients,3)

model1_11<-svyglm(tbdruind_bin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(HIV)+factor(Diab)+factor(HHTB)+factor(EverActiveTB)+factor(esrd)+factor(immu)+factor(PRED)+factor(tnf),family=quasibinomial,design=datdesign_11,na.action=na.omit)
mod1_11<-summary(model1_11)
round(mod1_11$coefficients,3)

model1_11<-svyglm(tbdruind_bin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(HIV)+factor(Diab)+factor(HHTB)+factor(EverActiveTB)+factor(esrd)+factor(immu)+factor(PRED)+factor(tnf),family=quasibinomial,design=datdesign_11,na.action=na.omit)
mod1_11<-summary(model1_11)
round(mod1_11$coefficients,3)

#datdesign_11
#no tnf
model1_99<-svyglm(tbdruind_bin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(HIV)+factor(Diab)+factor(HHTB)+factor(EverActiveTB)+factor(esrd)+factor(immu)+factor(PRED)+tnf,family=quasibinomial,design=datdesign_99,na.action=na.omit)
mod1_99<-summary(model1_99)
round(mod1_99$coefficients,3)

model2<-svyglm(tbdruind_bin~factor(SDDSRVYR)+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(HIV),family=quasibinomial,design=datdesign,na.action=na.omit)
mod2<-summary(model2)
round(mod2$coefficients,3)

model2<-svyglm(tbdruind_bin~factor(SDDSRVYR)+factor(AGECAT)+factor(Sex)+racenha+factor(USNUSB)+factor(HIV),family=quasibinomial,design=datdesign_11,na.action=na.omit)
mod2<-summary(model2)
round(mod2$coefficients,3)

model2_11<-svyglm(tbdruind_bin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(HIV),family=quasibinomial,design=datdesign_11,na.action=na.omit)
mod2_11<-summary(model2_11)
round(mod2_11$coefficients,3)

model2_99<-svyglm(tbdruind_bin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(HIV),family=quasibinomial,design=datdesign_99,na.action=na.omit)
mod2_99<-summary(model2_99)
round(mod2_99$coefficients,3)

model3<-svyglm(tbdruind_bin~factor(SDDSRVYR)+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(esrd),family=quasibinomial,design=datdesign,na.action=na.omit)
mod3<-summary(model3)
round(mod3$coefficients,3)

model3_11<-svyglm(tbdruind_bin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(esrd),family=quasibinomial,design=datdesign_11,na.action=na.omit)
mod3<-summary(model3_11)
round(mod3$coefficients,3)

model3_99<-svyglm(tbdruind_bin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(esrd),family=quasibinomial,design=datdesign_99,na.action=na.omit)
mod3<-summary(model3_99)
round(mod3$coefficients,3)

model4<-svyglm(tbdruind_bin~factor(SDDSRVYR)+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(Diab),family=quasibinomial,design=datdesign,na.action=na.omit)
mod4<-summary(model4)
round(mod4$coefficients,3)

model4_race<-svyglm(tbdruind_bin~factor(SDDSRVYR)+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(Diab),family=quasibinomial,design=datdesign_11,na.action=na.omit)
mod4<-summary(model4_race)
round(mod4$coefficients,3)

model4_11<-svyglm(tbdruind_bin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(Diab),family=quasibinomial,design=datdesign_11,na.action=na.omit)
mod4<-summary(model4_11)
round(mod4$coefficients,3)

model4_99<-svyglm(tbdruind_bin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(Diab),family=quasibinomial,design=datdesign_99,na.action=na.omit)
mod4<-summary(model4_99)
round(mod4$coefficients,3)

model5<-svyglm(tbdruind_bin~factor(SDDSRVYR)+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(HHTB),family=quasibinomial,design=datdesign,na.action=na.omit)
mod5<-summary(model5)
z<-round(mod5$coefficients,3)
write.table(z,file="zexport2.csv",sep=",")

model5_11<-svyglm(tbdruind_bin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(HHTB),family=quasibinomial,design=datdesign_11,na.action=na.omit)
mod5<-summary(model5_11)
z<-round(mod5$coefficients,3)
write.table(z,file="zzexport.csv",sep=",")

model5_99<-svyglm(tbdruind_bin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(HHTB),family=quasibinomial,design=datdesign_99,na.action=na.omit)
mod5<-summary(model5_99)
z<-round(mod5$coefficients,3)
write.table(z,file="zzexport.csv",sep=",")

model6<-svyglm(tbdruind_bin~factor(SDDSRVYR)+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(EverActiveTB),family=quasibinomial,design=datdesign,na.action=na.omit)
mod6<-summary(model6)
z<-round(mod6$coefficients,3)
write.table(z,file="zzexport4.csv",sep=",")

model6_11<-svyglm(tbdruind_bin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(EverActiveTB),family=quasibinomial,design=datdesign_11,na.action=na.omit)
mod6<-summary(model6_11)
z_11<-round(mod6$coefficients,3)
write.table(z_11,file="zzexport_11.csv",sep=",")

model6_99<-svyglm(tbdruind_bin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(EverActiveTB),family=quasibinomial,design=datdesign_99,na.action=na.omit)
mod6_99<-summary(model6_99)
z_99<-round(mod6_99$coefficients,3)
write.table(z_99,file="zzexport_99.csv",sep=",")

#IGRA
model7<-svyglm(lbxtbin~factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(HIV)+factor(Diab)+factor(HHTB)+factor(EverActiveTB)+factor(esrd)+factor(immu)+factor(PRED)+factor(tnf),family=quasibinomial,design=datdesign,na.action=na.omit)
mod7<-summary(model7)
z_7<-round(mod7$coefficients,3)
write.table(z_7,file="z7.csv",sep=",")

model8<-svyglm(lbxtbin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(HIV),family=quasibinomial,design=datdesign,na.action=na.omit)
mod8<-summary(model8)
z_8<-round(mod8$coefficients,3)
write.table(z_8,file="z8.csv",sep=",")

out.coef.mod8<-orfun(mod=model8)
z8<-round(out.coef.mod8,3)
write.table(z8,file="zexport8.csv",sep=",")

model9<-svyglm(lbxtbin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(esrd),family=quasibinomial,design=datdesign,na.action=na.omit)
mod9<-summary(model9)
z9<-round(mod9$coefficients,3)
write.table(z9,file="z9.csv",sep=",")

out.coef.mod9<-orfun(mod=model9)
z9<-round(out.coef.mod9,3)
write.table(z9,file="zexport9.csv",sep=",")

model10<-svyglm(lbxtbin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(Diab),family=quasibinomial,design=datdesign,na.action=na.omit)
mod10<-summary(model10)
z10<-round(mod10$coefficients,3)
write.table(z10,file="z10.csv",sep=",")

out.coef.mod10<-orfun(mod=model10)
z10<-round(out.coef.mod10,3)
write.table(z10,file="zexport10.csv",sep=",")

model11<-svyglm(lbxtbin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(HHTB),family=quasibinomial,design=datdesign,na.action=na.omit)
mod11<-summary(model11)
z11<-round(mod11$coefficients,3)
write.table(z11,file="z11.csv",sep=",")

out.coef.mod11<-orfun(mod=model11)
z11<-round(out.coef.mod11,3)
write.table(z11,file="zexport11.csv",sep=",")

model12<-svyglm(lbxtbin~SDDSRVYR+factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(EverActiveTB),family=quasibinomial,design=datdesign,na.action=na.omit)
mod12<-summary(model12)
z12<-round(mod12$coefficients,3)
write.table(z12,file="z12.csv",sep=",")

out.coef.mod12<-orfun(mod=model12)
z12<-round(out.coef.mod12,3)
write.table(z12,file="zexport12.csv",sep=",")

#install.packages("writexl")
#library("writexl")
#write.table(d,file="dddexport.csv",sep=",")

#Odds Ratio
#mod<-mod1
orfun<-function(mod){
  sum.mod<-summary(mod)
  mean=exp(sum.mod$coefficients[,1])
  cil=exp(sum.mod$coefficients[,1]+qnorm(c(0.025))*sum.mod$coefficients[,2])
  ciu=exp(sum.mod$coefficients[,1]+qnorm(c(0.975))*sum.mod$coefficients[,2])
  out=cbind(mean=mean,cil=cil,ciu=ciu)
  return(out=out)
}

#invlgt <- function(x) 1/(1+exp(-x))
out.coef.mod6_99<-orfun(mod=model6_99)
z6_99<-round(out.coef.mod6_99,3)
write.table(z6_99,file="zexport6_99.csv",sep=",")

orfun(mod=model4)
#round(out.coef.mod1[,2],3)

out.coef.mod2<-orfun(mod=model2)

save(newdata99_11,file="newdata99_11.rda")
save(newdata99,file="newdata99_riskfactors.rda")
save(newdata11,file="newdata11_riskfactors.rda")

# model3<-glm(tbdruind_bin~factor(AGECAT)*RIAGENDR*DMDBORN*SDDSRVYR+LBDHI+DIQ010+TBQ060+TBQ040+esrd+immu+tnf, data=newdata99_11, family = binomial)
# summary(model3)
# 
# model4<-glm(tbdruind_bin~factor(AGECAT)+factor(RIAGENDR)+DMDBORN+SDDSRVYR+LBDHI+DIQ010+TBQ060+TBQ040+esrd+immu+tnf, data=newdata99_11, family = quasibinomial)
# summary(model4)

#Logistic regression
#logit1 <- (svyglm(paq665~factor(hsd010)+ridageyr, family=quasibinomial, design=nhc, na.action = na.omit))
#summary(logit1)

#NHANES III, NHANES IV till March 2020 all using the same Prescription Medications- Drug Information (RXQ_DRUG) data file
# dat.drug.loc<-"https://wwwn.cdc.gov/Nchs/Nhanes/1999-2000/RXQ_DRUG.xpt"
# dat.drug.raw<-DownloadImport(dat.drug.loc)
# 
# #NHANES IV Prescription Medications (RXQ_RX_J)
# pm1112loc<-"https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/RXQ_RX_G.XPT"
# dat.pm1112<-DownloadImport(pm1112loc)
# sub.pm1112<-dat.pm1112[,c("SEQN","RXDUSE","RXDDRUG","RXDDRGID")]
# 
# pm9900loc<-"https://wwwn.cdc.gov/Nchs/Nhanes/1999-2000/RXQ_RX.XPT"
# dat.pm9900<-DownloadImport(pm9900loc)
# sub.pm9900<-dat.pm9900[,c("SEQN","RXD030","RXD240B","RXDDRGID")]
# colnames(sub.pm9900)[c(2,3)]<-c("RXDUSE","RXDDRUG")
# #sub.pmn4<-rbind(sub.pm9900,sub.pm0102,sub.pm0304,sub.pm0506,sub.pm0708,
# #                sub.pm0910,sub.pm1112,sub.pm1314,sub.pm1516,sub.pm1718)
# sub.pmn4<-rbind(sub.pm9900,sub.pm1112)
#                 
# #link prescription medicines dataset with drug information dataset
# dat.pm.druginf<-merge(sub.pmn4,dat.drug.raw,by="RXDDRGID")
# dim(sub.pmn4)
# dim(dat.drug.raw)
# dim(dat.pm.druginf)
# #NHANES IV demographic profile locations

# svyby(~take.tnfalpha,~SDDSRVYR,design = datdesign,svyciprop,vartype = "ci",na.rm=T)
# svyby(~take.immuagents,~SDDSRVYR,design=datdesign,svyciprop,vartype = "ci",na.rm=T)
# 
# svyby(~take.immuagents,~RIAGENDR,design=datdesign,svyciprop,vartype = "ci",na.rm=T)

