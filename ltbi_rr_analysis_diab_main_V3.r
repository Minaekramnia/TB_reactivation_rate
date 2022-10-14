#  setwd

# rm(list=ls())
##### SET UP  
# Libraries
library(survey)
library(plyr) 
library(MASS)
library(stringr)

# Functions
invlgt <- function(x) 1/(1+exp(-x))

mulohi <- function(x) as.numeric(c(mean(x),quantile(x,c(1,39)/40)))

pretty <- function(x0,r) { x <- gsub(" ","",format(round(x0,r),nsmall=r));
                          paste0(x[1],"\n(",x[2],", ",x[3],")") }

gammapar2 <- function(tgt) {
  tgt <- as.numeric(tgt)
  shape= tgt[1]^2/tgt[2]^2
  rate = tgt[1]/tgt[2]^2
  return(c(shape,rate))  }

betapar <- function(tgt) {
  tgt <- as.numeric(tgt)
  mn <- tgt[1]; cir <- (tgt[3]-tgt[2])
  xopt <- function(xx,mn=mn,cir=cir) {
    cir2 <- qbeta(c(1,39)/40,xx*(mn/(1-mn)),xx); 
    cir2 <- cir2[2]-cir2[1]
    sum((cir2-cir)^2) }
  zz <- optimize(xopt,c(0.2,100000),mn=mn,cir=cir)$minimum
  bp <-  c(zz*(mn/(1-mn)),zz)
  if(sum(abs(bp-1))<0.2) { c(1,1) } else { bp }  }

###### A. Estimate regression for IGRA positivity from NHANES DATA  
      # Regression used for both risk-group and overall

### Data import
load("datdesign_11.Rdata")
datdesign_11$variables$AGECAT[datdesign_11$variables$AGECAT=="age_75_80"] <- "age_75plus"

### Fit model
fit<-svyglm(lbxtbin~factor(AGECAT)+factor(Sex)+factor(USNUSB)+factor(racenha),
              family=quasibinomial,design=datdesign_11,na.action=na.omit)

### Create table for IGRA estimates by stratum, with/without RF
load("newdata11_sub.Rdata")

age_cats  <- c("age_15_24","age_25_34","age_35_44","age_45_54",
               "age_55_64","age_65_74","age_75plus")

new_dat <- as.data.frame(expand.grid(unique(newdata11_sub$Sex),age_cats,
            unique(newdata11_sub$racenha),unique(newdata11_sub$USNUSB)))
names(new_dat) <- c("Sex","AGECAT","racenha","USNUSB")

new_dat  <- new_dat[new_dat$USNUSB >=1 & new_dat$USNUSB <=2,]
new_dat$Diab <- unique(newdata11_sub$Diab)[2] # Has diabetes

# Duplicating for with/without RF
new_dat_rf <- new_dat_nr <- new_dat
new_dat_nr$Diab <- unique(newdata11_sub$Diab)[1] # No diabetes
new_dat <- rbind(new_dat_rf,new_dat_nr)

# Predict for new data from fitted regression model
pred <- predict(fit,vcov=T,se.fit=T,type="link",newdata=new_dat)
vcv  <- vcov(pred)
est  <- coef(pred)

# Draw random samples for igra values
n.sim <- 10000 # no samples
set.seed(1)
sims0 <- mvrnorm(n.sim, mu = est, Sigma = vcv)
igra_prev_sims  <- t(invlgt(sims0))

# Turning columns similar to Yunfei NHIS
new_dat$Sex <- replace(new_dat$Sex, new_dat$Sex==1, "Male")
new_dat$Sex <- replace(new_dat$Sex, new_dat$Sex==2, "Female")
new_dat$USNUSB <- replace(new_dat$USNUSB, new_dat$USNUSB==1, "US")
new_dat$USNUSB <- replace(new_dat$USNUSB, new_dat$USNUSB==2, "NUS")

# Concatenate columns to become similar to NHIS format: 
new_dat$concat <- gsub(" ", "", paste(new_dat$Sex,new_dat$AGECAT,new_dat$racenha,new_dat$USNUSB))

# Seperate out rf and nr
new_dat_nr  <- new_dat[new_dat$Diab==0,]
new_dat_rf  <- new_dat[new_dat$Diab==1,]

igra_prev_sims_nr  <- igra_prev_sims[new_dat$Diab==0,]
igra_prev_sims_rf  <- igra_prev_sims[new_dat$Diab==1,]

###### B. Attach population size and simulate pop uncertainty
### Data import
load("pop_dat_Oct-12-2022.rData") # pop_dat
estdiabetes <- read.csv("dat.estdiabetes.nhis.csv", header=TRUE, stringsAsFactors=FALSE)

estdiabetes$Sex    <- ifelse(str_detect(estdiabetes$X,"Male"),"Male","Female")
estdiabetes$USNUSB <- ifelse(str_detect(estdiabetes$X,"NUS"),"NUS","US")
estdiabetes$AGECAT <- NA
estdiabetes$AGECAT[str_detect(estdiabetes$X,"15-24")] <- "age_15_24"
estdiabetes$AGECAT[str_detect(estdiabetes$X,"25-34")] <- "age_25_34"
estdiabetes$AGECAT[str_detect(estdiabetes$X,"35-44")] <- "age_35_44"
estdiabetes$AGECAT[str_detect(estdiabetes$X,"45-54")] <- "age_45_54"
estdiabetes$AGECAT[str_detect(estdiabetes$X,"55-64")] <- "age_55_64"
estdiabetes$AGECAT[str_detect(estdiabetes$X,"65-74")] <- "age_65_74"
estdiabetes$AGECAT[str_detect(estdiabetes$X,"75-80")] <- "age_75plus"
estdiabetes$AGECAT[str_detect(estdiabetes$X,"80+"  )] <- "age_75plus"
estdiabetes$racenha <- NA
estdiabetes$racenha[str_detect(estdiabetes$X,"NHW")] <- "NHW"
estdiabetes$racenha[str_detect(estdiabetes$X,"NHA")] <- "NHA"
estdiabetes$racenha[str_detect(estdiabetes$X,"NHB")] <- "NHB"
estdiabetes$racenha[str_detect(estdiabetes$X,"others")] <- "others"
estdiabetes$racenha[str_detect(estdiabetes$X,"Hispanic")] <- "Hispanic"
# estdiabetes$variance <- estdiabetes$SE^2

names(estdiabetes)
estdiabetes2a <- aggregate(total~Sex+AGECAT+racenha+USNUSB,estdiabetes,sum)
estdiabetes2b <- aggregate(SE~Sex+AGECAT+racenha+USNUSB,estdiabetes,function(x) sqrt(sum(x^2)))

estdiabetes2 <- join(x = estdiabetes2a, y = estdiabetes2b,by=c("Sex","AGECAT","racenha","USNUSB"))

# Add pop data to new_dat2
new_dat2nr <- join(x = new_dat_nr, y = estdiabetes2,by=c("Sex","AGECAT","racenha","USNUSB"))
new_dat2rf <- join(x = new_dat_rf, y = estdiabetes2,by=c("Sex","AGECAT","racenha","USNUSB"))

# set missing pop categories to zero, for rf
new_dat2rf$total[is.na(new_dat2rf$total)] <- 0
new_dat2rf$SE[is.na(new_dat2rf$SE)] <- 0

# Replace pop data for overall
new_dat2nr <- join(x = new_dat2nr, y = pop_dat,by=c("Sex","AGECAT","racenha","USNUSB"))
new_dat2rf <- join(x = new_dat2rf, y = pop_dat,by=c("Sex","AGECAT","racenha","USNUSB"))

# Apply correction for different age groups
new_dat2nr$total2 <- new_dat2nr$total
new_dat2rf$total2 <- new_dat2rf$total
new_dat2nr$SE2 <- new_dat2nr$SE
new_dat2rf$SE2 <- new_dat2rf$SE

for(i in which(!is.na(new_dat2rf$pop_inflate))){
  new_dat2nr$total2[i] <- new_dat2nr$total[i]*new_dat2nr$pop_inflate[i]
  new_dat2rf$total2[i] <- new_dat2rf$total[i]*new_dat2rf$pop_inflate[i]
  new_dat2nr$SE2[i] <- new_dat2nr$SE[i]*new_dat2nr$pop_inflate[i]
  new_dat2rf$SE2[i] <- new_dat2rf$SE[i]*new_dat2rf$pop_inflate[i]
}

new_dat2nr$total2 <- new_dat2nr$pop2 - new_dat2rf$total2
new_dat2nr$SE2    <- sqrt(new_dat2nr$pop_se2^2 +  new_dat2rf$SE2^2)

new_dat2nr$SE2[new_dat2nr$total2<0]    <- 0
new_dat2nr$total2[new_dat2nr$total2<0] <- 0

# Simulate pop uncertainty
pop_sims_rf <- pop_sims_nr <- matrix(NA,nrow(new_dat2nr),n.sim)

set.seed(2)
for(i in 1:nrow(new_dat2nr)){ # i=15
  # rf
  if(new_dat2rf$total2[i]==0){
    pop_sims_rf[i,] <- 0
  } else {
    tmp <- gammapar2(c(new_dat2rf$total2[i], new_dat2rf$SE2[i]))
    pop_sims_rf[i,] <- rgamma(n.sim,tmp[1],tmp[2])
  }
  # nr
  if(new_dat2nr$total2[i]==0){
    pop_sims_nr[i,] <- 0
  } else {
    tmp <- gammapar2(c(new_dat2nr$total2[i], new_dat2nr$SE2[i]))
    pop_sims_nr[i,] <- rgamma(n.sim,tmp[1],tmp[2])
  }
}

###### C. Calc igra+ population size simulations (nr and rf)
igra_pop_sims_rf <- pop_sims_rf * igra_prev_sims_rf
igra_pop_sims_nr <- pop_sims_nr * igra_prev_sims_nr

###### D. Calc ltbi population size simulations
# Note that....
# igra = ltbi*sens + (1-ltbi) * (1-spec)
#      = ltbi*sens + 1-ltbi-spec+ltbi*spec 
#      = 1 + ltbi*(sens+spec-1) - spec 
# ltbi = (igra-1+spec)/(sens+spec-1)

# Stout : sens                   spec
# FB      78.9 (69.6 to 90.2)    98.5 (96.1 to 99.8)
# USB     78.0 (65.0 to 91.0)    97.9 (96.0 to 99.4)

# Simulate sens/spec of IGRA by US/NUS status
sens_par_nus <- betapar(c(78.9, 69.6, 90.2)/100 )
sens_par_usb <- betapar(c(78.0, 65.0,91.0)/100 )

spec_par_nus <- betapar(c(98.5, 96.1, 99.8)/100 )
spec_par_usb <- betapar(c(97.9, 96.0, 99.4)/100 )

set.seed(3)
sens_nus_sims <- rbeta(n.sim,sens_par_nus[1],sens_par_nus[2])
sens_usb_sims <- rbeta(n.sim,sens_par_usb[1],sens_par_usb[2])
spec_nus_sims <- rbeta(n.sim,spec_par_nus[1],spec_par_nus[2])
spec_usb_sims <- rbeta(n.sim,spec_par_usb[1],spec_par_usb[2])

sens_sims <- spec_sims <- pop_sims_rf
sens_sims[,] <- spec_sims[,] <- NA

sens_sims[new_dat2rf$USNUSB=="US", ] <- outer(rep(1,sum(new_dat2rf$USNUSB=="US" )),sens_usb_sims)
sens_sims[new_dat2rf$USNUSB=="NUS",] <- outer(rep(1,sum(new_dat2rf$USNUSB=="NUS")),sens_nus_sims)

spec_sims[new_dat2rf$USNUSB=="US", ] <- outer(rep(1,sum(new_dat2rf$USNUSB=="US" )),spec_usb_sims)
spec_sims[new_dat2rf$USNUSB=="NUS",] <- outer(rep(1,sum(new_dat2rf$USNUSB=="NUS")),spec_nus_sims)

# Calculate LTBI prev
ltbi_prev_sims_rf <- (igra_prev_sims_rf-1+spec_sims)/(sens_sims+spec_sims-1)
ltbi_prev_sims_nr <- (igra_prev_sims_nr-1+spec_sims)/(sens_sims+spec_sims-1)

ltbi_prev_sims_rf[ltbi_prev_sims_rf<0] <- 0  #  ??  Or bayesian estimation...
ltbi_prev_sims_nr[ltbi_prev_sims_nr<0] <- 0 

# Calc LTBI pop size
ltbi_pop_sims_rf <- pop_sims_rf * ltbi_prev_sims_rf
ltbi_pop_sims_nr <- pop_sims_nr * ltbi_prev_sims_nr

###### E. Create samples of cases
# Load data
load("reactivation_rates/ntss_dummy_50DC_2021-05-19.rData")

# Subset cases in 2011-12 age over 17
cases_all <- tbdummy50DC[tbdummy50DC$YEAR%in%(2011:2012) &
                       tbdummy50DC$AGE>14,]
cases_all <- cases_all[cases_all$ORIGIN!="UNK",]

cases_all$HIV <- "MISSING"
cases_all$HIV[cases_all$HIVSTAT=="NEG"] <- "NEG"
cases_all$HIV[cases_all$HIVSTAT=="POS"] <- "POS"

cases_all$PREVTB2 <- "MISSING"
cases_all$PREVTB2[cases_all$PREVTB=="N"] <- "N"
cases_all$PREVTB2[cases_all$PREVTB=="Y"] <- "Y"

# Fit imputation model
library(mgcv)

fit_z <- gam(rt~YEAR+s(AGE)+SEX+ORIGIN+RACEHISP+RISKDIAB+RISKRENAL+HIV+PREVTB2,data=cases_all[!is.na(cases_all$rt),],family=binomial())
pred_z <- predict.gam(fit_z,newdata=cases_all[is.na(cases_all$rt),],type="response")

cases_all$rt_imp <- cases_all$rt
cases_all$rt_imp[is.na(cases_all$rt)] <- pred_z

# Tabulate by nr and rf
cases_all$AGECAT <- NA
cases_all$AGECAT[cases_all$AGE%in%(15:24)]   <- age_cats[1]
cases_all$AGECAT[cases_all$AGE%in%(25:34)]   <- age_cats[2]
cases_all$AGECAT[cases_all$AGE%in%(35:44)]   <- age_cats[3]
cases_all$AGECAT[cases_all$AGE%in%(45:54)]   <- age_cats[4]
cases_all$AGECAT[cases_all$AGE%in%(55:64)]   <- age_cats[5]
cases_all$AGECAT[cases_all$AGE%in%(65:74)]   <- age_cats[6]
cases_all$AGECAT[cases_all$AGE%in%(75:100)]  <- age_cats[7]

race_cats <- unique(new_dat2nr$racenha)
cases_all$racenha <- NA
cases_all$racenha[cases_all$RACEHISP=="ASIAN"] <- "NHA"
cases_all$racenha[cases_all$RACEHISP=="WHITE"] <- "NHW"
cases_all$racenha[cases_all$RACEHISP=="BLACK"] <- "NHB"
cases_all$racenha[cases_all$RACEHISP=="HISP"]  <- "Hispanic"
cases_all$racenha[cases_all$RACEHISP=="NAHAW"] <- "others"
cases_all$racenha[cases_all$RACEHISP=="MULT"]  <- "others"
cases_all$racenha[cases_all$RACEHISP=="AMIND"] <- "others"
cases_all$racenha[cases_all$RACEHISP=="UNK"]   <- "others"

cases_all$USNUSB <- NA
cases_all$USNUSB[cases_all$ORIGIN=="NONUSB"] <- "NUS"
cases_all$USNUSB[cases_all$ORIGIN=="USBORN"] <- "US"

cases_all$Sex <- NA
cases_all$Sex[cases_all$SEX=="M"] <- "Male"
cases_all$Sex[cases_all$SEX=="F"] <- "Female"

# Multiple imputation
imp_cases <- matrix(NA,sum(is.na(cases_all$rt)),n.sim)
tmp       <- cases_all$rt_imp[is.na(cases_all$rt)]
set.seed(4)
for(i in 1:nrow(imp_cases)){
  imp_cases[i,] <- rbinom(n.sim,1,tmp[i])
}

# Tabulate cases
sum_cases_rf <- sum_cases_nr <- pop_sims_rf
sum_cases_rf[,] <- sum_cases_nr[,] <- NA
tmp2 <- cases_all$rt

for(i in 1:nrow(sum_cases_rf)){ # i=135
  id <- cases_all$AGECAT==new_dat2nr$AGECAT[i] &
        cases_all$racenha==new_dat2nr$racenha[i] &
        cases_all$USNUSB==new_dat2nr$USNUSB[i] &
        cases_all$Sex==new_dat2nr$Sex[i]
  id_rf <- id & cases_all$RISKDIAB=="Y"
  id_nr <- id & cases_all$RISKDIAB==""
  for(j in 1:ncol(sum_cases_rf)){ # j=1
    tmp2[is.na(cases_all$rt)] <- imp_cases[,j]  
    sum_cases_rf[i,j] <- sum((1-tmp2)[id_rf])
    sum_cases_nr[i,j] <- sum((1-tmp2)[id_nr])
  }
  cat('\r',i,"     "); flush.console()
}

# Simulate from Poisson
cases_sims_rf <- cases_sims_nr <- pop_sims_rf
cases_sims_rf[,] <- cases_sims_nr[,] <- NA

set.seed(5)
for(i in 1:nrow(cases_sims_rf)){
  cases_sims_rf[i,] <- rpois(n.sim,sum_cases_rf[i,])
  cases_sims_nr[i,] <- rpois(n.sim,sum_cases_nr[i,])
}

###### F. Calculate results
#  cases_sims_rf      sum_cases_nr
#  pop_sims_rf        pop_sims_nr
#  ltbi_pop_sims_rf   ltbi_pop_sims_nr
#  ltbi_prev_sims_rf  ltbi_prev_sims_nr
#  igra_pop_sims_rf   igra_pop_sims_nr
#  igra_prev_sims_rf  igra_prev_sims_nr

res_nam <- c("cases","igra_prev","ltbi_prev","risk_pop",
             "py_igra","py_ltbi","r_react_igra","r_react_ltbi",
             "rr_igra","rr_ltbi","rr_igra_adj","rr_ltbi_adj")
res_tab <- matrix(NA,3,length(res_nam))
rownames(res_tab) <- c("mean","ci_lo","ci_hi")
colnames(res_tab) <- res_nam

# cases
res_tab[,"cases"] <- mulohi(colSums(cases_sims_rf))

# igra_prev
res_tab[,"igra_prev"] <- mulohi(colSums(igra_pop_sims_rf)/colSums(pop_sims_rf))*100

# ltbi_prev
res_tab[,"ltbi_prev"] <- mulohi(colSums(ltbi_pop_sims_rf)/colSums(pop_sims_rf))*100

# risk_pop
res_tab[,"risk_pop"] <- mulohi(colSums(pop_sims_rf))/1e3

# py_igra
res_tab[,"py_igra"] <- mulohi(colSums(igra_pop_sims_rf*2))/1e3 # *2 for 2 year period

# py_ltbi
res_tab[,"py_ltbi"] <- mulohi(colSums(ltbi_pop_sims_rf*2))/1e3 # *2 for 2 year period

# r_react_igra
res_tab[,"r_react_igra"] <- mulohi(colSums(cases_sims_rf)/colSums(igra_pop_sims_rf*2))*100 

# r_react_ltbi
res_tab[,"r_react_ltbi"] <- mulohi(colSums(cases_sims_rf)/colSums(ltbi_pop_sims_rf*2))*100 

# rr_igra
rr_rf_i <- colSums(cases_sims_rf)/colSums(igra_pop_sims_rf*2)
rr_nr_i <- colSums(cases_sims_nr)/colSums(igra_pop_sims_nr*2)
res_tab[,"rr_igra"] <- mulohi(rr_rf_i/rr_nr_i) 

# rr_ltbi
rr_rf_i <- colSums(cases_sims_rf)/colSums(ltbi_pop_sims_rf*2)
rr_nr_i <- colSums(cases_sims_nr)/colSums(ltbi_pop_sims_nr*2)
res_tab[,"rr_ltbi"] <- mulohi(rr_rf_i/rr_nr_i) 

# rr_igra_adj
adj <- igra_pop_sims_rf/igra_pop_sims_nr
adj[adj>1] <- 1
rr_rf_i <- colSums(cases_sims_rf)/colSums(igra_pop_sims_rf*2)
rr_nr_i <- colSums(cases_sims_nr*adj)/colSums(igra_pop_sims_nr*adj*2)
res_tab[,"rr_igra_adj"] <- mulohi(rr_rf_i/rr_nr_i) 

# rr_ltbi_adj
adj <- ltbi_pop_sims_rf/ltbi_pop_sims_nr
adj[adj>1] <- 1
adj[is.nan(adj)] <- mean(adj[!is.nan(adj)])
rr_rf_i <- colSums(cases_sims_rf)/colSums(ltbi_pop_sims_rf*2)
rr_nr_i <- colSums(cases_sims_nr*adj)/colSums(ltbi_pop_sims_nr*adj*2)
res_tab[,"rr_ltbi_adj"] <- mulohi(rr_rf_i/rr_nr_i) 

# Save detailed results
write.csv(res_tab,file="diab_res_Oct-12-2022.csv")

# Round to 3 sig fig
res_tab2 <- res_tab
for(i in 1:ncol(res_tab)){ # i=1
  res_tab2[,i] <- signif(res_tab[,i],3)
}

# Seperate IGRA and LTBI tables
res_tab2_igra <- res_tab2[,c(1,2,4,5,7,9 ,11)]
res_tab2_ltbi <- res_tab2[,c(1,3,4,6,8,10,12)]

# Format for publication table
res_tab3_igra <- res_tab2_igra[1,]
res_tab3_ltbi <- res_tab2_ltbi[1,]


n_dec <- c(0,2,0,0,4,2,2)
for(i in 1:ncol(res_tab2_igra)){
  res_tab3_igra[i] <- pretty(res_tab2_igra[,i],n_dec[i])
  res_tab3_ltbi[i] <- pretty(res_tab2_ltbi[,i],n_dec[i])
}

# Save 
write.csv(res_tab3_igra,file="diab_res_igra_Oct-12-2022.csv")
write.csv(res_tab3_ltbi,file="diab_res_ltbi_Oct-12-2022.csv")


