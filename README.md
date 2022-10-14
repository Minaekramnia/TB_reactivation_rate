# TB_reactivation_rate

Calculating Latent TB reactivation rates for TB risk populations in the United States!


For individuals infected with Mycobacterium tuberculosis (Mtb), the risk of developing TB disease extends many years after the initial infection [1, 2]. As a consequence of this, in countries that have been able to reduce current Mtb transmission to low levels a substantial fraction of incident TB cases will result from historical infections. In the United States for 2019-2020, 87.5% of all TB cases could not be attributed to recent transmission, indicating the TB case was likely due to infection greater than 2 years prior to diagnosis [3].  These cases of ‘reactivation TB’ arise among the population of individuals with a latent TB infection (LTBI) acquired earlier in life. For countries like the United States, identifying individuals with LTBI and providing TB preventive therapy (TPT) is a major focus of TB prevention activities [4, 5].

Statistical analysis
We took a Monte Carlo simulation approach to propagate uncertainty in each analytic input through subsequent steps in the analysis. For sampling uncertainty in reactivation TB case totals derived from NTSS we constructed Poisson distributions centered at the observed TB case count for each analytic stratum and risk group, and incorporating the additional uncertainty introduced by the imputation of the recent transmission variable. For uncertainty in IGRA positivity estimates we extracted coefficient values and the variance-covariance matrix from the logistic regression models fit to NHANES data, and used these to parameterize a multivariate Normal distribution. We simulated from this distribution, calculated predicted values for each analytic stratum and inverse-logit transformed these values to obtain random IGRA positivity values. For IGRA sensitivity and specificity we constructed Gamma distributions matching the mean values and uncertainty intervals reported by Stout et al [12]. For population estimates derived from NHANES (prior TB), NHIS (diabetes, chronic renal disease), HIV Surveillance Reports (HIV) and ACS (general population) we constructed Gamma distributions matching the mean and uncertainty estimates available from each data source (Gamma distributions were used to avoid simulating negative population values), and assuming independence between the population totals in each stratum. We simulated 10,000 values from each of these distributions, and used these to calculate 10,000 estimates for each study outcome, representing the combined uncertainty in all model inputs. We calculated point estimates as the mean of the simulated values for each outcome, and equal-tailed 95% uncertainty intervals as the 2.5th and 97.5th percentiles of these distributions.


References: 
	Sutherland I. The ten-year incidence of clinical tuberculosis following “conversion” in 2550 individuals aged 14 to 19 years. TSRU Progress Report. The Hague: 1968.
2.	Ferebee SH, Mount FW. Tuberculosis morbidity in a controlled trial of the prophylactic use of isoniazid among household contacts. Am Rev Respir Dis. 1962;85:490-510.
3.	U.S. Centers for Disease Control and Prevention. Reported Tuberculosis in the United States, 2020 [retrieved from https://www.cdc.gov/tb/statistics/reports/2020/default.htm, Oct 6 2022]. Atlanta GA: U.S. Centers for Disease Control and Prevention, 2021.
4.	Bibbins-Domingo K, Grossman DC, Curry SJ, Bauman L, Davidson KW, Epling JW, et al. Screening for latent tuberculosis infection in adults: US Preventive Services Task Force recommendation statement. JAMA. 2016;316(9):962-9.
5.	US Centers for Disease Control and Prevention. Latent tuberculosis infection: A guide for primary health care providers [retrieved from https://www.cdc.gov/tb/publications/ltbi/pdf/targetedltbi.pdf, June 3 2019]. Atlanta GA: US Centers for Disease Control and Prevention, 2013.
