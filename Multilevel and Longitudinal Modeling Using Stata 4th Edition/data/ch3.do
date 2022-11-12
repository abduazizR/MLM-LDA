************************************************************
****Chapter 3: Random-intercept models with covariates *****
************************************************************


*** 3.2 Does smoking during pregnancy affect birthweight?

** 3.2.1 Data structure and descriptive statistics

use https://www.stata-press.com/data/mlmus4/smoking, clear

quietly xtset momid 
xtsum birwt smoke black

egen pickone = tag(momid)
summarize black if pickone==1

egen num = count(birwt), by(momid)
tabulate num if pickone==1

quietly xtset momid
xttab smoke


*** 3.4 Estimation using Stata

** 3.4.1 Using xtreg
quietly xtset momid
xtreg birwt smoke male mage hsgrad somecoll collgrad       ///
  married black kessner2 kessner3 novisit pretri2 pretri3, mle vce(robust) 

** 3.4.2 Using mixed

mixed birwt smoke male mage hsgrad somecoll collgrad married black ///
  kessner2 kessner3 novisit pretri2 pretri3 || momid:, mle vce(robust) ///
  stddeviations

estat icc
estat wcorrelation


*** 3.5 Coefficients of determination or variance explained

quietly xtset momid
xtreg birwt, mle vce(robust)

quietly xtset momid
xtreg birwt hsgrad somecoll collgrad married black, mle vce(robust)


*** 3.6 Hypothesis tests and confidence intervals

** 3.6.1 Hypothesis tests for individual regression coefficients

mixed birwt smoke male mage hsgrad somecoll collgrad      ///
  married black kessner2 kessner3 novisit pretri2 pretri3 ///
  || momid:, reml dfmethod(kroger) stddeviations
 
estat df

** 3.6.2 Joint hypothesis tests for several regression coefficients

quietly xtset momid
quietly xtreg birwt smoke male mage hsgrad somecoll collgrad  ///
  married black kessner2 kessner3 novisit pretri2 pretri3, mle vce(robust)
testparm kessner2 kessner3

* Likelihood-ratio test
quietly xtset momid
quietly xtreg birwt smoke male mage hsgrad somecoll collgrad  ///
  married black kessner2 kessner3 novisit pretri2 pretri3, mle
estimates store full
  
quietly xtreg birwt smoke male mage hsgrad somecoll collgrad   ///
  married black novisit pretri2 pretri3, mle
lrtest full .

** 3.6.3 Predicted means and confidence intervals
generate education = hsgrad*1 + somecoll*2 + collgrad*3
label define ed 0 "No HS Degree" 1 "HS Degree" 2 "Some Coll" 3 "College", replace
label values education ed

quietly xtset momid
quietly xtreg birwt i.smoke male mage i.education married black ///
  kessner2 kessner3 novisit pretri2 pretri3, mle vce(robust)
  
margins i.smoke#i.education 

marginsplot, xdimension(education)


*** 3.7 Between and within effects

** 3.7.1 Between-mother effects

quietly xtset momid
xtreg birwt smoke male mage hsgrad somecoll collgrad married black kessner2 ///
   kessner3 novisit pretri2 pretri3, be
  
** 3.7.2 Within-mother effects

quietly xtset momid
xtreg birwt smoke male mage kessner2 kessner3 novisit pretri2 pretri3, ///
   fe vce(robust)

** 3.7.5 Conventional Hausman test

quietly xtset momid
quietly xtreg birwt smoke male mage hsgrad somecoll collgrad married ///
   black kessner2 kessner3 novisit pretri2 pretri3, fe
estimates store fixed

quietly xtreg birwt smoke male mage hsgrad somecoll collgrad married ///
   black kessner2 kessner3 novisit pretri2 pretri3, mle
estimates store ml

hausman fixed ml, equations(1:1)

*** 3.7.6 Allowing for different within- and between effects

egen mn_smok = mean(smoke), by(momid)
generate dev_smok = smoke - mn_smok

quietly xtset momid
xtreg birwt dev_smok mn_smok male mage hsgrad somecoll collgrad married ///
  black kessner2 kessner3 novisit pretri2 pretri3, mle vce(robust)

lincom mn_smok - dev_smok

* Means for other time-varying covariates 
egen mn_male = mean(male), by(momid)
egen mn_mage = mean(mage), by(momid)
egen mn_kessner2 = mean(kessner2), by(momid)
egen mn_kessner3 = mean(kessner3), by(momid)
egen mn_novisit = mean(novisit), by(momid)
egen mn_pretri2 = mean(pretri2), by(momid)
egen mn_pretri3 = mean(pretri3), by(momid)

/*
foreach var of varlist male mage kessner* novisit pretri* {
   egen mn_`var' = mean(`var'), by(momid)
}
*/

*** 3.7.7 Robust Hausman test
 
quietly xtset momid
xtreg birwt smok male mage hsgrad somecoll collgrad married black kessner2 ///
  kessner3 novisit pretri2 pretri3 mn_smok mn_male mn_mage mn_kessner2  ///
  mn_kessner3 mn_novisit mn_pretri2 mn_pretri3, mle vce(robust)

testparm mn_*


*** 3.9 Assigning values to random effects: Residual diagnostics

quietly mixed birwt smoke male mage hsgrad somecoll collgrad married black ///
   kessner2 kessner3 novisit pretri2 pretri3, || momid:, reml 

predict lev2, reffects
predict comp_se, reses
generate diag_se = sqrt(exp(2*_b[lns1_1_1:_cons]) - comp_se^2)
replace lev2 = lev2/diag_se
predict lev1, rstandard

summarize gestat if lev1<-4

histogram lev1, normal xtitle(Standardized level-1 residuals)

histogram lev2 if idx==1, normal xtitle(Standardized level-2 residuals)


