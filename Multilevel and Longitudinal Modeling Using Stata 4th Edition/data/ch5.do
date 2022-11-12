*****************************************************************
****Chapter 5: Subject-specific effects and dynamic models *****
*****************************************************************


*** 5.1 Introduction

use https://www.stata-press.com/data/mlmus4/wagepan, clear

generate educt = educ - 12
generate yeart = year - 1980

xtset nr year


*** 5.2 Random effects approach: No endogeneity

xtreg lwage black hisp union married exper yeart educt, re vce(robust)


*** 5.3 Fixed effects: Level-2 endogeneity 

** 5.3.1 De-meaning and subject dummies

xtreg lwage union married exper, fe vce(robust)

*regress lwage ibn.nr union married exper, noconstant vce(robust)

** 5.3.2 Hausman test

generate exper_copy = exper
quietly xtreg lwage black hisp union married exper_copy yeart educt, re
estimates store ri
quietly xtreg lwage union married exper, fe 
estimates store fi
hausman fi ri

** 5.3.3 Mundlak

egen mn_union = mean(union), by(nr)
egen mn_married = mean(married), by(nr)

xtreg lwage black hisp union married exper yeart educt mn_union mn_married, ///
   re vce(robust)

testparm mn_union mn_married

** 5.3.4 First-differencing

* regress D.lwage D.union D.married, vce(robust)
regress D1.(lwage union married), vce(robust)


*** 5.4 Difference-in-difference

** 5.4.1 Does raising the minimum wage reduce employment?

use https://www.stata-press.com/data/mlmus4/minwage, clear

xtset store post

generate treated = state*post
xtreg fte post treated, fe vce(robust)

* xtdidregress (fte) (treated), group(store) time(post)

** 5.4.2 Repeated-measures ANOVA

anova fte state / store|state post state#post, dropemptycells

*display sqrt(e(F_4))


*** 5.5 Subject-specific coefficients

** 5.5.1 Random-coefficient model: No endogeneity

use https://www.stata-press.com/data/mlmus4/wagepan, clear

generate educt = educ - 12
generate yeart = year - 1980

mixed lwage black hisp union married exper yeart educt || nr: exper, ///
   covariance(unstructured) mle stddeviations vce(robust)
estimates store rc

quietly mixed lwage black hisp union married exper yeart educt || nr: , ///
   mle stddeviations vce(robust)
estimates store ri

lrtest ri rc, force

display (chi2tail(1,167.78)+chi2tail(2,167.78))/2

** 5.5.2. Fixed-coefficient model: Level-2 endogeneity

quietly xtset nr year 
xtreg D1.(lwage union married), fe vce(robust)

*regress D2.(lwage union married), vce(robust) noconstant


*** 5.6 Hausman-Taylor estimator: Level-2 endogeneity for level-1 and
* level-2 covariates

xthtaylor lwage black hisp union married exper educt, ///
    endog(union educt) vce(robust)

	
*** 5.7 Instrumental-variable methods: Level-1 and level-2 
* endogeneity

** 5.7.1 Do deterrents decrease crime rates?

use https://www.stata-press.com/data/mlmus4/crime, clear
xtset county year

xtreg lcrmrte lprbarr lprbconv lprbpris lavgsen lpolpc ldensity i.year, ///
   fe vce(robust)

** 5.7.3 Fixed-effects IV estimator
xtivreg lcrmrte lprbconv lprbpris lavgsen ldensity i.year ///
   (lprbarr lpolpc = ltaxpc lmix), fe vce(robust)

egen m_lprbarr = mean(lprbarr), by(county)
egen m_lpolpc = mean(lpolpc), by(county)
egen m_ltaxpc = mean(ltaxpc), by(county)
egen m_lmix = mean(lmix), by(county)
generate d_lprbarr = lprbarr - m_lprbarr
generate d_lpolpc = lpolpc - m_lpolpc 
generate d_ltaxpc = ltaxpc - m_ltaxpc
generate d_lmix = lmix - m_lmix

/*
foreach var of varlist lprbarr lpolpc ltaxpc lmix {
	egen m_`var' = mean(`var'), by(county)
	generate d_`var' = `var' - m_`var'
}
*/

correlate d_lprbarr d_lpolpc d_ltaxpc d_lmix

** 5.7.4 Random-effects IV estimator

xtivreg lcrmrte lprbconv lprbpris lavgsen ldensity lpctmin i.year ///
(lprbarr lpolpc = ltaxpc lmix), ec2sls vce(robust)

** 5.7.5 More Hausman tests

quietly xtivreg lcrmrte lprbconv lprbpris lavgsen ldensity  ///
   i.year (lprbarr lpolpc = ltaxpc lmix), fe 
quietly estimates store fe2sls

quietly xtivreg lcrmrte lprbconv lprbpris lavgsen ldensity lpctmin  ///
  i.year (lprbarr lpolpc = ltaxpc lmix), ec2sls 
quietly estimates store ec2sls

hausman fe2sls ec2sls


*** 5.8 Dynamic models

** 5.8.1 Dynamic model without subject-specific intercepts

use https://www.stata-press.com/data/mlmus4/wagepan, clear

generate educt = educ - 12
generate yeart = year - 1980

sort nr yeart
quietly xtset nr year

list nr yeart lwage L.lwage in 1/12, noobs sepby(nr)

regress lwage L.lwage black hisp union married exper yeart educt, vce(robust)


** 5.8.2 Dynamic models with subject-specific intercepts

quietly xtset nr year
ivregress 2sls D1.lwage D1.(union married) (L1D1.lwage = L2.lwage), vce(robust)

xtabond lwage union married exper, lags(1) twostep noconstant vce(robust)

estat abond

quietly xtabond lwage union married exper, lags(1) twostep noconstant 

estat sargan


*** 5.9 Missing data and dropout

** 5.9.1 Maximum likelihood estimation under MAR: A simulation

clear
set obs 1000000  /* 1 million */
set seed 123123123
generate zeta = rnormal(0,1)
generate y1 = 2 + zeta + rnormal(0,1)
generate y2 = 2 + 0 + zeta + rnormal(0,1)

* missingness indicator for time 2
generate missing_y2 =  y1>2 & runiform()<.9 

* reshape data
generate id = _n
reshape long y, i(id) j(occasion)

* complete data
regress y i.occasion, vce(cluster id)

* incomplete data
drop if missing_y2==1&occasion==2

regress y i.occasion, vce(cluster id)

quietly xtset id
xtreg y i.occasion, mle





