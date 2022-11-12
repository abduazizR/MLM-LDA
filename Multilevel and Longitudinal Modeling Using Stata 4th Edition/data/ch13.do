**************************
****Chapter 13: Counts****
**************************


*** 13.4 Did the German health-care reform reduce the number of doctor visits?

use https://www.stata-press.com/data/mlmus4/drvisits, clear


*** 13.5 Longitudinal data structure

quietly xtset id reform
xtdescribe if !missing(numvisit)


*** 13.6 Single-level Poisson regression

** 13.6.1 Model specification
* Estimation using poisson
poisson numvisit reform age educ married badh loginc summer, irr vce(cluster id)
estimates store ordinary 
* Estimation using glm
glm numvisit reform age educ married badh loginc summer, ///
 family(poisson) link(log) vce(cluster id) eform


*** 13.7 Random-intercept Poisson regression

** 13.7.3 Estimation 
* Using xtpoisson
quietly xtset id
xtpoisson numvisit reform age educ married badh loginc summer, normal irr 

display exp(sqrt(2)*.90407*invnormal(3/4))

* xtpoisson, coeflegend
display exp(sqrt(2*exp(_b[/lnsig2u]))*invnormal(3/4))

* Using mepoisson
mepoisson numvisit reform age educ married badh loginc summer || id:, irr ///
 intpoints(12)
estimates store meri

* Using gllamm
generate cons = 1
eq ri: cons
gllamm numvisit reform age educ married badh loginc summer, ///
 family(poisson) link(log) i(id) eqs(ri) eform adapt nip(12)
estimates store glri

display sqrt(.81735621)

* gllamm, robust eform

lrtest ordinary glri, force


*** 13.8 Random-coefficient Poisson regression 

** 13.8.1 Model specification
* Estimation using mepoisson
mepoisson numvisit reform age educ married badh loginc summer || ///
  id: reform, covariance(unstructured) irr
estimates store merc

lrtest meri merc

* Estimation using gllamm
estimates restore glri
matrix a = e(b)
eq rc: reform
gllamm numvisit reform age educ married badh loginc summer, ///
  family(poisson) link(log) i(id) nrf(2) eqs(ri rc) from(a) eform adapt 
estimates store glrc

* gllamm, robust eform


*** 13.9 Overdispersion in single-level models

** 13.9.1 Normally distributed random intercept 
* Estimation using xtpoisson
generate obs = _n
quietly xtset obs
xtpoisson numvisit reform age educ married badh loginc summer, normal irr 
   
*** 13.9.2 Negative binomial models
* Mean dispersion or NB2
nbreg numvisit reform age educ married badh loginc summer, ///
   irr dispersion(mean)
* Constant dispersion or NB1
nbreg numvisit reform age educ married badh loginc summer, ///
   irr dispersion(constant)

*** 13.9.3 Quasilikelihood 
* Estimation using glm 
glm numvisit reform age educ married badh loginc summer, ///
 irls family(poisson) link(log) eform scale(x2)


*** 13.10 Level-1 overdispersion in two-level models

** 13.10.1 Random-intercept model with robust standard errors
* Estimation using mepoisson
mepoisson numvisit reform age educ married badh loginc summer || id:, ///
 irr vce(robust)
 
** 13.10.3 Negative binomial models with random intercept 
* Estimation using menbreg
menbreg numvisit reform age educ married badh loginc summer || id:, ///
 irr dispersion(mean)

** 13.10.4 The HHG model 
*xtset id
*xtnbreg numvisit reform age educ married badh loginc summer, irr


*** 13.11 Other approaches to two-level count data

** 13.11.1 Conditional Poisson regression
* Estimation using xtpoisson, fe
quietly xtset id
xtpoisson numvisit reform educ married badh loginc summer, fe irr
  
*set matsize 2000

poisson numvisit reform educ married badh loginc summer i.id, irr

** 13.11.2 Conditional negative binomial regression
*xtnbreg numvisit reform age educ married badh loginc summer, fe irr

** 13.11.3 Generalized estimating equations

quietly xtset id
xtgee numvisit reform age educ married badh loginc summer, ///
  family(poisson) link(log) vce(robust) eform
  
estat wcorrelation
 
** 13.12 Estimating marginal and conditional effects when responses are missing at random

egen num = count(numvisit), by(id)
generate mult = 3 - num
expand mult
by id (reform), sort: generate occ = _n - 1
generate missing = (occ!=reform & num==1)
replace numvisit = . if missing==1
replace reform = occ if missing==1

sort id reform
set seed 1211

estimates restore glri
gllasim y, fsample

by id (reform), sort: generate drop = (y[1]>2 & runiform()<.9)
replace y = . if drop==1 & reform==1

quietly xtset id
xtpoisson y reform age educ married badh loginc summer, normal irr 
 
poisson y reform age educ married badh loginc summer, irr

quietly xtset id
xtgee y reform age educ married badh loginc summer, family(poisson) ///
 link(log) corr(exchangeable) vce(robust) eform



*** 13.13 Which Scottish counties have a high risk of lip cancer?

use https://www.stata-press.com/data/mlmus4/lips, clear


*** 13.14 Standardized mortality ratios

generate smr = 100*o/e


*** 13.15 Random-intercept Poisson regression

** 13.15.2 Estimation using gllamm

generate lne = ln(e)
gllamm o, i(county) offset(lne) family(poisson) link(log) adapt

** 13.15.3 Prediction of SMRs

gllapred mu, mu nooffset
generate thet = 100*mu
sort county
list county thet in 1/10, clean noobs

*copy https://www.stata-press.com/data/mlmus4/scotmaps.do scotmaps.do

twoway (scatter thet smr, msymbol(none) mlabpos(0) mlabel(county)) ///
 (function y=x, range(0 600)), xline(108) yline(108)               ///
 xtitle(Crude SMR) ytitle(Empirical Bayes SMR) legend(off)


*** 13.16 Nonparametric maximum likelihood estimation

** 13.16.1 Specification 
* Estimation using gllamm
gllamm o, i(county) offset(lne) family(poisson) link(log) ip(f) nip(2) 

matrix a = e(b)
local ll = e(ll)
local k = e(k)

gllamm o, i(county) offset(lne) family(poisson) link(log) ip(f) nip(3) ///
  gateaux(-5 5 100) from(a) lf0(`k' `ll')

matrix a = e(b)
local ll = e(ll)
local k = e(k)

gllamm o, i(county) offset(lne) family(poisson) link(log) ip(f) nip(4) ///
  gateaux(-5 5 100) from(a) lf0(`k' `ll')
  
matrix a = e(b)
local ll = e(ll)
local k = e(k)

gllamm o, i(county) offset(lne) family(poisson) link(log) ip(f) nip(5) ///
  gateaux(-5 5 100) from(a) lf0(`k' `ll')
  
matrix locs = e(zlc2)'
matrix lp = e(zps2)'
svmat locs
svmat lp

generate p = exp(lp1)
generate smrloc = 100*exp(_b[_cons] + locs1)
twoway (dropline p smrloc), xtitle(Location) ytitle(Probability)

** 13.16.2 Prediction
gllapred mu2, mu nooffset
generate thet2 = 100*mu2
sort county
list county thet thet2 in 1/10, clean noobs

