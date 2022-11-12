**********************************************
****Chapter 15: Continuous-time survival *****
**********************************************


*** 15.2 What makes marriages fail?

use https://www.stata-press.com/data/mlmus4/divorce, clear

generate mixrace = (heblack & !sheblack) | (sheblack & !heblack)
generate hedropout = hiseduc < 12
generate hecollege = hiseduc >= 16 if !missing(hiseduc)
generate heolder = agediff > 10 if !missing(agediff)
generate sheolder = agediff < -10


*** 15.3 Hazards and survival

list id dur divorce in 1/4, clean noobs

stset dur, failure(divorce) id(id)

sts graph

sts graph, hazard


*** 15.4 Proportional hazards models

** 15.4.1. Piecewise exponential model

stptime, at(0,1,5,9,15,25)

stsplit interval, at(1,5,9,15,25)
sort id _t
list id interval _t0 _t _d if id<34, sepby(id) noobs
 
generate pers_yr = _t-_t0
table interval, statistic(frequency) statistic(total pers_yr) ///
 statistic(total _d) nototal

* Estimation using streg
streg ibn.interval, distribution(exponential) noconstant

streg ibn.interval i.heblack i.mixrace i.hedropout i.hecollege i.heolder ///
   i.sheolder, distribution(exponential) noconstant

generate h0 = exp(_b[0.interval]*(0.interval) + _b[1.interval]*(1.interval)  ///
   + _b[5.interval]*(5.interval) + _b[9.interval]*(9.interval) ///
   + _b[15.interval]*(15.interval) + _b[25.interval]*(25.interval))
 
twoway line h0 _t, connect(J) sort xtitle(Time since wedding in years) ///
   ytitle(Baseline hazard) 


* Estimation using poisson
generate lexposure = ln(_t - _t0)
poisson _d ibn.interval i.heblack i.mixrace i.hedropout i.hecollege i.heolder ///
   i.sheolder, offset(lexposure) noconstant irr
   
drop interval pers_yr lexposure h0
stjoin


** 15.4.2 Cox regression model

* Estimation using stcox
stcox i.heblack i.mixrace i.hedropout i.hecollege i.heolder i.sheolder

stcurve, hazard xtitle(Time since wedding in years) at(sheolder=(0 1))

stcurve, hazard xtitle(Time since wedding in years) ///
   at(heblack=0 mixrace=0 hedropout=0 hecollege=0 heolder=0 sheolder=0) ///
   outfile(temp, replace)

** 15.4.3 Cox regression via Poisson regression for expanded data

stsplit, at(failures) riskset(click)
list id click _t0 _t _d _st if id==9 & (click<5 | click>624), separator(4) ///
 noobs

* Estimation using xtpoisson, fe
generate lexposure = ln(_t - _t0)
quietly xtset click
xtpoisson _d i.heblack i.mixrace i.hedropout i.hecollege i.heolder ///
   i.sheolder, fe offset(lexposure) irr

*clogit d i.heblack i.mixrace i.hedropout i.hecollege i.heolder i.sheolder, group (click)

** 15.4.4 Approximate Cox regression: Poisson regression with smooth baseline hazard

mkspline sp = _t, cubic knots(1,5,7,9,12,15,20,25)

* Estimation using poisson   
poisson _d sp* i.heblack i.mixrace i.hedropout i.hecollege i.heolder ///
   i.sheolder, offset(lexposure) irr

matrix a=e(b)
matrix coeff = a[1,1..7],a[1,20..20]
matrix score lhazard = coeff
generate h0=exp(lhazard)

append using temp

twoway (line h0 _t if !missing(divorce), sort)  (line haz1 _t, sort), ///
   xtitle(Time since wedding in years) ///
   ytitle(Baseline hazard) legend(order(1 "Cubic spline" 2 "Cox"))

drop if missing(_st)
drop sp* h0 haz1 lexposure lhazard click
stjoin


*** 15.5 Accelerated failure time models

** 15.5.1 Log-normal model

* Estimation using streg
streg i.heblack i.mixrace i.hedropout i.hecollege i.heolder i.sheolder,  ///
    distribution(lognormal) time
streg, tr

stcurve, hazard xtitle(Time since wedding in years) at(sheolder=(0 1))

* Estimation using stintreg

generate dur2 = dur if divorce==1
stintreg i.heblack i.mixrace i.hedropout i.hecollege i.heolder i.sheolder, ///
	interval(dur dur2) distribution(lognormal)


*** 15.6 Time-varying covariates

use https://www.stata-press.com/data/mlmus4/divorce2, clear

list id dur numkids divorce if id<34, sepby(id) noobs

stset dur, failure(divorce) id(id)

sort id dur
list id _t0 _t numkids divorce if id<34, sepby(id) noobs

generate mixrace = (heblack & !sheblack) | (sheblack & !heblack)
generate hedropout = hiseduc < 12
generate hecollege = hiseduc >= 16 if !missing(hiseduc)
generate heolder = agediff > 10 if !missing(agediff)
generate sheolder = agediff < -10

* Estimation using streg
streg i.heblack i.mixrace i.hedropout i.hecollege i.heolder i.sheolder  ///
    numkids, distribution(lognormal) time tr


	
*** 15.7 Does nitrate reduce the hazard of angina pectoris?

use https://www.stata-press.com/data/mlmus4/angina8, clear
list subj regime occ second uncen in 1/16, sepby(subj) noobs

generate treat = regime==1 & occ>1
generate id = _n
list id subj regime occ treat second uncen in 1/16, sepby(subj) noobs

stset second, failure(uncen) id(id)


*** 15.8 Marginal modeling

** 15.8.1 Cox regression with occasion-specific dummy variables

* Estimation using stcox
stcox i.occ i.treat, vce(cluster subj)

stcurve, hazard at(occ=2 treat=(1 0)) ///
  legend(order(1 "Treatment" 2 "Placebo")) ///
  xtitle(Time in seconds) ytitle(Hazard function) 

stcurve, hazard name(cox, replace) at(occ=1 treat=0) outfile(temp, replace)

** 15.8.2 Cox regression with occasion-specific baseline hazards
* Estimation using stcox
stcox i.treat, strata(occ) vce(cluster subj)


** 15.8.3 Approximate Cox regression

stsplit, at(failures) riskset(click)
 
generate lny = ln(_t - _t0)

* xtset click
* xtpoisson _d i.occ i.treat, fe offset(lny) irr

orthpoly _t, gen(t1-t4) degree(4)

* Estimation using poisson
poisson _d t1-t4 i.occ i.treat, offset(lny) irr vce(cluster subj)

matrix a = e(b)
matrix coeff = a[1,1..4], a[1,11..11]
matrix score lhazard = coeff
generate haz0 = exp(lhazard)

append using temp
twoway (line haz0 _t, sort) (line haz1 _t, sort), ///
   xtitle(Time in seconds) ytitle(Baseline hazard function) ///
   legend(order(1 "Polynomial" 2 "Cox"))


*** 15.9 Multilevel proportional hazards models

** 15.9.1 Cox regression with gamma shared frailty

drop click lny lhazard haz0 haz1 t1-t4
stjoin

* Estimation using stcox, shared
stcox i.occ i.treat, frailty(gamma) shared(subj)

stcurve, at(occ=2 treat=(1 0)) ///
   hazard xtitle(Time in seconds) ytitle(Hazard function) ///
   legend(order(1 "Treatment" 2 "Placebo")) range(87 430)

*stsplit, at(failures) riskset(click)
*generate lny = ln(_t - _t0)
*xtset subj
*xtpoisson _d i.click i.occ i.treat, re offset(lny) irr 


** 15.9.2 Approximate Cox regression with log-normal shared frailty

* Estimation using mepoisson
stsplit, at(failures) riskset(click)
generate lny = ln(_t - _t0)
orthpoly _t, gen(t1-t4) degree(4)

mepoisson _d t1-t4 i.occ i.treat, offset(lny) || subj:, irr

display exp(sqrt(2*_b[/var(_cons[subj])])*invnormal(3/4))

*mestreg t1-t4 i.occ i.treat || subj:, distribution(exponential)


** 15.9.3 Approximate Cox regression with normal random intercept and random coefficient

* Estimation using mepoisson
mepoisson _d t1-t4 i.occ i.treat, offset(lny) || subj: treat, ///
   covariance(unstructured) irr

/*
mestreg t1-t4 i.occ i.treat || subj: treat, ///
covariance(unstructured) distribution(exponential)
*/

drop lny click t1-t4
stjoin


*** 15.10 Multilevel accelerated failure time models

** 15.10.1 Lognormal model with gamma shared frailty

* Estimation using streg
streg i.occ i.treat, distribution(lognormal) shared(subj) frailty(gamma)

stcurve, at(occ=2 treat=(1 0)) hazard ///
   xtitle(Time in seconds) ytitle(Hazard function) ///
   legend(order(1 "Treatment" 2 "Placebo")) 

/*
stcurve, at(occ=2 treat=(1 0)) hazard ///
xtitle(Time in seconds) ytitle(Hazard function) ///
legend(order(1 "Treatment" 2 "Placebo")) unconditional
*/
   
** 15.10.2 Lognormal model with lognormal shared frailty

* Estimation using mestreg
mestreg i.occ i.treat || subj:, distribution(lognormal)

/*
generate lnsecond = ln(second)
generate dep_lo = lnsecond
generate dep_hi = lnsecond if uncen == 1
replace dep_hi = . if uncen == 0
meintreg dep_lo dep_hi i.occ i.treat || subj:
*/

** 15.10.3 Lognormal model with random intercept and random coefficient

* Estimation using mestreg
mestreg i.occ i.treat || subj: i.treat, covariance(unstructured) ///
   distribution(lognormal)


*** 15.11 A fixed-effects approach

* Estimation using stcox, strata
stcox i.occ i.treat, strata(subj)


*** 15.12 Different approaches to recurrent event data

use https://www.stata-press.com/data/mlmus4/recurrent, clear

list, sepby(id) noobs
egen idnum = group(id number)

** 15.12.1 Total time

stset stop, failure(event) id(idnum)
list id number stop event _t0 _t _d if _st==1, sepby(id) noobs

fillin id number
by id (number), sort: replace stop = stop[_n-1] if _fillin==1
replace event = 0 if _fillin==1

drop idnum
egen idnum = group(id number)

stset stop, failure(event) id(idnum)
list id number stop event _t0 _t _d if _st==1, sepby(id) noobs

*** 15.12.2 Counting process

by id (number), sort: generate start = stop[_n-1] if _n>1
replace start = 0 if start==.

stset stop, failure(event) enter(start) id(idnum)

list id number stop event _t0 _t _d if _st==1, sepby(id) noobs

*** 15.12.3 Gap time

stset stop, failure(event) origin(start) id(idnum)
list id number stop event _t0 _t _d if _st==1, sepby(id) noobs

