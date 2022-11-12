***********************************************
****Chapter 14: Discrete-time survival data****
***********************************************


*** 14.2 Single-level models for discrete-time survival data

** 14.2.1 Discrete-time hazard and discrete-time survival

use https://www.stata-press.com/data/mlmus4/promotion, clear

ltable dur event, noadjust
ltable dur event, hazard noadjust


** 14.2.2 Data expansion for disrete-time survival analysis
 
list id dur event if id<3, noobs

expand dur
by id, sort: generate year = _n
list id dur year event if id<3, sepby(id) noobs
generate y = 0
replace y = event if year==dur


** 14.2.3 Estimation via regression models for dichotomous responses

tabulate year y, row

* Estimation using logit
logit y i.year
predict haz, pr

list id year haz y event if id<3, sepby(id) noobs

twoway (line haz year if id==1, connect(stairstep)), legend(off) ///
 xtitle(Year) ytitle(Discrete-time hazard)

** 14.2.4 Including time-constant covariates
* Estimation using logit
logit y i.year undgrad i.phdmed phdprest

predict lo, xb

twoway (line lo year if id==1, connect(stairstep) lpatt(solid)) ///
       (line lo year if id==4, connect(stairstep) lpatt(dash)), ///
       legend(off) xtitle(Year) ytitle(Log odds)

generate ln_one_m_haz = ln(1-invlogit(lo))
by id (year), sort: generate ln_surv = sum(ln_one_m_haz)
generate surv = exp(ln_surv)

twoway (line surv year if id==1, connect(stairstep) lpatt(solid)) ///
       (line surv year if id==4, connect(stairstep) lpatt(dash)), ///
       legend(off) xtitle(Year) ytitle(Survival)

** 14.2.5 Including time-varying covariates

use https://www.stata-press.com/data/mlmus4/promotion, clear

reshape long art cit, i(id) j(year)
list id year art cit dur event if id<3, sepby(id) noobs

drop if year>dur
generate y = 0
replace y = event if year==dur
generate prestige = prest1
replace prestige = prest2 if year>=jobtime

* Estimation using logit
logit y i.year undgrad i.phdmed phdprest art cit prestige, or


** 14.2.5 Multiple absorbing events and competing risks

use https://www.stata-press.com/data/mlmus4/promotion, clear

generate event1 = 0
replace event1 = 1 if dur<jobtime & event==1

generate event2 = 0
replace event2 = 1 if dur>= jobtime

replace dur = jobtime if event2==1

reshape long art cit, i(id) j(year)

drop if year>dur

list id year jobtime dur event1 event2 if id==1|id==2|id==161, sepby(id) noobs

generate yy = 0
replace yy = 1 if year==dur & event1==1
replace yy = 2 if year==dur & event2==1

label define eventlabel 1 "promfirst" 2 "leave"
label values yy eventlabel
list id year jobtime dur event1 event2 yy if id==1|id==2|id==161, ///
 sepby(id) noobs

* Estimation using mlogit
mlogit yy ib6.year undgrad i.phdmed phdprest art cit prest1, rrr



*** 14.4 Data expansion

use https://www.stata-press.com/data/mlmus4/mortality, clear

egen discrete = cut(time), at(0 1 6 12 24 61) icodes
table discrete, statistic(min time) statistic(max time) nototal

replace discrete = discrete + 1

ltable discrete death, noadjust 

expand discrete
by kidid, sort: generate interval = _n
generate y = 0
replace y = death if interval == discrete
tabulate interval y, row 


*** 14.6 Complementary log-log models

** 14.6.1 Marginal baseline hazard
* Estimation using cloglog
cloglog y ibn.interval, noconstant 

predict haz, pr

sort kidid interval
list kidid interval haz if kidid==101, noobs
 
** 14.6.2 Including covariates 
* Estimation using cloglog
cloglog y ibn.interval mage c.mage#c.mage border i.prevbirth i.pdead ///
   1.nextbirth#4.interval 1.nextbirth#5.interval 2.nextbirth#5.interval, ///
   noconstant eform vce(cluster momid)


*** 14.7 Random-intercept complementary log-log model

** 14.7.1 Model specification
* Estimation using mecloglog
mecloglog y ibn.interval mage c.mage#c.mage border i.prevbirth i.pdead ///
   1.nextbirth#4.interval 1.nextbirth#5.interval 2.nextbirth#5.interval, ///
   noconstant || momid:, eform

*display _b[/var(_cons[momid])]/(_b[/var(_cons[momid])]+_pi^2/6)
estat icc

display exp(sqrt(2*_b[/var(_cons[momid])])*invnormal(3/4))



*** 14.8 Population-averaged or marginal vs. subject-specific or conditional ...

quietly tabulate interval, generate(int)
generate mage2 = mage^2
generate comp12 = 1.nextbirth*4.interval
generate comp24e = 1.nextbirth*5.interval
generate comp24l = 2.nextbirth*5.interval

gllamm y int1-int5 mage mage2 border p0014 p1523 p2435 p36up pdead ///
 comp12 comp24e comp24l, noconstant eform i(momid) link(cll) family(binomial) ///
 adapt

keep if kidid==101
replace pdead = 0
replace p2435 = 0
save junk, replace

replace kidid = 1
append using junk
replace kidid = 2 if kidid==101
append using junk
replace kidid = 3 if kidid==101
append using junk
replace kidid = 4 if kidid==101
append using junk
replace kidid = 5 if kidid==101

drop if interval>kidid

sort kidid interval
list kidid interval y, sepby(kidid) noobs

replace momid = kidid

gllapred loglik, ll fsample
generate msurv = exp(loglik)

generate zeta1 = 0
gllapred chaz, mu us(zeta) fsample

generate ln1mchaz = ln(1-chaz)
by kidid (interval), sort: generate sln1mchaz = sum(ln1mchaz)
generate csurv = exp(sln1mchaz)

twoway (line msurv interval if kidid==interval, sort lpatt(solid))  ///
       (line csurv interval if kidid==5, sort lpatt(dash)),         ///
       legend(order(1 "marginal" 2 "conditional"))                  ///
       ytitle(Survival probability)

