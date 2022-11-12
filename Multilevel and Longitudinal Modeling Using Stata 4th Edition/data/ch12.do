*************************************************************
****Chapter 12: Nominal responses and discrete choice  ******
*************************************************************


*** 12.2 Single-level models for nominal responses


** 12.2.1 Multinomial logit models

use https://www.stata-press.com/data/mlmus4/travel1, clear

list traveler alt Hinc Psize in 1/10, noobs clean

* Estimation using mlogit
mlogit alt Hinc Psize, baseoutcome(4)

label define modelabels 1 "Air" 2 "Train" 3 "Bus" 4 "Car"
label values alt modelabels

mlogit alt Hinc Psize, baseoutcome(4) rrr


** 12.2.2 Conditional logit models with alternative-specific covariates

use https://www.stata-press.com/data/mlmus4/travel2, clear

list traveler alt Choice GC Ttime in 1/8, sepby(traveler) noobs

* Estimation using clogit
clogit Choice GC Ttime, group(traveler) 
clogit Choice GC Ttime, group(traveler) or

* Estimation using cmclogit
cmset traveler alt
cmclogit Choice GC Ttime, noconstant or


** 12.2.2 Conditional logit models with alternative- and unit-specific covariates

* Estimation using clogit
label define modelabels 1 "Air" 2 "Train" 3 "Bus" 4 "Car"
label values alt modelabels

clogit Choice ib4.alt##(c.Hinc c.Psize) GC Ttime, group(traveler)

* Estimation using cmclogit
quietly cmset traveler alt
cmclogit Choice GC Ttime, casevars(Hinc Psize) basealternative(4)

*asclogit Choice GC Ttime, case(traveler) casevar(Hinc Psize) ///
*  alternatives(alt) basealternative(4)


*** 12.5 Does marketing affect choice of yogurt?

use https://www.stata-press.com/data/mlmus4/yogurt, clear

drop if brand==3
egen sum = total(choice), by(house occ)
drop if sum==0

recode brand(4=3)
label define b 1 "Yoplait" 2 "Dannon" 3 "WeightW", modify
label values brand b

egen set = group(house occ)

generate pricec = price*100
drop price

sort house occ
list if house==24, sepby(occ) noobs 


*** 12.6 Single-level conditional logit models

** 12.6.1 Conditional logit model with alternative-specific intercepts

* Estimation using clogit

clogit choice ib2.brand i.feature pricec, group(set) vce(cluster house)

* Estimation using cmclogit

cmset set brand
cmclogit choice i.feature pricec, basealternative(2) vce(cluster house)
* cmclogit choice ib2.brand i.feature pricec, noconstant vce(cluster house)


*** 12.7 Multilevel conditional logit models

** 12.7.1 Preference heterogeneity: Brand-specific random intercepts

* Estimation using cmxtmixlogit

cmset house occ brand
cmxtmixlogit choice i.feature pricec, noconstant random(ib2.brand, correlated)

cmxtmixlogit choice i.feature pricec, noconstant random(ib2.brand, correlated) ///
intpoints(3000)
cmxtmixlogit, or

estimates store rint
quietly cmclogit choice i.feature pricec
estimates store fixed
lrtest rint fixed, force

* Estimation using gllamm

tabulate brand, generate(br)
rename br1 Yoplait
rename br2 Dannon
rename br3 WeightW

gllamm brand Yoplait WeightW feature pricec, i(house) link(mlogit) ///
   expanded(set choice o) noconstant cluster(house) robust init
gllamm, eform

eq Yo: Yoplait
eq We: WeightW
gllamm brand Yoplait WeightW feature pricec, i(house) link(mlogit) /// 
 expanded(set choice o) noconstant nrf(2) eqs(Yo We) adapt 

gllamm brand Yoplait WeightW feature pricec, i(house) link(mlogit) ///
 expanded(set choice o) noconstant nrf(2) eqs(Yo We) adapt nip(16)

gllamm, eform


** 12.7.2 Response heterogeneity: Marketing variables with random coefficients

* Estimation using cmxtmixlogit
quietly cmset house occ brand
cmxtmixlogit choice, basealternative(2) random(i.feature pricec, correlated) ///
   intpoints(1000)
cmxtmixlogit, or

estimates store rcoeff
lrtest rcoeff fixed, force

* interpretation
display exp( .7866819 - 1.96*1.786191)
display exp( .7866819 + 1.96*1.786191)

display exp(-.5286953 - 1.96*.7113092)
display exp(-.5286953 + 1.96*.7113092)


* Estimation using gllamm
eq fea: feature
eq pr: pricec
gllamm brand Yoplait WeightW feature pricec,          ///
   i(house) link(mlogit) exp(set choice o) noconstant nrf(2) eqs(fea pr) ///
   adapt nip(12)
estimates store gll_rc
gllamm, eform



** 12.7.3 Preference and response heterogeneity

* Estimation using cmxtmixlogit
set processors 1
quietly cmset house occ brand
cmxtmixlogit choice, noconstant random(ib2.brand i.feature pricec, correlated)  ///
    intpoints(10000)
set processors `c(processors_lic)'

* Estimation using gllamm
eq Yo: Yoplait
eq We: WeightW
eq fea: feature
eq pr: pricec

gllamm brand Yoplait WeightW feature pricec, i(house) link(mlogit) ///
   expanded(set choice o) noconstant eqs(Yo We fea pr) nrf(4) skip adapt 


*** Prediction of marginal choice probabilities
   
estimates restore rcoeff
margins feature, alternative(1) noesample nose

* using gllamm
estimates restore gll_rc

gen feature_save=feature
replace feature=0 if brand==1
gllapred prob0, mu marg fsample
replace feature=1 if brand==1
gllapred prob1, mu marg fsample

table brand, statistic(mean prob0) statistic(mean prob1) nototal

replace feature=feature_save

  
** 12.9 Prediction of random effects and household-specific response probabilities 

estimates restore gll_rc

gllapred eb, u 
generate eb_pricec = ebm1 + _b[pricec]
generate eb_feature = ebm2 + _b[feature]

egen pickone = tag(house)
twoway scatter eb_pricec eb_feature if pickone==1, mlabel(house) ///
   xtitle(EB prediction of feature coefficient) ///
   ytitle(EB prediction of price coefficient) 

save temp, replace

drop _all
set obs 40
generate pricec =  _n*.5
generate set = 3000 + _n
expand 3
by set, sort: generate brand=_n
generate yopprice = pricec
replace pricec = 8 if brand>1
generate feature = 0
tabulate brand, generate(br)
rename br1 Yoplait
rename br2 Dannon
rename br3 WeightW
generate choice = Dannon
expand 6
by set brand, sort: generate house=_n
recode house 1=33 2=15 3=7 4=48 5=6 6=55
replace brand = .
append using temp

generate preddata = brand ==.

gllapred probs, mu fsample

egen select = anymatch(house), values(33 15 7 48 6 55)
twoway (line probs yopprice if Yoplait==1, sort) ///
  (line probs yopprice if WeightW==1, sort lpatt(longdash)) ///
  (line probs yopprice if Dannon==1, sort lpatt(shortdash)) ///
  if feature==0&select==1&preddata==1, ///
  by(house) legend(order(1 "Yoplait" 2 "WeightW" 3 "Dannon") rows(1)) ///
  xtitle(Price of Yoplait in cents/oz) ytitle(Probability of purchase)

