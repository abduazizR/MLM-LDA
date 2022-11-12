**************************************
****Chapter 11: Ordinal responses*****
**************************************


*** 11.3 Are antpsychotic drugs effective for patients with schizophrenia?

use https://www.stata-press.com/data/mlmus4/schiz, clear

generate impso = imps
recode impso -9=. 1/2.4=1 2.5/4.4=2 4.5/5.4=3 5.5/7=4


*** 11.4 Longitudinal data structure and graphs

** 11.4.1 Longitudinal data stucture

xtset id week
xtdescribe if !missing(impso)
table week treatment, statistic(count impso) nototal

** 11.4.2 Plotting cumulative proportions

egen propg1 = mean(impso>1), by(week treatment)
egen propg2 = mean(impso>2), by(week treatment)
egen propg3 = mean(impso>3), by(week treatment)

egen nonrare = anymatch(week), values(0,1,3,6)

label define t 0 "Control" 1 "Treatment", modify
label values treatment t
sort treatment id week
twoway (line propg1 week, sort)                                     ///
  (line propg2 week, sort lpatt(vshortdash))                        ///
  (line propg3 week, sort lpatt(dash)) if nonrare==1, by(treatment) ///
  legend(order(1 "Prop(y>1)" 2 "Prop(y>2)" 3 "Prop(y>3)") rows(1))  ///
  xtitle("Week")


** 11.4.3 Plotting cumulative logits and transforming the time scale

generate logodds1 = ln(propg1/(1-propg1))
generate logodds2 = ln(propg2/(1-propg2))
generate logodds3 = ln(propg3/(1-propg3))

twoway (line logodds1 week, sort)                                      ///
  (line logodds2 week, sort lpatt(vshortdash))                         ///
  (line logodds3 week, sort lpatt(dash)) if nonrare==1, by(treatment)  ///
  legend(order(1 "Log Odds(y>1)" 2 "Log Odds(y>2)" 3 "Log Odds(y>3)")  ///
  rows(1)) xtitle("Week")

generate weeksqrt = sqrt(week)
twoway (line logodds1 weeksqrt, sort)                                      ///
  (line logodds2 weeksqrt, sort lpatt(vshortdash))                         ///
  (line logodds3 weeksqrt, sort lpatt(dash)) if nonrare==1, by(treatment)  ///
  legend(order(1 "Log Odds(y>1)" 2 "Log Odds(y>2)" 3 "Log Odds(y>3)")      ///
  rows(1)) xtitle("Square root of week")


*** 11.5 Single-level proportional-odds model

** 11.5.2 Estimation using Stata

ologit impso c.weeksqrt##i.treatment, vce(cluster id) or

lincom c.weeksqrt+1.treatment#c.weeksqrt, or

predict prob1-prob4, pr
generate probg3 = prob4
generate probg2 = prob3 + probg3
generate probg1 = prob2 + probg2

twoway (line propg1 week if nonrare==1, sort lpatt(solid))          ///
  (line propg2 week if nonrare==1, sort lpatt(vshortdash))          ///
  (line propg3 week if nonrare==1, sort lpatt(dash))                ///
  (line probg1 week, sort lpatt(solid))                             ///
  (line probg2 week, sort lpatt(vshortdash))                        ///
  (line probg3 week, sort lpatt(dash)), by(treatment)               ///
  legend(order(1 "Prob(y>1)" 2 "Prob(y>2)" 3 "Prob(y>3)") rows(1))  ///
  xtitle("Week")


*** 11.6 A random-intercept proportional odds model

** 11.6.1 Model specification
* Estimation using meologit
meologit impso c.weeksqrt##i.treatment || id:, or 
estimates store model1
* Estimation using gllamm
generate interact = weeksqrt*treatment
gllamm impso weeksqrt treatment interact, i(id) link(ologit) adapt eform 
estimates store model1_gllamm

** 11.6.2 Measures of dependence and heterogeneity
* Residual intraclass correlation of latent responses
estimates restore model1
estat icc 

* Median odds ratio
display exp(sqrt(2*_b[/var(_cons[id])])*invnormal(3/4))

*** 11.7 Random-coefficient proportional odds model

** 11.7.1 Model specification
* Estimation using meologit
meologit impso c.weeksqrt##i.treatment || id: weeksqrt, /// 
 or covariance(unstructured)
estimates store model2

estat sd

lrtest model1 model2
display .5*chi2tail(1,r(chi2))+.5*chi2tail(2,r(chi2))

* Estimation using gllamm
generate cons = 1
eq inter: cons
eq slope: weeksqrt
gllamm impso weeksqrt treatment interact, i(id) nrf(2) eqs(inter slope) ///
  link(ologit) adapt eform
estimates store model2_gllamm


*** 11.8 Different kinds of predicted probabilities

** 11.8.1 Predicted population-averaged or marginal probabilities
estimates restore model2
predict pr1-pr4, pr marginal

twoway (line pr1 week, sort lpatt(solid))                                ///
  (line pr2 week, sort lpatt(vshortdash))                                ///
  (line pr3 week, sort lpatt(dash_dot))                                  ///
  (line pr4 week, sort lpatt(dash)), by(treatment)                       ///
  legend(order(1 "Prob(y=1)" 2 "Prob(y=2)" 3 "Prob(y=3)" 4 "Prob(y=4)")) /// 
  xtitle("Week")

generate mprobg1 = 1-pr1
generate mprobg2 = mprobg1-pr2
generate mprobg3 = pr4

twoway (line propg1 week if nonrare==1, sort lpatt(solid))        ///
 (line propg2 week if nonrare==1, sort lpatt(vshortdash))         ///
 (line propg3 week if nonrare==1, sort lpatt(dash))               ///
 (line mprobg1 week, sort lpatt(solid))                           ///
 (line mprobg2 week, sort lpatt(vshortdash))                      ///
 (line mprobg3 week, sort lpatt(dash)), by(treatment)             ///
 legend(order(1 "Prob(y>1)" 2 "Prob(y>2)" 3 "Prob(y>3)") rows(1)) xtitle("Week")

generate pr12 = 1-mprobg2
generate pr123 = 1-mprobg3
generate pr1234 = 1

twoway (area pr1 week, sort fintensity(inten10))                           ///
   (rarea pr12 pr1 week, sort fintensity(inten50))                         ///
   (rarea pr123 pr12 week, sort fintensity(inten70))                       ///
   (rarea pr1234 pr123 week, sort fintensity(inten90)), by(treatment)      ///
    legend(order(1 "Prob(y=1)" 2 "Prob(y=2)" 3 "Prob(y=3)" 4 "Prob(y=4)")) ///
	xtitle("Week")
 
 
** 11.8.2 Predicted subject-specific probabilities: Posterior mean

estimates restore model2_gllamm

fillin id week
replace weeksqrt = sqrt(week) if _fillin==1
egen trt = mean(treatment), by(id)
replace treatment = trt if _fillin==1
replace interact = weeksqrt*treatment if _fillin==1
replace cons = 1 if _fillin==1

gllapred pprobg1, mu above(1) fsample
gllapred pprobg2, mu above(2) fsample
gllapred pprobg3, mu above(3) fsample

generate week0 = week==0
sort treatment id week
by treatment: generate newid = sum(week0)
list id week week0 newid in 1/16, sepby(id) noobs

twoway (line pprobg1 week, sort lpatt(solid))                      ///
  (line pprobg2 week, sort lpatt(vshortdash))                      ///
  (line pprobg3 week, sort lpatt(dash))                            ///
  (scatter pprobg1 week if _fillin==1, msym(o))                    ///
  (scatter pprobg2 week if _fillin==1, msym(o))                    ///
  (scatter pprobg3 week if _fillin==1, msym(o) mcol(black))        ///
  if newid<13&treatment==0, by(newid)                              ///
  legend(order(1 "Prob(y>1)" 2 "Prob(y>2)" 3 "Prob(y>3)") rows(1)) xtitle("Week")

twoway (line pprobg1 week, sort lpatt(solid))                      ///
  (line pprobg2 week, sort lpatt(vshortdash))                      ///
  (line pprobg3 week, sort lpatt(dash))                            ///
  (scatter pprobg1 week if _fillin==1, msym(o))                    ///
  (scatter pprobg2 week if _fillin==1, msym(o))                    ///
  (scatter pprobg3 week if _fillin==1, msym(o) mcol(black))        ///
  if newid<13&treatment==1, by(newid)                              ///
  legend(order(1 "Prob(y>1)" 2 "Prob(y>2)" 3 "Prob(y>3)") rows(1)) xtitle("Week")


*** 11.9 Do experts differ in their grading of student essays?

use https://www.stata-press.com/data/mlmus4/essays, clear

tabulate grade


*** 11.10 A random-intercept probit model with grader bias

** 11.10.1 Model specification
* Estimation using gllamm
quietly tabulate grader, generate(grad)

gllamm grade grad2-grad5, i(essay) link(oprobit) adapt
estimates store model1

* meoprobit grade i.grader || essay:, intpoints(8)

*** 11.11 Including grader-specific measurement error variances

** 11.11.1 Model specification 
*Estimation using gllamm
eq het: grad2-grad5
matrix a = e(b)
gllamm grade grad2-grad5, i(essay) link(soprobit) s(het) from(a) adapt

estimates store model2
lrtest model1 model2

* hetoprobit grade i.grader, het(i.grader)

*** 11.12 Including grader-specific thresholds

** 11.12.2 Model specification
*Estimation using gllamm
eq thr: grad2-grad5
gllamm grade, i(essay) link(soprobit) s(het) thresh(thr) from(a) adapt skip

estimates store model3
lrtest model2 model3

test [_cut11=_cut12=_cut13=_cut14=_cut15=_cut16=_cut17]: grad2

