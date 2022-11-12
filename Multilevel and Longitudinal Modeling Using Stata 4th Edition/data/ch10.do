*****************************************************
****Chapter 10: Dichotomous or binary responses *****
*****************************************************


*** 10.2 Single level models for dichotomous responses

** 10.2.1 Generalized linear model formulation

use https://www.stata-press.com/data/mlmus4/womenlf, clear

recode workstat 2=1
logit workstat husbinc i.chilpres 
logit workstat husbinc i.chilpres, or

predict prob, p

twoway (line prob husbinc if chilpres==0, sort)         ///
 (line prob husbinc if chilpres==1, sort lpatt(dash)), ///
 legend(order(1 "No child" 2 "Child"))                  ///
 xtitle("Husband's income/$1000") ytitle("Probability that wife works")

twoway (function y=invlogit(_b[husbinc]*x+_b[_cons]), range(-100 100))    ///
 (function y=invlogit(_b[husbinc]*x+_b[1.chilpres]+_b[_cons]),            ///
 range(-100 100) lpatt(dash)),                                           ///
 xtitle("Husband's income/$1000") ytitle("Probability that wife works") ///
 legend(order(1 "No child" 2 "Child")) xline(1) xline(45)
 
quietly margins chilpres, at(husbinc=(-100(10)100)) 

marginsplot, noci recast(line) plot2opts(lpatt(dash))  ///
	  legend(order(1 "No child" 2 "Child"))                  ///
 	  ytitle("Probability that wife works") xline(1) xline(45)


glm workstat husbinc i.chilpres, link(logit) family(binomial)

** 10.2.2 Latent-response formulation
probit workstat husbinc i.chilpres

twoway (function y=invlogit(1.3358-0.0423*x), range(-100 100))           ///
 (function y=normal(0.7982-0.0242*x), range(-100 100) lpatt(dash)),        ///
 xtitle("Husband's income/$1000") ytitle("Probability that wife works") ///
 legend(order(1 "Logit link" 2 "Probit link")) xline(1) xline(45)

*** 10.3 Which treatment is best for toenail infection?

use https://www.stata-press.com/data/mlmus4/toenail, clear

*** 10.4 Longitudinal data structure

xtset patient visit
xtdescribe if !missing(outcome)


*** 10.5 Proportions and fitted population-averaged or marginal probabilities

label define tr 0 "Itraconazole" 1 "Terbinafine"
label values treatment tr

graph bar (mean) proportion = outcome, over(visit) by(treatment) ///
   ytitle(Proportion with onycholysis)

egen prop = mean(outcome), by(treatment visit)
egen mn_month = mean(month), by(treatment visit)

twoway line prop mn_month, by(treatment) sort ///
  xtitle(Time in months) ytitle(Probability with onycholysis)


logit outcome i.treatment##c.month, or vce(cluster patient)


predict prob, pr

twoway (line prop mn_month, sort) (line prob month, sort lpatt(dash)),          ///
   by(treatment) legend(order(1 "Observed proportions" 2 "Fitted probabilities")) ///
   xtitle(Time in months) ytitle(Probability of onycholysis)


*** 10.6 Random-intercept logistic regression

** 10.6.3 Estimation
* Using xtlogit
quietly xtset patient
xtlogit outcome i.treatment##c.month, intpoints(30)
xtlogit, or

* Using melogit
melogit outcome i.treatment##c.month || patient:, intpoints(30)
estimates store melogit

** Using gllamm
which gllamm
ssc install gllamm, replace

generate trt_month = treatment*month
gllamm outcome treatment month trt_month, i(patient) link(logit) /// 
 family(binomial) nip(30) adapt
estimates store gllamm
*gllamm, eform


*** 10.7 Subject-specific or conditional vs. population averaged or marginal relationships

display (invlogit(1) + invlogit(2))/2
display invlogit((1+2)/2)


*** 10.8 Measures of dependence and heterogeneity

** 10.8.1 Conditional or residual intraclass correlation of the latent responses

estimates restore melogit
estat icc

** 10.8.2 Median odds ratio

melogit, coeflegend
display exp(sqrt(2*_b[/var(_cons[patient])])*invnormal(3/4))

** 10.8.3 Measures of association for observed responses...

capture net install st0031, from(https://www.stata-journal.com/software/sj3-1)

* install xtrho: use findit xtrho, select st0031, then install

quietly xtset patient
quietly: xtlogit outcome i.treatment##c.month, re intpoints(30)
xtrho

*copy https://www.stata-press.com/data/mlmus4/ch10table.do ch10table.do


*** 10.9 Inference for random-intercept logistic models 

** 10.9.1 Tests and confidence intervals for odds ratios
lincom 1.treatment + 1.treatment#c.month*20
lincom 1.treatment + 1.treatment#c.month*20, or
*lincom -(1.treatment + 1.treatment#c.month*20), or
/* margins can also be used, but it doesn't save much effort */
  ** margins , dydx(treatment) at(month=20) predict(xb)

*** 10.10 Maximum likelihood estimation

** 10.10.2 Some speed and acuracy considerations
* Laplace approximation
melogit outcome i.treatment##c.month || patient:, intmethod(laplace)


*** 10.11 Assigning values to random effects

** 10.11.1 Maximum "likelihood" estimation

estimates restore melogit
predict offset, xb

statsby mlest=_b[_cons], by(patient) saving(ml, replace): logit outcome, ///
 offset(offset)

sort patient
merge m:1 patient using ml
drop _merge


** 10.11.2 Empirical Bayes prediction

estimates restore melogit
predict eb, reffects ebmeans reses(semean) 


** 10.12.3 Empirical Bayes modal prediction

predict ebmodal, reffects ebmodes reses(semode) 

twoway (rspike mlest ebmodal eb if visit==1)                          ///
  (scatter mlest  eb if visit==1, msize(small) msym(th) mcol(black))  ///
  (scatter ebmodal eb if visit==1, msize(small) msym(oh) mcol(black)) ///
  (function y=x, range(eb) lpatt(solid)),                             ///
  xtitle(Empirical Bayes prediction)                                    ///
  legend(order(2 "Maximum likelihood" 3 "Empirical Bayes modal"))

egen num0 = total(outcome==0), by(patient)
egen num1 = total(outcome==1), by(patient)
list patient num0 num1 eb ebmodal semean semode if visit==1&patient<=12, noobs


*** 10.12 Different kinds of predicted probabilities

** 10.12.1 Predicted population-averaged or marginal probabilities

predict margprob, pr marginal   

twoway (line prob month, sort) (line margprob month, sort lpatt(dash)),      ///
   by(treatment) legend(order(1 "Ordinary logit" 2 "Random intercept logit")) ///
   xtitle(Time in months) ytitle(Fitted marginal probabilities of onycholysis)

** 10.12.2 Predicted subject-specific probabilities

* Predictions for hypothetical subjects: Conditional probabilities

predict fixedpart, xb

gen condprobm4 = invlogit(fixedpart-4)
gen condprobm2 = invlogit(fixedpart-2)
gen condprob0  = invlogit(fixedpart)
gen condprob2  = invlogit(fixedpart+2)
gen condprob4  = invlogit(fixedpart+4)

twoway (line prop mn_month, sort)                               ///
 (line margprob month, sort lpatt(dash))                        ///
 (line condprob0 month, sort lpatt(shortdash_dot))              ///
 (line condprob4 month, sort lpatt(shortdash))                  ///
 (line condprobm4 month, sort lpatt(shortdash))                 ///
 (line condprob2 month, sort lpatt(shortdash))                  ///
 (line condprobm2 month, sort lpatt(shortdash)),               ///
 by(treatment)                                                  ///
 legend(order(1 "Observed proportion" 2 "Marginal probability"  ///
    3 "Median probability" 4 "Conditional probabilities"))      ///
 xtitle(Time in months) ytitle(Probabilities of onycholysis)

* Predictions for the subjects in the sample: Posterior mean

estimates restore gllamm
gllapred cmu, mu

sort patient visit
list patient visit if patient==2|patient==15, sepby(patient) noobs
fillin patient visit
list patient visit _fillin if patient==2|patient==15, sepby(patient) noobs

egen trt = mean(treatment), by(patient)
replace treatment = trt if _fillin==1
drop mn_month
egen mn_month = mean(month), by(treatment visit)
replace month = mn_month if _fillin==1
replace trt_month = treatment*month
drop cmu

gllapred cmu, mu fsample

list patient visit _fillin cmu if patient==2|patient==15, sepby(patient) noobs

set seed 1234421
sort patient
generate rand = runiform() if visit==1
by treatment (rand), sort: generate randid = _n if rand<.
egen randomid = mean(randid), by(patient)

twoway (line cmu month, sort) (scatter cmu month if _fillin==1, mcol(black)) ///
    if randomid<=16&treatment==0,  by(patient, compact legend(off) ///
	l1title("Posterior mean probabilities"))
	
twoway (line cmu month, sort) (scatter cmu month if _fillin==1, mcol(black)) ///
    if randomid<=16&treatment==1,  by(patient, compact legend(off) ///
	l1title("Posterior mean probabilities"))

 
*** 10.13 Other approaches to clustered dichotomous data

** 10.13.1 Conditional logistic regression

clogit outcome month i.treatment#c.month, group(patient) or

** 10.13.2 Generalized estimating equations (GEE)

quietly xtset patient
xtgee outcome i.treatment##c.month, link(logit) ///
  family(binomial) corr(exchangeable) vce(robust) eform
  
estat wcorrelation, format(%4.3f)
