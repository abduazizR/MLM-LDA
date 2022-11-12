***********************************************
****Chapter 2: Variance-components models *****
***********************************************


*** 2.2 How reliable are peak-expiratory-flow measurements?
use https://www.stata-press.com/data/mlmus4/pefr, clear

*** 2.3 Inspecting within-subject dependence
generate mean_wm = (wm1+wm2)/2
summarize mean_wm

twoway (scatter wm1 id, msymbol(circle))                                ///
       (scatter wm2 id, msymbol(circle_hollow)),                        ///
       xtitle(Subject id) xlabel(1/17) ytitle(Mini Wright measurements) ///
       legend( order(1 "Occasion 1" 2 "Occasion 2")) yline(453.9118)    

*** 2.5 Estimation using Stata

** 2.5.1 Data preparation: Reshaping to long form
list if id<6, clean noobs

* stack repeated observations for same method into single variable
reshape long wp wm, i(id) j(occasion)
list if id<6, clean noobs

** 2.5.2 Using xtreg
xtset id
xtreg wm, mle

** 2.5.3 Using mixed
mixed wm || id:, mle stddeviations
mixed wm || id:, reml stddeviations

estat icc


*** 2.6 Hypothesis tests and confidence intervals

mixed wm || id:, reml dfmethod(kroger) stddeviations

** 2.6.2 for between-cluster variance

* likelihood-ratio test
quietly mixed wm || id:, reml
estimates store ri
quietly mixed wm, reml
lrtest ri .

* score test 
quietly xtset id
quietly xtreg wm, re
xttest0

* F test
quietly xtset id
xtreg wm, fe


*** 2.9 Crossed versus nested effects

* random intercept model with occasion 2 dummy

generate occ2 = occasion==2
mixed wm occ2 || id:, reml stddeviations


*** 2.10 Parameter estimation

** 2.10.3 Inference for beta

* model-based and robust SE's for OLS

regress wm
regress wm, vce(cluster id)


*** 2.11 Assigning values to the random intercepts

** 2.11.1 Maximum "likelihood" estimation
* via OLS
quietly mixed wm || id:, reml 
predict pred, xb
generate res = wm - pred

regress res ibn.id, noconstant

quietly xtset id
quietly xtreg wm, mle
predict ml2, u

* via mean total residual
egen ml = mean(res), by(id)
sort id
display ml[1]

** 2.11.2 Empirical Bayes prediction

* using shrinkage factor
display 110.40^2/(110.40^2+(19.91^2)/2)

* find out how to access variances from stored estimates
mixed wm || id:, reml estmetric coeflegend
display exp(_b[lns1_1_1:_cons])
display exp(_b[lnsig_e:_cons])

display exp(_b[lns1_1_1:_cons])^2/(exp(_b[lns1_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2/2)

estimates store mixed
estat sd, post coeflegend

display _b[id:sd(_cons)]^2/(_b[id:sd(_cons)]^2 + ///
	_b[Residual:sd(e)]^2/2)

estimates restore mixed

generat eb1 = .98399606*ml

* using mixed
predict eb2, reffects
sort id
format eb1 eb2 %8.2f
list id eb1 eb2 if occasion==1, clean noobs

*** 2.11.3 Empirical Bayes standard errors

* comparative standard errors
generate comp_se = sqrt((1-.98399821)*exp(_b[lns1_1_1:_cons])^2)
display comp_se[1]


* mixed wm || id:, mle
* predict eb, reffects reses(comp_se)

* diagnostic standard errors
generate diag_se = sqrt(exp(_b[lns1_1_1:_cons])^2 - comp_se^2)
display diag_se[1]

* accounting for uncertainty in beta hat
quietly mixed wm || id:, reml
predict eb_reml, reffects reses(comp_seb)
display comp_seb[1]

generate diag_seb = sqrt(exp(_b[lns1_1_1:_cons])^2 - comp_seb^2)
display diag_seb[1]
