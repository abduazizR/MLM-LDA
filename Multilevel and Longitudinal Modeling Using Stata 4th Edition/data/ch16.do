*******************************************************************
****Chapter 16: Models with nested and crossed random effects *****
*******************************************************************


*** 16.2 Did the Guatemalan immunization campaign work?

use https://www.stata-press.com/data/mlmus4/guatemala, clear


*** 16.3 Three-level random-intercept logistic regression model

** 16.3.4 Estimation 

* Using melogit
melogit immun kid2p indNoSpa indSpa momEdPri momEdSec husEdPri husEdSec ///
 husEdDK rural pcInd81 || cluster: || mom:, intpoints(15)
 
melogit, or
 
display (_b[/var(_cons[cluster])]+_b[/var(_cons[cluster>mom])]) ///
    /(_b[/var(_cons[cluster])]+_b[/var(_cons[cluster>mom])]+_pi^2/3)
	
display (_b[/var(_cons[cluster])]) ///
   /(_b[/var(_cons[cluster])]+_b[/var(_cons[cluster>mom])]+_pi^2/3)

estat icc


display exp(sqrt(2*(_b[/var(_cons[cluster])] ///
+_b[/var(_cons[cluster>mom])]))*invnormal(3/4))
display exp(sqrt(2*(_b[/var(_cons[cluster>mom])]))*invnormal(3/4))

*display chi2tail(1,e(chi2_c))/2 + chi2tail(2,e(chi2_c))/4


* Using gllamm 

gllamm immun kid2p indNoSpa indSpa  momEdPri momEdSec husEdPri husEdSec  ///
   husEdDK rural pcInd81, family(binomial) link(logit) i(mom cluster) nip(5)

matrix a=e(b)
gllamm immun kid2p indNoSpa indSpa  momEdPri momEdSec husEdPri husEdSec  ///
   husEdDK rural pcInd81, family(binomial) link(logit) i(mom cluster) ///
   from(a) adapt
estimates store glri8_gllamm
gllamm, eform


*** 16.4 Three-level random-coefficient logistic regression model

** 16.4.1 Estimation

* Using melogit 

melogit immun kid2p rural pcInd81 || cluster: kid2p, covariance(unstructured) ///
 || mom:, intpoints(15) or
estimates store rc_0

estat sd

melogit immun kid2p rural pcInd81 || cluster: || mom:, intpoints(15) or
estimates store ri_0

lrtest rc_0 ri_0

* Using gllamm 
estimates restore glri8_gllamm
matrix a=e(b)
gllamm immun kid2p rural pcInd81, family(binomial) link(logit) ///
   i(mom cluster) from(a) skip adapt eform

generate cons=1
eq inter: cons
eq slope: kid2p
matrix a = e(b)
matrix a = (a,.2,0)
gllamm immun kid2p rural pcInd81, family(binomial) link(logit) i(mom cluster) ///
   nrf(1 2) eqs(inter inter slope) nip(8 4 4) from(a) copy adapt eform
estimates store rc0_gllamm



*** 16.5 Prediction of random effects

** 16.5.1 Empirical Bayes prediction

estimates restore rc_0
predict comms commi mother, reffects ebmeans
sort cluster mom
list cluster mom mother commi comms in 1/9, sepby(mom) noobs

egen pick_comm = tag(cluster)
twoway scatter comms commi if pick_com==1, xtitle(Intercept) ytitle(Slope)

** 16.5.2 Empirical Bayes modal prediction

predict comms2 commi2 mother2, reffects ebmodes
sort cluster mom
list cluster mom mother2 commi2 comms2 in 1/9, sepby(mom) noobs

*** 16.6 Different kinds of predicted probabilities

** 16.6.1 Predicted population-averaged or marginal probabilities: New clusters

predict margp, marginal

** 16.6.2 Predicted median or conditional probabilities

predict condp, conditional(fixedonly)	

label define r 0 "Urban" 1 "Rural"
label values rural r
twoway (line condp pcInd81 if kid2p==0, lpatt(solid) sort) ///
       (line condp pcInd81 if kid2p==1, lpatt(solid) sort) ///
       (line margp pcInd81 if kid2p==0, lpatt(dash) sort) ///
       (line margp pcInd81 if kid2p==1, lpatt(dash) sort), ///
        by(rural) legend(order(1 "Conditional" 3 "Marginal")) ///
        xtitle(Perentage Indigenous) ytitle(Probability)

** 16.6.3 Predicted posterior mean probabilities: Existing clusters

estimates restore rc0_gllamm

gllapred postp, mu
sort cluster mom
list cluster mom kid2p mother commi comms postp in 1/9, sepby(mom) noobs


*** 16.7 Do salamanders from different populations mate successfully?

use https://www.stata-press.com/data/mlmus4/salamander, clear

table female male if group==1, statistic(mean y) nototal


*** 16.8 Crossed random effects logistic regression

** 16.8.1 Setup for estimating crossed random-effects model using melogit

generate m = male - (group-1)*10
generate f = female - (group-1)*10

egen pickmale = tag(male)
sort m group
list group male m if pickmale==1&m<3, sepby(m) noobs

generate ww = wsf*wsm

** 16.8.2 Approximate maximum likelihood estimation
* Estimation using melogit
melogit y wsm wsf ww || group: R.m || f:

set processors 1
matrix a = e(b)
melogit y wsm wsf ww || group: R.m || f:, intpoints(2) from(a) 

matrix a=e(b)
melogit y wsm wsf ww || group: R.m || f:, intpoints(3) from(a)  
estimates store melogit_3p

set processors `c(processors_lic)'

display chi2tail(1,e(chi2_c))/2 + chi2tail(2,e(chi2_c))/4

lincom _cons, or
lincom _cons + wsf, or
lincom _cons + wsm, or
lincom _cons + wsm + wsf + ww, or

display exp(sqrt(2*(_b[/var(_cons[group>m])] ///
 +_b[/var(_cons[group>f])]))*invnormal(3/4))

 
** 16.8.3 Bayesian estimation

twoway function exp(lnigammaden(1,1,x)), range(0.001 25) xline(.5, lpatt(dash))

* Estimation using bayes: melogit
bayes, normalprior(5) igammaprior(1 1) dryrun: ///
  melogit y wsm wsf ww || group: R.m || f:

bayes, normalprior(5) igammaprior(1 1) initrandom rseed(232123): ///
  melogit y wsm wsf ww || group: R.m || f:

bayesstats ess

bayes, normalprior(5) igammaprior(1 1) initrandom rseed(232123) nchains(10) ///
  mcmcsize(2500) burnin(2500):                                            ///
  melogit y wsm wsf ww || group: R.m || f:
  
bayesstats grubin

bayes, normalprior(5) igammaprior(1 1) initrandom rseed(232123) ///
  mcmcsize(1000000) burnin(2500): ///
  melogit y wsm wsf ww || group: R.m || f:
 
bayes, saving(bigrun, replace)
estimates save bigrun, replace


** 16.8.5 Fully Bayesian versus empirical Bayesian inference for random effects

estimates use bigrun
bayesstats summary {UU0[1,1/3] VV0[1,1/3]} 

matrix Bayes_summary = r(summary)

bayesgraph histogram {UU0[1,1/3] VV0[1,1/3]}, skip(99) byparm(norescale) ///
   normal graphopts(fcolor(gs14) lcolor(gs6) lwidth(thin))

matrix list Bayes_summary

global B_mean = Bayes_summary["VV0[group>f]:1 2","Mean"]
global B_sd = Bayes_summary["VV0[group>f]:1 2","Std dev"]

bayesgraph kdensity {VV0[1,2]}, normal xline($B_mean) graphopts(lcol(gs8)) ///
 legend(off)

matrix Bayes_est = e(mean)

quietly melogit y wsm wsf ww || group: R.m || f:, noestimate
matrix colnames Bayes_est = `: colfullname e(b)'

melogit y wsm wsf ww || group: R.m || f:, from(Bayes_est) intpoints(3) ///
 noestimate

predict eb_mode*, reffects ebmodes reses(se_mode*)

set processors 1

predict eb_mean* if group==1, reffects ebmeans reses(se_mean*) intpoints(3)

set processors `c(processors_lic)'

tabstat eb_mean1 se_mean1 eb_mode1 se_mode1 if group==1&m<4, by(m) nototal
tabstat eb_mean2 se_mean2 eb_mode2 se_mode2 if group==1&f<4, by(f) nototal


quietly summ eb_mean2 if group==1&f==2
global eb_mean = r(mean)
quietly summ se_mean2 if group==1&f==2
global se_mean = r(mean)

quietly summ eb_mode2 if group==1&f==2
global eb_mode = r(mean)
quietly summ se_mode2 if group==1&f==2
global se_mode = r(mean)

   
twoway (function normalden(x, $B_mean, $B_sd), range(-3 6) lpatt(solid)        ///
   lcol(gs10) lwidth(thick))                                                   ///
   (function normalden(x, $eb_mean, $se_mean), range(-3 6) lpatt(dash))        ///
   (function normalden(x, $eb_mode, $se_mode), range(-3 6) lpatt(solid)),      ///
   legend(order(3 "EB Mode" 2 "EB Mean" 1 "Full Bayes") rows(1))               ///
   xtitle(Random effect) xline($B_mean, lpatt(solid) lcol(gs10) lwidth(thick)) ///
   xline($eb_mean, lpatt(dash))  xline($eb_mode, lpatt(solid))                 ///
   ytitle(Normal approximation)

