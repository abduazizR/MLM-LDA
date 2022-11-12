**********************************************
****Chapter 4: Random-coefficient models *****
**********************************************


*** 4.1 How effective are different schools?

use https://www.stata-press.com/data/mlmus4/gcse, clear


*** 4.3 Separate linear regressions for each school

regress gcse lrt if school==1
predict p_gcse, xb 
twoway (scatter gcse lrt) (line p_gcse lrt, sort) if school==1, ///
  xtitle(LRT) ytitle(GCSE)

twoway (scatter gcse lrt) (lfit gcse lrt, sort lpatt(solid)), ///
  by(school, compact legend(off) cols(5)) ///
  xtitle(LRT) ytitle(GCSE) ysize(3) xsize(2)

egen num = count(gcse), by(school)
statsby inter=_b[_cons] slope=_b[lrt], by(school) saving(ols): ///
  regress gcse lrt if num>4 
sort school
merge m:1 school using ols
drop _merge

twoway scatter slope inter, xtitle(Intercept) ytitle(Slope)

egen pickone = tag(school)
summarize inter slope if pickone==1
correlate inter slope if pickone==1, covariance

generate pred = inter + slope*lrt

sort school lrt
twoway (line pred lrt, connect(ascending)), xtitle(LRT) ///
   ytitle(Fitted regression lines)

   
*** 4.5 Estimation using mixed

** 4.5.1 Random-intercept model
mixed gcse lrt || school:, mle stddeviations vce(robust)
estimates store ri

** 4.5.2 Random-coefficient model
mixed gcse lrt || school: lrt, covariance(unstructured) mle stddeviations ///
 vce(robust)

* mixed, variance
estat recovariance
estimates store rc

sort school lrt
list school lrt if school==54, clean noobs
estat wcorrelation, at(school=54)


*** 4.6 Testing the slope variance

lrtest ri rc, force
display 0.5*chi2tail(1,40.37) + 0.5*chi2tail(2,40.37)   


*** 4.7 Interpretation of estimates

/*
display -.115 - 1.96*3.007
display -.115 + 1.96*3.007

display 0.557 - 1.96*0.121
display 0.557 + 1.96*0.121
*/

twoway function sqrt(9.0447+2*0.1804*x+0.0145*x^2+55.3653), range(-30 30) ///
 xtitle(LRT) ytitle(Estimated standard deviation of total residual)

 
*** 4.8 Assigning values to the random intercepts and slopes

** 4.8.1 Maximum "likelihood" estimation

estimates restore rc

predict fixed, xb
generate totres = gcse - fixed

statsby mli=_b[_cons] mls=_b[lrt], by(school) saving(ols, replace): ///
  regress totres lrt
sort school
merge m:1 school using ols
drop _merge

list lrt gcse mli mls if school==48, clean noobs

** 4.8.2 Empirical Bayes prediction

estimates restore rc
predict ebs ebi, reffects

list school mli ebi mls ebs if pickone==1 & (school<10 | school==48), noobs

twoway (scatter ebi mli if pickone==1 & school!=48, mlabel(school)) ///
  (function y=x, range(-10 10)), legend(off) xtitle(ML estimates of intercepts) ///
  ytitle(EB predictions of intercepts) legend(off) xline(0)

twoway (scatter ebs mls if pickone==1 & school!=48, mlabel(school)) ///
  (function y=x, range(-0.6 0.6)), xtitle(ML estimates of slopes)    ///
  ytitle(EB predictions of slopes) legend(off) xline(0)

** 4.8.3 Model visualization

predict murc, fitted

sort school lrt
twoway (line murc lrt, connect(ascending)), xtitle(LRT) ///
  ytitle(Empirical Bayes regression lines for model 2)

estimates restore ri
predict muri, fitted
sort school lrt
twoway (line muri lrt, connect(ascending)), xtitle(LRT) ///
  ytitle(Empirical Bayes regression lines for model 1)


** 4.8.4 Residual diagnostics

histogram ebi if pickone==1, normal xtitle(Predicted random intercepts)

histogram ebs if pickone==1, normal xtitle(Predicted random slopes)

scatter ebs ebi if pickone==1, saving(yx, replace)    ///
   xtitle("Random intercept") ytitle("Random slope")  ///
   ylabel(, nogrid)

histogram ebs if pickone==1, freq horizontal saving(hy, replace) normal ///
    yscale(alt) ytitle(" ") fxsize(35) ylabel(, nogrid)

histogram ebi if pickone==1, freq saving(hx, replace) normal ///
    xscale(alt) xtitle(" ") fysize(35) ylabel(, nogrid)

graph combine hx.gph yx.gph hy.gph, hole(2) imargin(0 0 0 0)

predict res1, residuals
histogram res1, normal xtitle(Predicted level-1 residuals)

** 4.8.5 Inferences for individual schools

quietly mixed gcse lrt || school: lrt, covariance(unstructured) reml
predict slope1 inter1, reffects reses(slope_se inter_se)

gsort + inter1 -pickone
generate rank = sum(pickone)

generate labpos = inter1 + 1.96*inter_se + .5

serrbar inter1 inter_se rank if pickone==1, addplot(scatter labpos rank,   ///
  mlabel(school) msymbol(none) mlabpos(0)) scale(1.96) xtitle(Rank)        ///
  ytitle(Prediction) legend(off)

generate lower = inter1 - 1.96*inter_se
generate upper = inter1 + 1.96*inter_se

twoway (rcap lower upper rank, blpatt(solid) lcol(black))       ///
  (scatter inter1 rank)                                         ///
  (scatter labpos rank, mlabel(school) msymbol(none) mlabpos(0) ///
           mlabcol(black) mlabsiz(medium)),                     ///
  xtitle(Rank) ytitle(Prediction) legend(off)                   ///
  xscale(range(1 65)) xlabel(1/65) ysize(1)

  
*** 4.9 Two-stage formulation 

mixed gcse i.schgend##c.lrt || school: lrt, covariance(unstructured) mle ///
   stddeviations vce(robust)

