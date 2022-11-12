*******************************************************************
****Chapter 8: Higher-level models with nested random effects *****
*******************************************************************


*** 8.2 Do peak-expiratory-flow measurements vary between methods within subjects?

use https://www.stata-press.com/data/mlmus4/pefr, clear

reshape long wm wp, i(id) j(occasion) 
generate i = _n 
reshape long w, i(i) j(meth) string 
sort id meth occasion
list id meth occasion w in 1/8, clean noobs

encode meth, generate(method) 
recode method 2=0
label define m 0 "Wright" 1 "Mini Wright"
label values method m


*** 8.3 Inspecting sources of variability

twoway (scatter w id if method==0, msymbol(circle))                 ///
       (scatter w id if method==1, msymbol(circle_hollow)),         ///
       xtitle(Subject id) ytitle(Peak-expiratory-flow measurements) ///
       legend( order(1 "Wright" 2 "Mini Wright")) xlabel(1/17) 

*** 8.6 Estimation using mixed

mixed w || id: || method:, reml dfmethod(kroger) stddeviations

estat icc
estat wcorrelation, at(id=1,method=1)
list id method occ if id==1
estat wcorrelation, at(id=1)

*** 8.7 Empirical Bayes prediction

predict subj instr, reffects
sort id method

list id method subj instr if id<8 & occasion==1, noobs sepby(id)


*** 8.8 Testing variance components

estimates store thrlev

quietly mixed w || id:, reml dfmethod(kroger) stddeviations
lrtest thrlev .


*** 8.9 Crossed versus nested random effects revisited

mixed w method || id: || method:, reml dfmethod(kroger) stddeviations


*** 8.10 Does nutrition affect cognitive development of Kenyan children?

use https://www.stata-press.com/data/mlmus4/kenya, clear


*** 8.11 Describing and plotting three-level data

** 8.11.1 Data structure and missing data

egen child = group(id schoolid)

quietly xtset child rn
xtdescribe if !missing(ravens)

* number of students and response rates for schools
egen numobs2 = count(ravens), by(child)
generate compl = numobs2==5
generate any = numobs2>0

egen pick_child = tag(child)
egen numcomp3 = total(complete*pick_child), by(schoolid)
egen numany3 = total(any*pick_child), by(schoolid)

generate prop = numcomp3/numany3

egen pick_school = tag(schoolid)
summarize numany3 numcomp3 prop if pick_school==1

** 8.11.2 Level-1 variables

quietly xtset child
xtsum ravens relyear

egen mn_raven = mean(ravens), by(child)
egen mn_relyr = mean(relyear), by(child)

** 8.11.3 Level-2 variables

quietly xtset schoolid
xtsum mn_raven mn_relyr boy age_at_time0 if pick_child==1

** 8.11.4 Level-3 variables

tabulate treatment if pick_school==1

tabulate treatment if pick_child==1

** 8.11.5 Plotting growth trajectories

sort schoolid id rn

twoway (line ravens relyear, connect(ascending)), by(schoolid, compact) ///
  xtitle(Time in years) ytitle(Raven's score)
  

*** 8.12 Three-level random-intercept model

** 8.12.3 Estimation using mixed

mixed ravens ib4.treatment##c.relyear age_at_time0 i.boy || ///
    schoolid: || id:, reml dfmethod(kroger) stddeviations
estimates store mod1


*** 8.13 Three-level random-coefficient models

** 8.13.1 Random coefficient at the child level   

mixed ravens ib4.treatment##c.relyear age_at_time0 i.boy || ///
   schoolid: || id: relyear, covariance(unstructured) reml ///
   dfmethod(kroger) stddeviations
estimates store mod2
lrtest mod1 mod2

** 8.13.2 Random coefficient at the child and school levels
   
mixed ravens ib4.treatment##c.relyear age_at_time0 i.boy ||  ///
   schoolid: relyear, covariance(unstructured) || ///
   id: relyear, covariance(unstructured) reml dfmethod(kroger) stddeviations
estimates restore mod2


*** 8.14 Residual diagnostics and predictions

predict ri3 rc2 ri2, reffects
predict res, residuals

replace ri3=. if pick_school!=1
replace ri2=. if pick_child!=1
graph box ri3 ri2 res, ascategory box(1, bstyle(outline)) ///
    yvaroptions(relabel(1 "School" 2 "Child" 3 "Occasion")) ///
	medline(lcolor(black))

scatter rc2 ri2 if pick_child==1,  saving(yx, replace) ///
   xtitle("Random intercept") ytitle("Random slope")
histogram rc2, freq horiz saving(hy, replace) ///
    yscale(alt) ytitle(" ") fxsize(35) normal
histogram ri2, freq saving(hx, replace) ///
    xscale(alt) xtitle(" ") fysize(35) normal
   
graph combine hx.gph yx.gph hy.gph, hole(2) imargin(0 0 0 0)

predict predtraj, fitted
sort schoolid id relyear

twoway (line predtraj relyear, connect(ascending)), by(schoolid, compact) ///
   xtitle(Time in years) ytitle(Raven's score)

summarize age_at_time0 if pick_child==1
replace age_at_time0 = r(mean)
replace boy = 1

predict means, xb

twoway (line means relyear if treatment==1, sort lpatt(solid)) ///
       (line means relyear if treatment==2, sort lpatt(dash)) ///
	   (line means relyear if treatment==3, sort lpatt(shortdash)) ///
	   (line means relyear if treatment==4, sort lpatt(longdash_dot)), ///
	   legend(order(1 "Meat" 2 "Milk" 3 "Calorie" 4 "Control")) ///
	   xtitle(Time in years) ytitle(Predicted mean Raven's score)
