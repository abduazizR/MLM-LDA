***********************************************
****Chapter 1: Review of Linear Regression*****
***********************************************


*** 1.2 Is there gender discrimination in faculty salaries?

use https://www.stata-press.com/data/mlmus4/faculty, clear


*** 1.3 Independent-samples t test

tabstat salary, by(male) statistics(mean sd n)

graph box salary, over(male) ytitle(Academic salary) asyvars

histogram salary, by(male, rows(2)) xtitle(Academic salary)

generate lsalary = log10(salary)

graph box lsalary, over(male) ytitle(Log academic salary) asyvars

histogram lsalary, by(male, rows(2)) xtitle(Log academic salary)

ttest salary, by(male)
ttest salary, by(male) unequal


*** 1.4 One-way analysis of variance

anova salary male

margins male


*** 1.5 Simple linear regression

twoway (scatter salary market, msize(small)) ///
       (lowess salary market, lwidth(medthick) lpatt(solid))

regress salary market

egen mn_market = mean(market)
generate marketc = market - mn_market

regress salary marketc

predict yhat, xb
twoway (scatter salary market, msize(small)) ///
       (line yhat market, sort lwidth(medthick) lpatt(solid)), ///
        ytitle(Academic salary) xtitle(Marketability)

*** 1.6 Dummy variables

regress salary male
regress salary male, vce(robust)


*** 1.7 Multiple linear regression

tabstat marketc, by(male) statistics(mean sd)

regress salary male marketc
predict yhat2, xb

twoway (scatter salary marketc if male==1, msymbol(o))    ///
       (line yhat2 marketc if male==1, sort lpatt(dash)) ///
       (scatter salary marketc if male==0, msymbol(oh))   ///
       (line yhat2 marketc if male==0, sort lpatt(solid)), ///
       ytitle(Academic salary) xtitle(Mean-centered marketability)    ///
       legend(order(1 " " 2 "Men" 3 " " 4 "Women"))

twoway (kdensity marketc if male==0) (kdensity marketc if male==1),  ///
       legend(order(2 "Men" 1 "Women")) xtitle(Mean-centered marketability) ///
	   ytitle(Estimated density)
	   
regress salary male marketc yearsdg

regress salary i.male marketc yearsdg

margins male, at(marketc=0 yearsdg=10)
margins male


*** 1.8 Interactions

generate male_years = male*yearsdg

regress salary male marketc yearsdg male_years

lincom male + male_years*10

twoway (function Women = 36773 + 763.19*x, range(0 41) lpatt(dash)) ///
 (function Men = 36773 + -593.31 + (763.19 + 227.15)*x,             ///
 range(0 41) lpatt(solid)), xtitle(Time since degree (years))       ///
 ytitle(Mean salary)

/*
twoway (function Women =_b[_cons] + _b[yearsdg]*x, range(0 41) lpatt(dash))  ///
 (function Men =_b[_cons] + _b[male] + (_b[yearsdg] + _b[male_years])*x,        ///
 range(0 41) lpatt(solid)), xtitle(Time since degree (years)) ytitle(Mean salary)
*/

regress salary male marketc yearsdg i.male#c.yearsdg

regress salary marketc i.male##c.yearsdg

lincom 1.male + 1.male#c.yearsdg*10

margins r.male, at(yearsdg=10) contrast(nowald effects)

quietly margins male, at(marketc==0 yearsdg=(0(10)40))

marginsplot, noci recast(line) plot1opts(lpatt(dash)) ///
	ytitle(Mean salary) xtitle(Time since degree (years)) title("")


*** 1.9 Dummies for more than two groups

generate associate = rank==2 if !missing(rank)
generate full = rank==3 if !missing(rank)

drop associate full
tabulate rank, generate(r)
rename r2 associate
rename r3 full

regress salary associate full

lincom _cons + associate
lincom full-associate

regress salary i.rank

lincom 3.rank - 2.rank
* contrast r3b2.rank

anova salary rank

regress salary i.male##c.yearsdg marketc i.rank
testparm i.rank

lincom 1.male + 1.male#c.yearsdg*10


*** 1.10 Other types of interactions

generate male_assoc = male*associate
generate male_full = male*full

regress salary male marketc yearsdg associate full male_assoc male_full
testparm male_assoc male_full

regress salary marketc yearsdg i.male##i.rank
testparm i.male#i.rank

generate market_yrs = marketc*yearsdg
regress salary male marketc yearsdg associate full market_yrs

regress salary i.male c.marketc##c.yearsdg i.rank


*** 1.11 Nonlinear effects

generate yearsdg2 = yearsdg^2

regress salary male marketc yearsdg male_years associate full ///
   market_yrs yearsdg2

twoway (function Women = _b[_cons] + _b[yearsdg]*x + _b[yearsdg2]*x^2,      ///
         range(0 41) lpatt(dash))                                           ///
       (function Men = _b[_cons] + _b[male] + (_b[yearsdg] + _b[male_years])*x ///
         + _b[yearsdg2]*x^2, range(0 41) lpatt(solid)),                     ///
        xtitle(Time since degree (years)) ytitle(Mean salary)

/*
regress salary i.male##c.yearsdg c.marketc c.marketc#c.yearsdg i.rank 
    c.yearsdg#c.yearsdg
*/
regress salary (i.male c.marketc c.yearsdg)##c.yearsdg i.rank


*** 1.12 Residual diagnostics

predict res, residuals
histogram res, normal

