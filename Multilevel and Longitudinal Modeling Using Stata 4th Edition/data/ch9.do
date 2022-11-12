*******************************************
****Chapter 9: Crossed random effects*****
*******************************************


*** 9.2 How does investment depend on expected profit and capital stock?

use https://www.stata-press.com/data/mlmus4/grunfeld, clear


*** 9.3 A two-way error-components model

** 9.3.3 Estimation using mixed

mixed I F C || _all: R.fn || yr:, reml stddeviations dfmethod(kroger)

* mixed I F C || _all: R.yr || fn:, reml stddeviations dfmethod(kroger)

estimates store crossed

estat sd, variance post coeflegend
display _b[_all:var(R_fn)]/(_b[_all:var(R_fn)]  ///
                   + _b[yr:var(_cons)] + _b[Residual:var(e)])
				   
estimates restore crossed

** 9.3.4 Prediction

predict firm_eff year_eff, reffects

sort fn yr
list fn firmname yr firm_eff year_eff if yr<1938&fn<5, sepby(fn) noobs

generate reffpart = firm_eff + year_eff
sort fn yr
twoway (line reffpart yr, connect(ascending))                                    ///
   (scatter reffpart yr if yr==1954, msymbol(none) mlabel(firmname) mlabpos(3)), ///
   xtitle(Year) ytitle(Predicted random effects of firm and year)                ///
   xscale(range(1935 1958)) legend(off)

twoway (line year_eff yr if fn==1, sort), xtitle(Year) ///
  ytitle(Predicted random effect for year)

  
*** 9.4 How much do primary and secondary schools affect attainment at age 16?

use https://www.stata-press.com/data/mlmus4/fife, clear


*** 9.5 Data structure

egen pick_comb= tag(pid sid)

egen numsid = sum(pick_comb), by(pid)

sort pid sid
list pid sid numsid if pick_comb==1 & pid<10, sepby(pid) noobs

egen pick_pid = tag(pid)
tabulate numsid if pick_pid==1

egen numpid = total(pick_comb), by(sid)
egen pick_sid = tag(sid)
tabulate numpid if pick_sid==1


*** 9.6 Additive crossed random-effects model

** 9.6.3 Estimation using mixed

mixed attain || _all: R.sid || pid:, reml stddeviations
estimates store additive

estat sd, variance post coeflegend

display _b[_all:var(R_sid)]/(_b[_all:var(R_sid)]  ///
	  + _b[pid:var(_cons)] + _b[Residual:var(e)])

display _b[pid:var(_cons)]/(_b[_all:var(R_sid)]  ///
         +  _b[pid:var(_cons)] + _b[Residual:var(e)])
 
display (_b[_all:var(R_sid)]+_b[pid:var(_cons)]) / ///
	  (_b[_all:var(R_sid)] + _b[pid:var(_cons)] + _b[Residual:var(e)])
	
estimates restore additive

sort pid sid
list pid sid if pid==4, noobs sepby(sid)

estat wcorrelation, at(pid==4)


*** 9.7 Crossed random-effects model with random interaction

** 9.7.3 Estimation using mixed

mixed attain || _all: R.sid || pid:, reml stddeviations

egen comb = group(sid pid)

mixed attain || _all: R.sid || pid: || sid:, reml stddeviations
estimates store interaction

estat sd, variance post coeflegend

display _b[_all:var(R_sid)]/(_b[_all:var(R_sid)]  ///
	  + _b[pid:var(_cons)] + _b[sid:var(_cons)] + _b[Residual:var(e)])
	  
display _b[pid:var(_cons)]/(_b[_all:var(R_sid)]  ///
	  + _b[pid:var(_cons)] + _b[sid:var(_cons)] + _b[Residual:var(e)])
 
display (_b[_all:var(R_sid)]+_b[pid:var(_cons)]  ///
	  + _b[sid:var(_cons)])/(_b[_all:var(R_sid)] ///
 	  + _b[pid:var(_cons)] + _b[sid:var(_cons)] + _b[Residual:var(e)])
	  
estimates restore interaction

sort pid sid
list pid sid if pid==4, noobs sepby(sid)

estat wcorrelation, at(pid==4)

** 9.7.4 Testing variance components

display 3/8*chi2tail(1,281.67) + 3/8*chi2tail(2,281.67) + 1/8*chi2tail(3,281.67) 

lrtest additive interaction

display 1/2*chi2tail(1,280.57) + 1/4*chi2tail(2,280.57)

quietly mixed attain || sid:, reml stddeviations
estimates store secondary
lrtest secondary additive

** 9.7.5 Some diagnostics

estimates restore additive

predict secondary primary, reffects

qnorm secondary if pick_sid, xtitle(Quantiles of normal distribution) ///
  ytitle(Quantiles of empirical Bayes predictions)

qnorm primary if pick_pid, xtitle(Quantiles of normal distribution) ///
  ytitle(Quantiles of empirical Bayes predictions)

  
*** 9.8 A trick requiring fewer random effects

ssc install supclust, replace

supclust pid sid, generate(region)

egen num = count(attain), by(sid pid)
drop if num<3
drop region
supclust sid pid, gen(region)

by region sid, sort: generate f = _n==1
by region: generate sec = sum(f)
table sec

mixed attain || region: R.sec || pid:, reml stddeviations


