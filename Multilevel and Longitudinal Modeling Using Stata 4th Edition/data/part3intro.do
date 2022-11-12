*****************************************************************
****Part III: Models for Longitudinal and panel data        *****
*****************************************************************


*** How and why do wages change over time?

use https://www.stata-press.com/data/mlmus4/wagepan, clear


*** Longitudinal data structure and descriptives

** long and wide form

keep nr lwage black hisp union married exper year educ

reshape wide lwage union married exper, i(nr) j(year)
format lwage* %5.3f
list nr lwage* in 1/5, clean noobs abbreviate(6)
reshape long lwage union married exper, i(nr) j(year)

list nr year lwage in 1/8, clean noobs

** Declaring data as longitudinal using stset

xtset nr year

** Missing data

quietly xtset nr year
xtdescribe if !missing(lwage)

** Time-varying and time-constant variables

quietly xtset nr year
xtsum lwage union educ year

quietly xtset nr year
xttab union


*** Graphical displays for longitudinal data

graph box lwage, over(year) intensity(0) medtype(line)    ///
   marker(1,mlabel(nr) mlabsize(vsmall) msym(i) mlabpos(0) ///
   mlabcol(black)) ytitle(Log hourly wage)

sort nr year
set seed 132152

egen pick_occ = tag(nr)
generate r = runiform() if pick_occ==1

egen num = rank(r) if r<.

egen number = mean(num), by(nr)

twoway line lwage year if number<=12, by(nr, compact)     ///
  ytitle(Log hourly wage) xtitle(Year) xlabel(,angle(45)) 

egen mn_lwage = mean(lwage), by(year)

sort nr year

twoway (scatter lwage year, jitter(2) msym(o) msize(tiny)) ///
  (line lwage year if number<=12, connect(ascending)       ///
  lwidth(vthin) lpatt(solid))                              ///
  (line mn_lwage year, sort lpatt(longdash)) if lwage>-2,  ///
  ytitle(Log hourly wage) xtitle(Year) ///
  legend(order(2 "Individual trajectories" 3 "Mean trajectory"))

  
*** Pooled ordinary least squares

generate educt = educ - 12
generate yeart = year - 1980

regress lwage black hisp union married exper yeart educt, vce(cluster nr)


*** Correlated rsiduals

predict res, residuals

* Variances and correlations (run 6 lines in one block)
preserve
keep nr res yeart
reshape wide res, i(nr) j(year)
tabstat res*, statistics(sd) format(%3.2f) 
correlate res*, wrap
restore

