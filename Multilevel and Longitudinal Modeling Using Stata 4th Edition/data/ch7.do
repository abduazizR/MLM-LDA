**********************************************
****Chapter 7: Growth curve models *****
**********************************************


*** 7.2 How do children grow?

use https://www.stata-press.com/data/mlmus4/asian, clear

** 7.2.1 Observed growth trajectories

label define g 1 "Boy" 2 "Girl"
label values gender g
sort id age

graph twoway (line weight age, connect(ascending)), ///
 by(gender) xtitle(Age in years) ytitle(Weight in Kg)

 
*** 7.3 Models for nonlinear growth

** 7.3.1 Polynomial models

*  Estimation using mixed

mixed weight i.gender age c.age#c.age || id: age, covariance(unstructured) ///
    mle stddeviations vce(robust)
estimates store rc

mixed weight i.gender c.age##c.age || id: c.age##c.age, ///
    covariance(unstructured) mle stddeviations vce(robust)
estimates restore rc

* predicting the mean trajectory

twoway (function Weight=_b[_cons]+_b[age]*x+_b[c.age#c.age]*x^2,  ///
   range(0.1 2.6)), xtitle(Age in years) ytitle(Weight in Kg)     ///
   yscale(range(0 20)) ylab(0(5)20)

estat sd, variance post coeflegend

scalar var_s = _b[id:var(age)]      
scalar var_i = _b[id:var(_cons)]    
scalar cov = _b[id:cov(age,_cons)]
estimates restore rc

twoway (function Weight=_b[_cons]+_b[age]*x+_b[c.age#c.age]*x^2,          /// 
        range(0.1 2.6) lwidth(medium))                                    /// 
        (function upper=_b[_cons]+_b[age]*x+_b[c.age#c.age]*x^2           /// 
        +1.96*sqrt(var_i+2*cov*x+var_s*x^2), range(0.1 2.6) lpatt(dash))  /// 
        (function lower=_b[_cons]+_b[age]*x+_b[c.age#c.age]*x^2           /// 
        -1.96*sqrt(var_i+2*cov*x+var_s*x^2), range(0.1 2.6) lpatt(dash)), /// 
        legend(order(1 "Mean" 2 "95% range")) xtitle(Age in years)        /// 
        ytitle(Weight in kg) yscale(range(0 20)) ylab(0(5)20)

* Predicting trajectories for individual children

predict traj, fitted

twoway (scatter weight age) (line traj age, sort) if gender==1, ///
 by(id, compact legend(off))

*** 7.3.2 Piecewise linear models

* Estimation using mixed

disp _N
set obs 201
replace age = 0.4 in 199
replace age = 0.9 in 200
replace age = 1.5 in 201
mkspline ages1 .4 ages2 .9 ages3 1.5 ages4 = age


twoway (line ages1 age, sort) (line ages2 age, sort)     ///
       (line ages3 age, sort) (line ages4 age, sort),    ///
       xscale(range(0 2.5)) xlabel(0 .4 .9 1.5)
   
mixed weight i.gender ages1 ages2 ages3 ages4 || id: age, ///
    covariance(unstructured) mle stddeviations vce(robust)


* Predicting the mean trajectory

replace gender = 1 in 199/201
predict mean, xb
estimates store piecewise
quietly estat sd, variance post

scalar var_s = _b[id:var(age)]
scalar var_i = _b[id:var(_cons)]
scalar cov = _b[id:cov(age,_cons)]
estimates restore piecewise

generate upper = mean + 1.96*sqrt(var_i+2*cov*age+var_s*age^2) 
generate lower = mean - 1.96*sqrt(var_i+2*cov*age+var_s*age^2) 

twoway (line mean age, sort lwidth(medium)) ///
       (line lower age, sort lpatt(dash)) ///
       (line upper age, sort lpatt(dash)) if gender==2, ///
		legend(order(1 "Mean" 2 "95% range")) xtitle(Age in years) ///
		ytitle(Weight in kg) xtick(.4 .9 1.5, grid) ylabel(,nogrid) ///
		yscale(range(0 20)) ylab(0(5)20)

		
*** 7.4 Two-stage formulation and cross-level interaction

mixed weight i.gender age c.age#i.gender c.age#c.age || ///
     id: age, covariance(unstructured) mle stddeviations vce(robust)

	 
*** 7.5 Heteroskedasticity

** 7.5.1 Heteroskedasticity at level 1

mixed weight i.gender age c.age#c.age || id: age, ///
   covariance(unstructured) mle residuals(independent, by(gender)) ///
   stddeviations vce(robust)
lrtest rc ., force

** 7.5.2 Heteroskedasticity at level 2

mixed weight i.gender age c.age#c.age || id: 1.gender#c.age 1.gender, ///
   noconstant covariance(unstructured) || id: 2.gender#c.age 2.gender, ///
   noconstant covariance(unstructured) mle stddeviations vce(robust)
lrtest rc ., force


*** 7.7 Growth-curve model as a structural equation model

** 7.7.1 Estimation using sem

use https://www.stata-press.com/data/mlmus4/reading, clear

misstable patterns read*

graph box read0 read1 read2 read3, ascategory intensity(0) medtype(line)

sem (read0 <- L1@1 L2@0 ) (read1 <- L1@1 L2@1 ) (read2 <- L1@1 L2@2 ) ///
    (read3 <- L1@1 L2@3 ) , means(L1 L2) noconstant method(mlmv) vce(robust)

/*
sem (L1 -> read0@1 read1@1 read2@1 read3@1) ///
    (L2 -> read1@1 read2@2 read3@3), means(L1 L2) noconstant method(mlmv) ///
	vce(robust)
*/

predict l1 l2, latent

** 7.7.2 Estimation using mixed

reshape long read math age, i(id) j(grade)

quietly xtset id grade
drop if missing(read)
xtdescribe 

egen mn_read = mean(read), by(grade)

twoway (connected mn_read grade, sort), xtitle(Grade) ///
   ytitle(Mean reading score)

/*
graph box read, over(grade) ytitle(Reading score) intensity(0) medtype(line)
*/
  
mixed read grade || id: grade, covariance(unstructured) mle ///
   residual(independent, by(grade)) vce(robust)

predict fixed, xb
twoway (connected mn_read grade, sort lpatt(solid)) ///
       (connected fixed grade, sort lpatt(dash)), xtitle(Grade) ///
      ytitle(Mean reading score) legend(order(1 "Raw mean" 2 "Fitted mean"))

predict zeta2 zeta1, reffects
generate eta1 = _b[_cons] + zeta1
generate eta2 = _b[grade] + zeta2


