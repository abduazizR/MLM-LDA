**********************************************
****Chapter 6: Marginal models *****
**********************************************


*** 6.3 Covariance structures

** 6.3.1 Unstructured covariance matrix

use https://www.stata-press.com/data/mlmus4/wagepan, clear

generate educt = educ - 12
generate yeart = year - 1980

mixed lwage black hisp union married exper yeart educt || nr:, ///
    noconstant residuals(unstructured, t(yeart)) nofetable nogroup mle ///
    stddeviations
estimates store un

estat wcorrelation

** 6.3.2. Random-intercept or compound-symmetric/exchangeable structure

quietly mixed lwage black hisp union married exper yeart educt || nr:, ///
   noconstant residuals(exchangeable) mle stddeviations
estimates store ri

** 6.3.3 Random-coefficient structure

mixed lwage black hisp union married exper yeart educt || nr: yeart, ///
   covariance(unstructured) nofetable nogroup mle stddeviations
estimates store rc

estat wcorrelation

** 6.3.4 Autoregressive and exponential structures  

mixed lwage black hisp union married exper yeart educt || nr:, noconstant ///
      residuals(ar 1, t(yeart)) nofetable nogroup mle stddeviations
estimates store ar1

estat wcorrelation

** 6.3.5. Moving-average residual structure

mixed lwage black hisp union married exper yeart educt || nr:, ///
   noconstant residuals(ma 1, t(yeart)) nofetable nogroup mle stddeviations
estimates store ma1

estat wcorrelation

** 6.3.6 Banded and Toeplitz structures

mixed lwage black hisp union married exper yeart educt || nr:, noconstant ///
     residuals(banded 1, t(yeart)) nofetable nogroup mle stddeviations
estimates store ba1

estat wcorrelation

mixed lwage black hisp union married exper yeart educt || nr:, noconstant ///
     residuals(toeplitz 2, t(yeart)) nofetable nogroup mle stddeviations
estimates store to2
	 
estat wcorrelation


*** 6.4 Hybrid and complex marginal models

** 6.4.1 Random effects and correlated level-1 residuals

mixed lwage black hisp union married exper yeart educt || nr:, ///
   residuals(ar 1, t(yeart)) nofetable nogroup mle stddeviations
estimates store ri_ar1

estat wcorrelation

** 6.4.2 Heteroskedastic level-1 residuals over occasions

mixed lwage black hisp union married exper yeart educt || nr: yeart, ///
    covariance(unstructured) residuals(independent, by(yeart)) nofetable nogroup  ///
    mle stddeviations
estimates store rc_het

estat wcorrelation

** 6.4.3 Heteroskedastic level-1 residuals over groups

generate ethnic = black*1 + hisp*2

mixed lwage black hisp union married exper yeart educt || nr:, ///
    residuals(independent, by(ethnic)) nofetable nogroup mle stddeviations
estimates store ri_het2

table ethnic, statistic(min nr)

* white
estat wcorrelation, at(nr=13)

* Black
estat wcorrelation, at(nr=383)

* Hispanic
estat wcorrelation, at(nr=1142)

** 6.4.4 Different covariance matrices over groups

mixed lwage black hisp union married exper yeart educt || nr:, noconstant ///
    residuals(ar 1, t(yeart) by(ethnic)) nofetable nogroup mle stddeviations

/*
estat wcorrelation, at(nr=13)
estat wcorrelation, at(nr=383)
estat wcorrelation, at(nr=1142)
*/


*** 6.5 Comparing the fit of marginal models

lrtest un ma1

estimates stats un ri rc ar1 ma1 ba1 to2 ri_ar1 rc_het

estimates table un ri rc ar1 ma1 ba1 to2 ri_ar1 rc_het, ///
  b(%4.3f) se(%4.3f) keep(lwage:)

  
*** 6.6 Generalized Estimating Equations (GEE)

quietly xtset nr yeart
xtgee lwage black hisp union married exper yeart educt, corr(ar 1) vce(robust)

/*
mixed lwage black hisp union married exper yeart educt || nr:, noconstant ///
   residuals(ar 1, t(yeart)) vce(robust) stddeviations
*/

estat wcorrelation, format(%4.3f)

/*
xtset nr yeart
xtreg lwage black hisp union married exper yeart educt, pa corr(ar 1) vce(robust)
matrix list e(R)
*/


*** 6.7 Marginal moceling with few units and many occasions

** 6.7.3 Fitting marginal models for long panels in Stata

use https://www.stata-press.com/data/mlmus4/garrett, clear
xtset country year
xtpcse gdp oild demand leftpow organ c.leftpow#c.organ
xtpcse gdp oild demand leftpow organ c.leftpow#c.organ, correlation(ar1)

