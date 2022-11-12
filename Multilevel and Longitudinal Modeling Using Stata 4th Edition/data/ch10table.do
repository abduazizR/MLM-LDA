
*** obtain xtlogit and gllamm estimates (as in Chapter 10)
use https://www.stata-press.com/data/mlmus/toenail, clear

generate trt_month = treatment*month
xtlogit outcome treatment month trt_month, i(patient) quad(30)
estimates store xtlogit

xtrho
xtrhoi -2.945253 4.006325 logit

gllamm outcome treatment month trt_month, i(patient) ///
  link(logit) family(binomial) nip(30) adapt
estimates store gllamm

summ month, detail

*** make prediction dataset
drop _all
set obs 2
generate treatment = _n-1 // treatment=0,1
generate patient = _n
expand 2
by patient, sort: generate month = 3*_n  // month=3,6
generate trt_month = treatment*month

* make all response patterns (to obtain likelihoods=probs of patterns)
generate y00 = 0
generate y11 = 1
generate y10 = 0
generate y01 = 0
by patient (month), sort: replace y10 = _n==1  // 1 at month=3, 0 at month=6
by patient (month), sort: replace y01 = _n==2  // 0 at month=3, 1 at month=6

* stack responses for the four patterns into one response variable y
generate obs = _n
reshape long y, i(obs) j(patt) string

* look at prediction dataset so far
sort patient patt treatment month 
list patient patt treatment month y, sepby(patt)

* treat each patient-pattern combination as one person (2 responses per person)
egen id = group(patt patient)
replace patient = id

* rename response variable to same name as estimation data
generate outcome = y

*** obtain log-likelihood contributions for the "patients"
gllapred ll, ll fsample

*** exponentiate log-likelihood contributions to get predicted probabilities 
* of pairs of responses
generate p = exp(ll)

*** make nice tables of joint probabilities

* create variables y3 for outcomes at month3 and y6 for outcomes at month 6
drop obs trt_month outcome
reshape wide y, i(patient) j(month)
tabulate y3 y6 [iweight=p] if treatment==0  // one obs per cell => weights will be shown
disp .702483654 *.14761442 /(.03541996 * .11450198 )  // OR

tabulate y3 y6 [iweight=p] if treatment==1
disp .75478452*.10620874/(.11535697*.02367046 ) // OR



/* 

compare with xtrhoi by setting month
to 3 for both occasions, and treatment
to 0
 
*/

drop _all
set obs 1
generate treatment = 0
generate patient = 1
expand 2
generate month = 3
generate trt_month = treatment*month
generate y00 = 0
generate y11 = 1
generate y10 = 0
generate y01 = 0
by patient (month), sort: replace y10 = _n==1
by patient (month), sort: replace y01 = _n==2

generate obs = _n
reshape long y, i(obs) j(patt) string

egen id = group(patt patient)
replace patient = id

generate outcome = y

gllapred ll, ll fsample
gen p = exp(ll)
list patt p if obs==1
disp .6637334 *.1879368 /(.0741796 ^2)


gllapred xb, xb fsample
disp xb[1]

xtrhoi -2.7936804 4.0104996 logit

