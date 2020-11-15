<<dd_version: 2>>
<<dd_include: header.txt >>


Tutorial: causal inference methods made easy for applied resarchers/epidemiologists/statisticians 
=================================================================================================

ICON-LSHTM, LONDON, 15th September 2020

Miguel Angel Luque Fernandez, PhD
	Assistant Professor of Epidemiology

Matthew J. Smith, MSc
	Research Degree Student in Biostatistics
	
Inequalities in Cancer Outcomes Network, LSHTM, London, UK

Copyright (c) 2017 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Bug reports: miguel-angel.luque@lshtm.ac.uk	     



## Contents

#### Introduction

#### 1) DATA: Nomenclature (Globals and Naive approach)

#### 2) Non-parametric G-Formula
###### 2.1) Introduction to the G-Formula
###### 2.2) Non-parametric G-Formula: one confounder
###### 2.3) Non-parametric implementation of the G-Formula using a fully saturated regression model

#### 3) Parametric G-Formula
###### 3.1) G-Computation or parametric G-Formula for one confounder
###### 3.2) G-Computation or parametric G-Formula for more than one confounder

#### 4) Inverse probability weighting (IPW) based on the propensity score
###### 4.1) IPW using Stata's *teffects*
###### 4.2) IPW by hand
###### 4.3) Exploring violations of the identifiability conditions

#### 5) Inverse Probability Weighting plus regression adjustment (IPWRA)
###### 5.1) First, let Stata calculate the IPWRA for the ATE
###### 5.2) Now, by hand

#### 6) Augmented Inverse Probability Weighting plus regression adjustment (AIPWRA)
###### 6.1) First, compute AIPWRA using Stata commands
###### 6.2) Now compute by hand

#### 7) Ensemble Learning Targeted Maximum Likelihood Estimation

#### 8) Simulation
###### 8.1) Data generation
###### 8.2) ATE estimation
###### 8.3) Results


# Introduction

In this tutorial we make complex causal inference methodology easily understood and applied.

Assumptions: 
1. Conditional exchageability
2. Positivity
3. Consistency



## 1) DATA: Nomenclature (Globals and Naive approach)

Our data is...

~~~~
<<dd_do>>
clear
set more off
use "~/Dropbox/ESTIMATORSCIproject/R_Stata_master_files/Data/rhc.dta", clear
describe
count
set seed 1234
<</dd_do>>
~~~~

We now define the outcome, exposure (or treatment), the main confounder, and other confounders

~~~~
<<dd_do>>
global Y death_d30 // Outcome
global A rhc // Exposure or treatment
global C gender // Confounder 
global W gender age edu race carcinoma // Confounders
<</dd_do>>
~~~~

A naive approach to assessing the main association can be performed using the *regress* command.

~~~~
<<dd_do>>
regress $Y $A $C // Naive approach
<</dd_do>>
~~~~

We can naively say there is strong evidence (p<0.001) that the risk of death within 30 days is 7% higher amongst those with RHC, adjusting for gender.
 
 
## 2) Non-parametric G-Formula

Before jumping straight into the black box of Stata's *teffects* command, we first estimate the Average Treatment Effect (ATE) by hand. 

### 2.1) Introduction to the G-Formula

The G-Formula is defined as

$$ E(Y^{a,c}) = \sum_{C} E(Y|A=a,C=c)P(C=c) $$

From the G-Formula, we estimate four possible combinations:

1. $ E(Y^{0,0}) $

2. $ E(Y^{0,1}) $

3. $ E(Y^{1,0}) $

4. $ E(Y^{1,1}) $ 

For example (and in terms of our data), $ E(Y^{0,0}) $ is obtained by summing the expected deaths amongst those without RHC and are female (weighted by the marginal probability of being female) and the expected deaths amongst those without RHC and are male (weighted by the marginal probability of being male). In terms of the G-Formula, this potential outcome is

$$ E(Y^{0,0}) = E(Y|A_{1}=0,C_{1}=0)P(C_{1}=0) + E(Y|A_{1}=0,C_{1}=1)P(C_{1}=1) $$

and, by the same process, the remaining combinations can be obtained:

2. $ E(Y^{0,1}) = E(Y|A_{1}=0,C_{1}=1)P(C_{1}=1) + E(Y|A_{1}=1,C_{1}=1)P(C_{1}=1) $

3. $ E(Y^{1,0}) = E(Y|A_{1}=1,C_{1}=0)P(C_{1}=0) + E(Y|A_{1}=0,C_{1}=0)P(C_{1}=0) $

4. $ E(Y^{1,1}) = E(Y|A_{1}=1,C_{1}=1)P(C_{1}=1) + E(Y|A_{1}=1,C_{1}=0)P(C_{1}=0) $ 

If the identifying assumptions hold, the potential outcomes can be estimated and the difference (i.e. $ Y^{a=1} - Y^{a=0} $ ) can be given a causal interpretation.


### 2.2) Non-parametric G-Formula: one confounder

To do this, first compute the marginal probability for confounder C (i.e. marginal probability of being married)

~~~~
<<dd_do>>
proportion $C
matrix m=e(b)
gen genderf = m[1,1]
sum genderf
gen genderm = m[1,2]
sum genderm
<</dd_do>>
~~~~

We see that, in this sample, the marginal probability of being male $ (C=1) $ is roughly 56%, thus 44% are female $ (C=0) $. 

Now we compute the G-Formula by hand (generalisation or standardisation)

~~~~
<<dd_do>>
ssc install sumup
sumup $Y, by($A $C)
matrix y00 = r(Stat1)
matrix y01 = r(Stat2)
matrix y10 = r(Stat3)
matrix y11 = r(Stat4)
gen ATE = ((y11[3,1]-y01[3,1]))*genderm + ((y10[3,1]-y00[3,1]))*genderf
qui: sum ATE
display "The ATE is: "  "`r(mean)'"
drop ATE 
<</dd_do>>
~~~~



In the code above we have estimated the four possible combinations:

1. $ E(Y^{0,0}) = y00 = 0.301772 $

2. $ E(Y^{0,1}) = y01 = 0.310345 $

3. $ E(Y^{1,0}) = y10 = 0.384106 $

4. $ E(Y^{1,1}) = y11 = 0.377152 $

We are now ready to estimate the average treatment effect of RHC on death within 30 days.  

Plugging the combinations into the formula for the average treatment effect

$$ ATE = [E(Y|A=1,C=1)-E(Y|A=0,C=1)]P(C=1) + [E(Y|A=1,C=0)-E(Y|A=0,C=0)]P(C=0) $$

we obtain $ ATE = 0.07369 $

Our naive regression analysis estimate (0.07352 grams) is a slight underestimate of the ATE.


Fortunately, Stata contains a command that calculates the ATE for us. Stata implements the results using the *teffects* command. We obtain the same value of the ATE from our by-hand approach.

~~~~
<<dd_do>>
teffects ra ($Y $C) ($A)
<</dd_do>>
~~~~            
            
We have calculated the ATE amongst those who had RHC or not, but often it is of public health interest to ask "what is the effect of having RHC for those who chose to do so?" In other words, we are interested in the average treatment effect amongst the treated (ATET), that is amongst only those who had RHC. This comes with a slightly weaker assumption of conditional exchangeability such that the potential outcome for the expected deaths amongst those who (possibly contrary to fact) did not have RHC is independent of the exposure given the confounder.

We can check the ATET, first by hand

~~~~
<<dd_do>>
proportion $C if $A==1
matrix m=e(b)
gen genderfatet = m[1,1]
gen gendermatet = m[1,2]
gen ATT = ((y11[3,1]-y01[3,1]))*gendermatet + ((y10[3,1]-y00[3,1]))*genderfatet
qui: sum ATT
display "The ATT is: "  "`r(mean)'"
drop ATT
<</dd_do>>
~~~~
            
We obtain $ ATET = 0.07325 $

We can check the Stata implementation results using the *teffects* command. We obtain the same value of the ATET from our by-hand approach

~~~~
<<dd_do>>
teffects ra ($Y $C) ($A), atet
<</dd_do>>
~~~~

again, we obtain $ ATET = 0.07325 $.


Let's now bootstrap the confidence intervals for the ATE

~~~~
<<dd_do>>
program define ATE, rclass
   capture drop y1
   capture drop y0
   capture drop ATE
   sumup $Y, by($A $C)
   matrix y00 = r(Stat1)
   matrix y01 = r(Stat2)
   matrix y10 = r(Stat3)
   matrix y11 = r(Stat4)
   gen ATE = ((y11[3,1]-y01[3,1]))*genderm + ((y10[3,1]-y00[3,1]))*genderf
   qui sum ATE
   return scalar ate = `r(mean)'
end
qui bootstrap r(ate), reps(1000): ATE
estat boot, all	
program drop ATE
drop ATE
<</dd_do>>
~~~~     

Let's now bootstrap the confidence intervals for the ATET

~~~~
<<dd_do>>
program define ATT, rclass
   capture drop y1
   capture drop y0
   capture drop ATT
   sumup $Y, by($A $C)
   matrix y00 = r(Stat1)
   matrix y01 = r(Stat2)
   matrix y10 = r(Stat3)
   matrix y11 = r(Stat4)
   gen ATT = ((y11[3,1]-y01[3,1]))*gendermatet + ((y10[3,1]-y00[3,1]))*genderfatet
   qui sum ATT
   return scalar atet = `r(mean)'
end
qui bootstrap r(atet), reps(1000): ATT
estat boot, all
program drop ATT
drop ATT
<</dd_do>>
~~~~

### 2.3) Non-parametric implementation of the G-Formula using a fully saturated regression model

~~~~
<<dd_do>>
regress $Y ibn.$A ibn.$A#c.($C) , noconstant vce(robust)
margins $A , vce(unconditional) // Marginal probability for C
margins r.$A , contrast(nowald) 
<</dd_do>>
~~~~

Note about *margins* Unconditional

~~~~
<<dd_do>>
regress $Y ibn.$A ibn.$A#c.($C), noconstant vce(robust) coeflegend
predictnl ATE = (_b[1bn.rhc] + _b[1bn.rhc#c.gender]*gender) - (_b[0bn.rhc] + _b[0bn.rhc#c.gender]*gender)
qui: sum ATE
display "The ATE is:  " "`r(mean)'"
drop ATE
<</dd_do>>
~~~~   	
  

## 3) Parametric G-Formula

### 3.1) G-Computation or parametric G-Formula for one confounder

~~~~
<<dd_do>>
teffects ra ($Y $C) ($A) //Parametric G-Formula implementation in Stata
regress $Y $C if $A==1
predict double y1hat
regress $Y $C if $A==0
predict double y0hat
mean y1hat y0hat
nlcom _b[y1hat] - _b[y0hat]
<</dd_do>>
~~~~
			
Bootstrapping to get the SE and 95%CI for the ATE 

~~~~
<<dd_do>>
capture program drop ATE
program define ATE, rclass
   capture drop y1
   capture drop y0
   reg $Y $C if $A==1 
   predict double y1, xb
   quiet sum y1
   reg $Y $C if $A==0 
   predict double y0, xb 
   quiet sum y0
   mean y1 y0 
   lincom _b[y1]-_b[y0]
   return scalar ace =`r(estimate)'
end
qui bootstrap r(ace), reps(1000): ATE
estat boot, all
<</dd_do>>
~~~~
			
### 3.2) G-Computation or parametric G-Formula for more than one confounder

Implement ATE 

~~~~
<<dd_do>>
teffects ra ($Y $W) ($A)
<</dd_do>>
~~~~

Perform non-parametric G-Formula with multiple variables

~~~~
<<dd_do>>
regress $Y ibn.$A ibn.$A#c.($W) , noconstant vce(robust)
margins $A, vce(unconditional)
margins r.$A, contrast(nowald)
<</dd_do>>
~~~~
	
G-Computation or parametric G-Formula computation

~~~~
<<dd_do>>
regress $Y $W if $A==1
predict double y1hat 
regress $Y $W if $A==0
predict double y0hat
mean y1hat y0hat
lincom _b[y1hat] - _b[y0hat]
<</dd_do>>
~~~~

Bootstrap the 95% confidence intervals for the ATE

~~~~
<<dd_do>>
capture program drop ATE
program define ATE, rclass
   capture drop y1
   capture drop y0
   reg $Y $W if $A==1 
   predict double y1, xb
   quiet sum y1
   reg $Y $W if $A==0 
   predict double y0, xb 
   quiet sum y0
   mean y1 y0 
   lincom _b[y1]-_b[y0]
   return scalar ace =`r(estimate)'
end
qui bootstrap r(ace), reps(1000): ATE dots
estat boot, all
<</dd_do>>
~~~~

Possibilities of G-Computation in Stata

Potential outcomes...

~~~~
<<dd_do>>
teffects ra ($Y $W) ($A), aequations
<</dd_do>>
~~~~

Relative risk

~~~~
<<dd_do>>
teffects ra ($Y $W) ($A), coeflegend
<</dd_do>>
~~~~

???

~~~~
<<dd_do>>
nlcom 100*_b[ATE:r1vs0.$A]/_b[POmean:0.$A]
<</dd_do>>
~~~~	


## 4) Inverse probability weighting (IPW) based on the propensity score
### 4.1) IPW using Stata's *teffects*

~~~~
<<dd_do>>
teffects ipw ($Y) ($A $W, logit), nolog vsquish
<</dd_do>>
~~~~

### 4.2) IPW by hand

Model for the propensity score (note here the possibility of potential misspecification)

~~~~
<<dd_do>>
logit $A $W, vce(robust) nolog
predict double ps
<</dd_do>>
~~~~

Sampling weights based on Horvitz-Thompson estimator

~~~~
<<dd_do>>
generate double ipw1 = ($A==1)/ps
regress $Y [pw=ipw1]
generate double ipw0 = ($A==0)/(1-ps)
regress $Y [pw=ipw0]
<</dd_do>>
~~~~

Bootstrap the SE and 95%CI for the ATE

~~~~
<<dd_do>>
program drop ATE
program define ATE, rclass
  capture drop y1
  capture drop y0
  regress $Y [pw=ipw1]
  matrix y1 = e(b)
  gen double y1 = y1[1,1]
  regress $Y [pw=ipw0]
  matrix y0 = e(b)
  gen double y0 = y0[1,1]
  mean y1 y0
  lincom _b[y1]-_b[y0]
  return scalar ace = `r(estimate)'
end
qui bootstrap r(ace), reps(1000): ATE
estat boot, all
<</dd_do>>
~~~~
    
Test the balance after adjustment

~~~~
<<dd_do>>
qui teffects ipw ($Y) ($A $W)   
tebalance summarize
<</dd_do>>
~~~~

### 4.3) Exploring violations of the identifiability conditions

Check for any violations from the positivity assumption by using overlap plots

~~~~
<<dd_do>>
qui: teffects ipw ($Y) ($A $W, logit), nolog vsquish
teffects overlap
<</dd_do>>
~~~~
<<dd_graph:sav(IPWoverlap1.png) replace>>
    
Check the overlap plots by hand

~~~~
<<dd_do>>
sort $A
by $A: summarize ps
kdensity ps if $A==1, generate(x1pointsa d1A) nograph n(10000)
kdensity ps if $A==0, generate(x0pointsa d0A) nograph n(10000)
label variable d1A "density for RHC=1"
label variable d0A "density for RHC=0"
twoway (line d0A x0pointsa , yaxis(1))(line d1A x1pointsa, yaxis(2))
<</dd_do>>
~~~~
<<dd_graph:sav(IPWoverlap2.png) replace>>



## 5) Inverse Probability Weighting plus regression adjustment (IPWRA)

### 5.1) First, let Stata calculate the IPWRA for the ATE

~~~~
<<dd_do>>
teffects ipwra ($Y $W) ($A $W)
<</dd_do>>
~~~~    
  
~~~~
<<dd_do>>
nlcom 100*_b[r1vs0.$A]/_b[POmean:0.$A]
<</dd_do>>
~~~~
  
  
### 5.2) Now, by hand

Step one (inverse probability weights):

~~~~
<<dd_do>>
logit $A $W
predict double dps, pr
gen ipw = .
replace ipw=($A==1)/dps if $A==1
replace ipw=($A==0)/(1-dps) if $A==0
sum ipw 
<</dd_do>>
~~~~
    
Stabilised weights

~~~~
<<dd_do>>
logit $A, vce(robust) nolog
predict double nps, pr
gen sws = .
replace sws = nps/dps if $A==1
replace sws = (1-nps)/(1-dps) if $A==0
sum sws
<</dd_do>>
~~~~

Step two (G-Computation):

~~~~
<<dd_do>>
reg $Y $W if $A==1 [pw=sws]
drop y1
predict double y1, xb
reg $Y $W if $A==0 [pw=sws]
drop y0
predict double y0, xb 
mean y1 y0
<</dd_do>>
~~~~
		

IPWRA (Problem SE) we need bootstrap

~~~~
<<dd_do>>
nlcom(_b[y1]-_b[y0])
<</dd_do>>
~~~~		


Bootstrap SE and 95% CI for the ATE (Stabilized weights)

~~~~
<<dd_do>>
capture program drop ACE
program define ACE, rclass
  capture drop y1
  capture drop y0
  reg $Y $W if $A==1 [pw=sws]
  predict double y1, xb
  quiet sum y1
  return scalar y1=`r(mean)'
  reg $Y $W if $A==0 [pw=sws]
  predict double y0, xb 
  quiet sum y0
  return scalar y0=`r(mean)'
  mean y1 y0 
  lincom _b[y1]-_b[y0]
  return scalar ace =`r(estimate)'
end
qui bootstrap r(ace), reps(1000): ACE
estat boot, all
<</dd_do>>
~~~~
	
Check that we get the same outcome as we did with our approach by hand

~~~~
<<dd_do>>
teffects ipwra ($Y $W) ($A $W, logit)
<</dd_do>>
~~~~
		

## 6) Augmented Inverse Probability Weighting plus regression adjustment (AIPWRA)

### 6.1) First, compute AIPWRA using Stata commands

~~~~
<<dd_do>>
teffects aipw ($Y $W) ($A $W, logit)
nlcom 100*_b[r1vs0.$A]/_b[POmean:0.$A]
<</dd_do>>
~~~~
		
### 6.2) Now compute by hand

Define the global variables

~~~~
<<dd_do>>
qui sum $Y
gen double Y = ($Y -`r(min)')/(`r(max)' - `r(min)')
qui sum $Y 
global cf = (`r(max)' - `r(min)')
<</dd_do>>
~~~~		
		
Step 1: prediction model for the outcome Q0 (g-computation)

~~~~
<<dd_do>>
qui glm Y $A $W, fam(bin) 
predict double QAW, xb
qui glm Y $W if $A==1, fam(bin) 
predict double Q1W, mu
qui glm Y $W if $A==0, fam(bin)
predict double Q0W, mu
<</dd_do>>
~~~~
		
Step 2: prediction model for the treatment g0 (IPTW): Inverse probability weights (Stabilized weights)

~~~~
<<dd_do>>
qui logit $A $W
drop dps
predict double dps, pr
qui logit $A
drop nps
predict double nps, pr
drop sws
gen sws = .
replace sws = nps/dps if $A==1
replace sws = (1-nps)/(1-dps) if $A==0
sum sws
<</dd_do>>
~~~~
		
Step 3: Estimation Equation

~~~~
<<dd_do>>
gen double ATE = (sws*(Y-QAW) + (Q1W - Q0W)*$cf)
sum ATE
<</dd_do>>
~~~~
		
Marginal Structural Model

~~~~
<<dd_do>>
glm $Y $A [pw=sws],nolog
<</dd_do>>
~~~~

Bootstrap the 95% confidence intervals for the AIPTW

~~~~
<<dd_do>>
drop ATE
capture program drop ATE
program define ATE, rclass
  capture drop ATE
  gen double ATE = (sws*(Y-QAW) + (Q1W - Q0W)*$cf) // Estimation Equation
  mean ATE
  lincom _b[ATE]
  return scalar ate =`r(estimate)'
end
qui bootstrap r(ate), reps(1000): ATE
estat boot, all
<</dd_do>>
~~~~


## 7) Ensemble Learning Targeted Maximum Likelihood Estimation

~~~~
<<dd_do>>
eltmle $Y $A $W, tmle   // install via 'help eltmle'
shell "C:\Program Files\R\R-3.6.3\bin\x64\R.exe" CMD BATCH SLS.R
<</dd_do>>
~~~~		
		

## 8) Simulation

### 8.1) Data generation

Let's generate some data...

Variables are defined as:
w1=Socieconomic status, w2=Age, w3=Cancer stage, w4=Comorbidities

Note: the interaction between age and comorbidities

~~~~
<<dd_do>>
clear
set obs 1000
set seed 777
gen w1    = round(runiform(1, 5)) //Quintiles of Socioeconomic Deprivation
gen w2    = rbinomial(1, 0.45) //Binary: probability age >65 = 0.45
gen w3    = round(runiform(0, 1) + 0.75*(w2) + 0.8*(w1)) //Stage 
recode w3 (5/6=1)       //Stage (TNM): categorical 4 levels
gen w4    = round(runiform(0, 1) + 0.75*(w2) + 0.2*(w1)) //Comorbidites: categorical four levels
gen A     = (rbinomial(1,invlogit(-1 -  0.15*(w4) + 1.5*(w2) + 0.75*(w3) + 0.25*(w1) + 0.8*(w2)*(w4)))) //Binary treatment
gen Y     = (rbinomial(1,invlogit(-3 + A + 0.25*(w4) + 0.75*(w3) + 0.8*(w2)*(w4) + 0.05*(w1)))) //Binary outcome
gen Yt1 = (invlogit(-3 + 1 + 0.25*(w4) + 0.75*(w3) + 0.8*(w2)*(w4) + 0.05*(w1))) // Potential outcome 1
gen Yt0 = (invlogit(-3 + 0 + 0.25*(w4) + 0.75*(w3) + 0.8*(w2)*(w4) + 0.05*(w1))) // Potential outcome 2
gen psi = Yt1-Yt0  // True ATE
<</dd_do>>
~~~~

Estimate the true simulated ATE (ATE = 17.8%)

~~~~
<<dd_do>>
mean psi
<</dd_do>>
~~~~

### 8.2) ATE estimation

Regression adjustment ATE = 24.1%

~~~~
<<dd_do>>
teffects ra (Y w1 w2 w3 w4) (A)
estimates store ra
<</dd_do>>
~~~~

Inverse probability weighting ATE = 37.9%

(Note: we will also store the estimates so that we can compare between models)

~~~~
<<dd_do>>
teffects ipw (Y) (A w1 w2 w3 w4)
estimates store ipw
<</dd_do>>
~~~~

IPTW-RA ATE = 31.9%

(Note: we will also store the estimates so that we can compare between models)

~~~~
<<dd_do>>
teffects ipwra (Y w1 w2 w3 w4) (A w1 w2 w3 w4)
estimates store ipwra
<</dd_do>>
~~~~

AIPTW ATE = 28.8%

(Note: we will also store the estimates so that we can compare between models)

~~~~
<<dd_do>>
teffects aipw (Y w1 w2 w3 w4) (A w1 w2 w3 w4)
estimates store aipw
<</dd_do>>
~~~~

### 8.3) Results

~~~~
<<dd_do>>
qui reg psi
estimates store psi
estout psi ra ipw ipwra aipw
<</dd_do>>
~~~~

Ensemble Learning Maximum Likelihood Estimation (ELTMLE) ATE = 22.1%

~~~~
<<dd_do>>
eltmle Y A w1 w2 w3 w4, tmle
shell "C:\Program Files\R\R-3.6.3\bin\x64\R.exe" CMD BATCH SLS.R 
<</dd_do>>
~~~~

Relative bias IPTW ( = 0.52631579)

~~~~
<<dd_do>>
display (0.38 - 0.18)/.38 
<</dd_do>>
~~~~

Relative bias AIPTW ( = 0.35714286)

~~~~
<<dd_do>>
display (0.28 - 0.18)/0.28
<</dd_do>>
~~~~

Relative bias ELTMLE ( = 0.18181818)

~~~~
<<dd_do>>
display (0.22 - 0.18)/0.22
<</dd_do>>
~~~~




dyndoc StataCode.do, replace

		
			
