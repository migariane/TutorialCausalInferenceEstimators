/* 
Tutorial: causal inference methods made easy for applied resarchers/epidemiologists/statisticians 
=================================================================================================

ICON-LSHTM, LONDON, 6th December 2017

Miguel Angel Luque Fernandez, PhD

Assistant Professor of Epidemiology

Inequalities in Cancer Outcomes Network, LSHTM, London, UK

Copyright (c) 2017 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Bug reports: miguel-angel.luque@lshtm.ac.uk	     
*/
  
/* Preliminaries */
   clear
    set more off
    cd "C:\Data" 
    use "C:\Data\rhc.dta", clear
    describe
    count
    
 /* Naive estimate of the ATE */
    // Define the outcome (Y), exposure (A) and confounders (C)
    global Y death_d30 // Outcome
    global A rhc // Exposure or treatment
    global C gender // Confounder 
    global W gender age edu race carcinoma // Confounders
    // Naive approach to estimate the causal effect
    regr $Y $A $C 	
    // The naive estimate of the causal effect is 0.07352
    
/* Non-parametric G-formula */
    // Non-parametric g-formula of the (1) ATE, and (2) ATET
        * 1) ATE
            * Estimate the marginal probabilities 
            proportion $C
            matrix m=e(b)
            gen genderf = m[1,1]
            sum genderf
            gen genderm = m[1,2]
            sum genderm
            * Compute the g-formula by hand
            sumup $Y, by($A $C)
            matrix y00 = r(Stat1)
            matrix y01 = r(Stat2)
            matrix y10 = r(Stat3)
            matrix y11 = r(Stat4)
            gen ATE = ((y11[3,1]-y01[3,1]))*genderm + ((y10[3,1]-y00[3,1]))*genderf
            qui: sum ATE
            display "The ATE is: "  "`r(mean)'"
            drop ATE 
            // The ATE from non-parametric estimator is:  0.073692

            * Check that Stata "teffects" command obtains the same estimate
            teffects ra ($Y $C) ($A)
            // The ATE from "teffects" implementation is:  0.073692
    	
        * 2) ATET
            * First, by hand...
            * Estimate the marginal probabilities
            proportion $C if $A==1
            matrix m=e(b)
            gen genderfatet = m[1,1]
            gen gendermatet = m[1,2]
            gen ATT = ((y11[3,1]-y01[3,1]))*gendermatet + ((y10[3,1]-y00[3,1]))*genderfatet
            qui: sum ATT
            display "The ATT is: "  "`r(mean)'"

            * Check using Stata "teffects" command
            teffects ra ($Y $C) ($A), atet	
            drop ATT
    
    // Bootstrap the confidence intervals
        * 1) For the ATE
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
               return scalar ate = "`r(mean)'"
            end
            
            qui bootstrap r(ate), reps(1000): ATE
            estat boot, all	
    	
        * 2) For the ATET
            program drop ATE
            program define ATE, rclass
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
               return scalar atet = "`r(mean)'"
            end
            
            qui bootstrap r(atet), reps(1000): ATE
            estat boot, all
            	
            drop ATE
    	
    // Non-parametric implementation of the g-formula using a fully saturated regression model
            regress $Y ibn.$A ibn.$A#c.($C) , noconstant vce(robust)
            margins $A , vce(unconditional) // Marginal probability for C
            margins r.$A , contrast(nowald) 
            regress $Y ibn.$A ibn.$A#c.($C), noconstant vce(robust) coeflegend
            predictnl ATE = (_b[1bn.rhc] + _b[1bn.rhc#c.gender]*gender) - (_b[0bn.rhc] + _b[0bn.rhc#c.gender]*gender)
            qui: sum ATE
            display "The ATE is:  " "`r(mean)'"
            drop ATE
			
			
/* Parametric G-formula */		
	// G-computation (or parametric G-formula) for one confounder
		teffects ra ($Y $C) ($A) //Parametric G-Formula implementation in Stata
		regress $Y $C if $A==1
		predict double y1hat
		predict double y0hat

	// Bootstrap to get the Standard errors and 95% confidence interval for the ATE
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

	// G-computation (or parametric G-formula) for more than one confounder
		teffects ra ($Y $W) ($A)
		regress $Y ibn.$A ibn.$A#c.($W) , noconstant vce(robust)
		margins $A, vce(unconditional)
		margins r.$A, contrast(nowald)
		regress $Y $W if $A==1
		predict double y1hat 
		regress $Y $W if $A==0
		predict double y0hat
		mean y1hat y0hat
		lincom _b[y1hat] - _b[y0hat]

	// Bootstrap the 95% confidence intervals for the ATE
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

	// Other potential outcomes
		teffects ra ($Y $W) ($A), aequations
		teffects ra ($Y $W) ($A), coeflegend
		nlcom 100*_b[ATE:r1vs0.mbsmoke]/_b[POmean:0.mbsmoke]			

		
/* Inverse probability weighting plus regression adjustment */		
	// IPTW based on the propensity score
		teffects ipw ($Y) ($A $W, logit), nolog vsquish

	// IPTW by hand
		logit $A $W, vce(robust) nolog
		predict double ps
		generate double ipw1 = ($A==1)/ps
		regress $Y [pw=ipw1]
		generate double ipw0 = ($A==0)/(1-ps)
		regress $Y [pw=ipw0]
    
    // Bootstrap the confidence intervals
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
    
    // Test the balance
		qui teffects ipw ($Y) ($A $W)   
		tebalance summarize
    
	// Check violations of identifiability conditions
		qui: teffects ipw ($Y) ($A $W, logit), nolog vsquish
		teffects overlap
		sort $A
		by $A: summarize ps
		kdensity ps if $A==1, generate(x1pointsa d1A) nograph n(10000) (n() set to 4642)
		kdensity ps if $A==0, generate(x0pointsa d0A) nograph n(10000) (n() set to 4642)
		label variable d1A “density for mbsmoke=1”
		label variable d0A “density for mbsmoke=0”
		twoway (line d0A x0pointsa , yaxis(1))(line d1A x1pointsa, yaxis(2))
    
    
	// IPTW with regression adjustment
		teffects ipwra ($Y $W) ($A $W)
		nlcom 100*_b[r1vs0.mbsmoke]/_b[POmean:0.mbsmoke]
    
    // IPTW-RA by hand
		
		* Step (i):  IPTW
		logit $A $W
		predict double dps, pr
		gen ipw = .
		replace ipw=($A==1)/dps if $A==1
		replace ipw=($A==0)/(1-dps) if $A==0
		sum ipw 
		
		* Stabilised weights
		logit $A, vce(robust) nolog
		predict double nps, pr
		gen sws = .
		replace sws = nps/dps if $A==1
		replace sws = (1-nps)/(1-dps) if $A==0
		sum sws
		
		* Step (ii) G-Computation
		reg $Y $W if $A==1 [pw=sws]
		drop y1
		predict double y1, xb
		reg $Y $W if $A==0 [pw=sws]
		drop y0
		predict double y0, xb 
		mean y1 y0
		
		* We need to bootstrap the confidence intervals 
		nlcom(_b[y1]-_b[y0])
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
		

/* Augmented inverse probability weighting */	
	// AIPTW 
		teffects aipw ($Y $W) ($A $W, logit)
		
	// By hand
		qui sum $Y
		gen double Y = ($Y -`r(min)')/(`r(max)' - `r(min)')
		qui sum $Y 
		global cf = (`r(max)' - `r(min)')
		
		* Step (i) prediction model for the outcome
		qui glm Y $A $W, fam(bin) 
		predict double QAW, xb
		qui glm Y $W if $A==1, fam(bin) 
		predict double Q1W, mu
		qui glm Y $W if $A==0, fam(bin)
		predict double Q0W, mu
		
		* Step (ii): prediction model for the treatment
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
		
		* Step (iii): Estimation equation
		gen double ATE = (sws*(Y-QAW) + (Q1W - Q0W)*$cf)
		sum ATE

		* Marginal structural model
		glm $Y $A [pw=sws],nolog

		* Bootstrap confidence intervals
		drop ATE
		capture program drop ACE
		program define ACE, rclass
		   capture drop ATE
		   gen double ATE = (sws*(Y-QAW) + (Q1W - Q0W)*$cf) // Estimation Equation
		   mean ATE
		   lincom _b[ATE]
		   return scalar ace =`r(estimate)'
		end
		qui bootstrap r(ace), reps(1000): ACE
		estat boot, all			
		
		
/* Ensemble learning targeted maximum likelihood estimation */
	// Run eltmle
		eltmle $Y $A $W, tmle   // install via 'help eltmle'
		shell "C:\Program Files\R\R-3.6.3\bin\x64\R.exe" CMD BATCH SLS.R

		
/* Simulation */
	// Data generation
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
	
	// Estimate the true simulated ATE
		mean psi
			
	//	ATE estimation
		* Regression adjustment
		teffects ra (Y w1 w2 w3 w4) (A)
		estimates store ra
			
		* IPTW
		teffects ipw (Y) (A w1 w2 w3 w4)
		estimates store ipw
			
		* IPTW-RA
		teffects ipwra (Y w1 w2 w3 w4) (A w1 w2 w3 w4)
		estimates store ipwra	
			
		* AIPTW
		teffects aipw (Y w1 w2 w3 w4) (A w1 w2 w3 w4)
		estimates store aipw
			
		* Results
		qui reg psi
		estimates store psi
		estout psi ra ipw ipwra aipw
	
	// Ensemble learning maximum likelihood estimation
		eltmle Y A w1 w2 w3 w4, tmle
		shell "C:\Program Files\R\R-3.6.3\bin\x64\R.exe" CMD BATCH SLS.R 
			
	// Relative bias of each ATE
		* IPTW
		display (0.38 - 0.18)/.38 
		
		* AIPTW
		display (0.28 - 0.18)/0.28
		
		* ELTMLE
		display (0.22 - 0.18)/0.22
		
			
