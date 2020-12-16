/*Name: Kareena Satia

Project: Replication of "Returns to Scale in Electricty Supply by Marc Nerlove (1960)"
and a partial replication of "Economies of Scale in U.S. Electric Power Generation
by Laurits R. Christensen; William H. Greene (1984)". I have used the data and model in Greene and 
I have further extended the model by computing the efficient scale using the Delta Method. 

Starting from a simple Cobb-Douglas production function, Nerlove derives the returns 
to scale by modifying the econometric specification based on results obtained in each steps. 

Christensen and Greene use a modification of Nerlove's model to compute the efficient scale
of output and its confidence interval. 
*/
clear
// This is the project working folder 
global projdir "C:\Users\18575\Desktop\Nerlove"

// Raw Data Folder
global raw "$projdir\Data"

// Folder where the Clean Data set is stored
global final "$projdir\Clean_Data"

// Folder for Output (tables and graphs) 
global outcomes "$projdir\Output"

log using "$outcomes\Log\logfile.log", replace

// Setting Working Directory
cd "$projdir"

import excel using "$raw\Nerlove_data.xlsx", firstrow

*Producing logs of variables* 
gen lnc = ln(COSTS)
gen lny = ln(KWH)
gen lnPL = ln(PL)
gen lnPF = ln(PF)
gen lnPK = ln(PK)

*Transforming variables for regression and labeling these variables appropriately* 
gen yvar = lnc - lnPF
gen p1p3 = lnPL - lnPF
gen p2p3 = lnPK - lnPF

label var lny "ln(Y)"
label var yvar "ln(C) - ln(PF)"
label var p1p3 "ln(PL) - ln(PF)"
label var p2p3 "ln(PK) - ln(PF)"

******************************MODEL A******************************************

*Table 3: Model A regression: Ln(C) - Ln(Pf) = K + 1/r(Y) + a1/r*(P_l - P_F) + a2/r*(P_K - P_F) + V 
reg yvar lny p1p3 p2p3
*outreg2 using "$outcomes\regNerlove", replace label bdec(2) sdec(2) rdec(2) noaster 

di "The elasticity of output with respect to the price of capital is" _b[p2p3]/_b[lny]

//The elasticity of output with respect to the price of capital is negative. Nerlove designed Model B to evade this difficulty. 

gen x = (KWH[_n] - KWH[_n-1])/KWH[_n-1]
label var x "change in Y"
reg yvar lny p1p3 p2p3 x
outreg2 using "$outcomes\regNerlove", append label tex bdec(2) sdec(2) rdec(2) noaster

*Plotting residuals of the regression against the logarithm of output
predict resid, residuals
twoway (qfitci resid lny)
graph export "$outcomes\Graphs\residuals.png", replace

//Visual evidence suggests that the outcome variable is not independent of the residuals 
// i.e the regression relationship defined in Model A is not linear in logs


// Since firms operate on their Short-run total cost curve, Nerlove divides the firms 
// into five quintiles based on output 

*Dividing firms in 5 quintiles based on output
sort KWH
xtile groups = KWH, nq(5)

*Running separate regressions for each of the quintiles that might be linear in logs
forvalues i=1/5{
    reg yvar lny p1p3 p2p3 if groups == `i'
	scalar ret_to_scale_A_`i' = (1/_b[lny])
	di "The returns to scale for firms in group `i' is", ret_to_scale_A_`i'
	outreg2 using "$outcomes\regNerlove", append label tex bdec(2) sdec(2) rdec(2) noaster
}

// Nerlove finds that the firms in the largest quintile are experiencing diminishing
// returns to scale

*Testing "neutral variations in returns to scale" by allowing the coefficients 
*of the various prices in the regression to be the same while allowing constant
*terms and the coefficients of output to differ

*Creating new variables of output for different groups to run the next set of regressions
forvalues i = 1/5{
    gen lny`i' = lny
	replace lny`i' = 0 if groups != `i'
	label var lny`i' "Y in Quantile `i' firms"
}

* This regresion allows the coefficient on ln(y) to differ across quintiles
reg yvar lny1 lny2 lny3 lny4 lny5 p1p3 p2p3
outreg2 using "$outcomes\regNerlove", append label tex bdec(2) sdec(2) rdec(2) noaster sortvar(lny p1p3 p2p3 x lny1 lny2 lny3 lny4 lny5)

// There is no discerning pattern in the coefficients on ln(y) for all the five groups 

******************************MODEL B****************************************

// In order to overcome the issue of negative elasticity of output with respect to the price of capital, 
// Nerlove is keeping the price of capital same for all the firms 
*Running Model B: Ln(C) = K' + 1/r(Y) + a1/r*(P_l) + a2/r*(P_F) + V 
reg lnc lny lnPL lnPF
outreg2 using "$outcomes\reg1Nerlove", replace label bdec(2) sdec(2) rdec(2) noaster

forvalues i=1/5{
    reg lnc lny lnPL lnPF if groups == `i'
	scalar ret_to_scale_B_`i' = (1/_b[lny])
	di "The returns to scale for firms in group `i' is", ret_to_scale_B_`i'
	outreg2 using "$outcomes\reg1Nerlove", append label tex bdec(2) sdec(2) rdec(2) noaster
}
// The results show increasing returns at a diminishing rate for all except the largest firms (group 5)


********************************MODEL C************************************

// The previous two models suggest that the returns to scale is heterogenous with respect to the output of the firms
// Hence, Nerlove treats the returns to scale as a continuous function of output 

*Transforming a variable 
gen lnysq = lny^2

//This model is derived from modifying Model A
*Model C regression: Ln(C) - Ln(Pf) = K + alpha(Y) + beta(Y)^2 + a1/r*(P_l - P_F) + a2/r*(P_K - P_F) + V
reg yvar lny lnysq p1p3 p2p3
outreg2 using "$outcomes\reg2Nerlove", replace tex label bdec(2) sdec(2) rdec(2) noaster

*Running Model C regression for all the five groups of firms 
forvalues i=1/5{
    reg yvar lny lnysq p1p3 p2p3 if groups == `i'
	scalar ret_to_scale_C_`i' = (1/_b[lny])
	di "The returns to scale for firms in group `i' is", ret_to_scale_`i'
	outreg2 using "$outcomes\reg4Nerlove", replace tex label bdec(2) sdec(2) rdec(2) noaster
}

*******************************MODEL D***************************************

// This model is derived from modifying Model B
*Model D regression: Ln(C) = K' + alpha(Y) + beta(Y)^2 + a1/r*(P_l) + a2/r*(P_F) + V
reg lnc lny lnysq lnPL lnPF 
outreg2 using "$outcomes\reg3Nerlove", replace tex label bdec(2) sdec(2) rdec(2) noaster

*Running Model D regression for all the five groups of firms
forvalues i=1/5{
    reg lnc lny lnysq lnPL lnPF if groups == `i'
	scalar ret_to_scale_D_`i' = (1/_b[lny])
	di "The returns to scale for firms in group `i' is", ret_to_scale_`i'
	outreg2 using "$outcomes\reg5Nerlove", replace tex label bdec(2) sdec(2) rdec(2) noaster
}

// Substantial increase in our estimate of the degree of returns to scale of the three largest size groups in Model C and Model D


******************************Greene's Paper*********************************
clear 

import delimited "$raw\Greene.csv"

*Generating transformations of variables for regression
gen lnq = ln(q)
gen lnsq = lnq^2
gen hlnsq = 0.5*lnsq
gen lnpk = ln(pk)
gen lnpf = ln(pf)
gen lnpl = ln(pl)
gen pkpf = pk/pf
gen lnpkpf = ln(pkpf)
gen plpf = pl/pf
gen lnplf = ln(plpf)
gen lnc = ln(cost)
gen cpf = cost/pf
gen lncpf = ln(cpf)

*Running Model E: ln(C/P_f) = alpha + beta*lnQ + gamma[1/2(lnQ)^2] + delta_k*ln(P_k/P_f) + delta_l*ln(P_l/P_f) + e 
*with the condition that the coefficients of lnpk+lnpf+lnpl = 1
reg lncpf lnq hlnsq lnpkpf lnplf

*Given the condition that the coefficient of lnpk+lnpf+lnpl = 1, we are finding the 
*coefficient and standard error of delta_f using the delta method
matlist e(V)
scalar delta_f = 1 - _b[lnpkpf] - _b[lnplf]

*Defining the jacobian matrix from delta method
matrix jacobian = (0,0,0,-1,-1)
matlist e(b)
matrix B = jacobian*e(V)*jacobian'
scalar C = B[1,1]
scalar standarderr = sqrt(C)
// The scalar standarderr contains the standard error of the coefficient on variable delta_f

// The efficient scale is the output at which the cost curve reaches its minimum. 
// I am attempting to find the efficient scale of output and its confidence interval.  
//Efficient scale of output = Q* = exp[(1-beta)/gamma]

scalar beta = _b[lnq]
scalar gamma = _b[hlnsq]

*Extracting the variance-covariance matrix of beta and gamma 
matrix N = e(V)[1..2,1..2]

*Defining the efficient scale of output
scalar qst = (1 - beta)/gamma
scalar qstar = exp(qst)

*Differentiating the efficient scale of output with respect to beta and gamma
scalar diff = qstar*(-1/gamma)
scalar diff2 = qstar*((-1+beta)/gamma^2)

*Defining the jacobian matrix from the delta method 
matrix jacobian = (diff, diff2)
matrix variance = L*N*L'

*Constructing confidence intervals of the efficient scale of output
scalar lowcon = qstar - 1.96*(variance[1,1])^0.5
scalar uppcon = qstar + 1.96*(variance[1,1])^0.5

di "The 95% confidence interval for the efficient scale of output is", lowcon, "to", uppcon 

log close 

