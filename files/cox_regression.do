/*===========================================================================
  STATA DO-FILE: Cox Proportional Hazards Model for DRKD Risk in T2D
  
  Study:    Development of Risk Prediction Models for DRKD Using AI
  Author:   Samuel Ato Gyasi Otabir
  Dataset:  SAIL Databank / Birmingham Insight (replace synthetic data below)
  Date:     April 2026
  
  WORKFLOW:
    1.  Load and prepare data
    2.  Clinical range filters
    3.  Feature engineering
    4.  Descriptive statistics
    5.  Kaplan-Meier curves
    6.  Cox PH model (main effects)
    7.  Cox PH model (with interactions)
    8.  Proportional hazards assumption check
    9.  Model performance (C-statistic + bootstrap CI)
    10. Export results to Excel
===========================================================================*/

version 17
clear all
set more off
capture log close
log using "cox_regression_log.smcl", replace

/*---------------------------------------------------------------------------
  SECTION 1: LOAD DATA
  Replace this section with your real SAIL / Birmingham dataset.
  Expected variables:
    AGE            - Age in years
    SEX            - 1=Male, 2=Female
    ETHNICITY      - Ethnicity code
    HBA1C_MMOL     - HbA1c in mmol/mol
    SBP            - Systolic blood pressure (mmHg)
    BMI            - Body mass index (kg/m²)
    UACR           - Urine albumin-creatinine ratio (mg/mmol)
    FOLLOW_UP_YEARS - Follow-up time in years
    EVENT          - 1=DRKD event occurred, 0=censored
---------------------------------------------------------------------------*/

* ── Synthetic data generation (remove when using real data) ───────────────
set seed 42
set obs 5000

gen SEQN           = _n + 100000
gen AGE            = round(rnormal(55, 15))
replace AGE        = max(18, min(85, AGE))
gen SEX            = runiformint(1, 2)
gen ETHNICITY      = runiformint(1, 6)
gen HBA1C          = rnormal(6.5, 2.0)
replace HBA1C      = max(4.0, min(15.0, HBA1C))
gen SBP            = rnormal(130, 20)
replace SBP        = max(90, min(200, SBP))
gen BMI            = rnormal(28, 6)
replace BMI        = max(15, min(60, BMI))
gen UACR_RAW       = exp(rnormal(2, 1.5))
replace UACR_RAW   = max(0.1, min(1000, UACR_RAW))
gen EVENT          = rbinomial(1, 0.15)
gen FOLLOW_UP_MONTHS = round(rexponential(120))
replace FOLLOW_UP_MONTHS = max(1, min(240, FOLLOW_UP_MONTHS))
gen T2D            = (runiform() < 0.15)        // ~15% T2D patients

* Keep T2D patients only
keep if T2D == 1

/*---------------------------------------------------------------------------
  SECTION 2: FEATURE ENGINEERING
---------------------------------------------------------------------------*/

* HbA1c: DCCT % to IFCC mmol/mol
gen HBA1C_MMOL = (HBA1C - 2.15) * 10.929
label variable HBA1C_MMOL "HbA1c (mmol/mol)"

* Log-transform UACR (right-skewed)
gen LOG_UACR = log10(UACR_RAW + 1)
label variable LOG_UACR "Log10(UACR + 1)"

* Follow-up in years
gen FOLLOW_UP_YEARS = FOLLOW_UP_MONTHS / 12
label variable FOLLOW_UP_YEARS "Follow-up time (years)"

* Variable labels
label variable AGE        "Age (years)"
label variable SEX        "Sex (1=Male, 2=Female)"
label variable BMI        "BMI (kg/m²)"
label variable SBP        "Systolic BP (mmHg)"
label variable EVENT      "DRKD event (1=yes, 0=censored)"

/*---------------------------------------------------------------------------
  SECTION 3: CLINICAL RANGE FILTERS
---------------------------------------------------------------------------*/

keep if HBA1C_MMOL >= 20 & HBA1C_MMOL <= 195
keep if BMI        >= 14 & BMI        <= 70
keep if SBP        >= 60 & SBP        <= 280
keep if AGE        >= 18

* Drop missing values in key variables
foreach v of varlist AGE SEX BMI SBP HBA1C_MMOL LOG_UACR EVENT FOLLOW_UP_YEARS ETHNICITY {
    drop if missing(`v')
}

di "Final analytic sample: " _N " patients"

/*---------------------------------------------------------------------------
  SECTION 4: DESCRIPTIVE STATISTICS
---------------------------------------------------------------------------*/

di _newline "=== DESCRIPTIVE STATISTICS ==="
summarize AGE SEX BMI SBP HBA1C_MMOL LOG_UACR FOLLOW_UP_YEARS EVENT, detail

* Descriptive statistics by event status
di _newline "=== DESCRIPTIVE STATISTICS BY EVENT STATUS ==="
table EVENT, statistic(mean AGE BMI SBP HBA1C_MMOL LOG_UACR) ///
             statistic(sd   AGE BMI SBP HBA1C_MMOL LOG_UACR) nformat(%6.2f)

* Frequency table for sex and ethnicity
tabulate SEX EVENT, col chi2
tabulate ETHNICITY EVENT, col chi2

/*---------------------------------------------------------------------------
  SECTION 5: SURVIVAL TIME SETUP
  stset declares the survival time and event indicator for Stata's
  survival analysis commands (st* family).
---------------------------------------------------------------------------*/

stset FOLLOW_UP_YEARS, failure(EVENT==1) id(SEQN)
stdes

/*---------------------------------------------------------------------------
  SECTION 6: KAPLAN-MEIER CURVES
---------------------------------------------------------------------------*/

* KM curve overall
sts graph, ///
    title("Kaplan-Meier Survival Curve — Overall") ///
    ytitle("Survival probability") xtitle("Time (years)") ///
    ylabel(0(0.2)1) ci
graph export "km_overall.png", replace

* KM by HbA1c tertile
xtile HBA1C_TERTILE = HBA1C_MMOL, nquantiles(3)
label define tertile_lbl 1 "Low HbA1c" 2 "Medium HbA1c" 3 "High HbA1c"
label values HBA1C_TERTILE tertile_lbl

sts graph, by(HBA1C_TERTILE) ///
    title("Kaplan-Meier by HbA1c Tertile") ///
    ytitle("Survival probability") xtitle("Time (years)") ///
    ylabel(0(0.2)1) ci legend(order(1 "Low" 2 "Medium" 3 "High") title("HbA1c tertile"))
graph export "km_hba1c_tertile.png", replace

* Log-rank test for HbA1c tertile
sts test HBA1C_TERTILE
di "Log-rank test p-value for HbA1c tertile above"

* KM by sex
sts graph, by(SEX) ///
    title("Kaplan-Meier by Sex") ///
    ytitle("Survival probability") xtitle("Time (years)") ///
    ylabel(0(0.2)1) ci legend(order(1 "Male" 2 "Female") title("Sex"))
graph export "km_sex.png", replace

sts test SEX
di "Log-rank test p-value for sex above"

/*---------------------------------------------------------------------------
  SECTION 7: COX PH MODEL — MAIN EFFECTS
---------------------------------------------------------------------------*/

di _newline "=== COX MODEL — MAIN EFFECTS ==="

* Ethnicity dummy variables (reference = group 1)
tab ETHNICITY, gen(ETH_)

stcox AGE SEX BMI SBP HBA1C_MMOL LOG_UACR ETH_2 ETH_3 ETH_4 ETH_5 ETH_6, ///
    vce(robust) efron

* Store estimates
estimates store cox_main

* Results table: hazard ratios + 95% CI + p-values
stcox, hr
di _newline "Concordance (C-statistic): " e(C)

/*---------------------------------------------------------------------------
  SECTION 8: COX PH MODEL — WITH INTERACTION TERMS
  Gavin's suggestion: include clinically motivated interactions.
  Three interactions:
    1. HbA1c × Age      — glycaemic impact amplifies with age
    2. HbA1c × LOG_UACR — both T2D/CKD markers — likely synergistic
    3. SBP × LOG_UACR   — hypertension + albuminuria known risk multiplier
---------------------------------------------------------------------------*/

di _newline "=== COX MODEL — WITH INTERACTIONS ==="

gen HBA1C_x_AGE  = HBA1C_MMOL * AGE
gen HBA1C_x_UACR = HBA1C_MMOL * LOG_UACR
gen SBP_x_UACR   = SBP        * LOG_UACR

label variable HBA1C_x_AGE  "HbA1c × Age interaction"
label variable HBA1C_x_UACR "HbA1c × Log(UACR) interaction"
label variable SBP_x_UACR   "SBP × Log(UACR) interaction"

stcox AGE SEX BMI SBP HBA1C_MMOL LOG_UACR ///
      HBA1C_x_AGE HBA1C_x_UACR SBP_x_UACR ///
      ETH_2 ETH_3 ETH_4 ETH_5 ETH_6, ///
      vce(robust) efron

estimates store cox_interactions

stcox, hr
di _newline "Concordance (C-statistic): " e(C)

* Likelihood ratio test: interactions model vs main effects model
lrtest cox_interactions cox_main
di "LR test above: do interactions significantly improve fit?"

/*---------------------------------------------------------------------------
  SECTION 9: PROPORTIONAL HAZARDS ASSUMPTION CHECK
  Schoenfeld residuals test — p > 0.05 suggests PH assumption holds
  for that variable. Global test p > 0.05 suggests overall PH holds.
---------------------------------------------------------------------------*/

di _newline "=== PROPORTIONAL HAZARDS ASSUMPTION CHECK ==="

* Use the interactions model as final model
quietly stcox AGE SEX BMI SBP HBA1C_MMOL LOG_UACR ///
              HBA1C_x_AGE HBA1C_x_UACR SBP_x_UACR ///
              ETH_2 ETH_3 ETH_4 ETH_5 ETH_6, ///
              vce(robust) efron schoenfeld(sch*) scaledsch(sca*)

estat phtest, detail
di "Global PH test p-value above — p>0.05 means PH assumption holds"

* Plot scaled Schoenfeld residuals for key variables
foreach v in AGE HBA1C_MMOL SBP LOG_UACR BMI {
    capture stphplot, by(`v') ///
        title("Log-log plot: `v'") ///
        ytitle("log(-log(Survival))") xtitle("log(Time)")
    capture graph export "phplot_`v'.png", replace
}

/*---------------------------------------------------------------------------
  SECTION 10: C-STATISTIC WITH BOOTSTRAP 95% CI
---------------------------------------------------------------------------*/

di _newline "=== C-STATISTIC WITH BOOTSTRAP 95% CI ==="

* Store observed C-statistic
quietly stcox AGE SEX BMI SBP HBA1C_MMOL LOG_UACR ///
              HBA1C_x_AGE HBA1C_x_UACR SBP_x_UACR ///
              ETH_2 ETH_3 ETH_4 ETH_5 ETH_6, vce(robust) efron
scalar c_obs = e(C)
di "Observed C-statistic: " c_obs

* Bootstrap 95% CI
bootstrap c_stat=e(C), reps(200) seed(42) nodots: ///
    stcox AGE SEX BMI SBP HBA1C_MMOL LOG_UACR ///
          HBA1C_x_AGE HBA1C_x_UACR SBP_x_UACR ///
          ETH_2 ETH_3 ETH_4 ETH_5 ETH_6, vce(robust) efron

estat bootstrap, percentile
di "Bootstrap 95% CI for C-statistic above"

/*---------------------------------------------------------------------------
  SECTION 11: EXPORT RESULTS TO EXCEL
  Uses putexcel to write all key results into a downloadable .xlsx file.
---------------------------------------------------------------------------*/

di _newline "=== EXPORTING RESULTS TO EXCEL ==="

putexcel set "cox_regression_results_stata.xlsx", replace sheet("Model_Results")

* Header
putexcel A1 = "Cox Regression Results — DRKD Risk Prediction in T2D", bold
putexcel A2 = "Author: Samuel Ato Gyasi Otabir | Date: April 2026"

* Column headers for results table
putexcel A4 = "Variable"       , bold
putexcel B4 = "Hazard Ratio"   , bold
putexcel C4 = "95% CI Lower"   , bold
putexcel D4 = "95% CI Upper"   , bold
putexcel E4 = "p-value"        , bold
putexcel F4 = "Interpretation" , bold

* Run model and collect results into a matrix
quietly stcox AGE SEX BMI SBP HBA1C_MMOL LOG_UACR ///
              HBA1C_x_AGE HBA1C_x_UACR SBP_x_UACR ///
              ETH_2 ETH_3 ETH_4 ETH_5 ETH_6, vce(robust) efron

* Write HR, CI and p-values row by row
local vars "AGE SEX BMI SBP HBA1C_MMOL LOG_UACR HBA1C_x_AGE HBA1C_x_UACR SBP_x_UACR ETH_2 ETH_3 ETH_4 ETH_5 ETH_6"
local row = 5
foreach v of local vars {
    local hr   = exp(_b[`v'])
    local ll   = exp(_b[`v'] - 1.96*_se[`v'])
    local ul   = exp(_b[`v'] + 1.96*_se[`v'])
    local pval = 2*(1 - normal(abs(_b[`v']/_se[`v'])))
    
    putexcel A`row' = "`v'"
    putexcel B`row' = `hr',   nformat("0.000")
    putexcel C`row' = `ll',   nformat("0.000")
    putexcel D`row' = `ul',   nformat("0.000")
    putexcel E`row' = `pval', nformat("0.000")
    
    local row = `row' + 1
}

* Write C-statistic
local bottom = `row' + 1
putexcel A`bottom' = "C-statistic (concordance)", bold
putexcel B`bottom' = e(C), nformat("0.000")

* Add descriptive stats sheet
putexcel set "cox_regression_results_stata.xlsx", modify sheet("Descriptive_Stats")
putexcel A1 = "Variable", bold
putexcel B1 = "N",        bold
putexcel C1 = "Mean",     bold
putexcel D1 = "SD",       bold
putexcel E1 = "Min",      bold
putexcel F1 = "Max",      bold

local vars2 "AGE BMI SBP HBA1C_MMOL LOG_UACR FOLLOW_UP_YEARS"
local row2 = 2
foreach v of local vars2 {
    quietly summarize `v'
    putexcel A`row2' = "`v'"
    putexcel B`row2' = r(N),      nformat("0")
    putexcel C`row2' = r(mean),   nformat("0.00")
    putexcel D`row2' = r(sd),     nformat("0.00")
    putexcel E`row2' = r(min),    nformat("0.00")
    putexcel F`row2' = r(max),    nformat("0.00")
    local row2 = `row2' + 1
}

di "Results exported to: cox_regression_results_stata.xlsx"
di "Sheets: Model_Results, Descriptive_Stats"

log close
di _newline "=== DO-FILE COMPLETE ==="
