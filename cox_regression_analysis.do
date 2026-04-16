/* ============================================================
   Cox Regression Analysis — CKD Risk in Type 2 Diabetes
   Equivalent of cox-regression.ipynb (Python)
   ============================================================ */

clear all
set more off
capture log close
log using "cox_regression_analysis.log", replace


/* ============================================================
   SECTION 1: LOAD & MERGE DATA
   In Python we generated synthetic data. Replace the use
   statements below with your actual NHANES/SAIL file paths.
   ============================================================ */

* Load and merge component files (adjust paths as needed)
* use "demo.dta",  clear
* merge 1:1 SEQN using "biol.dta",  nogen
* merge 1:1 SEQN using "bp.dta",    nogen
* merge 1:1 SEQN using "bmx.dta",   nogen
* merge 1:1 SEQN using "alb.dta",   nogen
* merge 1:1 SEQN using "mort.dta",  nogen
* merge 1:1 SEQN using "diab.dta",  nogen

* ── Rename variables to match analysis names ──────────────────
rename RIDAGEYR   AGE
rename RIAGENDR   SEX
rename RIDRETH3   ETHNICITY
rename LBXGH      HBA1C
rename BPXSY1     SBP
rename BMXBMI     BMI
rename URDACT     UACR
rename MORTSTAT   EVENT
rename PERMTH_EXM FOLLOW_UP_MONTHS
rename DIQ010     DIABETES


/* ============================================================
   SECTION 2: FILTER TO DIABETIC PATIENTS ONLY
   ============================================================ */

keep if DIABETES == 1


/* ============================================================
   SECTION 3: VARIABLE TRANSFORMATIONS & RANGE CHECKS
   ============================================================ */

* Convert HbA1c from DCCT % to IFCC mmol/mol
gen HBA1C_MMOL = (HBA1C - 2.15) * 10.929
label variable HBA1C_MMOL "HbA1c (mmol/mol)"

* Log-transform UACR (right-skewed)
gen LOG_UACR = log10(UACR + 1)
label variable LOG_UACR "log10(UACR + 1)"

* Convert follow-up to years
gen FOLLOW_UP_YEARS = FOLLOW_UP_MONTHS / 12
label variable FOLLOW_UP_YEARS "Follow-up (years)"

* Apply clinical range filters
keep if HBA1C_MMOL >= 20  & HBA1C_MMOL <= 195
keep if BMI        >= 14  & BMI        <= 70
keep if SBP        >= 60  & SBP        <= 280
keep if AGE        >= 18

* Summary after filtering
sum AGE SEX BMI SBP HBA1C_MMOL LOG_UACR FOLLOW_UP_YEARS EVENT


/* ============================================================
   SECTION 4: MULTIPLE IMPUTATION (MICE)
   Replaces listwise deletion. Outcome variables (EVENT,
   FOLLOW_UP_YEARS) are never imputed.
   ============================================================ */

* Register variables for imputation
mi set wide
mi register imputed AGE SEX BMI SBP HBA1C_MMOL LOG_UACR ETHNICITY
mi register regular EVENT FOLLOW_UP_YEARS

* Report missingness
misstable summarize AGE SEX BMI SBP HBA1C_MMOL LOG_UACR ETHNICITY

* MICE — 10 iterations, 5 imputed datasets
* Continuous vars: predictive mean matching (pmm)
* Binary/categorical: logistic / multinomial logistic
mi impute chained ///
    (pmm, knn(5))  AGE BMI SBP HBA1C_MMOL LOG_UACR ///
    (logit)        SEX ///
    (mlogit)       ETHNICITY ///
    = EVENT FOLLOW_UP_YEARS, ///
    add(5) rseed(42) dots


/* ============================================================
   SECTION 5: INTERACTION TERMS (clinically motivated)
   1. HbA1c × Age      — glycaemic damage amplifies with age
   2. HbA1c × LOG_UACR — synergistic CKD/T2DM markers
   3. SBP   × LOG_UACR — hypertension + albuminuria multiplier
   ============================================================ */

mi passive: gen HBA1C_x_AGE  = HBA1C_MMOL * AGE
mi passive: gen HBA1C_x_UACR = HBA1C_MMOL * LOG_UACR
mi passive: gen SBP_x_UACR   = SBP        * LOG_UACR

label variable HBA1C_x_AGE  "HbA1c × Age interaction"
label variable HBA1C_x_UACR "HbA1c × log(UACR) interaction"
label variable SBP_x_UACR   "SBP × log(UACR) interaction"


/* ============================================================
   SECTION 6: KAPLAN-MEIER CURVES
   ============================================================ */

* Declare survival data (use _1 = first imputed dataset for KM)
mi extract 1, clear

stset FOLLOW_UP_YEARS, failure(EVENT==1)

* KM by HbA1c tertile
xtile HBA1C_TERTILE = HBA1C_MMOL, nq(3)
label define tertile 1 "Low" 2 "Medium" 3 "High"
label values HBA1C_TERTILE tertile

sts graph, by(HBA1C_TERTILE) ///
    title("Kaplan-Meier: By HbA1c Tertile") ///
    xtitle("Time (years)") ytitle("Survival probability") ///
    legend(order(1 "Low" 2 "Medium" 3 "High") title("HbA1c tertile"))
graph export "km_hba1c_tertile.png", replace

* Log-rank test
sts test HBA1C_TERTILE
* Note p-value for write-up

* KM by sex
sts graph, by(SEX) ///
    title("Kaplan-Meier: By Sex") ///
    xtitle("Time (years)") ytitle("Survival probability") ///
    legend(order(1 "Male" 2 "Female") title("Sex"))
graph export "km_sex.png", replace

sts test SEX


/* ============================================================
   SECTION 7: COX PROPORTIONAL HAZARDS MODEL
   Run on MI data using Rubin's rules via -mi estimate-
   ============================================================ */

* Re-load MI dataset
use "your_mi_dataset.dta", clear   // or continue from mi impute above

stset FOLLOW_UP_YEARS, failure(EVENT==1)

* Fit Cox model with main effects + interactions + ethnicity dummies
mi estimate, post: ///
    stcox AGE SEX BMI SBP HBA1C_MMOL LOG_UACR ///
          HBA1C_x_AGE HBA1C_x_UACR SBP_x_UACR ///
          i.ETHNICITY

* Store estimates
estimates store cox_mi

* Display results with hazard ratios
mi estimate, post hr: ///
    stcox AGE SEX BMI SBP HBA1C_MMOL LOG_UACR ///
          HBA1C_x_AGE HBA1C_x_UACR SBP_x_UACR ///
          i.ETHNICITY


/* ============================================================
   SECTION 8: PROPORTIONAL HAZARDS ASSUMPTION CHECK
   Schoenfeld residuals test (equivalent to R/Python cph.check_assumptions)
   ============================================================ */

* Run on a single complete/imputed dataset for diagnostics
mi extract 1, clear
stset FOLLOW_UP_YEARS, failure(EVENT==1)

stcox AGE SEX BMI SBP HBA1C_MMOL LOG_UACR ///
      HBA1C_x_AGE HBA1C_x_UACR SBP_x_UACR ///
      i.ETHNICITY, schoenfeld(sch*) scaledsch(sca*)

estat phtest, detail
* p > 0.05 for each covariate = PH assumption holds


/* ============================================================
   SECTION 9: C-STATISTIC WITH BOOTSTRAP 95% CI
   ============================================================ */

* Harrell's C from a single imputed dataset
stcox AGE SEX BMI SBP HBA1C_MMOL LOG_UACR ///
      HBA1C_x_AGE HBA1C_x_UACR SBP_x_UACR ///
      i.ETHNICITY
estat concordance
* Note: Stata reports Somers' D; C = (D + 1) / 2

* Bootstrap 95% CI (200 reps)
bootstrap c=r(C), reps(200) seed(42): ///
    stcox AGE SEX BMI SBP HBA1C_MMOL LOG_UACR ///
          HBA1C_x_AGE HBA1C_x_UACR SBP_x_UACR ///
          i.ETHNICITY
estat concordance


/* ============================================================
   SECTION 10: EXPORT RESULTS TO EXCEL
   Uses -putexcel- (built-in) and -esttab- (ssc install estout)
   ============================================================ */

* Install estout if not already installed
* ssc install estout, replace

* ── Sheet 1: Cox model coefficients (HR, CI, p) ──────────────
mi estimate, post hr: ///
    stcox AGE SEX BMI SBP HBA1C_MMOL LOG_UACR ///
          HBA1C_x_AGE HBA1C_x_UACR SBP_x_UACR ///
          i.ETHNICITY

esttab using "cox_regression_results.xlsx", ///
    eform ci(2) b(3) p(4) ///
    title("Cox PH Model — Hazard Ratios") ///
    label star(* 0.05 ** 0.01 *** 0.001) ///
    replace sheet("2_Cox_Model_Results")

* ── Sheet 2: Descriptive statistics ──────────────────────────
putexcel set "cox_regression_results.xlsx", sheet("1_Descriptive_Stats") modify
tabstat AGE BMI SBP HBA1C_MMOL LOG_UACR FOLLOW_UP_YEARS EVENT, ///
    stats(n mean sd min p25 p50 p75 max) col(stats) save
putexcel A1 = matrix(r(StatTotal)), names

* ── Sheet 3: Descriptive stats by event status ───────────────
putexcel set "cox_regression_results.xlsx", sheet("6_Stats_by_Event") modify
tabstat AGE BMI SBP HBA1C_MMOL LOG_UACR FOLLOW_UP_YEARS, ///
    by(EVENT) stats(mean sd) col(stats) save
putexcel A1 = "Event status"
putexcel A2 = "No event (0)"
putexcel A3 = "Event (1)"


/* ============================================================
   NOTE ON MACHINE LEARNING COMPARISON (GBSA/RSF)
   Gradient Boosting Survival Analysis is not available natively
   in Stata. Options:
     • Use R (sksurv / randomForestSRC) and import results
     • Use Stata's -stpm2- for flexible parametric survival models
       as an alternative non-parametric benchmark
   ============================================================ */

log close
