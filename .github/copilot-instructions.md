# Copilot Instructions for Cox Regression Notebook Project

## 1. Project snapshot
- Single artifact: `cox-regression.ipynb`.
- Focus area: survival analysis and Cox proportional hazards modeling for chronic kidney disease in T2DM.
- No traditional package/module layout; treat this as a data analysis notebook repo.
- Uses synthetic NHANES-like data due to download URL issues with real CDC data.

## 2. Primary AI goals
- Enhance reproducibility: prefer explicit steps in notebook cells (load data, preprocess, fit model, evaluate, visualize).
- Keep analysis flow linear: dataset input -> data cleaning -> feature engineering -> Cox model fit -> validation metrics.

## 3. Environment & execution
- Most commands:
  - `jupyter notebook cox-regression.ipynb`
  - `jupyter lab`
- Install common libs likely used by this domain: `pandas numpy scipy lifelines statsmodels matplotlib seaborn`

## 4. Project-specific conventions
- Keep cell outputs compact and annotate with markdown summaries.
- Use `lifelines` style API for CoxPH (e.g., `CoxPHFitter().fit(df, duration_col='duration', event_col='event')`).
- Synthetic data generation uses realistic distributions matching NHANES variable structures.
- Random seed (42) ensures reproducible synthetic data generation.

## 5. What to avoid
- Avoid adding deep application architecture assumptions; this repo is analysis-focused.
- Avoid generating separate `src/` frameworks unless user asks for packaging.
- Don't attempt real NHANES downloads - they fail due to URL authentication issues.

## 6. Extension points
- Add one or two regression comparisons (e.g., Cox vs RandomForest survival from `scikit-survival`) only after verifying data and problem statement.
- Recommend adding a versioned data ingestion cell like:
  - `data = pd.read_csv('data/ckd_t2dm.csv')`

## 7. Quick sanity checks for automated patches
- Does the notebook include target columns? (`"duration"`, `"event"` are common in survival data).
- Do output cells include model diagnostics (e.g., log-likelihood, concordance index)?

## 8. Feedback request
- If the user has additional source files outside the current tree, ask to expose them before major automation.
- Confirm whether we should add typical unit-test scaffolding (pytest) or keep it notebook-only.
