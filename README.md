
# The disconnect between development and intended use of clinical prediction models for Covid-19: a systematic review and real-world data illustration


  R code repository for the paper "The disconnect between development and intended use of clinical prediction models for Covid-19: a systematic review and real-world data illustration".

  The repository contains the following:

* [model_strategies.R](model_strategies.R) : this R script reproduces all analysis strategies and figures of the manuscript

* [Covid-19 2020 RIVM data](Data/datahosp.rda): the data used for the analysis. For privacy reasons, 27 patients were removed. The analysis results provided by these data only differ very slightly from the results presented in the paper.


  Data description:

* **Age_cat**: Patient age, categorized ($\small\leq$ 50, 50 - 59, 60 - 69, 70 - 79, 80 - 89, $\small>$ 90)
* **Sex_m**: Patient gender (0 = Male, 1 = Female)
* **N_comobidities**: Number of medical conditions, capped at 3
* **DaysToDeath**: Days from hospitalization and positive test to death, censored at 28 days
* **Death**: Death status indicator related to DaysToDeath (0 = alive, 1 = dead)
* **DaysToICU**: Days from hospitalization and positive test to death, censored at 28 days
* **Death**: ICU treatment status indicator related to DaysToICU (0 = untreated, 1 = treated)
