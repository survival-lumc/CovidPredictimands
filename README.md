# The disconnect between development and intended use of clinical prediction models for Covid-19: a systematic review and real-world data illustration

  R code repository for the paper "The disconnect between development and intended use of clinical prediction models for Covid-19: a systematic review and real-world data illustration".

  The repository contains the following:

* [model_strategies.R](model_strategies.R) : this R script reproduces all analysis strategies and figures of the manuscript

* [Covid-19 2020 RIVM data](Data/datahosp.rda): the data used for the analysis. For privacy reasons, 27 patients were removed. The analysis results provided by these data only differ very slightly from the results presented in the paper.

You can download a zip file containing the directory by clicking Code -> Download ZIP at the top-right of this Github page.

### Data description:

The rights to the original dataset used in this paper belong to the Dutch National Institute for Public Health and the Environment (RIVM). Requests to access the original data should be directed to COVID-19EPI@rivm.nl. To provide insight into the analyses that were conducted for this study, we provide here a subset of the data, to comply with privacy regulations. 

The original data consist of 22324 hospitalized patients that tested positive for SARS-CoV-2 on a PCR  test before December 31st, 2020, with follow-up until January 31st, 2021. We excluded patients with a positive test result obtained after death (n = 64) as well as patients with missing information on age and sex  (n = 9). The final dataset used in the analysis set consisted of 22251 cases.

### Authors:

| Author                   | Affiliation                                                                                                                                                |
| ------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Ilaria Prosepe**       | Department of Biomedical Data Sciences, Leiden University Medical Center, Leiden, Netherlands                                                              |
| **Rolf H. H. Groenwold** | Department of Clinical Epidemiology, Leiden University Medical Center, Leiden, Netherlands                                                                 |
| **Rachel Knevel**        | Department of Rheumatology, Leiden University Medical Center, Leiden, Netherlands                                                                          |
| **Romin Pajouheshnia**   | Division of Pharmacoepidemiology and Clinical Pharmacology, Utrecht Institute for Pharmaceutical Sciences (UIPS), Utrecht University, Utrecht, Netherlands |
| **Nan van Geloven**      | Department of Biomedical Data Sciences, Leiden University Medical Center, Leiden, Netherlands                                                              |
