# BayesTTE

These folders contain all the components required to simulate and analyse TTE endpoints for the case studies in my thesis. The data section includes IPD files in CSV format representing historical PFS and OS reconstructed as pseudo-IPD from published Kaplan–Meier survival curves.

The project is organised using RStudio project files (Simulacion_Survival_X.Rproj) to maintain a structured workflow and an R Markdown report (Informe_X.rmd) that compiles results, figures and interpretations in a single document.

For data generation and simulation, scripts such as gen_surv_data.r create synthetic survival datasets, while rand.r handles patient randomisation in simulated trials. Simulation scripts (sim_bayes.r, sim_freq.r, and sim_trials.r) allow for running Bayesian and frequentist survival simulations and for replicating full clinical trial scenarios with various features.

Analysis scripts are provided to process the simulated datasets. analysis_surv_data_bayes.r and analysis_surv_data_freq.r perform Bayesian and frequentist analyses. The scripts surv_analysis_bayes.r and surv_analysis_freq.r handle post-processing and extract results while get_boundaries_IA.r calculates IAs boundaries for applying stopping rules.

Finally, the repository includes a collection of Bayesian models written in Stan which are used throughout the analyses performed in my thesis.
