# Understanding the Impact of Correlation within Pair-Bonds on Cormack-Jolly-Seber Models

The `R` code used to conduct and analyze the simulation study from our article "Understanding the Impact of Correlation within Pair-Bonds on Cormack-Jolly-Seber Models" is available here. 

Details about our methods and approach can be found in the article which is available at the link below. 

## Scripts


- `01_Fn_SimPropModel.R` and `02_Fn_ModelCode.R` contain custom functions written by the authors that are used to make conducting the study multiple times for different settings convenient and fast. Making changes to these files may result in issues running the code. To run the code, these files should be stored in your working directory (R code to check working directory: `getwd()`) in a folder called `scripts`. 

- `11_RunSimulationStudy.R` can be used to run the simulation study in the main text and in Appendix B.2. Parameter settings can be changed in lines 58-70 and lines 86-100. Warning this code can be slow, furthermore, it is recommended to run this code from a linux terminal  if possible using `Rscript 11_RunSimulationStudy.R`.

- `12_Run_CompareChatEstimators.R` allows the user to conduct the study completed in Appendix B.3, in which the deviance, Pearson's, and Fletcher's \hat{c} estimators are compared. Settings are changed in lines 90-103. Again, we recommend running this from a linux terminal. 

- `13_Run_FurtherInvestigation.R` allows users to investigate the behavior of the CJS model on a more granular level than is feasible to present within our manuscript. Some useful functions that appear in the script (stored in `02_Fn_ModelCode.R`) are discussed in the following line. 
  - `sim_cjs_dat()` will simulate datasets from the extension discussed in the manuscript (the amount is controlled by the iterations argument) based on the parameters that appear in the `parameter_list` object. 
  - `double_observed()` will take every dataset in `sim_cjs_dat()` and add an exact set of duplicate recapture histories to the data. This allows for investigation of the impact of replicated data on the likelihood ratio test and c-hat statistic when groups are unaccounted for. 
  - `gender_randomize()` will assign a new random gender to each entry in the each dataset generated by `sim_cjs_dat()`. Essentially removes correlation structure from the gender data. Data from this set and from `double_observed()` can be investigated. 
  - `model_cjs_data()` fits the standard CJS model to the data generated from the extended model. The parameter settings are determined by the `grouping` argument, where
    - grouping = B: (\phi^G, P^G)
    - grouping = S: (\phi^G, P)
    - grouping = R: (\phi, P^G)
    - grouping = N: (\phi, P)
  
  - `lrt_plots()`, `chat_plots`, `qlrt_plots` produce comparison plots between a simpler standrd CJS model and a more complex one (for instance (\phi^G, P) vs (\phi, P)) for the likelihood ratio test, c-hat statistic, and qausi-likelihood ratio test (argument of c-hat is needed here). 
  - `model_plots()` plots the estimate and confidence intervals of each run of the CJS model. For a given parameter (survival rate or recapture rate) and gender.
 - `21_PlotMSResults.R` can be used to produce the same plots as the ones available in the manuscript. Assuming the user has conducting the simulation studies from the paper using `11_RunSimulationStudy.R` and `12_Run_CompareChatEstimators.R`. The results are fairly stable across different runs of sufficient size so setting a seed is not necessary here. 
  
Article Link: TBA 

Reference: TBA
