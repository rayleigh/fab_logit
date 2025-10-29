This contains code needed to compute the Frequentist, Assisted by Bayes intervals for logistic regression. In code/stan, all Stan files can be found. Meanwhile, in code/R/, the code is organized as following: 

* **gen_functions:** Contains helper functions used by all experiments. 
* **sim_1:** There is no folder because the FAB intervals are generated with the standard FAB functions with a normal(0, 1) prior and n = 100.
* **sim_2:** Holds code to fit the models with the different priors to the simulated data set, to generate the FAB intervals, and then derive statistical processes related to the model
* **real_data_analysis:** Contains code to fit the model to real data and generate the FAB intervals

Note that the data_files/zip_code_db.Rdata is pulled from the R package: zipcodeR. 
