# Predicting Task and Subject Identity using Functional Connectivity and BOLD Variability 

This repository provides code for the analyses used in the manuscript: Predicting Task and Subject Identity using Functional Connectivity and BOLD Variability submitted to Brain Connectivity. 

## Preprocessing
The main functions for preprocessing the data and computing features are: 

`preprocessROIdata_v1.m`: formats the time course for further analysis and calls `preparedatastruct2.m`    
`compute_features.m`: computes and saves machine learning features so that we 
don't have to recompute features every modeling run     


## Analysis 

The main functions for running the analyses are: 

`runclassifier_fully_parallelized.m`: runs the logistic regression analyses    
`runclassifier_1NN.m`: runs the single nearest neighbor analyses    

These functions are called by `run_classifiers_in_parallel.m`

## Results 

The main functions for result aggregation are: 

`finn_et_al_analysis.m`: performs the analyses in Finn et al. 2015    
`finn_plot_all_accuracies.m`: produces bar plots comparing accuracies for each preprocessing variation     
`plot_all.m`: generates plots based on the type of analysis ran (task classification, 
subject classification, data exploration, preprocessing robustness) 
<!-- `plot_1NN_heatmaps.m`: plots the Finn et al. 2015 heatmap comparison     
`plot_FV_heatmap.m`: visualizes the raw BOLD variability features     
`print_result_table_regularization_optimized`: prints out the logistic regression accuracies    
`plot_accuracies_all_datasets.m`: plots the logistic regression accuracies for each preprocessing variation      
-->


