function [mad] = compute_mad(ROI_NOW)
% function to compute the median absolue deviation 

m = mean(ROI_NOW);
abs_dev = abs(ROI_NOW - m) ;
mad = quantile( abs_dev, 0.5); 

