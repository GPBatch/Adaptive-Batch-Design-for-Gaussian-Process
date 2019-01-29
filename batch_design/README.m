% These two files "seq_design.m" and "updategp" are supplementary documents 
% for the paper "Adaptive Batching for Gaussian Process Surrogates with 
% Application in Noisy Level Set Estimation?. 

% "seq_design.m" is the main function to implement the design algorithms MCU
% with FB, MCU with RB, MCU with MLB and ABSUR. "updategp" is the wrapper
% function for the whole sequential design process in synthetic
% experiments, including initialization, sequential design, updating
% models and measuring performance.

% The GP is fitted with functions "gp_optim" and "gp_pred" in GPstuff
% package written by Vanhatalo, J., Riihimäki, J., Hartikainen, J., 
% Jylänki, P., Tolvanen, V., & Vehtari, A. (2013). GPstuff: Bayesian 
% modeling with Gaussian processes. Journal of Machine Learning Research, 
% 14(Apr), 1175-1179. To run the codes smoothly, you also need to download
% the GPstuff package.

% The complete implementation codes can be found in GitHub: 

% https://github.com/GPBatch/Gaussian-Process-with-Batch-Design/