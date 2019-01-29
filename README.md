# Adaptive-Batch-Design-for-Gaussian-Process
This is the supplementary codes for the paper "Adaptive Batching for Gaussian Process Surrogates with Application in Noisy Level Set Estimation" 

## GPstuff
The whole package is built on GPstuff package, which is initially written by Jarno Vanhatalo, Jaakko Riihimäki, Jouni Hartikainen, Pasi Jylänki, Ville Tolvanen, and Aki Vehtari (2013). GPstuff: Bayesian Modeling with Gaussian Processes. Journal of Machine Learning Research, 14(Apr):1175-1179. The GPstuff toolbox is a collection of Gaussian Process models and computational tools required for inference. The tools include variance inference methods and sparse approximation. The toolbox works with MATLAB and can be called in R with Octave. Detailed installing instructions are included in the manual paper.

## Batch Design
We implemented the adaptive batch design mentioned in the paper: FB, RB, MLB and ABSUR. All adaptive batch design implementations are included in folder 'batch_design'. We also coded the synthetic experiments including 1D exponential function, 2D Brannin-Hoo function and 6D Hartman function in 'syntheticExperiments.m'. It also serves as a demo file, where users can change the true functions combined with different noise and use different model and batch sequential design algorithms. The 'updategp' function is a wrapper function to generate the whole designs with observations and batch size.
