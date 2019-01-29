function [x_seq, y_seq, r_seq, ee, er, bias, metric, t, Ef, Varf, l, sigma2, sigman, steps, t_optim, t_gen, overhead, gamma, post_var] = updategppar(I, k, m0, d, budget, fun, noisestructure, noisevar, model, design, batch, r, t0)

% Generates and records the synthetic experiments results in parallel
% computing. 
%
% Description:
%   [x_seq, y_seq, r_seq,...] = UPDATEGPPAR(I, k, m0, d, budget, fun, 
%   noisestructure, noisevar, model, design, batch, r, t0) returns the
%   entire designs x_seq with observations y_seq and batch size r_seq for a
%   given true function and a give noise structure.

%   PARAMS:
%       I - the initial design
%       k - number of runs of experiments
%       m0 - size of test set 
%       d - dimension
%       budget - total budget (N)
%       fun - the true function
%       noisestructure - the distribution of noise: 'normal', 't_constdf',
%       't_heterodf', 'mixed'
%       noisevar - variance of noise: 'small' (0.1), 'large'(1), 'mixed', 
%       'hetero'
%       model - GP or t-GP
%       design - 'MCU', 'cSUR', 'ABSUR'
%       batch - 'FB', 'MLB', 'RB', 'ABSUR'
%       r - the candidate set for r levels
%       t0 - hyperparameter in ABSUR

%   Outputs:
%       x_seq - designs 
%       y_seq - noisy observations
%       r_seq - batch size
%       ee - empirical error 
%       er - error rate 
%       bias - bias 
%       metric - optimized acquisition function value 
%       t - total time 
%       Ef - posterior mean function 
%       Varf - posterior function variance 
%       l - lengthscale 
%       sigma2 - sigma2 in covariance function 
%       sigman - sigma2 in likelihood 
%       steps - design size
%       t_optim - time to optimize acquisition function 
%       t_gen - time to generate observation
%       overhead - c_over in ABSUR/ the estimated overhead to optimize the 
%       acquisition function 
%       gamma - hyperparameter in MLB and RB 
%       post_var - the one-step-ahead posterior variance of selected design

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intialize batch and r %
if (nargin < 13)
    if (strcmp(model, 'gauss'))
        if (d < 6)
            t0 = 0.01;
        else
            t0 = 0.1;
        end
    else
        if (d < 6)
            t0 = 0.005;
        else
            t0 = 0.05;
        end
    end
end

if (strcmp(batch, 'null'))
    r = 1;
end

% initialize variables %
[x_seq, y_seq, r_seq, ee, er, bias, metric, t, Ef, Varf, xt, xtt, l, sigma2, sigman, steps, t_optim, t_gen, overhead, gamma] = intgp(k, m0, d, budget, fun, r);


ft = fun(xtt);

% critical probability pcr in eq. 5.2 %
pcr = 0.1;
lambda = 0.8;

% do experiments in parellel %
parfor j = 1:k
        
    rng(j+2);

    Xint = lhsdesign(I,d);

    % update gp %
    [x_seq_tmp, y_seq_tmp, r_seq_tmp, metric_tmp, t_tmp, ee_tmp, er_tmp, bias_tmp, Ef_tmp, Varf_tmp, l_tmp, sigma2_tmp, sigman_tmp, steps(j), t_optim_tmp, t_gen_tmp, overhead_tmp, gamma_tmp] = updategp(fun, noisestructure, noisevar, Xint, xt, model, design, budget, xtt, r, ft, pcr, lambda, batch, t0);

    
    % record performance %
    l(:,:,j) = l_tmp;
    sigma2(:,j) = sigma2_tmp;
    sigman(:,j) = sigman_tmp;
    

    Ef(:,:,j) = Ef_tmp;
    Varf(:,:,j) = Varf_tmp;
    x_seq(:,:,j) = x_seq_tmp;
    y_seq(:,j) = y_seq_tmp;
    r_seq(:,j) = r_seq_tmp;
    metric(:,j) = metric_tmp;
    ee(:,j) = ee_tmp;
    er(:,j) = er_tmp;
    bias(:,j) = bias_tmp;
    t(:,j) = t_tmp;
    t_optim(:,j) = t_optim_tmp;
    t_gen(:,j) = t_gen_tmp;
    overhead(:, j) = overhead_tmp;
    
    gamma(:, j) = gamma_tmp;
end
end