% AM_PUT_DEMO the main file for American option examples. The parameters in
% model list can be changed by user. Main function to calculate the
% stopping criteria is osp_seq_design.

% Description:
%   Calculates the stoppint criteria for the whole path in American
%   put/call option using Gaussian Process (Gaussian or t likelihood) with
%   adaptive batch design.
%   This is an example for 2D Average basket American put option example
%   with Gaussian GP and ABSUR batch design algorithm.

clear all; close all; clc;

startup;

currentFolder = pwd;
parentFolder = fileparts(currentFolder);
addpath(parentFolder);
addpath(fullfile(parentFolder, 'batch_design'));

% Initializes the parameters for Gaussian Processes and adaptive batch
% design
model.batch = 'ABSUR';
model.budget = 2000;
model.r_cand = [20 30 40 50 60 80]; 
model.adaptive_grid_loop = 100;
model.look_ahead = 1;
model.init_size = 20;
model.final_runs = 0;
model.design = 'ABSUR';
model.method = 'gauss';
model.km_batch = 20;

% Initializes the parameters for American put option (Can be changed
% accordingly for call option)
model.K = 40;
model.x0 = [40,40];
model.sigma = [0.2,0.2];
model.r = 0.06;
model.div = 0;
model.T = 1;
model.dt = 0.04;
model.dim = 2;
model.sim_func = @sim_gbm;
model.option_payoff = @put_payoff; % Change to maxCall for a call option

% Initializes the initial design
rng(8);

lhs20 = lhsdesign(40,2);
lhs20 = 25 + 30*lhs20;
lhs20 = lhs20(lhs20(:,1) + lhs20(:,2) <= 80, :);
model.init_grid = lhs20;
model.lhs_rect = repmat([25,55], 2, 1);

% Calculates the stopping criteria for the whole path
[put2d, timeElapsed, nsims, empLoss] = osp_seq_design(model);

% Calculates the payoff
NN = 200000;
MM = 25;
mygr = zeros(NN,2,MM+1); % test set
mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
for i = 2:(MM+1)
   mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
end
oos = forward_sim_policy( mygr, MM, put2d, model);
payoff_gp = oos.payoff;

% Validate plot
plot_2d_amput(put2d(15))