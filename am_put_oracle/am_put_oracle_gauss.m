clear all; close all; clc;

currentFolder = pwd;
parentFolder = fileparts(currentFolder);
addpath(parentFolder);
addpath(fullfile(parentFolder, 'others'));
addpath(fullfile(parentFolder, 'batch_design'));

startup;

% % budget = 2000; r_cand = [20 30 40 50 60 80]; adaptive_grid_loop = 100; km_batch = 20;
% %%%%%% 2D am_put %%%%%%
% 
% n = 20;
% 
% model.batch = 'MFA';
% model.budget = 2000;
% model.r_cand = [20 30 40 50 60 80]; 
% 
% adaptive_grid_loop = 100;
% model.look_ahead = 1;
% model.init_size = 20;
% model.final_runs = 0;
% model.design = 'MCU';
% model.method = 'gauss';
% 
% model.cand_len = 1000;
% model.km_batch = 20;
% model.K = 40;
% model.x0 = [40,40];
% model.sigma = [0.2,0.2];
% model.r = 0.06;
% model.div = 0;
% model.T = 1;
% model.dt = 0.04;
% model.dim = 2;
% model.sim_func = @sim_gbm;
% model.option_payoff = @put_payoff;
% model.search = false;
% model.el_thresh = 0.0001;
% 
% rng(8);
% 
% lhs20 = lhsdesign(40,2);
% lhs20 = 25 + 30*lhs20;
% lhs20 = lhs20(lhs20(:,1) + lhs20(:,2) <= 80, :);
% 
% model.init_grid = lhs20;
% 
% model.lhs_rect = repmat([25,55], 2, 1);
% 
% rng(10);
% 
% NN = 200000;
% MM = 25;
% 
% mygr = zeros(NN,2,MM+1);
% 
% mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
% 
% for i = 2:(MM+1)
%    mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
% end
% 
% M = model.T / model.dt;
% 
% put2d = repmat(struct('gp',[],'x',[],'y',[], 'r', [], 't', []), M, n);
% timeElapsed = zeros(n,1);
% nsims = zeros(n,1);
% empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
% payoff_gp = zeros(NN,n);
% 
% parfor j = 1:n
%     [put2d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
%     oos = forward_sim_policy( mygr, MM, put2d(:,j), model);
%     payoff_gp(:,j) = oos.payoff;
% end
% 
% payoff_gp_gauss_mfa = payoff_gp;
% [mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 40)), mean(mean(payoff_gp))]
% 
% save(fullfile(parentFolder,'results/am_put_2d/', strcat(model.method, '_', model.design, '_', model.batch)), 'put2d', 'timeElapsed', 'nsims', 'empLoss', 'payoff_gp')
% 
% plot_2d_amput(put2d(15))
% saveas(gcf, strcat(parentFolder,'/results/am_put_2d/plots/', model.method, '_', model.design, '_', model.batch),'png')
% % % 
% % % % mean_payoff: gp_lhs: 1.2270 gp_mfa: 1.4419
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%% 2D am_put %%%%%%
% 
% n = 20;
% 
% model.batch = 'plain';
% model.budget = 2000;
% model.r_cand = [20 30 40 50 60 80]; 
% 
% adaptive_grid_loop = 100;
% model.look_ahead = 1;
% model.init_size = 20;
% model.final_runs = 0;
% model.design = 'MCU';
% model.method = 'gauss';
% 
% model.cand_len = 1000;
% model.km_batch = 20; 
% model.K = 40;
% model.x0 = [40,40];
% model.sigma = [0.2,0.2];
% model.r = 0.06;
% model.div = 0;
% model.T = 1;
% model.dt = 0.04;
% model.dim = 2;
% model.sim_func = @sim_gbm;
% model.option_payoff = @put_payoff;
% model.search = false;
% model.el_thresh = 0.0001;
% 
% rng(8);
% 
% lhs20 = lhsdesign(40,2);
% lhs20 = 25 + 30*lhs20;
% lhs20 = lhs20(lhs20(:,1) + lhs20(:,2) <= 80, :);
% 
% model.init_grid = lhs20;
% 
% model.lhs_rect = repmat([25,55], 2, 1);
% 
% rng(10);
% 
% NN = 200000;
% MM = 25;
% 
% mygr = zeros(NN,2,MM+1);
% 
% mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
% 
% for i = 2:(MM+1)
%    mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
% end
% 
% M = model.T / model.dt;
% 
% put2d = repmat(struct('gp',[],'x',[],'y',[], 'r', [], 't', []), M, n);
% timeElapsed = zeros(n,1);
% nsims = zeros(n,1);
% empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
% payoff_gp = zeros(NN,n);
% 
% parfor j = 1:n
%     [put2d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
%     oos = forward_sim_policy( mygr, MM, put2d(:,j), model);
%     payoff_gp(:,j) = oos.payoff;
% end
% 
% payoff_gp_gauss_plain = payoff_gp;
% [mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 40)), mean(mean(payoff_gp))]
% 
% save(fullfile(parentFolder,'results/am_put_2d/', strcat(model.method, '_', model.design, '_', model.batch)), 'put2d', 'timeElapsed', 'nsims', 'empLoss', 'payoff_gp')
% 
% plot_2d_amput(put2d(15))
% saveas(gcf, strcat(parentFolder,'/results/am_put_2d/', model.method, '_', model.design, '_', model.batch),'png')
% % % 1.2270    1.4475 %
% 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% 2D am_put %%%%%%
% 
% n = 20;
% 
% model.batch = 'aSUR';
% model.budget = 2000;
% model.r_cand = [20 30 40 50 60 80]; 
% 
% adaptive_grid_loop = 100;
% model.look_ahead = 1;
% model.init_size = 20;
% model.final_runs = 0;
% model.design = 'aSUR';
% model.method = 'gauss';
% 
% model.cand_len = 1000;
% model.km_batch = 20; 
% model.K = 40;
% model.x0 = [40,40];
% model.sigma = [0.2,0.2];
% model.r = 0.06;
% model.div = 0;
% model.T = 1;
% model.dt = 0.04;
% model.dim = 2;
% model.sim_func = @sim_gbm;
% model.option_payoff = @put_payoff;
% model.search = false;
% model.el_thresh = 0.0001;
% 
% rng(8);
% 
% lhs20 = lhsdesign(40,2);
% lhs20 = 25 + 30*lhs20;
% lhs20 = lhs20(lhs20(:,1) + lhs20(:,2) <= 80, :);
% 
% model.init_grid = lhs20;
% 
% model.lhs_rect = repmat([25,55], 2, 1);
% 
% rng(10);
% 
% NN = 200000;
% MM = 25;
% 
% mygr = zeros(NN,2,MM+1);
% 
% mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
% 
% for i = 2:(MM+1)
%    mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
% end
% 
% M = model.T / model.dt;
% 
% put2d = repmat(struct('gp',[],'x',[],'y',[], 'r', [], 't', []), M, n);
% timeElapsed = zeros(n,1);
% nsims = zeros(n,1);
% empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
% payoff_gp = zeros(NN,n);
% 
% parfor j = 1:n
%     [put2d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
%     oos = forward_sim_policy( mygr, MM, put2d(:,j), model);
%     payoff_gp(:,j) = oos.payoff;
% end
% 
% payoff_gp_gauss_asur = payoff_gp;
% [mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 40)), mean(mean(payoff_gp))]
% 
% save(fullfile(parentFolder,'results/am_put_2d/', strcat(model.method, '_', model.design, '_', model.batch)), 'put2d', 'timeElapsed', 'nsims', 'empLoss', 'payoff_gp')
% 
% plot_2d_amput(put2d(15))
% saveas(gcf, strcat(parentFolder,'/results/am_put_2d/', model.method, '_', model.design, '_', model.batch),'png')
% % 
% % % mean_payoff: gp_lhs: 1.2270 gp_asur: 1.4509
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% 2D am_put %%%%%%
% 
% n = 20;
% 
% model.batch = 'MFA_NoConstraint';
% model.budget = 2000;
% model.r_cand = [20 30 40 50 60 80]; 
% 
% model.adaptive_grid_loop = 100;
% model.look_ahead = 1;
% model.init_size = 20;
% model.final_runs = 0;
% model.design = 'MCU';
% model.method = 'gauss';
% 
% model.cand_len = 1000;
% model.km_batch = 20; 
% model.K = 40;
% model.x0 = [40,40];
% model.sigma = [0.2,0.2];
% model.r = 0.06;
% model.div = 0;
% model.T = 1;
% model.dt = 0.04;
% model.dim = 2;
% model.sim_func = @sim_gbm;
% model.option_payoff = @put_payoff;
% model.search = false;
% model.el_thresh = 0.0001;
% 
% % rng(1);
% % 
% % lhs20 = lhsdesign(20,2);
% % lhs20 = 25 + 30*lhs20;
% % lhs20 = lhs20(lhs20(:,1) + lhs20(:,2) <= 80, :);
% 
% rng(8);
% 
% lhs20 = lhsdesign(40,2);
% lhs20 = 25 + 30*lhs20;
% lhs20 = lhs20(lhs20(:,1) + lhs20(:,2) <= 80, :);
% 
% 
% model.init_grid = lhs20;
% 
% model.lhs_rect = repmat([25,55], 2, 1);
% 
% rng(10);
% 
% NN = 200000;
% MM = 25;
% 
% mygr = zeros(NN,2,MM+1);
% 
% mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
% 
% for i = 2:(MM+1)
%    mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
% end
% 
% M = model.T / model.dt;
% 
% put2d = repmat(struct('gp',[],'x',[],'y',[], 'r', [], 't', []), M, n);
% timeElapsed = zeros(n,1);
% nsims = zeros(n,1);
% empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
% payoff_gp = zeros(NN,n);
% 
% parfor j = 1:n
%     [put2d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
%     oos = forward_sim_policy( mygr, MM, put2d(:,j), model);
%     payoff_gp(:,j) = oos.payoff;
% end
% 
% payoff_gp_gauss_freemfa = payoff_gp;
% [mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 40)), mean(mean(payoff_gp))]
% 
% save(fullfile(parentFolder,'results/am_put_2d/', strcat(model.method, '_', model.design, '_', model.batch)), 'put2d', 'timeElapsed', 'nsims', 'empLoss', 'payoff_gp')
% 
% plot_2d_amput(put2d(15))
% saveas(gcf, strcat(parentFolder,'/results/am_put_2d/', model.method, '_', model.design, '_', model.batch),'png')
% % 
% % % % mean_payoff: gp_lhs: 1.2270 gp_freemfa: 1.4411
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% 2D am_put %%%%%%
% 
% n = 20;
% 
% model.batch = 'plain';
% model.budget = 2000;
% model.r_cand = [20 30 40 50 60 80]; 
% 
% model.adaptive_grid_loop = 100;
% model.look_ahead = 1;
% model.init_size = 10;
% model.final_runs = 0;
% model.design = 'lhs';
% model.method = 'gauss';
% 
% model.cand_len = 1000;
% model.km_batch = 20; 
% model.K = 40;
% model.x0 = [40,40];
% model.sigma = [0.2,0.2];
% model.r = 0.06;
% model.div = 0;
% model.T = 1;
% model.dt = 0.04;
% model.dim = 2;
% model.sim_func = @sim_gbm;
% model.option_payoff = @put_payoff;
% model.search = false;
% model.el_thresh = 0.0001;
% 
% rng(1);
% 
% lhs20 = lhsdesign(20,2);
% lhs20 = 25 + 30*lhs20;
% lhs20 = lhs20(lhs20(:,1) + lhs20(:,2) <= 80, :);
% 
% model.init_grid = lhs20;
% 
% model.lhs_rect = repmat([25,55], 2, 1);
% 
% rng(10);
% 
% NN = 200000;
% MM = 25;
% 
% mygr = zeros(NN,2,MM+1);
% 
% mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
% 
% for i = 2:(MM+1)
%    mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
% end
% 
% M = model.T / model.dt;
% 
% put2d = repmat(struct('gp',[],'x',[],'y',[], 'r', [], 't', []), M, n);
% timeElapsed = zeros(n,1);
% nsims = zeros(n,1);
% empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
% payoff_gp = zeros(NN,n);
% 
% parfor j = 1:n
%     [put2d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
%     oos = forward_sim_policy( mygr, MM, put2d(:,j), model);
%     payoff_gp(:,j) = oos.payoff;
% end
% 
% payoff_gp_gauss_lhs = payoff_gp;
% [mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 40)), mean(mean(payoff_gp))]
% 
% save(fullfile(parentFolder,'results/am_put_2d/', strcat(model.method, '_', model.design, '_', model.batch)), 'put2d', 'timeElapsed', 'nsims', 'empLoss', 'payoff_gp')
% 
% plot_2d_amput(put2d(15))
% saveas(gcf, strcat(parentFolder,'/results/am_put_2d/', model.method, '_', model.design, '_', model.batch),'png')
% 
% % mean_payoff: gp_lhs: 1.2270 gp_lhs: 1.4397


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3d %%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

currentFolder = pwd;
parentFolder = fileparts(currentFolder);
addpath(parentFolder);
addpath(fullfile(parentFolder, 'others'));
addpath(fullfile(parentFolder, 'batch_design'));

startup;

% %%%%%% 2D am_put %%%%%%
% 
% n = 2;
% 
% model.batch = 'plain';
% model.budget = 12000;
% model.r_cand = [5 10 20 30 40 50];
% 
% model.adaptive_grid_loop = 150;
% model.look_ahead = 1;
% model.init_size = 50;
% model.final_runs = 0;
% model.design = 'lhs';
% model.method = 'gauss';
% 
% model.cand_len = 1000;
% model.km_batch = 80;
% model.K = 100;
% model.x0 = [90,90];
% model.sigma = [0.2,0.2];
% model.r = 0.05;
% model.div = 0.1;
% model.T = 3;
% model.dt = 1/3;
% model.dim = 2;
% model.sim_func = @sim_gbm;
% model.option_payoff = @maxCall;
% model.search = false;
% model.el_thresh = 0.0001;
% 
% rng(1);
% 
% lhs_int = lhsdesign(model.init_size,model.dim);
% lhs_int = 50 + 100*lhs_int;
% 
% model.init_grid = lhs_int;
% 
% model.lhs_rect = repmat([50,150], model.dim, 1);
% 
% 
% NN = 200000;
% MM = 9;
% 
% rng(10);
%     
% mygr = zeros(NN,model.dim,MM);
% 
% mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
% 
% for i = 2:MM
%    mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
% end
%     
% M = model.T / model.dt;
% 
% put3d = repmat(struct('gp',[],'x',[],'y',[], 'r', [], 't', []), M, n);
% timeElapsed = zeros(n,1);
% nsims = zeros(n,1);
% empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
% payoff_gp = zeros(NN,n);
% 
% parfor j = 1:n
%     [put3d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
%     oos = forward_sim_policy( mygr, MM, put3d(:,j), model);
%     payoff_gp(:,j) = oos.payoff;
% end
% 
% plot_2d_amcall(put3d(5,1))
% saveas(gcf, strcat(parentFolder,'/results/am_put_3d/', model.method, '_', model.design, '_', model.batch),'png')
% 
% payoff_gp_gauss_lhs = payoff_gp;
% [mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 40)), mean(mean(payoff_gp))]
% % 
% save(fullfile(parentFolder,'results/am_put_3d/', strcat(model.method, '_', model.design, '_', model.batch)), 'put3d', 'timeElapsed', 'nsims', 'empLoss', 'payoff_gp')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%% 2D am_put %%%%%%
% 
% n = 2;
% 
% model.batch = 'plain';
% model.budget = 12000;
% model.r_cand = [5 10 20 30 40 50];
% 
% model.adaptive_grid_loop = 150;
% model.look_ahead = 1;
% model.init_size = 50;
% model.final_runs = 0;
% model.design = 'cSUR';
% model.method = 'gauss';
% 
% model.cand_len = 1000;
% model.km_batch = 80;
% model.K = 100;
% model.x0 = [90,90];
% model.sigma = [0.2,0.2];
% model.r = 0.05;
% model.div = 0.1;
% model.T = 3;
% model.dt = 1/3;
% model.dim = 2;
% model.sim_func = @sim_gbm;
% model.option_payoff = @maxCall;
% model.search = false;
% model.el_thresh = 0.0001;
% 
% rng(1);
% 
% lhs_int = lhsdesign(model.init_size,model.dim);
% lhs_int = 50 + 100*lhs_int;
% 
% model.init_grid = lhs_int;
% 
% model.lhs_rect = repmat([50,150], model.dim, 1);
% 
% 
% NN = 200000;
% MM = 9;
% 
% rng(10);
%     
% mygr = zeros(NN,model.dim,MM);
% 
% mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
% 
% for i = 2:MM
%    mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
% end
%     
% M = model.T / model.dt;
% 
% put3d = repmat(struct('gp',[],'x',[],'y',[], 'r', [], 't', []), M, n);
% timeElapsed = zeros(n,1);
% nsims = zeros(n,1);
% empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
% payoff_gp = zeros(NN,n);
% 
% parfor j = 1:n
%     [put3d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
%     oos = forward_sim_policy( mygr, MM, put3d(:,j), model);
%     payoff_gp(:,j) = oos.payoff;
% end
% 
% plot_2d_amcall(put3d(5,1))
% saveas(gcf, strcat(parentFolder,'/results/am_put_3d/', model.method, '_', model.design, '_', model.batch),'png')
% 
% payoff_gp_gauss_plain = payoff_gp;
% [mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 40)), mean(mean(payoff_gp))]
% 
% save(fullfile(parentFolder,'results/am_put_3d/', strcat(model.method, '_', model.design, '_', model.batch)), 'put3d', 'timeElapsed', 'nsims', 'empLoss', 'payoff_gp')
% 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%% 2D am_put %%%%%%
% 
% n = 2;
% 
% model.batch = 'MFA_NoConstraint';
% model.budget = 12000;
% model.r_cand = [40 50 60 70 80 100 120];
% 
% model.adaptive_grid_loop = 150;
% model.look_ahead = 1;
% model.init_size = 50;
% model.final_runs = 0;
% model.design = 'MCU';
% model.method = 'gauss';
% 
% model.cand_len = 1000;
% model.km_batch = 80;
% model.K = 100;
% model.x0 = [90,90];
% model.sigma = [0.2,0.2];
% model.r = 0.05;
% model.div = 0.1;
% model.T = 3;
% model.dt = 1/3;
% model.dim = 2;
% model.sim_func = @sim_gbm;
% model.option_payoff = @maxCall;
% model.search = false;
% model.el_thresh = 0.0001;
% 
% rng(1);
% 
% lhs_int = lhsdesign(model.init_size,model.dim);
% lhs_int = 50 + 100*lhs_int;
% 
% model.init_grid = lhs_int;
% 
% model.lhs_rect = repmat([50,150], model.dim, 1);
% 
% 
% NN = 200000;
% MM = 9;
% 
% rng(10);
%     
% mygr = zeros(NN,model.dim,MM);
% 
% mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
% 
% for i = 2:MM
%    mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
% end
%     
% M = model.T / model.dt;
% 
% put3d = repmat(struct('gp',[],'x',[],'y',[], 'r', [], 't', []), M, n);
% timeElapsed = zeros(n,1);
% nsims = zeros(n,1);
% empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
% payoff_gp = zeros(NN,n);
% 
% parfor j = 1:n
%     [put3d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
%     oos = forward_sim_policy( mygr, MM, put3d(:,j), model);
%     payoff_gp(:,j) = oos.payoff;
% end
% 
% plot_2d_amcall(put3d(5,1))
% saveas(gcf, strcat(parentFolder,'/results/am_put_3d/', model.method, '_', model.design, '_', model.batch),'png')
% 
% payoff_gp_gauss_freemfa = payoff_gp;
% [mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 40)), mean(mean(payoff_gp))]
% 
% save(fullfile(parentFolder,'results/am_put_3d/', strcat(model.method, '_', model.design, '_', model.batch)), 'put3d', 'timeElapsed', 'nsims', 'empLoss', 'payoff_gp')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%% 2D am_put %%%%%%
% 
% n = 2;
% 
% model.batch = 'aSUR';
% model.budget = 12000;
% model.r_cand = [40 50 60 70 80 100 120];
% 
% model.adaptive_grid_loop = 150;
% model.look_ahead = 1;
% model.init_size = 50;
% model.final_runs = 0;
% model.design = 'aSUR';
% model.method = 'gauss';
% 
% model.cand_len = 1000;
% model.km_batch = 80;
% model.K = 100;
% model.x0 = [90,90];
% model.sigma = [0.2,0.2];
% model.r = 0.05;
% model.div = 0.1;
% model.T = 3;
% model.dt = 1/3;
% model.dim = 2;
% model.sim_func = @sim_gbm;
% model.option_payoff = @maxCall;
% model.search = false;
% model.el_thresh = 0.0001;
% 
% rng(1);
% 
% lhs_int = lhsdesign(model.init_size,model.dim);
% lhs_int = 50 + 100*lhs_int;
% 
% model.init_grid = lhs_int;
% 
% model.lhs_rect = repmat([50,150], model.dim, 1);
% 
% 
% NN = 200000;
% MM = 9;
% 
% rng(10);
%     
% mygr = zeros(NN,model.dim,MM);
% 
% mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
% 
% for i = 2:MM
%    mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
% end
%     
% M = model.T / model.dt;
% 
% put3d = repmat(struct('gp',[],'x',[],'y',[], 'r', [], 't', []), M, n);
% timeElapsed = zeros(n,1);
% nsims = zeros(n,1);
% empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
% payoff_gp = zeros(NN,n);
% 
% parfor j = 1:n
%     [put3d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
%     oos = forward_sim_policy( mygr, MM, put3d(:,j), model);
%     payoff_gp(:,j) = oos.payoff;
% end
% 
% plot_2d_amcall(put3d(5,1))
% saveas(gcf, strcat(parentFolder,'/results/am_put_3d/', model.method, '_', model.design, '_', model.batch),'png')
% 
% payoff_gp_gauss_asur = payoff_gp;
% [mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 40)), mean(mean(payoff_gp))]
% 
% save(fullfile(parentFolder,'results/am_put_3d/', strcat(model.method, '_', model.design, '_', model.batch)), 'put3d', 'timeElapsed', 'nsims', 'empLoss', 'payoff_gp')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 2D am_put %%%%%%

% % v1: N = 20000; R = 80; n = 200; n_0 = 100; model.r_cand = [40 50 60 70 80 100 120 140];
% % v2: N = 20000; R = 40; n = 400; n_0 = 150; model.r_cand = [20 30 40 50 60 70 80 100 120];
% %%%%%% 2D am_put %%%%%%
% 
% n = 20;
% 
% model.batch = 'plain';
% model.budget = 30000;
% model.r_cand = [30 40 50 80 120 160 200];
% 
% model.adaptive_grid_loop = 1000;
% model.look_ahead = 1;
% model.init_size = 600;
% model.final_runs = 0;
% model.design = 'lhs';
% model.method = 'gauss';
% 
% model.cand_len = 1000;
% model.km_batch = 30;
% model.K = 100;
% model.x0 = [90,90,90];
% model.sigma = [0.2,0.2,0.2];
% model.r = 0.05;
% model.div = 0.1;
% model.T = 3;
% model.dt = 1/3;
% model.dim = 3;
% model.sim_func = @sim_gbm;
% model.option_payoff = @maxCall;
% model.search = false;
% model.el_thresh = 0.0001;
% 
% rng(1);
% 
% lhs_int = lhsdesign(model.init_size,3);
% lhs_int = 50 + 100*lhs_int;
% 
% model.init_grid = lhs_int;
% 
% model.lhs_rect = repmat([50,150], model.dim, 1);
% 
% rng(1);
% 
% lhs_int = lhsdesign(model.init_size,3);
% lhs_int = 50 + 100*lhs_int;
% 
% model.init_grid = lhs_int;
% 
% model.lhs_rect = repmat([50,150], 3, 1);
% 
% 
% NN = 200000;
% MM = 9;
% 
% rng(10);
%     
% mygr = zeros(NN,3,MM);
% 
% mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
% 
% for i = 2:MM
%    mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
% end
%     
% M = model.T / model.dt;
% 
% put3d = repmat(struct('gp',[],'x',[],'y',[], 'r', [], 't', []), M, n);
% timeElapsed = zeros(n,1);
% nsims = zeros(n,1);
% empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
% payoff_gp = zeros(NN,n);
% 
% parfor j = 1:n
%     [put3d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
%     oos = forward_sim_policy( mygr, MM, put3d(:,j), model);
%     payoff_gp(:,j) = oos.payoff;
% end
% 
% payoff_gp_gauss_lhs = payoff_gp;
% [mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 40)), mean(mean(payoff_gp))]
% 
% save(fullfile(parentFolder,'results/am_put_3d/', strcat(model.method, '_', model.design, '_', model.batch)), 'put3d', 'timeElapsed', 'nsims', 'empLoss', 'payoff_gp')
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%% 2D am_put %%%%%%
% 
% n = 20;
% 
% model.batch = 'plain';
% model.budget = 30000;
% % model.r_cand = [40 50 60 70 80 100 120];
% model.r_cand = [30 40 50 80 120 160 200];
% 
% model.adaptive_grid_loop = 1000;
% model.look_ahead = 1;
% model.init_size = 600;
% model.final_runs = 0;
% model.design = 'MCU';
% model.method = 'gauss';
% 
% model.cand_len = 1000;
% model.km_batch = 30;
% model.K = 100;
% model.x0 = [90,90,90];
% model.sigma = [0.2,0.2,0.2];
% model.r = 0.05;
% model.div = 0.1;
% model.T = 3;
% model.dt = 1/3;
% model.dim = 3;
% model.sim_func = @sim_gbm;
% model.option_payoff = @maxCall;
% model.search = false;
% model.el_thresh = 0.0001;
% 
% rng(1);
% 
% lhs_int = lhsdesign(model.init_size,3);
% lhs_int = 50 + 100*lhs_int;
% 
% model.init_grid = lhs_int;
% 
% model.lhs_rect = repmat([50,150], model.dim, 1);
% 
% rng(1);
% 
% lhs_int = lhsdesign(model.init_size,3);
% lhs_int = 50 + 100*lhs_int;
% 
% model.init_grid = lhs_int;
% 
% model.lhs_rect = repmat([50,150], 3, 1);
% 
% 
% NN = 200000;
% MM = 9;
% 
% rng(10);
%     
% mygr = zeros(NN,3,MM);
% 
% mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
% 
% for i = 2:MM
%    mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
% end
%     
% M = model.T / model.dt;
% 
% put3d = repmat(struct('gp',[],'x',[],'y',[], 'r', [], 't', []), M, n);
% timeElapsed = zeros(n,1);
% nsims = zeros(n,1);
% empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
% payoff_gp = zeros(NN,n);
% 
% parfor j = 1:n
%     [put3d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
%     oos = forward_sim_policy( mygr, MM, put3d(:,j), model);
%     payoff_gp(:,j) = oos.payoff;
% end
% 
% payoff_gp_gauss_plain = payoff_gp;
% [mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 40)), mean(mean(payoff_gp))]
% 
% save(fullfile(parentFolder,'results/am_put_3d/', strcat(model.method, '_', model.design, '_', model.batch)), 'put3d', 'timeElapsed', 'nsims', 'empLoss', 'payoff_gp')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%% 2D am_put %%%%%%

n = 20;

model.batch = 'aSUR';
model.budget = 30000;
% model.r_cand = [20 30 40 50 60 70 80 100];
% model.r_cand = [30 40 50 80 120 160 200];
model.r_cand = [30 50 80 120 160 200 280 400];

model.adaptive_grid_loop = 1000;
model.look_ahead = 1;
model.init_size = 600;
model.final_runs = 0;
model.design = 'aSUR';
model.method = 'gauss';

model.cand_len = 1000;
model.km_batch = 30;
model.K = 100;
model.x0 = [90,90,90];
model.sigma = [0.2,0.2,0.2];
model.r = 0.05;
model.div = 0.1;
model.T = 3;
model.dt = 1/3;
model.dim = 3;
model.sim_func = @sim_gbm;
model.option_payoff = @maxCall;
model.search = false;
model.el_thresh = 0.0001;

rng(1);

lhs_int = lhsdesign(model.init_size,3);
lhs_int = 50 + 100*lhs_int;

model.init_grid = lhs_int;

model.lhs_rect = repmat([50,150], model.dim, 1);

rng(1);

lhs_int = lhsdesign(model.init_size,3);
lhs_int = 50 + 100*lhs_int;

model.init_grid = lhs_int;

model.lhs_rect = repmat([50,150], 3, 1);


NN = 200000;
MM = 9;

rng(10);
    
mygr = zeros(NN,3,MM);

mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);

for i = 2:MM
   mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
end
    
M = model.T / model.dt;

put3d = repmat(struct('gp',[],'x',[],'y',[], 'r', [], 't', []), M, n);
timeElapsed = zeros(n,1);
nsims = zeros(n,1);
empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
payoff_gp = zeros(NN,n);

parfor j = 1:n
    [put3d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
    oos = forward_sim_policy( mygr, MM, put3d(:,j), model);
    payoff_gp(:,j) = oos.payoff;
end

payoff_gp_gauss_asur = payoff_gp;
[mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 40)), mean(mean(payoff_gp))]

save(fullfile(parentFolder,'results/am_put_3d/', strcat(model.method, '_', model.design, '_', model.batch)), 'put3d', 'timeElapsed', 'nsims', 'empLoss', 'payoff_gp')
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%% 2D am_put %%%%%%
% 
% n = 20;
% 
% model.batch = 'MFA';
% model.budget = 30000;
% % model.r_cand = [20 30 40 50 60 70 80 100];
% model.r_cand = [30 50 80 120 160 200 280 400];
% 
% model.adaptive_grid_loop = 1000;
% model.look_ahead = 1;
% model.init_size = 600;
% model.final_runs = 0;
% model.design = 'MCU';
% model.method = 'gauss';
% 
% model.cand_len = 1000;
% model.km_batch = 30;
% model.K = 100;
% model.x0 = [90,90,90];
% model.sigma = [0.2,0.2,0.2];
% model.r = 0.05;
% model.div = 0.1;
% model.T = 3;
% model.dt = 1/3;
% model.dim = 3;
% model.sim_func = @sim_gbm;
% model.option_payoff = @maxCall;
% model.search = false;
% model.el_thresh = 0.0001;
% 
% rng(1);
% 
% lhs_int = lhsdesign(model.init_size,3);
% lhs_int = 50 + 100*lhs_int;
% 
% model.init_grid = lhs_int;
% 
% model.lhs_rect = repmat([50,150], model.dim, 1);
% 
% rng(1);
% 
% lhs_int = lhsdesign(model.init_size,3);
% lhs_int = 50 + 100*lhs_int;
% 
% model.init_grid = lhs_int;
% 
% model.lhs_rect = repmat([50,150], 3, 1);
% 
% 
% NN = 200000;
% MM = 9;
% 
% rng(10);
%     
% mygr = zeros(NN,3,MM);
% 
% mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
% 
% for i = 2:MM
%    mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
% end
%     
% M = model.T / model.dt;
% 
% put3d = repmat(struct('gp',[],'x',[],'y',[], 'r', [], 't', []), M, n);
% timeElapsed = zeros(n,1);
% nsims = zeros(n,1);
% empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
% payoff_gp = zeros(NN,n);
% 
% parfor j = 1:n
%     [put3d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
%     oos = forward_sim_policy( mygr, MM, put3d(:,j), model);
%     payoff_gp(:,j) = oos.payoff;
% end
% 
% payoff_gp_gauss_mfa = payoff_gp;
% [mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 40)), mean(mean(payoff_gp))]
% 
% save(fullfile(parentFolder,'results/am_put_3d/', strcat(model.method, '_', model.design, '_', model.batch)), 'put3d', 'timeElapsed', 'nsims', 'empLoss', 'payoff_gp')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%% 2D am_put %%%%%%
% 
% n = 20;
% 
% model.batch = 'MFA_NoConstraint';
% model.budget = 30000;
% % model.r_cand = [20 30 40 50 60 70 80 100];
% model.r_cand = [30 50 80 120 160 200 280 400];
% 
% model.adaptive_grid_loop = 1000;
% model.look_ahead = 1;
% model.init_size = 600;
% model.final_runs = 0;
% model.design = 'MCU';
% model.method = 'gauss';
% 
% model.cand_len = 1000;
% model.km_batch = 30;
% model.K = 100;
% model.x0 = [90,90,90];
% model.sigma = [0.2,0.2,0.2];
% model.r = 0.05;
% model.div = 0.1;
% model.T = 3;
% model.dt = 1/3;
% model.dim = 3;
% model.sim_func = @sim_gbm;
% model.option_payoff = @maxCall;
% model.search = false;
% model.el_thresh = 0.0001;
% 
% rng(1);
% 
% lhs_int = lhsdesign(model.init_size,3);
% lhs_int = 50 + 100*lhs_int;
% 
% model.init_grid = lhs_int;
% 
% model.lhs_rect = repmat([50,150], model.dim, 1);
% 
% rng(1);
% 
% lhs_int = lhsdesign(model.init_size,3);
% lhs_int = 50 + 100*lhs_int;
% 
% model.init_grid = lhs_int;
% 
% model.lhs_rect = repmat([50,150], 3, 1);
% 
% 
% NN = 200000;
% MM = 9;
% 
% rng(10);
%     
% mygr = zeros(NN,3,MM);
% 
% mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
% 
% for i = 2:MM
%    mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
% end
%     
% M = model.T / model.dt;
% 
% put3d = repmat(struct('gp',[],'x',[],'y',[], 'r', [], 't', []), M, n);
% timeElapsed = zeros(n,1);
% nsims = zeros(n,1);
% empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
% payoff_gp = zeros(NN,n);
% 
% parfor j = 1:n
%     [put3d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
%     oos = forward_sim_policy( mygr, MM, put3d(:,j), model);
%     payoff_gp(:,j) = oos.payoff;
% end
% 
% payoff_gp_gauss_freemfa = payoff_gp;
% [mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 40)), mean(mean(payoff_gp))]
% 
% save(fullfile(parentFolder,'results/am_put_3d/', strcat(model.method, '_', model.design, '_', model.batch)), 'put3d', 'timeElapsed', 'nsims', 'empLoss', 'payoff_gp')
% % % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
