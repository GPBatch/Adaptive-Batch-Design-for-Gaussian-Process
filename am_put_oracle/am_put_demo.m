clear all; close all; clc;

currentFolder = pwd;
parentFolder = fileparts(currentFolder);
addpath(parentFolder);
addpath(fullfile(parentFolder, 'others'));
addpath(fullfile(parentFolder, 'batch_design'));

% % 
% startup;
% % 
% 
% %%%%%% 1D am_put %%%%%%
% 
% n = 1;
% 
% model.option_payoff = @put_payoff;
% model.adaptive_grid_loop = 100;
% model.km_batch = 10;
% model.init_size = 10; 
% 
% model.design = 'ICU';
% model.method = 't';
% model.T = 1;
% model.dt = 0.04;
% model.dim = 1;
% model.sim_func = @sim_gbm;
% 
% model.look_ahead = 1;
% model.K = 40;
% model.x0 = 40;
% model.sigma = 0.2;
% model.r = 0.06;
% model.div = 0;
% model.lhs_rect = [25,40];
% model.cand_len = 1000;
% model.el_thresh = 0.0001;
% model.search = true;
% 
% rng(1)
% model.init_grid = lhsdesign(10,1)*15+25;
%     
% M = model.T / model.dt;
% 
% NN = 10000;
% MM = 25;
% 
% rng(1);
% 
% mygr = zeros(NN,1,MM+1);
% 
% mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
% 
% for i = 2:(MM+1)
%    mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
% end
% 
% put1d = repmat(struct('gp',[],'x',[],'y',[]), M, n);
% timeElapsed = zeros(n,1);
% nsims = zeros(n,1);
% empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
% payoff_gp = zeros(NN,n);
% 
% parfor j = 1:n
%     [put1d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
%     oos = forward_sim_policy( mygr, MM, put1d(:,j), model);
%     payoff_gp(:,j) = oos.payoff;
% end
% 
% 
% [mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 40)), mean(mean(payoff_gp))]
% 
% xt = linspace(25,40,500)';
% [m, var] = gp_pred(put1d(10).gp, put1d(10).x, put1d(10).y, xt);
% 
% figure(1), hold on;
% plot(xt,m);
% plot(xt,m+2*sqrt(var),'--');
% plot(xt,m-2*sqrt(var),'--');
% plot(put1d(10).x, put1d(10).y,'o')
% plot([25,40],[0,0],'--')
% set(gca, 'FontSize', 18)
% xlabel('x')
% ylabel('timing value')
% title(strcat(model.method, model.design))
% 
% % plot(put1d(10).x,'o');
% % set(gca, 'FontSize', 18)
% % xlabel('step')
% % ylabel('x')
% % title('MCl-MCU')
% 
% save(fullfile(parentFolder,'results/am_put_1d_r200n100/', strcat(model.method, model.design)), 'put1d', 'timeElapsed', 'nsims', 'empLoss')
% 
% saveas(gcf, strcat(parentFolder,'/results/am_put_1d_r200n100/', model.method, model.design),'png')

%%%%%% 2D am_put %%%%%%

n = 2;

model.batch = 'MFA';
model.budget = 1500;
model.r_cand = [10 20 40 80 160 320 640 1280 2560 5120 10240];

model.adaptive_grid_loop = 150;
model.look_ahead = 1;
model.init_size = 10;
model.final_runs = 0;
model.design = 'MCU';
model.method = 'gauss';

model.cand_len = 1000;
model.km_batch = 10;
model.K = 40;
model.x0 = [40,40];
model.sigma = [0.2,0.2];
model.r = 0.06;
model.div = 0;
model.T = 1;
model.dt = 0.04;
model.dim = 2;
model.sim_func = @sim_gbm;
model.option_payoff = @put_payoff;
model.search = false;
model.el_thresh = 0.0001;

% rng(4);
% lhs200 = lhsdesign(408,2);
% lhs200 = lhs200(lhs200(:,1) + lhs200(:,2) <= 1, :);
% size(lhs200)
% 
% % model.domain = lhs200;
% % 
rng(1);

lhs20 = lhsdesign(20,2);
lhs20 = 25 + 30*lhs20;
lhs20 = lhs20(lhs20(:,1) + lhs20(:,2) <= 80, :);

model.init_grid = lhs20;

% rng(1);
% 
% lhs20 = lhsdesign(20,2);
% lhs20 = 25 + 30*lhs20;
% 
% model.init_grid = lhs20;

model.lhs_rect = repmat([25,55], 2, 1);

rng(1);

NN = 160000;
MM = 25;

mygr = zeros(NN,2,MM+1);

mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);

for i = 2:(MM+1)
   mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
end

M = model.T / model.dt;

put2d = repmat(struct('gp',[],'x',[],'y',[], 'r', []), M, n);
timeElapsed = zeros(n,1);
nsims = zeros(n,1);
empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
payoff_gp = zeros(NN,n);

for j = 1:n
    [put2d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
    oos = forward_sim_policy( mygr, MM, put2d(:,j), model);
    payoff_gp(:,j) = oos.payoff;
end


[mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 40)), mean(mean(payoff_gp))]

save(fullfile(parentFolder,'results/am_put_2d/', strcat(model.method, '_', model.design, '_', model.batch)), 'put2d', 'timeElapsed', 'nsims', 'empLoss')

plot_2d_amput(put2d(15))
saveas(gcf, strcat(parentFolder,'/results/am_put_2d/', model.method, '_', model.design, '_', model.batch),'png')

% 
% 
% % % boundary matrix 1D %
% % m = 500;
% % xt = linspace(25,40,m)';
% % y = zeros(m,24);
% % 
% % for i = 1:24
% %   y(:,i) = gp_pred(put1d(i,1).gp, put1d(i,1).x, put1d(i,1).y, xt);
% % end
% % 
% % figure (2), hold on;
% % h1=pcolor(xt,(1:1:24)',y');
% % set(h1, 'edgealpha', 0), set(h1, 'facecolor', 'interp');
% % contour(xt,(1:1:24)',y', [0,0], 'r-', 'linewidth', 3);
% % xlabel('x');
% % ylabel('time step');
% % title('boundary matrix')
% % colorbar
% % set(gca, 'FontSize', 18);
% % 
% % saveas(gcf, strcat(parentFolder,'/results/am_put_1d_r200n100/boundaryMatrix_', model.method, model.design),'png')
% % 
% 
% %%%%%% 3D am_put %%%%%%
% 
% n = 10;
% 
% model.adaptive_grid_loop = 100;
% model.look_ahead = 1;
% model.init_size = 20;
% model.final_runs = 0;
% model.design = 'lhs';
% model.method = 'gauss';
% 
% model.cand_len = 1000;
% model.km_batch = 3;
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
% lhs_int = lhsdesign(20,3);
% lhs_int = 50 + 100*lhs_int;
% 
% model.init_grid = lhs_int;
% 
% model.lhs_rect = repmat([50,150], 3, 1);
% 
% 
% NN = 160000;
% MM = 9;
% 
% M = model.T / model.dt;
% 
% put3d = repmat(struct('gp',[],'x',[],'y',[]), M, n);
% timeElapsed = zeros(n,1);
% nsims = zeros(n,1);
% empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
% payoff_gp = zeros(NN,n);
% 
% parfor j = 1:n
%     [put3d(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
%     
%     rng(j);
%     
%     mygr = zeros(NN,3,MM);
% 
%     mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
% 
%     for i = 2:MM
%        mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
%     end
%     oos = forward_sim_policy( mygr, MM, put3d(:,j), model);
%     payoff_gp(:,j) = oos.payoff;
% end
% 
% tcsur = mean(payoff_gp);
% % [mean( exp(-model.r*model.T)*model.option_payoff(mygr(:,:,MM), 100)), mean(payoff_gp)]
% 
% save(fullfile(parentFolder,'results/am_put_3d/', strcat(model.method, model.design)), 'put3d', 'timeElapsed', 'nsims', 'empLoss', 'payoff_gp')
% % 
