function [xoptim, roptim, metric, t_optim, gamma] = gaoptim(gprocess, x, y, r_seq, design, r, batch, gamma, r_lower, r_upper, overhead, t_0)

% Optimizes the aquisition function and select the next design x and batch 
% size r with constrained genetic optimization algorithm (200 generations 
% and tolerance value 1e-3) and grid search

d = size(x,2);
optim_t = tic;
tau2 = 1;
% tau2 = gprocess.lik.sigma2; for estimated tau in gp
% tau2 = gprocess.lik.sigma2 * (gprocess.lik.nu + 1) / (gprocess.lik.nu -
% 1);

switch batch
    case 'null'
        [xoptim, metric] = optimMCU(gprocess, x, y);        
        roptim = r(1);
    case 'plain'
        % chooses the next design location %
        [xoptim, metric] = optimMCU(gprocess, x, y);
        % fixed batch size %
        roptim = r(1);
    case 'RB'
        eta = 0.8;
        [xoptim, metric] = optimMCU(gprocess, x, y);
        rindex = find(r == r_seq(end));
        if (rindex == size(r, 2))
            % if r reaches the upper bound, just keep it at the same level
            roptim = r(rindex);
        else
            % compare whether to stay at the current r or move to the next
            % level
            rcand = [r(rindex) r(rindex + 1)];
            sigman2_est = tau2;
            % one-step-ahead variance for different levels
            sigman2_est_by_level = sigman2_est./(rcand - 1);
            [Ef, Varf] = gp_pred(gprocess, x, y, xoptim);
            post_var = Varf.*sigman2_est_by_level./(sigman2_est_by_level + Varf);
            sigman_est_optim = sqrt(post_var);
            % if no solution, update value of gamma
            while (sigman_est_optim(1) < gamma)
                gamma = gamma*eta;
            end
            r_index = rindex - 1 + sum(sigman_est_optim > gamma);
            roptim = r(r_index);
        end
    case 'MLB'
        eta = 0.5;
        % chooses the next design location with MCU %
        [xoptim, metric] = optimMCU(gprocess, x, y);
        sigman2_est = tau2;
        sigman2_est_by_level = sigman2_est./(r - 1);
        [Ef, Varf] = gp_pred(gprocess, x, y, xoptim);
        % one-step-ahead variance for all batch levels %
        post_var = Varf.*sigman2_est_by_level./(sigman2_est_by_level + Varf);
        sigman_est_optim = sqrt(post_var);
        % updates gamma and chooses r %
        r_index = sum(sigman_est_optim > gamma);
        while (r_index == 0)
            gamma = gamma*eta;
            r_index = sum(sigman_est_optim > gamma);
        end
        roptim = r(r_index);
    case 'ABSUR'
        % chooses the next design location with ABSUR %
        opts = gaoptimset('Display','off','TolFun', 1e-6, 'Generations', 500, 'Vectorized','on','UseParallel',true);
        [xroptim, metric]= ga(@(xt) acquiFun(gprocess, x, y, xt, design, overhead, t_0), d + 1, [],[], [], [], [zeros(d,1); r_lower], [ones(d,1); r_upper],[], opts);
        xoptim = xroptim(1:(end - 1));
        roptim = ceil(xroptim(end));
    otherwise
end
metric = -metric;
t_optim = toc(optim_t);
end

function EI = acquiFun(gprocess, x, y, xt, design, overhead, t_0)
%%%%%%%%%% The acquisition function optimization with ga %%%%%%%%%%
tau2 = 1;
switch design
    case 'MCU'
        [Ef, Varf] = gp_pred(gprocess, x, y, xt);
        Ef = Ef(:,1);
        % calculate the local empirical error %
        EI = -normcdf(-abs(Ef)./sqrt(abs(Varf)));
    case 'ABSUR'
        [Ef, Varf] = gp_pred(gprocess, x, y, xt(:,1:(end - 1)));
        Ef = Ef(:,1);
        sigman2_est = tau2;
        sigmanoise = sigman2_est./xt(:,end);
        % approximate the one-step-ahead variance %
        Varfnew = sigmanoise .* Varf./(sigmanoise + Varf);
        % calculate the reduced empirical error %
        EI = ( -normcdf( -abs(Ef)./sqrt(abs(Varf)) ) + normcdf( -abs(Ef)./sqrt(abs(Varfnew)) ) )./ (xt(:,end)*t_0 + overhead);
    case 'cSUR'
        [Ef, Varf] = gp_pred(gprocess, x, y, xt);
        Ef = Ef(:,1);
        sigman2_est = tau2;
        sigmanoise = sigman2_est./xt(:,end);
        % approximate the one-step-ahead variance %
        Varfnew = sigmanoise .* Varf./(sigmanoise + Varf);
        % calculate the reduced empirical error %
        EI = -normcdf( -abs(Ef)./sqrt(abs(Varf)) ) + normcdf( -abs(Ef)./sqrt(abs(Varfnew)) );
end
end


function [xoptim, metric] = optimMCU(gprocess, x, y)
% MCU optimization with greedy search %
d = size(x, 2);
if (d == 2)
    xt = lhs(200, d); % candidate set
else
    xt = lhs(1000, d); % candidate set
end
[Ef, Varf] = gp_pred(gprocess, x, y, xt);
Ef = Ef(:,1);
% calculate the local empirical error %
EI = -normcdf(-abs(Ef)./sqrt(abs(Varf)));
% choose the location of sample %
[metric, index] = max(EI);
xoptim = xt(index,:);
end

