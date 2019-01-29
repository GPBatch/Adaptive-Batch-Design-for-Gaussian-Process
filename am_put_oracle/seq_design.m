function [x_optim, r_optim, gamma, metric] = seq_design(fit, x, y, var, r, model, lhs_rect, i, gamma)

% Optimizes the acquisition function in sequential design
% Description:
%   SEQ_DESIGN optimizes the acquisition function in sequential design for
%   American put/call option
%   [X_OPTIM, R_OPTIM,...] = SEQ_DESIGN(FIT, X, Y, VAR, R, MODEL, LHS_RECT,
%   I, GAMMA) returns the optimized sampling location and batch size in
%   sequential design
%       FIT - the GP fit used to calculate the acquisition function
%       X - inputs
%       Y - observations
%       VAR - variance of observations
%       R - batch size
%       MODEL - the GP likelihood (now supports Gaussian and t)
%       LHS_RECT - the edge of design space
%       I - the time stamp in American put/call option
%       GAMMA - parameter in MLB and RB

gprocess = fit.gp; % fitted GP
design = model.design; % sequantial design algorithm
r_cand = model.r_cand; % candidate batch size in RB and MLB
d = size(x,2); % dimension
batch = model.batch; % adaptive batch design algorithm
tau2 = gprocess.lik.sigma2; % current noise varianceﬂ

switch batch
    case 'FB'
        [x_optim, metric] = optimMCU(gprocess, x, y, model, lhs_rect);                
        r_optim = model.km_batch;
    case 'RB'
        eta = 0.8;
        [x_optim, metric] = optimMCU(gprocess, x, y, model, lhs_rect);
        rindex = find(r_cand == r(end));
        if (rindex == size(r_cand, 2))
            % if r reaches the upper bound, just keep it at the same level
            r_optim = r_cand(rindex);
        else
            % compare whether to stay at the current r or move to the next
            % level
            rcand = [r_cand(rindex) r_cand(rindex + 1)];
            sigman2_est = tau2;
            % one-step-ahead variance for different levels
            sigman2_est_by_level = sigman2_est./(rcand - 1);
            [Ef, Varf] = gp_pred(gprocess, x, y, x_optim);
            post_var = Varf.*sigman2_est_by_level./(sigman2_est_by_level + Varf);
            sigman_est_optim = sqrt(post_var);
            % if no solution, update value of gamma
            while (sigman_est_optim(1) < gamma)
                gamma = gamma*eta;
            end
            r_index = rindex - 1 + sum(sigman_est_optim > gamma);
            r_optim = r_cand(r_index);
        end
    case 'MLB'
        eta = 0.5;
        % chooses the next design location with MCU %
        [x_optim, metric] = optimMCU(gprocess, x, y, model, lhs_rect);
        sigman2_est = tau2;
        sigman2_est_by_level = sigman2_est./(r_cand - 1);
        [Ef, Varf] = gp_pred(gprocess, x, y, x_optim);
        % one-step-ahead variance for all batch levels %
        post_var = Varf.*sigman2_est_by_level./(sigman2_est_by_level + Varf);
        sigman_est_optim = sqrt(post_var);
        % updates gamma and chooses r %
        r_index = sum(sigman_est_optim > gamma);
        while (r_index == 0)
            gamma = gamma*eta;
            r_index = sum(sigman_est_optim > gamma);
        end
        r_optim = r_cand(r_index);
    case 'ABSUR'
        opts = gaoptimset('Display','off','TolFun', 1e-6, 'Generations', 500, 'Vectorized','on','UseParallel',true);
        [xroptim, metric]= ga(@(xt) acquiFun(gprocess, x, y, r, xt, model, design, i), d + 1, [],[], [], [], [lhs_rect(:,1); min(r_cand)], [lhs_rect(:,2); max(r_cand)],[], opts);
        x_optim = xroptim(1:(end - 1));
        r_optim = ceil(xroptim(end));
    otherwise
end

metric = -metric;
end

function EI = acquiFun(gprocess, x, y, r, xt, model, design, i)
% Acquisition function for ABSUR and cSUR

% Description:
%   ACQUIFUN calculates the weighted acquisition function in sequential
%   design
%   EI = ACQUIFUN(GPROCESS, X, Y, R, XT, MODEL, DESIGN, i) returns the
%   weighted expected improvement at ith time stamp in put/call option
%       GPROCESS - fitted Gaussian Process to calculate the EI
%       X - inputs
%       Y - observations
%       R - batch size
%       XT - candidate samples
%       MODEL - GP model(now supports Gaussian likelihood and t likelihood)

d = size(x,2);
if (strcmp(model.method, 'gauss'))
    t_0 = 0.01;
else
    t_0 = 0.005;
end
theta_for_optim = [0.1371 0.000815 1.9871E-6];
overhead = 3*d*CalcOver(theta_for_optim(1), theta_for_optim(2), theta_for_optim(3), size(x, 1));

switch design
    case 'ABSUR'
        [Ef, Varf] = gp_pred(gprocess, x, y, xt(:,1:(end - 1)));
        Ef = Ef(:,1);
        sigmanoise = gprocess.lik.sigma2;
        sigmanoise = sigmanoise./xt(:,end);
        % approximates the new variance %
        Varfnew = sigmanoise .* Varf./(sigmanoise + Varf);
        % calculates the reduced empirical error %
        EI = ( -normcdf( -abs(Ef)./sqrt(abs(Varf)) ) + normcdf( -abs(Ef)./sqrt(abs(Varfnew)) ) )./ (xt(:,end)*t_0 + overhead);
        % density of x %
        x_dens = lognpdf(xt(:,1), log(model.x0(1))+(model.r - model.div - model.sigma(1)^2/2)*i*model.dt, model.sigma(1)*sqrt(i*model.dt));

        if ( d >= 2 )
            for j = 2:d
                x_dens = x_dens.*lognpdf(xt(:,j), log(model.x0(j))+(model.r - model.div - model.sigma(j)^2/2)*i*model.dt, model.sigma(j)*sqrt(i*model.dt));
            end
        end
        
        % excludes deep in-the-money designs %
        if (d == 2  && model.K == 40)
            x_dens(xt(:,1) + xt(:,2) > 80) = 0;
        end
        if (d == 2  && model.K == 100)
            x_dens(xt(:,1) < 100 & xt(:,2) < 100) = 0;
        end
        if (d == 3  && model.K == 100)
            x_dens(xt(:,1) < 100 & xt(:,2) < 100 & xt(:,3) < 100) = 0;
        end
        EI = EI.*x_dens;
    case 'cSUR'
        [Ef, Varf] = gp_pred(gprocess, x, y, xt);
        Ef = Ef(:,1);
        sigmanoise = gprocess.lik.sigma2;
        
        % approximates the new variance %
        Varfnew = sigmanoise .* Varf./(sigmanoise + Varf);
        % calculates the reduced empirical error %
        EI = -normcdf( -abs(Ef)./sqrt(abs(Varf)) ) + normcdf( -abs(Ef)./sqrt(abs(Varfnew)));
        % density of x %
        x_dens = lognpdf(xt(:,1), log(model.x0(1))+(model.r - model.div - model.sigma(1)^2/2)*i*model.dt, model.sigma(1)*sqrt(i*model.dt));

        if ( d >= 2 )
            for j = 2:d
                x_dens = x_dens.*lognpdf(xt(:,j), log(model.x0(j))+(model.r - model.div - model.sigma(j)^2/2)*i*model.dt, model.sigma(j)*sqrt(i*model.dt));
            end
        end
        
        % excludes deep in-the-money designs %
        if (d == 2 && model.K == 40)
            x_dens(xt(:,1) + xt(:,2) > 80) = 0;
        end
        if (d == 2  && model.K == 100)
            x_dens(xt(:,1) < 100 & xt(:,2) < 100) = 0;
        end
        if (d == 3  && model.K == 100)
            x_dens(xt(:,1) < 100 & xt(:,2) < 100 & xt(:,3) < 100) = 0;
        end
        EI = EI.*x_dens;
end
end

function [xoptim, metric] = optimMCU(gprocess, x, y, model, lhs_rect)
% MCU optimization with greedy search %
d = size(x, 2);
if (d == 2)
    xt = lhsCons(400, lhs_rect); % candidate set
else
    xt = lhsCons(1000, lhs_rect); % candidate set
end
[Ef, Varf] = gp_pred(gprocess, x, y, xt);

% density of xt calculated with Geometric Brownian Motion %
x_dens = lognpdf(xt(:,1), log(model.x0(1))+(model.r - model.div - model.sigma(1)^2/2)*i*model.dt, model.sigma(1)*sqrt(i*model.dt));

if ( d >= 2 )
    for j = 2:d
        x_dens = x_dens.*lognpdf(xt(:,j), log(model.x0(j))+(model.r - model.div - model.sigma(j)^2/2)*i*model.dt, model.sigma(j)*sqrt(i*model.dt));
    end
end
% excludes deep in-the-money designs %
if (d == 2  && model.K == 40)
    x_dens(xt(:,1) + xt(:,2) > 80) = 0;
end
if (d == 2  && model.K == 100)
    x_dens(xt(:,1) < 100 & xt(:,2) < 100) = 0;
end
if (d == 3  && model.K == 100)
    x_dens(xt(:,1) < 100 & xt(:,2) < 100 & xt(:,3) < 100) = 0;
end
        
% calculate the local empirical error %
EI = normcdf(-abs(Ef)./sqrt(abs(Varf)));
EI = EI.*x_dens;
% choose the location of sample %
[metric, index] = max(EI);
xoptim = xt(index,:);
end