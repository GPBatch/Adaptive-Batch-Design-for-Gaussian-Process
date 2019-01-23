function [x_optim, r_optim, gamma, metric] = seq_design(fit, x, y, var, r, model, lhs_rect, i, xtest, xt_dens, gamma)

% optimize the performance metric and select the next design with constrained genetic
% optimization algorithm (200 generations and tolerance value 1e-3)
gprocess = fit.gp;
design = model.design;
r_cand = model.r_cand;
d = size(x,2);
opts = gaoptimset('Display','off','TolFun', 1e-3, 'Generations', 200, 'Vectorized','on','UseParallel',true);
batch = model.batch;

switch batch
    case 'null'
        opts = gaoptimset('Display','off','TolFun', 1e-5, 'Generations', 1000, 'Vectorized','on','UseParallel',true);
        [x_optim, metric]= ga(@(xt) acquiFun(gprocess, x, y, r, xt, model, design, i, xtest, xt_dens), d,[],[], [], [], lhs_rect(:,1), lhs_rect(:,2),[], opts);
        r_optim = -1;
    case 'plain'
        opts = gaoptimset('Display','off','TolFun', 1e-5, 'Generations', 1000, 'Vectorized','on','UseParallel',true);
        [x_optim, metric]= ga(@(xt) acquiFun(gprocess, x, y, r, xt, model, design, i, xtest, xt_dens), d,[],[], [], [], lhs_rect(:,1), lhs_rect(:,2),[], opts);
        r_optim = model.km_batch;
    case 'MFA'
        opts = gaoptimset('Display','off','TolFun', 1e-5, 'Generations', 1000, 'Vectorized','on','UseParallel',true);
        [x_optim, metric]= ga(@(xt) acquiFun(gprocess, x, y, r, xt, model, design, i, xtest, xt_dens), d,[],[], [], [], lhs_rect(:,1), lhs_rect(:,2),[], opts);
%         sigman2_est = calculate_overall_sigman2(y, r, var);
        sigman2_est = gprocess.lik.sigma2;
        sigman2_est_by_level = sigman2_est./(r_cand - 1);
        [Ef, Varf] = gp_pred(gprocess, x, y, x_optim);
        post_var = Varf.*sigman2_est_by_level./(sigman2_est_by_level + Varf);
        sigman_est_optim = sqrt(post_var);
        r_index = sum(sigman_est_optim > gamma);
        while (r_index == 0 || r_cand(r_index) < r(end))
            gamma = gamma*0.95;
            r_index = sum(sigman_est_optim > gamma);
        end
        r_optim = r_cand(r_index);
    case 'MFA_NoConstraint'
        opts = gaoptimset('Display','off','TolFun', 1e-5, 'Generations', 1000, 'Vectorized','on','UseParallel',true);
        [x_optim, metric]= ga(@(xt) acquiFun(gprocess, x, y, r, xt, model, design, i, xtest, xt_dens), d,[],[], [], [], lhs_rect(:,1), lhs_rect(:,2),[], opts);
%         sigman2_est = calculate_overall_sigman2(y, r, gprocess.lik.sigma2);
        sigman2_est = gprocess.lik.sigma2;
        sigman2_est_by_level = sigman2_est./(r_cand - 1);
        [Ef, Varf] = gp_pred(gprocess, x, y, x_optim);
        post_var = Varf.*sigman2_est_by_level./(sigman2_est_by_level + Varf);
        sigman_est_optim = sqrt(post_var);
        r_index = sum(sigman_est_optim > gamma);
        while (r_index == 0)
            gamma = gamma*0.6;
            r_index = sum(sigman_est_optim > gamma);
        end
        r_optim = r_cand(r_index);
    case 'aSUR'
        opts = gaoptimset('Display','off','TolFun', 1e-6, 'Generations', 500, 'Vectorized','on','UseParallel',true);
        [xroptim, metric]= ga(@(xt) acquiFun(gprocess, x, y, r, xt, model, design, i, xtest, xt_dens), d + 1, [],[], [], [], [lhs_rect(:,1); min(r_cand)], [lhs_rect(:,2); max(r_cand)],[], opts);
        x_optim = xroptim(1:(end - 1));
        r_optim = ceil(xroptim(end));
    otherwise
end

metric = -metric;
end

% function [c, ceq] = simple_constraint(x)
% c = x(:,1) + x(:,2) - 80;
% ceq = [];
% end

function EI = acquiFun(gprocess, x, y, r, xt, model, design, i, xtest, xt_dens)
d = size(x,2);
if (strcmp(model, 'gauss'))
%     t_0 = 0.01;
    t_0 = 1;
else
    t_0 = 0.005;
end
theta_for_optim = [0.1371 0.000815 1.9871E-6];
overhead = 3*d*CalcOver(theta_for_optim(1), theta_for_optim(2), theta_for_optim(3), size(x, 1));

switch design
    case 'MCU'
        [Ef, Varf] = gp_pred(gprocess, x, y, xt);
        Ef = Ef(:,1);
        % density of x %
        x_dens = lognpdf(xt(:,1), log(model.x0(1))+(model.r - model.div - model.sigma(1)^2/2)*i*model.dt, model.sigma(1)*sqrt(i*model.dt));

        if ( d >= 2 )
            for j = 2:d
                x_dens = x_dens.*lognpdf(xt(:,j), log(model.x0(j))+(model.r - model.div - model.sigma(j)^2/2)*i*model.dt, model.sigma(j)*sqrt(i*model.dt));
            end
        end
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
        EI = -normcdf(-abs(Ef)./sqrt(abs(Varf)));
        EI = EI.*x_dens;
    case 'tMSE'
        [Ef, Varf] = gp_pred(gprocess, x, y, xt);
%         [Ef, Varf] = gp_pred(gprocess, standardize(x), y, standardize(xt));
        Ef = Ef(:,1);
        sigmaepsilon = 0.2;
        w = exp(-0.5.*(Ef./sqrt(Varf + sigmaepsilon^2)).^2)./sqrt(2*pi*(Varf + sigmaepsilon^2));
        % eq. 3.5 %
        EI = -Varf.*w;
        % density of x %
        x_dens = lognpdf(xt(:,1), log(model.x0(1))+(model.r - model.div - model.sigma(1)^2/2)*i*model.dt, model.sigma(1)*sqrt(i*model.dt));

        if ( d >= 2 )
            for j = 2:d
                x_dens = x_dens.*lognpdf(xt(:,j), log(model.x0(j))+(model.r - model.div - model.sigma(j)^2/2)*i*model.dt, model.sigma(j)*sqrt(i*model.dt));
            end
        end
        
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
    case 'aSUR'
        [Ef, Varf] = gp_pred(gprocess, x, y, xt(:,1:(end - 1)));
%         [Ef, Varf] = gp_pred(gprocess, standardize(x), y, standardize(xt));
        Ef = Ef(:,1);
%         switch model.method
%             case {'gauss', 'mgauss'}
%                 sigmanoise = gprocess.lik.sigma2;
%             case 't'
%                 sigmanoise = gprocess.lik.sigma2;
%             case {'probit', 'mprobit'}
%                 sigmanoise = inverseterm_probit(Ef, Varf);
%         end
%         sigmanoise = calculate_overall_sigman2(y, r, gprocess.lik.sigma2);
        sigmanoise = gprocess.lik.sigma2;
        sigmanoise = sigmanoise./xt(:,end);
        % approximate the new variance in eq. 4.3 %
        Varfnew = sigmanoise .* Varf./(sigmanoise + Varf);
        % calculate the reduced empirical error in eq. 3.2 %
        EI = ( -normcdf( -abs(Ef)./sqrt(abs(Varf)) ) + normcdf( -abs(Ef)./sqrt(abs(Varfnew)) ) )./ (xt(:,end)*t_0 + overhead);
        % density of x %
        x_dens = lognpdf(xt(:,1), log(model.x0(1))+(model.r - model.div - model.sigma(1)^2/2)*i*model.dt, model.sigma(1)*sqrt(i*model.dt));

        if ( d >= 2 )
            for j = 2:d
                x_dens = x_dens.*lognpdf(xt(:,j), log(model.x0(j))+(model.r - model.div - model.sigma(j)^2/2)*i*model.dt, model.sigma(j)*sqrt(i*model.dt));
            end
        end
        
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
%         EI( option_payoff( xt,model.K) <= 0) = 0;
    case 'cSUR'
        [Ef, Varf] = gp_pred(gprocess, x, y, xt);
%         [Ef, Varf] = gp_pred(gprocess, standardize(x), y, standardize(xt));
        Ef = Ef(:,1);
%         switch model.method
%             case {'gauss', 'mgauss'}
%                 sigmanoise = gprocess.lik.sigma2;
%             case 't'
%                 sigmanoise = gprocess.lik.sigma2;
%             case {'probit', 'mprobit'}
%                 sigmanoise = inverseterm_probit(Ef, Varf);
%         end
%         sigmanoise = calculate_overall_sigman2(y, r, gprocess.lik.sigma2);
        sigmanoise = gprocess.lik.sigma2;
        
        % approximate the new variance in eq. 4.3 %
        Varfnew = sigmanoise .* Varf./(sigmanoise + Varf);
        % calculate the reduced empirical error in eq. 3.2 %
        EI = -normcdf( -abs(Ef)./sqrt(abs(Varf)) ) + normcdf( -abs(Ef)./sqrt(abs(Varfnew)));
        % density of x %
        x_dens = lognpdf(xt(:,1), log(model.x0(1))+(model.r - model.div - model.sigma(1)^2/2)*i*model.dt, model.sigma(1)*sqrt(i*model.dt));

        if ( d >= 2 )
            for j = 2:d
                x_dens = x_dens.*lognpdf(xt(:,j), log(model.x0(j))+(model.r - model.div - model.sigma(j)^2/2)*i*model.dt, model.sigma(j)*sqrt(i*model.dt));
            end
        end
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
%         EI( option_payoff( xt,model.K) <= 0) = 0;
    otherwise
        n = size(xtest,1);
%         n = 100^d;
            
%         xtest = lhsCons(n, lhs_rect);
%         xt_dens = lognpdf(xtest(:,1), log(model.x0(1))+(model.r - model.div - model.sigma(1)^2/2)*i*model.dt, model.sigma(1)*sqrt(i*model.dt));
% 
%         if ( d >= 2 )
%             for j = 2:d
%                 xt_dens = xt_dens.*lognpdf(xtest(:,j), log(model.x0(j))+(model.r - model.div - model.sigma(j)^2/2)*i*model.dt, model.sigma(j)*sqrt(i*model.dt));
%             end
%         end

%         [Efall, Covfall] = gp_jpred(gprocess, standardize(x), y, standardize([xtest;xt]));
        [Efall, Covfall] = gp_jpred(gprocess, x, y, [xtest;xt]);
        Varfall = diag(Covfall);
        m = size(xt,1);
        Varftest = repmat(Varfall(1:n),1,m);
        Varf = Varfall((n+1):end);
        Eftest = repmat(Efall(1:n,1),1,m);
        Ef = Efall((n+1):end,1);
        switch model.method
            case {'gauss', 'mgauss', 't'}
                sigmanoise = gprocess.lik.sigma2/model.km_batch;
            case {'probit', 'mprobit'}
                sigmanoise = inverseterm_probit(Ef, Varf);
        end
        % approximate the new variance in eq. 4.5 %
        Varfupdate = Varftest - Covfall(1:n,(n+1):end).^2 ./repmat((sigmanoise + Varf)',n,1);
        ilocal = normcdf(-abs(Eftest)./sqrt(abs(Varfupdate)));
        ilocal = ilocal.*repmat(xt_dens,1,m);
        % calculate the expected empirical error in eq. 3.4 %
        EI = sum(ilocal)';
%         EI( option_payoff( xt,model.K) <= 0) = 0;
end
end
