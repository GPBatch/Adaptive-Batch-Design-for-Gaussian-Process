function [fit, l, sigma2, sigman, nu] = gp_setup_amput(model, x, y, r, opt)

% This function is used to construct GP simulator for American Put/Call Option
% function.

% Params:
% @model: the likelihood function used in GP: GP or t-GP
% @x: current designs
% @y: current observations
% @r: current batch

d = size(x,2);

x0 = x;

if (nargin == 4)
    opt=optimset('TolFun',1e-3,'TolX',1e-3);
end

switch model
            case 'gauss'
                gpcf = gpcf_sexp('lengthScale',repmat(0.5,1,d), 'magnSigma2', 1);
                gprocess = gp_set('lik', lik_gaussian(), 'cf', gpcf, 'jitterSigma2', 1e-9);
                gprocess = gp_optim_amput(gprocess, x, y, 'opt', opt);
                w = gp_pak(gprocess);
                sigman = exp(w(end))^2;
                sigma2 = exp(w(1))^2;
                l = exp(w(2:(end-1)));
                nu = 0;
            case 't'
                gpcf = gpcf_sexp('lengthScale',repmat(0.5,1,d), 'magnSigma2', 1);
                gprocess = gp_set('lik', lik_t(), 'cf', gpcf, 'jitterSigma2', 1e-9);
                gprocess = gp_optim_amput(gprocess, x, y, 'opt', opt);
                w = gp_pak(gprocess);
                sigman = exp(w(end-1))^2;
                sigma2 = exp(w(1))^2;
                l = exp(w(2:(end-2)));
                nu = exp(w(end));
end


fit.gp = gprocess;
fit.x = x0;
fit.y = y;
fit.r = r;
fit.t = 0;
end