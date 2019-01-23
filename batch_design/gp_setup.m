function [gprocess, l, sigma2] = gp_setup(model, x, y, sigman, opt)

%%%%%%%%%% This function is used to set up GP for further use %%%%%%%%%%%
d = size(x,2);

if (nargin == 4)
    opt=optimset('TolFun',1e-3,'TolX',1e-3);
end

gpcf = gpcf_sexp('lengthScale',repmat(0.5,1,d), 'magnSigma2', 1);
switch model
    case 'gauss'
        gprocess = gp_set('lik', lik_gaussian('sigma2', sigman, 'sigma2_prior', []), 'cf', gpcf, 'jitterSigma2', 1e-9);
%         gprocess = gp_set('lik', lik_gaussian(), 'cf', gpcf, 'jitterSigma2', 1e-9);
    case 't'
        gprocess = gp_set('lik', lik_t('sigma2', sigman, 'sigma2_prior', []), 'cf', gpcf, 'jitterSigma2', 1e-9);
%         gprocess = gp_set('lik', lik_t(), 'cf', gpcf, 'jitterSigma2', 1e-9);
end
gprocess = gp_optim(gprocess, x, y, 'opt', opt);
w = gp_pak(gprocess);
sigma2 = exp(w(1))^2;
l = exp(w(2:(2+d-1)));
  
end