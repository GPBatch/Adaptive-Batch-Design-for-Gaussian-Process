function fit = updateFit(fit, add_grid, add_mean, add_var, r_cur)
fit.x = [fit.x; add_grid];
fit.y = [fit.y; add_mean];
fit.r = [fit.r; r_cur];
% fit.gp.lik.sigma2 = [fit.gp.lik.sigma2; add_var];