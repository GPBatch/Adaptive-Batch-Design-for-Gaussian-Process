function sigman2_overall = calculate_overall_sigman2(y, r, sigman2)
% Calculates the overall mean noise variance, especially for heteroskedastic
% cases

% Description:
%   SIGMAN2_OVERALL = CALCULATE_OVERALL_SIGMAN2(Y, R, SIGMAN2) calculates the
%   overall mean noise variance given the observations, batch size and
%   estimated noise variance.
%       Y - current observations
%       R - current batch size
%       SIGMAN2 - estimated noise variance

y2_sum = sum(sigman2.*(r-1) + r.*y.^2);
n = sum(r);
n_ybar2 = sum(r.*y)^2/n;
sigman2_overall = (y2_sum - n_ybar2) / (n-1);
end
