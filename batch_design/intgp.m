function [x_seq, y_seq, r_seq, ee, er, bias, metric, t, Ef, Varf, xt, xtt, l, sigma2, sigman, steps, t_optim, t_gen, overhead, gamma, post_var] = intgp(k, m0, d, budget, fun, r)

%%%%%%%%%% This function is used to initialize the variable space for
%%%%%%%%%% updatapar.

%%% construct testing points %%%

if ( d == 1 )
    xt = linspace(0,1,m0)';
    xtt = xt;
end

if ( d == 2 )
    x1=repmat(linspace(0,1,m0)',1,m0);
    x2=repmat(linspace(0,1,m0)',1,m0)';
    xt=[x1(:) x2(:)];
    % Approximate empirical error rate %
    rng default  % For reproducibility
    p = sobolset(2,'Skip',1e3,'Leap',1e2);
    p = scramble(p,'MatousekAffineOwen');
    xtt = net(p,100000);
    f = fun(xtt);
    xt1 = xtt(abs(f)<0.12,:);
    xt1 = xt1(1:4000,:);
    xt2 = xtt(abs(f)>=0.12,:);
    xt2 = xt2(1:1000,:);
    xtt = [xt1;xt2];
end

if ( d == 6 )
    % Approximate empirical error rate %
    rng default  % For reproducibility
    p = sobolset(6,'Skip',1e3,'Leap',1e2);
    p = scramble(p,'MatousekAffineOwen');
    xtt = net(p,10000000);
    f = fun(xtt);
    xt1 = xtt(abs(f) < 0.23,:);
    xt1 = xt1(1:4000,:);
    xt2 = xtt(abs(f) >= 0.23,:);
    xt2 = xt2(1:1000,:);
    xtt = [xt1;xt2];
    xt = xtt;
end

%%% initialize variables %%%

m = size(xt,1);
num_diff_samples = min([budget/r(1),5000]); % the variable size is set up to 5000
x_seq = zeros(num_diff_samples,d,k);
y_seq = zeros(num_diff_samples,k);
r_seq = zeros(num_diff_samples, k);

er = zeros(num_diff_samples,k);
ee = zeros(num_diff_samples,k);
bias = zeros(num_diff_samples,k);
metric = zeros(num_diff_samples,k);
t = zeros(num_diff_samples,k);
t_optim = zeros(num_diff_samples,k);
t_gen = zeros(num_diff_samples,k);
overhead = zeros(num_diff_samples,k);

steps = zeros(k,1);

Ef = zeros(m,num_diff_samples,k);
Varf = zeros(m,num_diff_samples,k);

l = zeros(num_diff_samples,d,k);
sigma2 = zeros(num_diff_samples,k);
sigman = zeros(num_diff_samples,k);
gamma = zeros(num_diff_samples,k);
post_var = zeros(num_diff_samples,size(r, 2), k);
end