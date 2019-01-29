clear all; close all; clc;
startup;

k = 20;
I = 20;
d = 2;

fun = @braninsc2; 

noisestructure = 't_constdf';
noisevar = 'large';

candidatesize = 1000;

budget = 150;

x1=repmat(linspace(0,1,50)',1,50);
x2=repmat(linspace(0,1,50)',1,50)';
xt=[x1(:) x2(:)];
m = size(xt,1);

rng default  % For reproducibility
p = sobolset(2,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
xtt = net(p,100000);
f = fun(xtt);
xt1 = xtt(abs(f)<0.5,:);
xt1 = xt1(1:400,:);
xt2 = xtt(abs(f)>=0.5,:);
xt2 = xt2(1:100,:);
xtt = [xt1;xt2];

x_seq_t_mee = zeros(budget,d,k);
x_seq_gauss_mee = zeros(budget,d,k);
x_seq_probit_mee = zeros(budget,d,k);
x_seq_t_eee = zeros(budget,d,k);
x_seq_gauss_eee = zeros(budget,d,k);
x_seq_probit_eee = zeros(budget,d,k);
x_seq_t_dee = zeros(budget,d,k);
x_seq_gauss_dee = zeros(budget,d,k);
x_seq_probit_dee = zeros(budget,d,k);
y_seq_t_mee = zeros(budget,k);
y_seq_gauss_mee = zeros(budget,k);
y_seq_probit_mee = zeros(budget,k);
y_seq_t_eee = zeros(budget,k);
y_seq_gauss_eee = zeros(budget,k);
y_seq_probit_eee = zeros(budget,k);
y_seq_t_dee = zeros(budget,k);
y_seq_gauss_dee = zeros(budget,k);
y_seq_probit_dee = zeros(budget,k);
x_random_t = zeros(budget,d,k);
x_random_gauss = zeros(budget,d,k);
x_random_probit = zeros(budget,d,k);
y_random_t = zeros(budget,k);
y_random_gauss = zeros(budget,k);
y_random_probit = zeros(budget,k);

x_seq_mgauss_mee = zeros(budget,d,k);
x_seq_mprobit_mee = zeros(budget,d,k);

x_seq_mgauss_eee = zeros(budget,d,k);
x_seq_mprobit_eee = zeros(budget,d,k);

x_seq_mgauss_dee = zeros(budget,d,k);
x_seq_mprobit_dee = zeros(budget,d,k);

y_seq_mgauss_mee = zeros(budget,k);
y_seq_mprobit_mee = zeros(budget,k);

y_seq_mgauss_eee = zeros(budget,k);
y_seq_mprobit_eee = zeros(budget,k);

y_seq_mgauss_dee = zeros(budget,k);
y_seq_mprobit_dee = zeros(budget,k);

x_random_mgauss = zeros(budget,d,k);
x_random_mprobit = zeros(budget,d,k);

y_random_mgauss = zeros(budget,k);
y_random_mprobit = zeros(budget,k);

ermgauss_mee = zeros(budget-I,k);
ermprobit_mee = zeros(budget-I,k);

vmgauss_mee = zeros(budget-I,k);
vmprobit_mee = zeros(budget-I,k);

biasmgauss_mee = zeros(budget-I,k);
biasmprobit_mee = zeros(budget-I,k);

ermgauss_eee = zeros(budget-I,k);
ermprobit_eee = zeros(budget-I,k);

vmgauss_eee = zeros(budget-I,k);
vmprobit_eee = zeros(budget-I,k);

biasmgauss_eee = zeros(budget-I,k);
biasmprobit_eee = zeros(budget-I,k);

ermgauss_dee = zeros(budget-I,k);
ermprobit_dee = zeros(budget-I,k);

vmgauss_dee = zeros(budget-I,k);
vmprobit_dee = zeros(budget-I,k);

biasmgauss_dee = zeros(budget-I,k);
biasmprobit_dee = zeros(budget-I,k);

ermgaussrandom = zeros(budget-I,k);
ermprobitrandom = zeros(budget-I,k);

vmgaussrandom = zeros(budget-I,k);
vmprobitrandom = zeros(budget-I,k);

biasmgaussrandom = zeros(budget-I,k);
biasmprobitrandom = zeros(budget-I,k);

eemgauss_mee = zeros(budget-I,k);
eemprobit_mee = zeros(budget-I,k);

eemgauss_dee = zeros(budget-I,k);
eemprobit_dee = zeros(budget-I,k);

eemgauss_eee = zeros(budget-I,k);
eemprobit_eee = zeros(budget-I,k);

eemgaussrandom = zeros(budget-I,k);
eemprobitrandom = zeros(budget-I,k);


meemgauss = zeros(budget-I,k);
meemprobit = zeros(budget-I,k);

eeemgauss = zeros(budget-I,k);
eeemprobit = zeros(budget-I,k);

deemgauss = zeros(budget-I,k);
deemprobit = zeros(budget-I,k);

ert_mee = zeros(budget-I,k);
ergauss_mee = zeros(budget-I,k);
erprobit_mee = zeros(budget-I,k);
vt_mee = zeros(budget-I,k);
vgauss_mee = zeros(budget-I,k);
vprobit_mee = zeros(budget-I,k);
biast_mee = zeros(budget-I,k);
biasgauss_mee = zeros(budget-I,k);
biasprobit_mee = zeros(budget-I,k);
ert_eee = zeros(budget-I,k);
ergauss_eee = zeros(budget-I,k);
erprobit_eee = zeros(budget-I,k);
vt_eee = zeros(budget-I,k);
vgauss_eee = zeros(budget-I,k);
vprobit_eee = zeros(budget-I,k);
biast_eee = zeros(budget-I,k);
biasgauss_eee = zeros(budget-I,k);
biasprobit_eee = zeros(budget-I,k);
ert_dee = zeros(budget-I,k);
ergauss_dee = zeros(budget-I,k);
erprobit_dee = zeros(budget-I,k);
vt_dee = zeros(budget-I,k);
vgauss_dee = zeros(budget-I,k);
vprobit_dee = zeros(budget-I,k);
biast_dee = zeros(budget-I,k);
biasgauss_dee = zeros(budget-I,k);
biasprobit_dee = zeros(budget-I,k);
ertrandom = zeros(budget-I,k);
ergaussrandom = zeros(budget-I,k);
erprobitrandom = zeros(budget-I,k);
vtrandom = zeros(budget-I,k);
vgaussrandom = zeros(budget-I,k);
vprobitrandom = zeros(budget-I,k);
biastrandom = zeros(budget-I,k);
biasgaussrandom = zeros(budget-I,k);
biasprobitrandom = zeros(budget-I,k);
eet_mee = zeros(budget-I,k);
eegauss_mee = zeros(budget-I,k);
eeprobit_mee = zeros(budget-I,k);
eet_dee = zeros(budget-I,k);
eegauss_dee = zeros(budget-I,k);
eeprobit_dee = zeros(budget-I,k);
eet_eee = zeros(budget-I,k);
eegauss_eee = zeros(budget-I,k);
eeprobit_eee = zeros(budget-I,k);
eetrandom = zeros(budget-I,k);
eegaussrandom = zeros(budget-I,k);
eeprobitrandom = zeros(budget-I,k);

meet = zeros(budget-I,k);
meegauss = zeros(budget-I,k);
meeprobit = zeros(budget-I,k);
eeet = zeros(budget-I,k);
eeegauss = zeros(budget-I,k);
eeeprobit = zeros(budget-I,k);
deet = zeros(budget-I,k);
deegauss = zeros(budget-I,k);
deeprobit = zeros(budget-I,k);



tmeeprobit = zeros(k,1);
teeeprobit = zeros(k,1);
tdeeprobit = zeros(k,1);
trandomprobit = zeros(k,1);
tmeegauss = zeros(k,1);
teeegauss = zeros(k,1);
tdeegauss = zeros(k,1);
trandomgauss = zeros(k,1);
tmeet = zeros(k,1);
teeet = zeros(k,1);
tdeet = zeros(k,1);
trandomt = zeros(k,1);
tmeemprobit = zeros(k,1);
teeemprobit = zeros(k,1);
tdeemprobit = zeros(k,1);
trandommprobit = zeros(k,1);
tmeemgauss = zeros(k,1);
teeemgauss = zeros(k,1);
tdeemgauss = zeros(k,1);
trandommgauss = zeros(k,1);

ttmsegauss = zeros(k,1);

x_seq_gauss_tmse = zeros(budget, d,k);
y_seq_gauss_tmse = zeros(budget, k);
er_gauss_tmse = zeros(budget-I,k);
ee_gauss_tmse = zeros(budget-I,k);
bias_gauss_tmse = zeros(budget-I,k);
tmsegauss = zeros(budget-I,k);

ttmset = zeros(k,1);
x_seq_t_tmse = zeros(budget, d,k);
y_seq_t_tmse = zeros(budget, k);
er_t_tmse = zeros(budget-I,k);
ee_t_tmse = zeros(budget-I,k);
bias_t_tmse = zeros(budget-I,k);
tmset = zeros(budget-I,k);

ttmsegauss = zeros(k,1);

x_seq_mgauss_tmse = zeros(budget,d,k);
y_seq_mgauss_tmse = zeros(budget, k);
er_mgauss_tmse = zeros(budget-I,k);
ee_mgauss_tmse = zeros(budget-I,k);
bias_mgauss_tmse = zeros(budget-I,k);
tmsemgauss = zeros(budget-I,k);

ttmseprobit = zeros(k,1);
x_seq_probit_tmse = zeros(budget,d,k);
y_seq_probit_tmse = zeros(budget, k);
er_probit_tmse = zeros(budget-I,k);
ee_probit_tmse = zeros(budget-I,k);
bias_probit_tmse = zeros(budget-I,k);
tmseprobit = zeros(budget-I,k);

ttmsemprobit = zeros(k,1);
x_seq_mprobit_tmse = zeros(budget,d,k);
y_seq_mprobit_tmse = zeros(budget, k);
er_mprobit_tmse = zeros(budget-I,k);
ee_mprobit_tmse = zeros(budget-I,k);
bias_mprobit_tmse = zeros(budget-I,k);
tmsemprobit = zeros(budget-I,k);

erlastgauss_mee = zeros(k,1);
erlastgauss_tmse = zeros(k,1);
erlastgauss_dee = zeros(k,1);
erlastgauss_eee = zeros(k,1);

erlastmgauss_mee = zeros(k,1);
erlastmgauss_tmse = zeros(k,1);
erlastmgauss_dee = zeros(k,1);
erlastmgauss_eee = zeros(k,1);

erlastprobit_mee = zeros(k,1);
erlastprobit_tmse = zeros(k,1);
erlastprobit_dee = zeros(k,1);
erlastprobit_eee = zeros(k,1);

erlastmprobit_mee = zeros(k,1);
erlastmprobit_tmse = zeros(k,1);
erlastmprobit_dee = zeros(k,1);
erlastmprobit_eee = zeros(k,1);

erlastt_mee = zeros(k,1);
erlastt_tmse = zeros(k,1);
erlastt_dee = zeros(k,1);
erlastt_eee = zeros(k,1);

% Ef_gauss_mee = zeros(candidatesize+m,k);
% Varf_gauss_mee = zeros(candidatesize+m,k);
% 
% Ef_mgauss_mee = zeros(candidatesize+m,k);
% Varf_mgauss_mee = zeros(candidatesize+m,k);
% 
% Ef_probit_mee = zeros(candidatesize+m,k);
% Varf_probit_mee = zeros(candidatesize+m,k);
% 
% Ef_mprobit_mee = zeros(candidatesize+m,k);
% Varf_mprobit_mee = zeros(candidatesize+m,k);
% 
% Ef_t_mee = zeros(candidatesize+m,k);
% Varf_t_mee = zeros(candidatesize+m,k);
% 
% Ef_gauss_tmse = zeros(candidatesize+m,k);
% Varf_gauss_tmse = zeros(candidatesize+m,k);
% 
% Ef_mgauss_tmse = zeros(candidatesize+m,k);
% Varf_mgauss_tmse = zeros(candidatesize+m,k);
% 
% Ef_probit_tmse = zeros(candidatesize+m,k);
% Varf_probit_tmse = zeros(candidatesize+m,k);
% 
% Ef_mprobit_tmse = zeros(candidatesize+m,k);
% Varf_mprobit_tmse = zeros(candidatesize+m,k);
% 
% Ef_t_tmse = zeros(candidatesize+m,k);
% Varf_t_tmse = zeros(candidatesize+m,k);
% 
% Ef_gauss_dee = zeros(candidatesize+m,k);
% Varf_gauss_dee = zeros(candidatesize+m,k);
% 
% Ef_mgauss_dee = zeros(candidatesize+m,k);
% Varf_mgauss_dee = zeros(candidatesize+m,k);
% 
% Ef_probit_dee = zeros(candidatesize+m,k);
% Varf_probit_dee = zeros(candidatesize+m,k);
% 
% Ef_mprobit_dee = zeros(candidatesize+m,k);
% Varf_mprobit_dee = zeros(candidatesize+m,k);
% 
% Ef_t_dee = zeros(candidatesize+m,k);
% Varf_t_dee = zeros(candidatesize+m,k);
% 
% Ef_gauss_eee = zeros(candidatesize+m,k);
% Varf_gauss_eee = zeros(candidatesize+m,k);
% 
% Ef_mgauss_eee = zeros(candidatesize+m,k);
% Varf_mgauss_eee = zeros(candidatesize+m,k);
% 
% Ef_probit_eee = zeros(candidatesize+m,k);
% Varf_probit_eee = zeros(candidatesize+m,k);
% 
% Ef_mprobit_eee = zeros(candidatesize+m,k);
% Varf_mprobit_eee = zeros(candidatesize+m,k);
% 
% Ef_t_eee = zeros(candidatesize+m,k);
% Varf_t_eee = zeros(candidatesize+m,k);
% 
% Ef_gauss_random = zeros(candidatesize+m,k);
% Varf_gauss_random = zeros(candidatesize+m,k);
% 
% Ef_mgauss_random = zeros(candidatesize+m,k);
% Varf_mgauss_random = zeros(candidatesize+m,k);
% 
% Ef_probit_random = zeros(candidatesize+m,k);
% Varf_probit_random = zeros(candidatesize+m,k);
% 
% Ef_mprobit_random = zeros(candidatesize+m,k);
% Varf_mprobit_random = zeros(candidatesize+m,k);
% 
% Ef_t_random = zeros(candidatesize+m,k);
% Varf_t_random = zeros(candidatesize+m,k);

Ef_gauss_mee = zeros(m,budget-I,k);
Varf_gauss_mee = zeros(m,budget-I,k);

Ef_mgauss_mee = zeros(m,budget-I,k);
Varf_mgauss_mee = zeros(m,budget-I,k);

Ef_probit_mee = zeros(m,budget-I,k);
Varf_probit_mee = zeros(m,budget-I,k);

Ef_mprobit_mee = zeros(m,budget-I,k);
Varf_mprobit_mee = zeros(m,budget-I,k);

Ef_t_mee = zeros(m,budget-I,k);
Varf_t_mee = zeros(m,budget-I,k);

Ef_gauss_tmse = zeros(m,budget-I,k);
Varf_gauss_tmse = zeros(m,budget-I,k);

Ef_mgauss_tmse = zeros(m,budget-I,k);
Varf_mgauss_tmse = zeros(m,budget-I,k);

Ef_probit_tmse = zeros(m,budget-I,k);
Varf_probit_tmse = zeros(m,budget-I,k);

Ef_mprobit_tmse = zeros(m,budget-I,k);
Varf_mprobit_tmse = zeros(m,budget-I,k);

Ef_t_tmse = zeros(m,budget-I,k);
Varf_t_tmse = zeros(m,budget-I,k);

Ef_gauss_dee = zeros(m,budget-I,k);
Varf_gauss_dee = zeros(m,budget-I,k);

Ef_mgauss_dee = zeros(m,budget-I,k);
Varf_mgauss_dee = zeros(m,budget-I,k);

Ef_probit_dee = zeros(m,budget-I,k);
Varf_probit_dee = zeros(m,budget-I,k);

Ef_mprobit_dee = zeros(m,budget-I,k);
Varf_mprobit_dee = zeros(m,budget-I,k);

Ef_t_dee = zeros(m,budget-I,k);
Varf_t_dee = zeros(m,budget-I,k);

Ef_gauss_eee = zeros(m,budget-I,k);
Varf_gauss_eee = zeros(m,budget-I,k);

Ef_mgauss_eee = zeros(m,budget-I,k);
Varf_mgauss_eee = zeros(m,budget-I,k);

Ef_probit_eee = zeros(m,budget-I,k);
Varf_probit_eee = zeros(m,budget-I,k);

Ef_mprobit_eee = zeros(m,budget-I,k);
Varf_mprobit_eee = zeros(m,budget-I,k);

Ef_t_eee = zeros(m,budget-I,k);
Varf_t_eee = zeros(m,budget-I,k);

Ef_gauss_random = zeros(m,budget-I,k);
Varf_gauss_random = zeros(m,budget-I,k);

Ef_mgauss_random = zeros(m,budget-I,k);
Varf_mgauss_random = zeros(m,budget-I,k);

Ef_probit_random = zeros(m,budget-I,k);
Varf_probit_random = zeros(m,budget-I,k);

Ef_mprobit_random = zeros(m,budget-I,k);
Varf_mprobit_random = zeros(m,budget-I,k);

Ef_t_random = zeros(m,budget-I,k);
Varf_t_random = zeros(m,budget-I,k);

nut_mee = zeros(budget-I,k);
nut_eee = zeros(budget-I,k);
nut_dee = zeros(budget-I,k);
nut_tmse = zeros(budget-I,k);

lt_mee = zeros(budget-I,d,k);
sigma2t_mee = zeros(budget-I,k);
sigmant_mee = zeros(budget-I,k);

lt_eee = zeros(budget-I,d,k);
sigma2t_eee = zeros(budget-I,k);
sigmant_eee = zeros(budget-I,k);

lt_dee = zeros(budget-I,d,k);
sigma2t_dee = zeros(budget-I,k);
sigmant_dee = zeros(budget-I,k);

lt_tmse = zeros(budget-I,d,k);
sigma2t_tmse = zeros(budget-I,k);
sigmant_tmse = zeros(budget-I,k);


lgauss_mee = zeros(budget-I,d,k);
sigma2gauss_mee = zeros(budget-I,k);
sigmangauss_mee = zeros(budget-I,k);

lgauss_eee = zeros(budget-I,d,k);
sigma2gauss_eee = zeros(budget-I,k);
sigmangauss_eee = zeros(budget-I,k);

lgauss_dee = zeros(budget-I,d,k);
sigma2gauss_dee = zeros(budget-I,k);
sigmangauss_dee = zeros(budget-I,k);

lgauss_tmse = zeros(budget-I,d,k);
sigma2gauss_tmse = zeros(budget-I,k);
sigmangauss_tmse = zeros(budget-I,k);

lmgauss_mee = zeros(budget-I,d,k);
sigma2mgauss_mee = zeros(budget-I,k);
sigmanmgauss_mee = zeros(budget-I,k);

lmgauss_eee = zeros(budget-I,d,k);
sigma2mgauss_eee = zeros(budget-I,k);
sigmanmgauss_eee = zeros(budget-I,k);

lmgauss_dee = zeros(budget-I,d,k);
sigma2mgauss_dee = zeros(budget-I,k);
sigmanmgauss_dee = zeros(budget-I,k);

lmgauss_tmse = zeros(budget-I,d,k);
sigma2mgauss_tmse = zeros(budget-I,k);
sigmanmgauss_tmse = zeros(budget-I,k);

lprobit_mee = zeros(budget-I,d,k);
sigma2probit_mee = zeros(budget-I,k);
% sigmanprobit_mee = zeros(budget-I,k);

lprobit_eee = zeros(budget-I,d,k);
sigma2probit_eee = zeros(budget-I,k);
% sigmanprobit_eee = zeros(budget-I,k);

lprobit_dee = zeros(budget-I,d,k);
sigma2probit_dee = zeros(budget-I,k);
% sigmanprobit_dee = zeros(budget-I,k);

lprobit_tmse = zeros(budget-I,d,k);
sigma2probit_tmse = zeros(budget-I,k);
% sigmanprobit_tmse = zeros(budget-I,k);

lmprobit_mee = zeros(budget-I,d,k);
sigma2mprobit_mee = zeros(budget-I,k);
% sigmanmprobit_mee = zeros(budget-I,k);

lmprobit_eee = zeros(budget-I,d,k);
sigma2mprobit_eee = zeros(budget-I,k);
% sigmanmprobit_eee = zeros(budget-I,k);

lmprobit_dee = zeros(budget-I,d,k);
sigma2mprobit_dee = zeros(budget-I,k);
% sigmanmprobit_dee = zeros(budget-I,k);

lmprobit_tmse = zeros(budget-I,d,k);
sigma2mprobit_tmse = zeros(budget-I,k);
% sigmanmprobit_tmse = zeros(budget-I,k);



% gauss_mee %

parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);
    
    [x_seq, y_seq, metric, tmeegauss(j), ee, er, bias, Ef, Varf, l, sigma2, sigman] = updategp6(fun, noisestructure, noisevar, Xint, xt, 'gauss', 'MEE', budget, xtt);
    
    lgauss_mee(:,:,j) = l;
    sigma2gauss_mee(:,j) = sigma2;
    sigmangauss_mee(:,j) = sigman;
    
%     [Ef, Varf] = gp_pred(gprocess, x_seq, y_seq, xt);
    Ef_gauss_mee(:,:,j) = Ef;
    Varf_gauss_mee(:,:,j) = Varf;
    x_seq_gauss_mee(:,:,j) = x_seq;
    y_seq_gauss_mee(:,j) = y_seq;
    meegauss(:,j) = metric;
    eegauss_mee(:,j) = ee;
    ergauss_mee(:,j) = er;
    biasgauss_mee(:,j) = bias;
%     erlastgauss_mee(j) = erlast;
end

disp('gaussmee finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/gaussmee', 'x_seq_gauss_mee', 'y_seq_gauss_mee', 'meegauss', 'eegauss_mee', 'ergauss_mee', 'biasgauss_mee', 'tmeegauss', 'Ef_gauss_mee','Varf_gauss_mee','lgauss_mee','sigma2gauss_mee','sigmangauss_mee')

% gauss_tmse %

parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);
    
    [x_seq, y_seq, metric, ttmsegauss(j), ee, er, bias, Ef, Varf, l, sigma2, sigman] = updategp6(fun, noisestructure, noisevar, Xint, xt, 'gauss', 'tMSE', budget,xtt);
    
    lgauss_tmse(:,:,j) = l;
    sigma2gauss_tmse(:,j) = sigma2;
    sigmangauss_tmse(:,j) = sigman;
    
    Ef_gauss_tmse(:,:,j) = Ef;
    Varf_gauss_tmse(:,:,j) = Varf;
    x_seq_gauss_tmse(:,:,j) = x_seq;
    y_seq_gauss_tmse(:,j) = y_seq;
    tmsegauss(:,j) = metric;
    eegauss_tmse(:,j) = ee;
    ergauss_tmse(:,j) = er;
    biasgauss_tmse(:,j) = bias;
%     erlastgauss_tmse(j) = erlast;
end

disp('gausstmse finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/gausstmse', 'x_seq_gauss_tmse', 'y_seq_gauss_tmse', 'tmsegauss','eegauss_tmse', 'ergauss_tmse', 'biasgauss_tmse','ttmsegauss','Ef_gauss_tmse', 'Varf_gauss_tmse','lgauss_tmse','sigma2gauss_tmse','sigmangauss_tmse')

% gauss_meesur %

parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);

    [x_seq, y_seq, metric, tdeegauss(j), ee, er, bias, Ef, Varf, l, sigma2, sigman] = updategp6(fun, noisestructure, noisevar, Xint, xt, 'gauss', 'MEESUR', budget,xtt);
    
    lgauss_dee(:,:,j) = l;
    sigma2gauss_dee(:,j) = sigma2;
    sigmangauss_dee(:,j) = sigman;
    
    Ef_gauss_dee(:,:,j) = Ef;
    Varf_gauss_dee(:,:,j) = Varf;
    x_seq_gauss_dee(:,:,j) = x_seq;
    y_seq_gauss_dee(:,j) = y_seq;
    deegauss(:,j) = metric;
    eegauss_dee(:,j) = ee;
    ergauss_dee(:,j) = er;
    biasgauss_dee(:,j) = bias;
%     erlastgauss_dee(j) = erlast;
end


disp('gaussmeesur finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/gaussmeesur', 'x_seq_gauss_dee', 'y_seq_gauss_dee', 'deegauss', 'eegauss_dee', 'ergauss_dee', 'biasgauss_dee','tdeegauss','Ef_gauss_dee','Varf_gauss_dee','lgauss_dee','sigma2gauss_dee','sigmangauss_dee')

% % %% gauss_eee %
% 
parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);

    [x_seq, y_seq, metric, teeegauss(j), ee, er, bias, Ef, Varf, l, sigma2, sigman] = updategp6(fun, noisestructure, noisevar, Xint, xt, 'gauss', 'EEE', budget,xtt);
    
    lgauss_eee(:,:,j) = l;
    sigma2gauss_eee(:,j) = sigma2;
    sigmangauss_eee(:,j) = sigman;
    
    Ef_gauss_eee(:,:,j) = Ef;
    Varf_gauss_eee(:,:,j) = Varf;
    x_seq_gauss_eee(:,:,j) = x_seq;
    y_seq_gauss_eee(:,j) = y_seq;
    eeegauss(:,j) = metric;
    eegauss_eee(:,j) = ee;
    ergauss_eee(:,j) = er;
    biasgauss_eee(:,j) = bias;
%     erlastgauss_eee(j) = erlast;
end

disp('gausseee finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/gausseee', 'x_seq_gauss_eee', 'y_seq_gauss_eee', 'eeegauss', 'eegauss_eee', 'ergauss_eee', 'biasgauss_eee', 'teeegauss','Ef_gauss_eee','Varf_gauss_eee','lgauss_eee','sigma2gauss_eee','sigmangauss_eee')
 
% % % % % % t_mee %
% 
parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);
    
    [x_seq, y_seq, metric, tmeet(j), ee, er, bias, Ef, Varf, l, sigma2, sigman, nu] = updategp6(fun, noisestructure, noisevar, Xint, xt, 't', 'MEE', budget,xtt);
    
    x_seq_t_mee(:,:,j) = x_seq;
    y_seq_t_mee(:,j) = y_seq;
    meet(:,j) = metric;
    eet_mee(:,j) = ee;
    ert_mee(:,j) = er;
    nut_mee(:,j) = nu;
    
    biast_mee(:,j) = bias;
    lt_mee(:,:,j) = l;
    sigma2t_mee(:,j) = sigma2;
    sigmant_mee(:,j) = sigman;
    
    Ef_t_mee(:,:,j) = Ef;
    Varf_t_mee(:,:,j) = Varf;
%     erlastt_mee(j) = erlast;
end

disp('tmee finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/tmee', 'x_seq_t_mee', 'y_seq_t_mee', 'meet', 'eet_mee', 'ert_mee', 'biast_mee','tmeet','Ef_t_mee', 'Varf_t_mee','lt_mee','sigma2t_mee', 'sigmant_mee', 'nut_mee')
% 
% % % % t_tmse %

parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);
    

    [x_seq, y_seq, metric, ttmset(j), ee, er, bias, Ef, Varf, l, sigma2, sigman, nu] = updategp6(fun, noisestructure, noisevar, Xint, xt, 't', 'tMSE', budget,xtt);
    
    lt_tmse(:,:,j) = l;
    sigma2t_tmse(:,j) = sigma2;
    sigmant_tmse(:,j) = sigman;
    nut_tmse(:,j) = nu;
    
    x_seq_t_tmse(:,:,j) = x_seq;
    y_seq_t_tmse(:,j) = y_seq;
    tmset(:,j) = metric;
    eet_tmse(:,j) = ee;
    ert_tmse(:,j) = er;
    biast_tmse(:,j) = bias;
    
    Ef_t_tmse(:,:,j) = Ef;
    Varf_t_tmse(:,:,j) = Varf;
%     erlastt_tmse(j) = erlast;
end

disp('ttmse finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/ttmse', 'x_seq_t_tmse', 'y_seq_t_tmse', 'tmset', 'eet_tmse', 'ert_tmse', 'biast_tmse','ttmset','Ef_t_tmse','Varf_t_tmse','lt_tmse','sigma2t_tmse', 'sigmant_tmse', 'nut_tmse')

% % t_eee %

parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);
    
    [x_seq, y_seq, metric, teeet(j), ee, er, bias, Ef, Varf, l, sigma2, sigman, nu] = updategp6(fun, noisestructure, noisevar, Xint, xt, 't', 'EEE', budget,xtt);
    
    lt_eee(:,:,j) = l;
    sigma2t_eee(:,j) = sigma2;
    sigmant_eee(:,j) = sigman;
    nut_eee(:,j) = nu;
    
    Ef_t_eee(:,:,j) = Ef;
    Varf_t_eee(:,:,j) = Varf;
    x_seq_t_eee(:,:,j) = x_seq;
    y_seq_t_eee(:,j) = y_seq;
    eeet(:,j) = metric;
    eet_eee(:,j) = ee;
    ert_eee(:,j) = er;
    biast_eee(:,j) = bias;
%     erlastt_eee(j) = erlast;
end

disp('teee finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/teee', 'x_seq_t_eee', 'y_seq_t_eee', 'eeet', 'eet_eee', 'ert_eee', 'biast_eee','teeet','Ef_t_eee','Varf_t_eee','lt_eee','sigma2t_eee', 'sigmant_eee', 'nut_eee')

% % t_dee %
% 
parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);
    
    [x_seq, y_seq, metric, tdeet(j), ee, er, bias, Ef, Varf, l, sigma2, sigman, nu] = updategp6(fun, noisestructure, noisevar, Xint, xt, 't', 'MEESUR', budget,xtt);
    
    lt_dee(:,:,j) = l;
    sigma2t_dee(:,j) = sigma2;
    sigmant_dee(:,j) = sigman;
    nut_dee(:,j) = nu;
    
    Ef_t_dee(:,:,j) = Ef;
    Varf_t_dee(:,:,j) = Varf;
    x_seq_t_dee(:,:,j) = x_seq;
    y_seq_t_dee(:,j) = y_seq;
    deet(:,j) = metric;
    eet_dee(:,j) = ee;
    ert_dee(:,j) = er;
    biast_dee(:,j) = bias;
%     erlastt_dee(j) = erlast;
end

disp('tmeesur finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/tdee', 'x_seq_t_dee', 'y_seq_t_dee', 'deet', 'eet_dee', 'ert_dee', 'biast_dee','tdeet','Ef_t_dee','Varf_t_dee','lt_dee','sigma2t_dee', 'sigmant_dee', 'nut_dee')

% % probit_mee %
% 
parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);
    
    [x_seq, y_seq, metric, tmeeprobit(j), ee, er, bias,Ef, Varf, l, sigma2, sigman] = updategp6(fun, noisestructure, noisevar, Xint, xt, 'probit', 'MEE', budget,xtt);
    
    lprobit_mee(:,:,j) = l;
    sigma2probit_mee(:,j) = sigma2;

    Ef_probit_mee(:,:,j) = Ef;
    Varf_probit_mee(:,:,j) = Varf;
    x_seq_probit_mee(:,:,j) = x_seq;
    y_seq_probit_mee(:,j) = y_seq;
    meeprobit(:,j) = metric;
    eeprobit_mee(:,j) = ee;
    erprobit_mee(:,j) = er;
    biasprobit_mee(:,j) = bias;
%     erlastprobit_mee(j) = erlast;
end

disp('probitmee finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/probitmee', 'x_seq_probit_mee', 'y_seq_probit_mee', 'meeprobit', 'eeprobit_mee', 'erprobit_mee', 'biasprobit_mee','tmeeprobit','Ef_probit_mee', 'Varf_probit_mee','lprobit_mee','sigma2probit_mee')

% probit_tmse %

parfor j = 1:k
      
    rng(j)
    
    Xint = lhsdesign(I,2);
    

    [x_seq, y_seq, metric, ttmseprobit(j), ee, er, bias, Ef, Varf, l, sigma2, sigman] = updategp6(fun, noisestructure, noisevar, Xint, xt, 'probit', 'tMSE', budget,xtt);
    
    lprobit_tmse(:,:,j) = l;
    sigma2probit_tmse(:,j) = sigma2;
    
    Ef_probit_tmse(:,:,j) = Ef;
    Varf_probit_tmse(:,:,j) = Varf;
    x_seq_probit_tmse(:,:,j) = x_seq;
    y_seq_probit_tmse(:,j) = y_seq;
    tmseprobit(:,j) = metric;
    eeprobit_tmse(:,j) = ee;
    erprobit_tmse(:,j) = er;
    biasprobit_tmse(:,j) = bias;
%     erlastprobit_tmse(j) = erlast;
end

disp('probittmse finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/probittmse', 'x_seq_probit_tmse', 'y_seq_probit_tmse', 'tmseprobit', 'eeprobit_tmse', 'erprobit_tmse', 'biasprobit_tmse','ttmseprobit', 'Ef_probit_tmse', 'Varf_probit_tmse','lprobit_tmse','sigma2probit_tmse')

% % % % % probit_eee %

parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);
    
    [x_seq, y_seq, metric, teeeprobit(j), ee, er, bias, Ef, Varf, l, sigma2, sigman] = updategp6(fun, noisestructure, noisevar, Xint, xt, 'probit', 'EEE', budget,xtt);
    
    lprobit_eee(:,:,j) = l;
    sigma2probit_eee(:,j) = sigma2;
   

    Ef_probit_eee(:,:,j) = Ef;
    Varf_probit_eee(:,:,j) = Varf;
    x_seq_probit_eee(:,:,j) = x_seq;
    y_seq_probit_eee(:,j) = y_seq;
    eeeprobit(:,j) = metric;
    eeprobit_eee(:,j) = ee;
    erprobit_eee(:,j) = er;
    biasprobit_eee(:,j) = bias;
%     erlastprobit_eee(j) = erlast;
end

disp('probiteee finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/probiteee', 'x_seq_probit_eee', 'y_seq_probit_eee', 'eeeprobit', 'eeprobit_eee', 'erprobit_eee', 'biasprobit_eee','teeeprobit', 'Ef_probit_eee', 'Varf_probit_eee','lprobit_eee','sigma2probit_eee')

% % % % % probit_dee %

parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);
    
    [x_seq, y_seq, metric, tdeeprobit(j), ee, er, bias, Ef, Varf, l, sigma2, sigman] = updategp6(fun, noisestructure, noisevar, Xint, xt, 'probit', 'MEESUR', budget,xtt);
    
    lprobit_dee(:,:,j) = l;
    sigma2probit_dee(:,j) = sigma2;
    

    Ef_probit_dee(:,:,j) = Ef;
    Varf_probit_dee(:,:,j) = Varf;
    x_seq_probit_dee(:,:,j) = x_seq;
    y_seq_probit_dee(:,j) = y_seq;
    deeprobit(:,j) = metric;
    eeprobit_dee(:,j) = ee;
    erprobit_dee(:,j) = er;
    biasprobit_dee(:,j) = bias;
%     erlastprobit_dee(j) = erlast;
end

disp('probitmeesur finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/probitdee', 'x_seq_probit_dee', 'y_seq_probit_dee', 'deeprobit', 'eeprobit_dee', 'erprobit_dee', 'biasprobit_dee','tdeeprobit', 'Ef_probit_dee', 'Varf_probit_dee','lprobit_dee','sigma2probit_dee')

% % % % % mgauss_mee %
% 
parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);
                                                                                      
    [x_seq, y_seq, metric, tmeemgauss(j), ee, er, bias, Ef, Varf, l, sigma2, sigman] = updategp6(fun, noisestructure, noisevar, Xint, xt, 'mgauss', 'MEE', budget,xtt);
    
    lmgauss_mee(:,:,j) = l;
    sigma2mgauss_mee(:,j) = sigma2;
    sigmanmgauss_mee(:,j) = sigman;
    
    x_seq_mgauss_mee(:,:,j) = x_seq;
    y_seq_mgauss_mee(:,j) = y_seq;
    meemgauss(:,j) = metric;
    eemgauss_mee(:,j) = ee;
    ermgauss_mee(:,j) = er;
    biasmgauss_mee(:,j) = bias;
    Ef_mgauss_mee(:,:,j) = Ef;
    Varf_mgauss_mee(:,:,j) = Varf;
%     erlastmgauss_mee(j) = erlast;
end

disp('mgaussmee finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/mgaussmee', 'x_seq_mgauss_mee', 'y_seq_mgauss_mee', 'meemgauss', 'eemgauss_mee', 'ermgauss_mee', 'biasmgauss_mee', 'tmeemgauss', 'Ef_mgauss_mee', 'Varf_mgauss_mee','lmgauss_mee','sigma2mgauss_mee','sigmanmgauss_mee')

% % mgauss_tmse %

parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);

    [x_seq, y_seq, metric, ttmsemgauss(j), ee, er, bias, Ef, Varf, l, sigma2, sigman] = updategp6(fun, noisestructure, noisevar, Xint, xt, 'mgauss', 'tMSE', budget,xtt);
    
    lmgauss_tmse(:,:,j) = l;
    sigma2mgauss_tmse(:,j) = sigma2;
    sigmanmgauss_tmse(:,j) = sigman;
    
    x_seq_mgauss_tmse(:,:,j) = x_seq;
    y_seq_mgauss_tmse(:,j) = y_seq;
    tmsemgauss(:,j) = metric;
    eemgauss_tmse(:,j) = ee;
    ermgauss_tmse(:,j) = er;
    biasmgauss_tmse(:,j) = bias;
    Ef_mgauss_tmse(:,:,j) = Ef;
    Varf_mgauss_tmse(:,:,j) = Varf;
%     erlastmgauss_tmse(j) = erlast;
    
end

disp('mgausstmse finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/mgausstmse', 'x_seq_mgauss_tmse', 'y_seq_mgauss_tmse', 'tmsemgauss','eemgauss_tmse', 'ermgauss_tmse', 'biasmgauss_tmse','ttmsemgauss','Ef_mgauss_tmse','Varf_mgauss_tmse','lmgauss_tmse','sigma2mgauss_tmse','sigmanmgauss_tmse')

% % mgauss_meesur %

parfor j = 1:k
    
    rng(j+20)
    
    Xint = lhsdesign(I,2);
    
    [x_seq, y_seq, metric, tdeemgauss(j), ee, er, bias, Ef, Varf, l, sigma2, sigman] = updategp6(fun, noisestructure, noisevar, Xint, xt, 'mgauss', 'MEESUR', budget,xtt);
    
    lmgauss_dee(:,:,j) = l;
    sigma2mgauss_dee(:,j) = sigma2;
    sigmanmgauss_dee(:,j) = sigman;
    
    x_seq_mgauss_dee(:,:,j) = x_seq;
    y_seq_mgauss_dee(:,j) = y_seq;
    deemgauss(:,j) = metric;
    eemgauss_dee(:,j) = ee;
    ermgauss_dee(:,j) = er;
    biasmgauss_dee(:,j) = bias;
    Ef_mgauss_dee(:,:,j) = Ef;
    Varf_mgauss_dee(:,:,j) = Varf;
%     erlastmgauss_dee(j) = erlast;
end

disp('mgaussmeesur finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/mgaussmeesur', 'x_seq_mgauss_dee', 'y_seq_mgauss_dee', 'deemgauss', 'eemgauss_dee', 'ermgauss_dee', 'biasmgauss_dee','tdeemgauss', 'Ef_mgauss_dee', 'Varf_mgauss_dee','lmgauss_dee','sigma2mgauss_dee','sigmanmgauss_dee')

% mgauss_eee %

parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);

    [x_seq, y_seq, metric, teeemgauss(j), ee, er, bias, Ef, Varf, l, sigma2, sigman] = updategp6(fun, noisestructure, noisevar, Xint, xt, 'mgauss', 'EEE', budget,xtt);
    
    lmgauss_eee(:,:,j) = l;
    sigma2mgauss_eee(:,j) = sigma2;
    sigmanmgauss_eee(:,j) = sigman;
    
    x_seq_mgauss_eee(:,:,j) = x_seq;
    y_seq_mgauss_eee(:,j) = y_seq;
    eeemgauss(:,j) = metric;
    eemgauss_eee(:,j) = ee;
    ermgauss_eee(:,j) = er;
    biasmgauss_eee(:,j) = bias;
    
    Ef_mgauss_eee(:,:,j) = Ef;
    Varf_mgauss_eee(:,:,j) = Varf;
%     erlastmgauss_eee(j) = erlast;
end
% 
disp('mgausseee finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/mgausseee', 'x_seq_mgauss_eee', 'y_seq_mgauss_eee', 'eeemgauss', 'eemgauss_eee', 'ermgauss_eee', 'biasmgauss_eee', 'teeemgauss', 'Ef_mgauss_eee', 'Varf_mgauss_eee','lmgauss_eee','sigma2mgauss_eee','sigmanmgauss_eee')


% % mprobit_mee %

parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);
    
    [x_seq, y_seq, metric, tmeemprobit(j), ee, er, bias, Ef, Varf, l, sigma2, sigman] = updategp6(fun, noisestructure, noisevar, Xint, xt, 'mprobit', 'MEE', budget,xtt);
    
    lmprobit_mee(:,:,j) = l;
    sigma2mprobit_mee(:,j) = sigma2;
    

    Ef_mprobit_mee(:,:,j) = Ef;
    Varf_mprobit_mee(:,:,j) = Varf;
    x_seq_mprobit_mee(:,:,j) = x_seq;
    y_seq_mprobit_mee(:,j) = y_seq;
    meemprobit(:,j) = metric;
    eemprobit_mee(:,j) = ee;
    ermprobit_mee(:,j) = er;
    biasmprobit_mee(:,j) = bias;
%     erlastmprobit_mee(j) = erlast;
end

disp('mprobitmee finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/mprobitmee', 'x_seq_mprobit_mee', 'y_seq_mprobit_mee', 'meemprobit', 'eemprobit_mee', 'ermprobit_mee', 'biasmprobit_mee','tmeemprobit','Ef_mprobit_mee','Varf_mprobit_mee','lmprobit_mee','sigma2mprobit_mee')

% mprobit_tmse %

parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);
    

    [x_seq, y_seq, metric, ttmsemprobit(j), ee, er, bias, Ef, Varf, l, sigma2, sigman] = updategp6(fun, noisestructure, noisevar, Xint, xt, 'mprobit', 'tMSE', budget,xtt);
    
    lmprobit_tmse(:,:,j) = l;
    sigma2mprobit_tmse(:,j) = sigma2;
    

    x_seq_mprobit_tmse(:,:,j) = x_seq;
    y_seq_mprobit_tmse(:,j) = y_seq;
    tmsemprobit(:,j) = metric;
    eemprobit_tmse(:,j) = ee;
    ermprobit_tmse(:,j) = er;
    biasmprobit_tmse(:,j) = bias;
    Ef_mprobit_tmse(:,:,j) = Ef;
    Varf_mprobit_tmse(:,:,j) = Varf;
%     erlastmprobit_tmse(j) = erlast;
end

disp('mprobittmse finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/mprobittmse', 'x_seq_mprobit_tmse', 'y_seq_mprobit_tmse', 'tmsemprobit', 'eemprobit_tmse', 'ermprobit_tmse', 'biasmprobit_tmse','ttmsemprobit','Ef_mprobit_tmse','Varf_mprobit_tmse','lmprobit_tmse','sigma2mprobit_tmse')

% mprobit_eee %

parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);
    
    [x_seq, y_seq, metric, teeemprobit(j), ee, er, bias, Ef, Varf, l, sigma2, sigman] = updategp6(fun, noisestructure, noisevar, Xint, xt, 'mprobit', 'EEE', budget,xtt);
    
    lmprobit_eee(:,:,j) = l;
    sigma2mprobit_eee(:,j) = sigma2;
    
    x_seq_mprobit_eee(:,:,j) = x_seq;
    y_seq_mprobit_eee(:,j) = y_seq;
    eeemprobit(:,j) = metric;
    eemprobit_eee(:,j) = ee;
    ermprobit_eee(:,j) = er;
    biasmprobit_eee(:,j) = bias;
    Ef_mprobit_eee(:,:,j) = Ef;
    Varf_mprobit_eee(:,:,j) = Varf;
%     erlastmprobit_eee(j) = erlast;
end

disp('mprobiteee finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/mprobiteee', 'x_seq_mprobit_eee', 'y_seq_mprobit_eee', 'eeemprobit', 'eemprobit_eee', 'ermprobit_eee', 'biasmprobit_eee','teeemprobit','Ef_mprobit_eee','Varf_mprobit_eee','lmprobit_eee','sigma2mprobit_eee')

% mprobit_dee %

parfor j = 1:k
    
    rng(j)
    
    Xint = lhsdesign(I,2);
    
    [x_seq, y_seq, metric, tdeemprobit(j), ee, er, bias, Ef, Varf, l, sigma2, sigman] = updategp6(fun, noisestructure, noisevar, Xint, xt, 'mprobit', 'MEESUR', budget,xtt);
    
    lmprobit_dee(:,:,j) = l;
    sigma2mprobit_dee(:,j) = sigma2;
    
    x_seq_mprobit_dee(:,:,j) = x_seq;
    y_seq_mprobit_dee(:,j) = y_seq;
    deemprobit(:,j) = metric;
    eemprobit_dee(:,j) = ee;
    ermprobit_dee(:,j) = er;
    biasmprobit_dee(:,j) = bias;
    Ef_mprobit_dee(:,:,j) = Ef;
    Varf_mprobit_dee(:,:,j) = Varf;
%     erlastmprobit_dee(j) = erlast;
end

disp('mprobitmeesur finished')
save('/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/results/2d_t_hetero_var7/mprobitdee', 'x_seq_mprobit_dee', 'y_seq_mprobit_dee', 'deemprobit', 'eemprobit_dee', 'ermprobit_dee', 'biasmprobit_dee','tdeemprobit','Ef_mprobit_dee','Varf_mprobit_dee','lmprobit_dee','sigma2mprobit_dee')


% 
% surf(linspace(0,1,50)',linspace(0,1,50)',vec2mat(f2,50))
% xlabel('x1')
% ylabel('x2')
% zlabel('f')
% title('Branin-Hoo Function')
% set(gca,'FontSize',18)
% saveas(gcf,'/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/plots/bran3d.png')
% 
% surf(linspace(0,1,50)',linspace(0,1,50)',vec2mat(f,50))
% xlabel('x1')
% ylabel('x2')
% zlabel('f')
% title('Modified Branin-Hoo Function')
% set(gca,'FontSize',18)
% saveas(gcf,'/Users/LittleBear/Desktop/Gaussian process report/Fall 2017/plots/bran3d_modif.png')