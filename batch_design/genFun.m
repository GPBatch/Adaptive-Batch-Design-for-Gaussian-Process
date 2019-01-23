function [y, r_init, t, sigman]= genFun(X, fun, noisestructure, noisevar, r)

%%%%%%%%%% Generates noisy observation for the true function: the noise is
%%%%%%%%%% normal or t with constant degree or t with heteroskedastic
%%%%%%%%%% degree or mixed Gaussian distributed.

gen_t = tic;

f = fun(X);
n = size(f,1);

switch noisevar
    case 'small'
        sd = 0.1;
    case 'middle'
        sd = 0.2;
    case 'large'
        sd = 1;
    case 'hetero'
        sd = 0.4.*(4.*X(:,1)+1);
    case 'mixed'
        sd = [0.5 1];
end

if (isscalar(r) && r >= 1)
    % handles constant r %
    switch noisestructure
        case 'normal'
            ytmp = random('normal', repmat(f,1,r), sd, n, r);
        case 't_constdf'
            ytmp = random('tlocationscale', repmat(f,1,r), sd, 3, n, r);
        case 't_heterodf'
            ytmp = random('tlocationscale', repmat(f,1,r), repmat(sd, 1, r), repmat(6-4*X(:,1), 1, r), n, r);
        case 'mixed'
            ytmp = zeros(n,r);
            for i = 1:size(X,1)
                for j = 1:r
                    ytmp(i,j) = random(gmdistribution([f(i);f(i)], cat(3,sd(1)^2,sd(2)^2), [0.5, 0.5]));
                end
            end
        case 'null'
            ytmp = zeros(n,r);
            for j = 1:r
                ytmp(:,j) = fun(X);
            end
    end
    y = mean(ytmp,2);
    sigman = var(ytmp,0,2);
    r_init = repmat(r, size(y,1), 1);
else
    % handles r as a vector % 
    l = size(r, 2);
    r = fliplr(r);
    n_seq = [repmat(floor(n/l), l-1, 1); n - floor(n/l)*(l-1)];
    y = zeros(n, 1);
    r_init = zeros(n, 1);
    sigman = zeros(n, 1);
    ind = 0;
    r_seq = r;
    for k = 1:l
            switch noisestructure
                case 'normal'
                    ytmp = random('normal', repmat(f((ind + 1):(ind + n_seq(k))),1,r_seq(k)), sd, n_seq(k), r_seq(k));
                case 't_constdf'
                    ytmp = random('tlocationscale', repmat(f((ind + 1):(ind + n_seq(k))),1,r_seq(k)), sd, 3, n_seq(k), r_seq(k));
                case 't_heterodf'
                    ytmp = random('tlocationscale', repmat(f((ind + 1):(ind + n_seq(k))),1,r_seq(k)), sd, 6-4*X(:,1), n_seq(k), r_seq(k));
                case 'mixed'
                    ytmp = zeros(n_seq(k),r_seq(k));
                    for i = (ind + 1):(ind + n_seq(k))
                        for j = 1:r_seq(k)
                            ytmp(i,j) = random(gmdistribution([f(i);f(i)], cat(3,sd(1)^2,sd(2)^2), [0.5, 0.5]));
                        end
                    end
                case 'null'
                  ytmp = repmat(f((ind + 1):(ind + n_seq(k))), 1, r_seq(k));
            end
        y_tmp = mean(ytmp,2);
        sigman_tmp = var(ytmp,0,2);
        sigman((ind + 1):(ind + n_seq(k))) = sigman_tmp;
        r_init_tmp = repmat(r_seq(k), n_seq(k), 1);
        y((ind + 1):(ind + n_seq(k))) = y_tmp;
        r_init((ind + 1):(ind + n_seq(k))) = r_init_tmp;
        ind = ind + n_seq(k);
    end
end
t = toc(gen_t);
end

