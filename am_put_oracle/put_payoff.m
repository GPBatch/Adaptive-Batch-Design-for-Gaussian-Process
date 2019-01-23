%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%'  Arithmetic basket Put for d-dim x
%'  @ title American Put payoff
%'
%' @param K is the strike
%' @param x is a vector of asset prices
%' @details in more than 1D, the prices are averaged and maxed with zero.
%' @export
function y = put_payoff(x,K)

if (nargin == 1)
    K = 40;
end

y = K-mean(x,2);
y(y<0) = 0;

end

