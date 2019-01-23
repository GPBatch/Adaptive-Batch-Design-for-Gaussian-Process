%%%%%%%%%%%%%%%%%%%%%%
%' @title Max Call payoff
%'
%' @param x: asset prices
%' @export
%' @inheritParams put.payoff
function y = maxCall(x, K)

if (nargin == 1)
    K = 100;
end

y = max(x,[],2)-K;

y(y<0) = 0;

end