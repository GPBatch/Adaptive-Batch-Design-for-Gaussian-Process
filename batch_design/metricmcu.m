function Imcu = metricmcu(Ef, Varf)
% Calculates the local empirical error

% Description:
%   Imcu = METRICMCU(EF, VARF) returns the local empirical error at
%   samples.
%       EF - posterior mean at sample x (can be a vector)
%       VARF - posterior variance at sample x (length same as Ef)

Imcu = normcdf(-abs(Ef)./sqrt(abs(Varf)));
end