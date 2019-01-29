function overhead = CalcOver(theta0, theta1, theta2, n)

% Calculates c_over in ABSUR
%   theta0, theta1 and theta2 - parameters in linear regression
%   n - current design size
overhead = theta0 + theta1*n + theta2.*n.^2;
