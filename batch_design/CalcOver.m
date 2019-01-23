function overhead = CalcOver(theta0, theta1, theta2, n)

%%%%%%%%% This function is used to calculated c_over in ABSUR %%%%%%
overhead = theta0 + theta1*n + theta2.*n.^2;
