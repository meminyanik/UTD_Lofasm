[x,k] = secant(@fringeRate,1,5);

function y = fringeRate(x)
% y = x * sin(x) - 1;
y = x^2 - 4*sin(x);
end

function [x,k] = secant(f,x_0,x_1)
% f = nonlinear function
% x_0 = first initial guess
% x_1 = second initial guess
 
% x = root of the function
% k = number of iterations
 
tol = 10^-10;   % tolerance value
x_k_1 = x_0;    % x(k-1)
x_k = x_1;      % x(k)
k = 0;
while abs(f(x_k)) > tol   % termination criterion
    x = x_k - f(x_k) * (x_k - x_k_1) / (f(x_k) - f(x_k_1));
    x_k_1 = x_k;
    x_k = x;
    k = k + 1;
end
end