function [x, k] = cheapNewton(L, x0)

e = sqrt(eps);
Lp = @(x) (L(x + e) - L(x))/(e);

if nargin < 2
    x0 = 0;
end
xk = x0;
x = [xk];
k = [];
tol = 1e-2;
err = 2*tol;
while err > tol
    ki = Lp(xk);
    xkp1 = xk - (ki)^(-1)*L(xk);
    err = norm(L(xkp1));
    
    k = [k ki];
    x = [x xkp1];
    xk = xkp1;
end

end