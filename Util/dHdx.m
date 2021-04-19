function yp = dHdx(x, b, eta)
yp = b/((tanh(b*eta) + tanh(b*(1 - eta))).*cosh(b*(x - eta)).^2);
end