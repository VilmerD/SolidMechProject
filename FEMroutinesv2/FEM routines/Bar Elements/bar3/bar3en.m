function IntE = bar3en(ec, ed, k, r)
x0 = ec(:, 1) - ec(:, 2);   l0 = norm(x0);
dx = ed(1:3) - ed(4:6);     l = norm(x0 + dx);
lam = l/l0;                 lamc = r / l0;

IntE = 0;
if lam < lamc
    IntE = l0*k/(2*lamc)*(lam - lamc)^2;
end
end