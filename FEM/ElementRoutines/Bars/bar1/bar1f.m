function F = bar1f(k, ed)
% bar1f(k, dx) computes the internal force in the linear spring given
% the stiffness k and displacement dx
du = ed(1) - ed(2);
F = k * [du; -du];