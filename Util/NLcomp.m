function c = NLcomp(model, solver, x, xp, nn)
solver.forceFactorization();
K = @(u) model.Kk(x, u);
r = @(u, f) model.fintk(x, u) - f;
[P, u, u0, R] = NRDispCont(solver, K, r, xp, nn);
c = P'*u;
end