vq = 0.3;

geomfile = 'beamFullCoarseNew.mat';
load(geomfile);

eltype = '2D4t';
t = 1e-3;
mpara = [1, 0.3];

model = NLCont2D([ex, ey], edof, ndof, mpara, t, eltype, bc, 2);
x0 = ones(nelm, 1)*vq;

fr = 10e-3;
model.fr = fr;
c0 = 32;

solver = LinearSolver()
[objective, sol] = setupNLCompOpt(model, vq, x0, c0);
clear eltype t mpara rho coord dof edof enod F ndof nelm bc

mmainit;
mmamain;

figure()
fill(ex', ey', xk, 'LineStyle', 'none');
axis image;

stats = sol.statistics;
if length(stats.factorizations) < stats.ncalls
    stats.factorizations(stats.ncalls) = 0;
end