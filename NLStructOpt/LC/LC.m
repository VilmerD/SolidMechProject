vq = 0.3;

geomfile = 'beamFullCoarseNew.mat';
load(geomfile);

eltype = '2D4t';
t = 1e-3;
mpara = [1, 0.3];

model = init2D([ex, ey], edof, ndof, mpara, t, eltype, bc, 2);
x0 = ones(nelm, 1)*vq;

fr = 10e-3;
model.fr = fr;
c0 = 32;

maxits = 6;
nb = [2, 4, 12, 16];
lnb = length(nb);

solver = LinearSolver(maxits, nbasis);

objective = SetupLC(model, solver, vq, x0, c0);

mmainit;
mmamain;