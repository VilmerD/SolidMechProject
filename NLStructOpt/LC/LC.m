% FEM Model
geomfile = 'beamFullCoarseNew.mat';
load(geomfile);

eltype = '2D4t';
t = 1e-3;
mpara = [1, 0.3];

model = init2D([ex, ey], edof, ndof, mpara, t, eltype, bc, 2);

%% Filtering and other
fr = 10e-3;
model.fr = fr;
c0 = 311;

%% Linear Solver and setting up probelm
maxits = 6;
nbasis = 8;

solver = LinearSolver(maxits, nbasis);

vq = 0.3;
x0 = ones(nelm, 1)*vq;

F = zeros(ndof, 1);
Fmax = -1e4;
F(dofs_disp) = Fmax;
objective = SetupLC(model, solver, F, vq, x0, c0);

%% Solving problem using MMA
mmainit;
mmamain;