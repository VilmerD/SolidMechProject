%% Load data 
h = 100e-3;
resolution = 'Coarse';
F = SymmetricStructureFactory(resolution, h);
geomfile = F.makeStructure(1, [1, 1]);
load(geomfile);

eltype = '2D4t';
t = 1e-3;
mpara = [210e9, 0.3];
materialModel = 2;

model = NLCont2D([ex, ey], edof, ndof, mpara, t, eltype, bc, materialModel);
solver = LinearSolver(0, 0);
%% Setup displacement controlled data
dmax = -0.3;
xp = bc;
xp(:, 2) = xp(:, 2)*dmax*h;

%% Test displacement controlled algo
K = @(u) model.K(u);
r = @(u, f) model.fint(u) - f;
nmax = 10;
m = model.ndof;
options = struct('solver', @solver.solveq);
[P, u] = NRDC(K, r, xp, nmax, m, options);
uend = u(:, end);