%% Load data 
load beamFullCoarseDC.mat;

%% Make model instance
matmod = 2;
eltype = '2D4t';
t = 1e-3;
mpara = [210e9, 0.3];
model = NLCont2D([ex ey], edof, ndof, mpara, t, eltype, bc, matmod);

%% Setup displacement controlled data
dmax = -0.3;
h = 0.1;
bc(:, 2) = bc(:, 2)*dmax*h;

%% Test displacement controlled algo
K = @model.K;
r = @(u, f) model.fint(u) - f;
nmax = 3;
m = model.ndof;
options = struct();
[P, u] = NRDC(K, r, bc, nmax, m, options);

%% Plot results
ue = u(model.edof(:, 2:end));
eldisp2(ex, ey, ue, [1 4 1], 1);

%% End compliance
np = bc(:, 1);
c = u(np)'*P(np);