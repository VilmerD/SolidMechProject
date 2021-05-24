%% Load data 
h = 0.1;
resolution = 'Fine';
height = 100e-3;
load('solutions.mat');
z = solutions{end}.stats.designs(:, end);

%% Make model instance
matmod = 2;
eltype = '2D4t';
t = 1e-3;
mpara = [210e9, 0.3];
model = NLCont2D([ex ey], edof, ndof, mpara, t, eltype, bc, matmod);

%% Setup displacement controlled data
dmax = -0.5;
xp = bc;
xp(:, 2) = xp(:, 2)*dmax*h;

%% Test displacement controlled algo
K = @(u) model.Kk(z, u);
r = @(u, f) model.fintk(z, u) - f;
nmax = 12;
m = model.ndof;
options = struct();
[P, u] = NRDC(K, r, xp, nmax, m, options);
uend = u(:, end);
%% Plot results
ue = uend(model.edof(:, 2:end));
fill((ex + ue(:, 1:2:end))', (ey + ue(:, 2:2:end))', z);
axis image