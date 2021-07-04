%% Load data 
h = 100e-3;
load(geomfile);
z = xk;

%% Make model instance
mpara = [210e9, 0.3];
material = NHCont(1, mpara);
t = 1e-3;
model = NLCont2D([ex ey], edof, ndof, t, element, material);
solver = LinearSolver(0, 0);
%% Setup displacement controlled data
dmax = -0.5;
xp = bc;
xp(:, 2) = xp(:, 2)*dmax*h;

%% Test displacement controlled algo
K = @() model.K(z);
r = @(u, f) model.fint(u, z) - f;
nmax = 12;
m = model.ndof;
options = struct();
options.solver = @solver.solveq;
[P, u] = NRDCFAST(K, r, xp, nmax, m, options);
uend = u(:, end);
%% Plot results
ue = uend(model.edof(:, 2:end));
edx = ue(:, 1:2:end);
edy = ue(:, 2:2:end);

%% Stress
S = model.stresses(uend);
mask = z < 1e-2;
S(mask, :) = 0;
Snod = zeros(ndof, 1);
for i = 1:2:ndof
    [r, c] = find(edof(:, 2:end) == i);
    Snod(i) = mean(S(r));
end
mask = z > 1e-2;
Selm = Snod(edof(:, 2:2:end));
x = ex + edx;
y = ey + edy;
fill(x(mask, :)', y(mask, :)', Selm(mask, :)', 'Linestyle', 'None');
axis image;