function simulateDesign(geomfile, dmax)
h = 100e-3;
load(geomfile, 'bc', 'edof', 'ex', 'ey', 'eltype', 'ndof');

% Make model instance
matmod = 2;
t = 1e-3;
mpara = [210e9, 0.3];
model = NLCont2DFast([ex ey], edof, ndof, mpara, t, eltype, matmod);
solver = LinearSolver(0, 0);

% Setup displacement controlled data
if nargin < 2
    dmax = -0.3;
end
xp = bc;
xp(:, 2) = xp(:, 2)*dmax*h;

% Run displacement simulation
K = @() model.K(z);
r = @(u, f) model.fint(u, z) - f;
nmax = 12;
m = model.ndof;
options = struct();
options.solver = @solver.solveq;
[~, u] = NRDCFAST(K, r, xp, nmax, m, options);
uend = u(:, end);

% Plot results
ue = uend(model.edof(:, 2:end));
myeldisp2(ex, ey, ue, xk);
end