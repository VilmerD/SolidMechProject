%% Load data 
element = cont2D4();

width = 300e-3; height = width/3;    
xres = 30;      yres = xres/3;
F = StructureFactory([xres, yres], [width, height]);
F.addBoundaryCondition(@(x, y) abs(x - 0) < 1e-6, 1, 0);
F.addBoundaryCondition(@(x, y) abs(x - width) < 1e-6, [1, 2], 0);
F.addBoundaryCondition(@(x, y) logical((abs(x - 0) < width/xres + 1e-6).*...
    (abs(y - height) < 1e-6)), 2, 1);

geomfile = F.make(element);
load(geomfile);

mpara = [210e9, 0.3];
material = NHCont(1, mpara);

t = 1e-3;
model = NLCont2D([ex, ey], edof, ndof, t, element, material);
solver = LinearSolver(0, 0);
%% Setup displacement controlled data
dmax = -0.3;
xp = bc;
xp(:, 2) = xp(:, 2)*dmax*height;

%% Test displacement controlled algo
K = @model.K;
r = @(u, f) model.fint(u) - f;
nmax = 20;
m = model.ndof;
options = struct('solver', @solver.solveq);
[P, u] = NRDCFAST(K, r, xp, nmax, m, options);
uend = u(:, end);
ed = uend(edof(:, 2:end));

%% Plot
myeldisp2(ex, ey, ed);

%% Stress
S = model.stresses(uend);
Snod = zeros(ndof, 1);
for i = 1:2:ndof
    [r, c] = find(edof(:, 2:end) == i);
    Snod(i) = mean(S(r));
end
fill(ex', ey', Snod(edof(:, 2:2:end))', 'Linestyle', 'None');
axis image;