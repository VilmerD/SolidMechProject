%% Test solution
load('solutions.mat');
data = solutions{end};
stats = data.stats;
load(data.geomfile);

%% Make model instance
matmod = 2;
eltype = '2D4t';
t = 1e-3;
mpara = [1, 0.3];
model = NLCont2D([ex ey], edof, ndof, mpara, t, eltype, bc, matmod);
model.fr = 10e-3;
%%
Emin = 210e4;    Emax = 210e9;
q = 3;
E =         @(rho) Emin + (Emax - Emin)*rho.^q;

Mf = model.dfilter();  
z = stats.designs(:, end);

K = @(u) model.Kk(E(Mf*z), u);
r = @(u, f) model.fintk(E(Mf*z), u) - f;
bc(:, 2) = bc(:, 2)*(-0.05);
options = struct();
[P, u] = NRDC(K, r, bc, 6, model.ndof, options);
uend = u(:, end);
%%
figure(1);
ue = uend(model.edof(:, 2:end));
ux = ue(:, 1:2:end);
uy = ue(:, 2:2:end);
fill((ex + ux)', (ey + uy)', Mf*z);
axis image;