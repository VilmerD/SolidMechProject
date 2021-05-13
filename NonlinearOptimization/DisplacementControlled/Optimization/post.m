%% Test solution
load('solutions.mat');
data = solutions{3};
stats = data.stats;
load(data.geomfile);
nits = length(stats.factorizations);

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
bc(:, 2) = bc(:, 2)*(-0.03);
options = struct();
[P, u] = NRDC(K, r, bc, 5, model.ndof, options);
uend = u(:, end);
%%
figure(1);
ue = uend(model.edof(:, 2:end));
ux = ue(:, 1:2:end);
uy = ue(:, 2:2:end);
fill((ex + ux)', (ey + uy)', Mf*z);
axis image;

%% Plot the objective function
figure(2)
g0 = stats.g0;  g1 = stats.g1;
ax = subplot(2, 1, 1);
plot(ax, g0, 'Displayname', 'g0')
ymax = max(g0)*1.1; ymin = min(g0)*1.1;
axis(ax, [0, nits, ymin, ymax]);
xlabel('Iteration');
legend();

ax = subplot(2, 1, 2);
plot(ax, g1, 'Displayname', 'g1')
ymax = max(g1)+1e-2; ymin = min(g1)*1.3;
axis(ax, [0, nits, ymin, ymax]);
xlabel('Iteration');
legend();

sgtitle('Function value and kkt condition')
%% Number of factorizations
figure(3);
n = stats.lineqs;
f = stats.factorizations;
coth = computeCoth(stats.designs);
ymax = max(max(n), max(f));
ymin = 0;
xmax = nits;
xmin = 1;
subplot(211)
axis([xmin, xmax, ymin, ymax]);
hold on;
plot(f, 'g', 'Displayname', '# of factorizations');
plot(n, 'r', 'Displayname', '# of lineqs. solved');
xlabel("Iteration")
legend();
hold off;

subplot(212);
semilogy(coth, 'b', 'Displayname', 'Similarity in design');
axis([1, xmax, min(coth), 1]);
xlabel("Iteration")
legend();

%% DZ plot
figure(4)
designs = stats.designs;
[nelm, nd] = size(designs);
dz = zeros(nd-1, 1);
for k = 1:nd-1
    dz(k) = norm(designs(:, k+1) - designs(:, k));
end
semilogy(abs(dz));