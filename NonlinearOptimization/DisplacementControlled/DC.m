%% FEM Model
% Load data
dims = [3000e-3 1000e-3];
res = [90 30];
geomfile = generateMBB(res, dims);
load(geomfile);

t = 10e-3;
material = NHCont(1, [210e9, 0.3]);
element = cont2D4();

% k0, bmax, kstep, kincr, eta
threshpara = {150, 50, 5, 1, 0.5};

% about 2 elements wide
radius = dims(1)/res(1)*2;
filterpara = {radius};

model = NLContModel_opt.makeModel(geomfile, t, element, material, ...
    threshpara, filterpara);
%% Linear Solver and setting up probelm
vq = 0.3;
x0 = ones(model.nelm, 1)*vq;

alph = -0.5;
nsteps = 10;
objective = SetupDC(model, vq, alph, nsteps);

%% Solving problem using MMA
mmainitCODC;
mmamainCODC;
