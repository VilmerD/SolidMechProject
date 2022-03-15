%% FEM Model
% Load data
dims = [3000e-3 1000e-3];
res = [60 20];
geomfile = generateMBB(res, dims);
load(geomfile);

mpara = [210e9, 0.3];
material = NHCont(1, mpara);

t = 10e-3;
element = cont2D4();
model = NLCont2D(ex, ey, edof, ndof, t, element, material, bc, F);
%% Linear Solver and setting up probelm
vq = 0.3;
x0 = ones(model.nelm, 1)*vq;

alph = -0.5;
nsteps = 10;
radius = dims(1)/res(1)*2;
objective = SetupDC(model, vq, radius, alph, nsteps);

%% Solving problem using MMA
logfile_optID = fopen('Data/Runs/Logs/log.txt', 'a');
mmainitCODC;
mmamainCODC;
