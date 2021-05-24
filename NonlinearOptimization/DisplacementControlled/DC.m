% FEM Model
height = 100e-3;
geomfile = 'NonlinearOptimization//Mats//beamSymCoarseNew.mat';
load(geomfile);

eltype = '2D4t';
t = 1e-3;
mpara = [1, 0.3];
materialModel = 2;

model = NLCont2D([ex, ey], edof, ndof, mpara, t, eltype, bc, materialModel);

% Filtering
filterRadius = 15e-3;
model.fr = filterRadius;

%% Linear Solver and setting up probelm
vq = 0.3;
x0 = ones(nelm, 1)*vq;

amountDisplaced = -0.5;
xp = bc;
xp(:, 2) = xp(:, 2)*height*amountDisplaced;

if ~exist('solver') || ~isa(solver, 'LinearSolver')
    solver = LinearSolver(16, 8);
end

[objective, Listener] = SetupDC(model, solver, vq, xp, x0);

%% Solving problem using MMA
mmainit;
mmamain;

%% Save Data
stats = Listener.statistics;
d = datevec(now); date = sprintf('%0.2i%0.2i%0.2i%0.2i', d(2:5));
data = struct('stats', stats, ...
              'maxits', maxits, ...
              'nbasis', nbasis, ...
              'vq', vq, ...
              'filterRadius', filterRadius, ...
              'amountDisplaced', amountDisplaced, ...
              'date', date, ...
              'geomfile', geomfile);

addpath(genpath('NonlinearOptimization'))
try
    load('solutions.mat');
    solutions{end + 1} = data;
catch ME
    if strcmpi(ME.identifier, 'MATLAB:load:couldNotReadFile')
        solutions{1} = data;
    else
        disp(ME.identifier)
        rethrow(ME)
    end
end
save('solutions.mat', 'solutions')
