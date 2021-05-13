% FEM Model
resolution = 'Coarse';
height = 0.1;
Factory = StructureFactory(resolution, height);
prescribeForce = 1;
geomfile = Factory.makeStructure(prescribeForce);
load(geomfile);

eltype = '2D4t';
t = 1e-3;
mpara = [1, 0.3];
materialModel = 2;

model = NLCont2D([ex, ey], edof, ndof, mpara, t, eltype, bc, materialModel);

%% Filtering
filterRadius = 10e-3;
model.fr = filterRadius;

%% Linear Solver and setting up probelm
% maxits = 4;     % Set to 0 to always force factorization
nbasis = 10;    
solver = LinearSolver(maxits, nbasis);

vq = 0.3;
x0 = ones(nelm, 1)*vq;

amountDisplaced = -0.5;
xp = bc;
xp(:, 2) = xp(:, 2)*height*amountDisplaced;
if amountDisplaced == -0.3
    c0 = 1e2;
elseif amountDisplaced == -0.4
    c0 = 5e2;
elseif amountDisplaced == -0.5
    c0 = 8e2;
end
[objective, Listener] = ...
    SetupDC(model, solver, vq, xp, x0, c0, dzGuess);

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
              'dzGuess', dzGuess, ...
              'amountDisplaced', amountDisplaced, ...
              'date', date, ...
              'geomfile', geomfile);

addpath(genpath('NonlinearOptimization'))
try
    load('NonlinearOptimization\Mats\solutions');
    solutions{end + 1} = data;
catch ME
    if strcmpi(ME.identifier, 'MATLAB:load:couldNotReadFile')
        solutions{1} = data;
    else
        disp(ME.identifier)
        rethrow(ME)
    end
end
save('NonlinearOptimization\Mats\solutions', 'solutions')
