% FEM Model
resolution = 'Coarse';
height = 0.1;
Factory = StructureFactory(resolution, height);
geomfile = Factory.makeStructure(1, 1, 1);
load(geomfile);

eltype = '2D4t';
t = 1e-3;
mpara = [1, 0.3];
materialModel = 2;

model = NLCont2D([ex, ey], edof, ndof, mpara, t, eltype, bc, materialModel);
%% Filtering and other
FilterRadius = 10e-3;
model.fr = FilterRadius;
% c0 = 1e4; % if disp = 30%
c0 = 800;   % if disp = 50%
%% Linear Solver and setting up probelm
maxits = 4;
nbasis = 10;
solver = LinearSolver(maxits, nbasis);

vq = 0.3;
x0 = ones(nelm, 1)*vq;

height = 0.1;
amount_displaced = -0.5;
xp = bc;
xp(:, 2) = xp(:, 2)*height*amount_displaced;

objective = SetupDC(model, solver, vq, xp, x0, c0);
%% Solving problem using MMA
mmainit;
mmamain;

%% Save Data
stats = solver.statistics;
d = datevec(now); date = sprintf('%0.2i%0.2i%0.2i%0.2i', d(2:5));
data = struct('stats', stats, 'maxits', maxits, 'nbasis', nbasis, ...
    'vq', vq, 'FilterRadius', FilterRadius, 'date', date, ...
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
