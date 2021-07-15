%% FEM Model
% Load data 
element = cont2D4();

q = 3;
width = 300e-3;
height = width/q;    
xres = 45;
yres = xres/q;
F = StructureFactory([xres, yres], [width, height]);
F.addBoundaryCondition(@(x, y) abs(x - 0) < 1e-6, 1, 0);
F.addBoundaryCondition(@(x, y) abs(x - width) < 1e-6, [1, 2], 0);
F.addBoundaryCondition(@(x, y) logical((abs(x - 0) < 5e-3 + 1e-6).*...
    (abs(y - height) < 2e-3 + 1e-6)), 2, 1);

geomfile = F.make(element);
load(geomfile);

mpara = [210e9, 0.3];
material = NHCont(1, mpara);

t = 1e-3;
model = NLCont2D([ex, ey], edof, ndof, t, element, material);
%% Linear Solver and setting up probelm
amountDisplaced = -0.5;

xp = bc;
xp(:, 2) = xp(:, 2)*height*amountDisplaced;

if ~exist('solver', 'var') || ~isa(solver, 'LinearSolver')
    solver = LinearSolver(2, 10);
end

[objective, Listener, x0] = SetupDC(model, solver, xp);

%% Solving problem using MMA
mmainit;
mmamain;

%% Save Data
stats = Listener.statistics;
d = datevec(now); date = sprintf('%0.2i%0.2i%0.2i%0.2i', d(2:5));
data = struct('stats', stats, ...
              'maxits', solver.maxits, ...
              'nbasis', solver.nbasis, ...
              'amountDisplaced', amountDisplaced, ...
              'date', date, ...
              'geomfile', geomfile, ...
              'element', element);
              

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
