function file = generateMBB(res, dims)
% Dimentions and resolution
width = dims(1);
height = dims(2);    
F = StructureFactory(res, dims);

% Symmetry condition
bcsym = @(x, y) abs(x - 0) < 1e-6;
F.addBoundaryCondition(bcsym, 1, 0);

% Bottom right
bcbot = @(x, y) abs(x - width) < 1e-6;
F.addBoundaryCondition(bcbot, [1, 2], 0);

% Load
bcload = @(x, y) logical((abs(x - 0) < 5e-3 + 1e-6).*...
    (abs(y - height) < 2e-3 + 1e-6));
F.addBoundaryCondition(bcload, 2, 1);

% Make material
file = sprintf('Data/Mesh/MBB/%ix%i.mat', res(:));
F.make(4, file);
end