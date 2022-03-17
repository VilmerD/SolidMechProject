function file = generateMBB(res, dims)
% Check if geom exists
[exists, file] = findGeometry(res);
if ~exists
    fprintf('Generating geometry\n');
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
    loadw = width*0.05;
    loadh = height*0.05;
    bcload = @(x, y) logical((x < loadw + 1e-6).*...
        (height - y < loadh + 1e-6));
    F.addBoundaryCondition(bcload, 2, 1);
    
    % Make material
    F.make(4, file);
end
end

function [found, filename] = findGeometry(res)
filename = sprintf('%ix%i.mat', res(:));
files_spec = dir(fullfile('**/MBB/', filename));
if isempty(files_spec)
    found = false;
else
    filename = fullfile(files_spec(1).folder, files_spec(1).name);
    found = true;
end
end