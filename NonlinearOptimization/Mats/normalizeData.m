function normalizeData()
% Normalizes the dimensions to 1
files = ["NonlinearOptimization\Mats\beamFullCoarse.mat", ...
    "NonlinearOptimization\Mats\beamFullFine.mat"];

for f = files
    load(f, 'F', 'bc', 'coord', 'dof', 'edof', 'enod')
    width = max(coord(:, 1));
    height = max(coord(:, 2));
    coord(:, 1) = coord(:, 1)/width;
    coord(:, 2) = coord(:, 2)/height;
    save(f, 'F', 'bc', 'coord', 'dof', 'edof', 'enod');
end
end