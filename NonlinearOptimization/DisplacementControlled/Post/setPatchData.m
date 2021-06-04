function setPatchData(P, feild, data)
%SETPATCHDATA sets the CData property of patches in P to data
%   Speeds the animation as new patch objects don't need to be created each
%   time the image is rendered
[np, ~] = size(P);
[nd, ~] = size(data);
if np ~= nd
    errorStruct.message = ...
        'Mismatch in number of patches and number of data points';
    error(errorStruct);
end
for i = 1:np
    set(P(i), feild, data(i, :));
end
end

