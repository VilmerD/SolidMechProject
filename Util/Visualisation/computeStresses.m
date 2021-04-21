function S = computeStresses(sys, u, coord, edof)
% Computes the stresses in the system
%
% Inputs:
%   sys:    init2D object
%   u:      displacements
%   coord:  coordinates for the elements
%   edof:   element-degree-of-freedom matrix

[m, n] = size(u);

S = zeros(m, n);
for k = 1:n
    Selm = sys.stresses(u(:, k));
    
    St = zeros(m, 1);
    for i = 1:2:2*size(coord, 1)
        [c0, ~] = find(edof(:, 2:7) == i);
        St(i) = sum(Selm(c0)) / size(c0, 1);
    end
    S(:, k) = St;
    
end