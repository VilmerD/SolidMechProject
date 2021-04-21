function sg = stress1D(E, stretch)
% Computes the stress
sg = E/5*(stretch ^ 2 - stretch ^ (-3));