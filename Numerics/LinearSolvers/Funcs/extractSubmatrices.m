function [A, B, C, D] = extractSubmatrices(K, bc, dd)
%EXTRACTSUBMATRICES extracts the submatrices of K correspoding to the
%indeces bc and it's complement dd
%   Detailed explanation goes here
            
A = sparse(K(dd, dd));
B = sparse(K(dd, bc));
C = sparse(K(bc, dd));
D = K(bc, bc);
end

