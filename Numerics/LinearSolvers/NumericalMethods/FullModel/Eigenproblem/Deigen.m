function dwdz = Deigen(psi, w, edof, K0, M0)
% Deigen computes the derivative of w with respect to the design z
% 
% INPUT:
%       psi:
%       w: 
%       edof:
%       K0:
%       M0:
%
% OUTPUT:
%       dwdz:
%

nelm = size(edof, 1);
ndof = size(psi, 1);
ne = size(w, 1);
dl = zeros(nelm, ne);
index_matrix_full = reshape(0:ndof:(ne-1)*ndof, 1, 1, ne) + ...
    repmat(edof(:, 2:end)', 1, 1, ne);
psielm = psi(index_matrix_full);
for elm = 1:nelm
    % Preparing matrices
    me = M0{elm};
    ke = K0{elm};
    we = w;
    
    % Computing derivatives
    psie = psielm(:, elm, :);
    for k = 1:ne
        dl(elm, k) = psie(:, 1, k)' * (ke - we(k)^2 * me) * psie(:, 1, k);
    end
end
dwdz = dl./(2*w');
