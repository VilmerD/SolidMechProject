function dc = Dcomp(u, edof, K0)
% Dcomp computes the derivatives of the compliance with respect to the
% design z
%
% INPUT:
%
% OUTPUT:
nelm = size(edof, 1);
dc = zeros(nelm, 1);
uelm = u(edof(:, 2:end)');
for elm = 1:nelm
    ue = uelm(:, elm);
    ke = K0{elm};
    dc(elm) = -(ue'*ke*ue);
end
end