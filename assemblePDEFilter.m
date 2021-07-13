function [KpMpMs, T] = assemblePDEFilter(res, dims, l0, ls)
% Discretizing the space
[ex, ey, coord, ~, edof] = makeMatQ4Conv(res);
coord = [coord(:, 1)*dims(1) coord(:, 2)*dims(2)];
Ae = dims(1)*dims(2)/(res(1)*res(2));

% Finding boundary dofs for the surface integral
cx = coord(:, 1);   cy = coord(:, 2);
W = max(ex(:));     H = max(ey(:));
So = find(cy == 0);
Ea = find(cx == W);
No = find(cy == H);
We = find(cx == 0);
pairs = [So(1:end-1) So(2:end)
         Ea(1:end-1) Ea(2:end)
         No(1:end-1) No(2:end)
         We(1:end-1) We(2:end)];
     
% Creating base matrices and allocating memory for the sparse matrix
[Kb, Mb, Msb, Tb] = heatFlowMats();
nelm = size(edof, 1);
np = size(edof, 2) - 1;
nne = np^2;
I = zeros(nne*nelm, 1);
J = I;
X = I;
It = zeros(np*nelm, 1);
Jt = reshape(repmat(1:nelm, np, 1), [], 1);
Xt = It;
for e = 1:nelm
    dof = edof(e, 2:end);
    
    % Making the indices of the next entries
    k0 = (e - 1)*nne + 1;
    ke = e*nne;
    ij = meshgrid(dof);
    I(k0:ke) = reshape(ij', 1, nne);
    J(k0:ke) = reshape(ij, 1, nne);
    
    % Entries for big left matrix
    X(k0:ke) = reshape(Ae*(Kb*l0^2 + Mb), [], 1);
    
    % Entries for T-matrix
    k0t = (e - 1)*np + 1;
    ket = e*np;
    It(k0t:ket) = dof;
    Xt(k0t:ket) = Tb*Ae;
end
KpM = sparse(I, J, X);
T = sparse(It, Jt, Xt);

% Surface integral is treated specially
npp = 2;
nnpp = npp^2;
Is = zeros(nnpp*size(pairs, 1), 1);
Js = Is;
Xs = Is;
for e = 1:size(pairs, 1)
    k0 = (e - 1)*nnpp + 1;
    ke = e*nnpp;
    p = pairs(e, :);
    ij = meshgrid(p);
    Is(k0:ke) = reshape(ij', 1, nnpp);
    Js(k0:ke) = reshape(ij, 1, nnpp);
    cp = coord(p, :);
    dl = sqrt(sum((cp(1, :) - cp(2, :)).^2));
    Xs(k0:ke) = Msb*ls*dl;
end
Ms = sparse(Is, Js, Xs);
KpMpMs = KpM + Ms;
end