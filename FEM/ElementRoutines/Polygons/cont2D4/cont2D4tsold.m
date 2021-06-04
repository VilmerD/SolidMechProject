function defgrad = cont2D4ts(ec, ed, dN)
% Computes the deformation gradient
dNdxi = dN(1:2, :);
dNdeta = dN(3:4, :);

defgrad = cell(4, 1);
ind = [1 1
       2 1
       2 2
       1 2];
for k = 1:4
    H = zeros(4, 8);
    dNdxik = dNdxi(ind(k, 2), :);
    dNdetak = dNdeta(ind(k, 1), :);
    dNk = [dNdxik * ec'
           dNdetak * ec'] \ [dNdxik; dNdetak];
    H(1:2, 1:2:7) = dNk;
    H(3:4, 2:2:8) = dNk;
    defgrad{k} = H * ed + [1 0 0 1]';
end

end