function defgrad = cont2D4ts(H, ed)
% Computes the deformation gradient

defgrad = cell(4, 1);
for k = 1:4
    defgrad{k} = H{k} * ed + [1 0 0 1]';
end
end