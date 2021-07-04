function defgrad = defgrad(H, ed)
% Computes the deformation gradient

defgrad = cell(4, 1);
d = [1 0 0 1]';
for k = 1:4
    defgrad{k} = H{k} * ed + d;
end
end