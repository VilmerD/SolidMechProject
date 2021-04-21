function defgrad = cont2D4ts(ec, ed)
% Computes the deformation gradient

[dN, ~] = cont2D4N(ec);

function defgrad = F(p1, p2)
    defm = zeros(4, 8);
    defm(1:2, 1:2:7) = dN(p1, p2);
    defm(3:4, 2:2:8) = dN(p1, p2);
    defgrad = defm * ed + [1 0 0 1]';
end
           
x1 = -0.577350269189626;
x2 = -x1;
defgrad = {F(x1, x1), F(x2, x1), F(x2, x2), F(x1, x2)};

end