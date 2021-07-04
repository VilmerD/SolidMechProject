function Ke = stiffness(H, B, dJ, t, D, defgrad, stress)
% Computes the stiffness matrix

Ke = zeros(8, 8);
Nl = zeros(3, 8);
for k = 1:4
    % Non-linear part
    Hk = H{k};
    H1 = Hk(1, :);
    H2 = Hk(2, :);
    H3 = Hk(3, :);
    H4 = Hk(4, :);
    F = defgrad{k};

    % Linear part
    % A*H in a faster form
    Nl(1, :) = (F(1) - 1)*H1 + F(3)*H3;
    Nl(2, :) = F(2)*H2 + (F(4) - 1)*H4;
    Nl(3, :) = F(2)*H1 + (F(1) - 1)*H2 + (F(4) - 1)*H3 + F(3)*H4;
    Bk = B{k} + Nl;
    
    s = stress{k};
    N = s(1)*(H1'*H1 + H3'*H3) + s(2)*(H2'*H2 + H4'*H4) + s(3)*(H1'*H2 + H2'*H1 + H3'*H4 + H4'*H3);
    Ke = Ke + (Bk'*D{k}*Bk + N)*dJ{k}*t;
end
end