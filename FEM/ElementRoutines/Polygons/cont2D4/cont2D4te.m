function Ke = cont2D4te(H, B, dJ, t, D, defgrad, stress)
% Computes the stiffness matrix

Ke = zeros(8, 8);
for k = 1:4
    % Non-linear part
    Hk = H{k};
    H1 = Hk(1, :);
    H2 = Hk(2, :);
    H3 = Hk(3, :);
    H4 = Hk(4, :);
    F = defgrad{k};
    F = F - [1 0 0 1]';
    Nl = zeros(3, 8);

    % Linear part
    Nl(1, :) = F(1)*H1 + F(3)*H3;
    Nl(2, :) = F(2)*H2 + F(4)*H4;
    Nl(3, :) = F(2)*H1 + F(1)*H2 + F(4)*H3 + F(3)*H4;
    Bk = B{k} + Nl;
    
    s = stress{k};
    H1tH2 = H1'*H2;
    H3tH4 = H3'*H4;
    N = s(1)*H1'*H1 + s(3)*(H1tH2 + H1tH2') + s(2)*H2'*H2 +...
        s(1)*H3'*H3 + s(3)*(H3tH4 + H3tH4') + s(2)*H4'*H4;
    Ke = Ke + (Bk'*D{k}*Bk + N)*dJ{k}*t;
end
end