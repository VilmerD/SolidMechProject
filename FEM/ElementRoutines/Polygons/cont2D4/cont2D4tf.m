function fe = cont2D4tf(H, B, dJ, t, defgrad, stress)
% Computes the element forces

fe = zeros(8, 1);
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
    
    Nl(1, :) = F(1)*H1 + F(3)*H3;
    Nl(2, :) = F(2)*H2 + F(4)*H4;
    Nl(3, :) = F(2)*H1 + F(1)*H2 + F(4)*H3 + F(3)*H4;
    Bk = B{k} + Nl;
    
    s = stress{k};
    fe = fe + Bk'*s*dJ{k}*t;
end
end