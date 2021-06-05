function fe = cont2D4tf(H, B, dJ, t, defgrad, stress)
% Computes the element forces

fe = zeros(8, 1);
Nl = zeros(3, 8);
for k = 1:4
    Hk = H{k};
    F = defgrad{k};
    Nl(1, :) = (F(1)-1)*Hk(1,:)+F(3)*Hk(3,:);
    Nl(2, :) = F(2)*Hk(2,:)+(F(4)-1)*Hk(4,:);
    Nl(3, :) = F(2)*Hk(1,:)+(F(1)-1)*Hk(2,:)+(F(4)-1)*Hk(3,:)+F(3)*Hk(4,:);
    Bk = B{k} + Nl;
    
    s = stress{k};
    fe = fe + Bk'*s*dJ{k}*t;
end
end