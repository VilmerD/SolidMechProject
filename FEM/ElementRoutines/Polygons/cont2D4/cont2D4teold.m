function Ke = cont2D4te(ec, t, D, defgrad, stress, dN)
% Computes the stiffness matrix

function [B, H] = BHmat(dNcurr, F)
    Bl = zeros(3, 8);
    Bl(1, 1:2:7) = dNcurr(1, :);
    Bl(2, 2:2:8) = dNcurr(2, :);
    Bl(3, 2:2:8) = dNcurr(1, :);
    Bl(3, 1:2:7) = dNcurr(2, :);
    
    H = zeros(4, 8);
    H(1:2, 1:2:7) = dNcurr;
    H(3:4, 2:2:8) = dNcurr;
    
    F = F - [1 0 0 1]';
    A = [F(1) 0      F(3) 0;                    
         0    F(2)   0    F(4);
         F(2) F(1)   F(4) F(3)];
    B = Bl + A*H;
end

dNdxi = dN(1:2, :);
dNdeta = dN(3:4, :);

ind = [1 1
       2 1
       2 2
       1 2];      

Ke = zeros(8, 8);
for k = 1:4
    s = stress{k};
    Rk = [s(1) s(3); s(3) s(2)];
    R0 = [Rk        zeros(2)
          zeros(2)  Rk      ];
    
    dNdxik = dNdxi(ind(k, 2), :);
    dNdetak = dNdeta(ind(k, 1), :);
    Jk = [dNdxik * ec'
          dNdetak * ec'];
    dNk = Jk \ [dNdxik; dNdetak];
    [Bk, Hk] = BHmat(dNk, defgrad{k});
    
    Ke = Ke + (Bk'*D{k}*Bk + Hk'*R0*Hk)*det(Jk)*t;
end
end