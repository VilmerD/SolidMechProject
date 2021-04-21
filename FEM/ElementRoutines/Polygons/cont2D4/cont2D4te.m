function Ke = cont2D4te(ec, t, D, ed, stress)
% Computes the stiffness matrix
[dN, J] = cont2D4N(ec);

% Linear term
function Bl = Blmat(p1, p2)
    Bl = zeros(3, 8);
    dNcurr = dN(p1, p2);
    Bl(1, 1:2:7) = dNcurr(1, :);
    Bl(2, 2:2:8) = dNcurr(2, :);
    Bl(3, 2:2:8) = dNcurr(1, :);
    Bl(3, 1:2:7) = dNcurr(2, :);
end

function H = Hmat(p1, p2)
    H = zeros(4, 8);
    dNcurr = dN(p1, p2);
    H(1:2, 1:2:7) = dNcurr;
    H(3:4, 2:2:8) = dNcurr;
end

% Non-linear term
function A = Amat(F)
    F = F - [1 0 0 1]';
    A = [F(1) 0      F(3) 0;                    
         0    F(2)   0    F(4);
         F(2) F(1)   F(4) F(3)];
end

function B = Bmat(p1, p2, F)
    B = Blmat(p1, p2) + Amat(F)*Hmat(p1, p2);
end

x1 = -0.577350269189626;
x2 = -x1;
points = [x1, x1;
          x2, x1;
          x2, x2;
          x1, x2];
defgrad = cont2D4ts(ec, ed);


Ke = zeros(8, 8);
for k = 1:4
    s = stress{k};
    R0 = [s(1) s(3); s(3) s(2)];
    R0 = [R0 zeros(2, 2); zeros(2, 2) R0];
    
    p1 = points(k, 1);
    p2 = points(k, 2);
    Bk = Bmat(p1, p2, defgrad{k});
    Hk = Hmat(p1, p2);
    Jk = J(p1, p2);
    
    Ke = Ke + (Bk'*D{k}*Bk + Hk'*R0*Hk)*det(Jk)*t;
end
end