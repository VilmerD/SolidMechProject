function Ke = cont2D6ue(ec, t, D, stress)
% Computes the stiffness matrix
[dN, J] = cont2D6N(ec);

% Linear term
function B = Bmat(p1, p2, p3)
    B = zeros(3, 12);
    dNcurr = dN(p1, p2, p3);
    B(1, 1:2:11) = dNcurr(1, :);
    B(2, 2:2:12) = dNcurr(2, :);
    B(3, 2:2:12) = dNcurr(1, :);
    B(3, 1:2:11) = dNcurr(2, :);
end

function H = Hmat(p1, p2, p3)
    H = zeros(4, 12);
    dNcurr = dN(p1, p2, p3);
    H(1:2, 1:2:11) = dNcurr;
    H(3:4, 2:2:12) = dNcurr;
end

points = [2, 1, 1
          1, 2, 1
          1, 1, 2]/3;
Hi = 1/3;
Ke = zeros(12, 12);
for k = 1:3
    s = stress{k};
    R0 = [s(1) s(3); s(3) s(2)];
    R0 = [R0 zeros(2, 2); zeros(2, 2) R0];
    
    p1 = points(k, 1); p2 = points(k, 2); p3 = points(k, 3);
    Bk = Bmat(p1, p2, p3);
    Hk = Hmat(p1, p2, p3);
    Jk = J(p1, p2, p3);
    
    Ke = Ke + (Bk'*D{k}*Bk + Hk'*R0*Hk)*det(Jk)*t * Hi;
end
end