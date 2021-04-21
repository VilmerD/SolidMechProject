function [dN, J] = cont2D4N(ec)
dNdxi = @(xi, eta) 1/4*[-(1 - eta), (1 - eta), (1 + eta) , -(1 + eta)];
dNdeta = @(xi, eta) 1/4*[-(1 - xi), - (1 + xi), (1 + xi) , (1 - xi)];

J = @(p1, p2) [dNdxi(p1, p2) * ec(1, :)' dNdxi(p1, p2) * ec(2, :)'; ...
               dNdeta(p1, p2) * ec(1, :)' dNdeta(p1, p2) * ec(2, :)'];
         
dN = @(p1, p2) J(p1, p2) \ [dNdxi(p1, p2); dNdeta(p1, p2)];
end