function [dN, J] = cont2D6N(ec)
dNdxi1 = @(xi1, xi2, xi3) [4*xi1 - 1, 0, 0, 4*xi2, 0, 4*xi3];
dNdxi2 = @(xi1, xi2, xi3) [0, 4*xi2 - 1, 0, 4*xi1, 4*xi3, 0];
dNdxi3 = @(xi1, xi2, xi3) [0, 0, 4*xi3 - 1, 0, 4*xi2, 4*xi1];
dNdxi = @(p1, p2, p3) [dNdxi1(p1, p2, p3)
                       dNdxi2(p1, p2, p3)
                       dNdxi3(p1, p2, p3)];

J = @(p1, p2, p3) ...
    [1 dNdxi1(p1, p2, p3) * ec(1, :)' dNdxi1(p1, p2, p3) * ec(2, :)'
     1 dNdxi2(p1, p2, p3) * ec(1, :)' dNdxi2(p1, p2, p3) * ec(2, :)'
     1 dNdxi3(p1, p2, p3) * ec(1, :)' dNdxi3(p1, p2, p3) * ec(2, :)'];

dN = @(p1, p2, p3) [0 1 0; 0 0 1] * (J(p1, p2, p3) \ dNdxi(p1, p2, p3));
end