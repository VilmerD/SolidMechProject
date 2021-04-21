function f = cont2D3uf(ec, t, stress)
% Computes the internval force vector for the 2 dimentional three-node
% element given the thickness t and ec and ed and stress
% 
% ec = [x1 x2 x3
%       y1 y2 y3];
%
% ed = [a1 ... a6];
%
% stress = [Sxx;
%           Syy;
%           Sxy];

% f = B0T*S*A0|t

A2 = det([1 1 1; ec]);

x32 = ec(1, 3) - ec(1, 2);
x13 = ec(1, 1) - ec(1, 3);
x21 = ec(1, 2) - ec(1, 1); 

y23 = ec(2, 2) - ec(2, 3);
y31 = ec(2, 3) - ec(2, 1);
y12 = ec(2, 1) - ec(2, 2);

B = [y23 0     y31 0     y12 0;
     0   x32   0   x13   0   x21;
     x32 y23   x13 y31   x21 y12] / A2;
      
f = B'*stress*t*A2/2;
end