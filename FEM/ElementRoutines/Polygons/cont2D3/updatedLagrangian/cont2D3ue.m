function Ke = cont2D3ue(ec, t, D, stress)
% Computes the stiffness matrix for a 2-dimentional three-node element with
% thickness t and
% 
% ec = [x1 x2 x3
%       y1 y2 y3];
%
% stress = [sxx;
%           syy;
%           sxy];
%
% D = [D11 D12 D13;
%      D21 D22 D23;
%      D31 D32 D33];

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

% Nonlinear part
H = [y23 0     y31 0     y12 0;
     x32 0     x13 0     x21 0;
     0   y23   0   y31   0   y12;
     0   x32   0   x13   0   x21] / A2;

% R0 matrix
R = [stress(1), stress(3);
     stress(3), stress(2)];
R = [R         zeros(2, 2); 
     zeros(2, 2) R];

Ke = (B'*D*B + H'*R*H)*A2*t/2;
end