function Ke = cont2D3te(ec, t, D, ed, stress)
% Computes the stiffness matrix for a 2-dimentional three-node element with
% thickness t and
% 
% ec = [x1 x2 x3
%       y1 y2 y3];
%
% ed = [a1 ... a6];
%
% stress = [Sxx;
%           Syy;
%           Sxy];
%
% D = [D11 D12 D13;
%      D21 D22 D23;
%      D31 D32 D33];

A02 = det([1 1 1; ec]);

x32 = ec(1, 3) - ec(1, 2);
x13 = ec(1, 1) - ec(1, 3);
x21 = ec(1, 2) - ec(1, 1); 

y23 = ec(2, 2) - ec(2, 3);
y31 = ec(2, 3) - ec(2, 1);
y12 = ec(2, 1) - ec(2, 2);

% Deformation gradient
F = cont2D3ts(ec, ed) - [1 0 0 1]'; 

% A(u)
A = [F(1) 0      F(3) 0;                    
     0    F(2)   0    F(4);
     F(2) F(1)   F(4) F(3)];

% Linear part
Bl0 = [y23 0     y31 0     y12 0;
       0   x32   0   x13   0   x21;
       x32 y23   x13 y31   x21 y12] / A02;

% Nonlinear part
H0 = [y23 0     y31 0     y12 0;
      x32 0     x13 0     x21 0;
      0   y23   0   y31   0   y12;
      0   x32   0   x13   0   x21] / A02;

B0 = Bl0 + A*H0;

% R0 matrix
R0 = [stress(1), stress(3);
      stress(3), stress(2)];
R0 = [R0         zeros(2, 2); 
     zeros(2, 2) R0];

Ke = (B0'*D*B0 + H0'*R0*H0)*A02*t/2;
end