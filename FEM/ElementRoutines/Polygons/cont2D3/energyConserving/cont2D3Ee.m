function Ke = cont2D3Ee(ec, t, D, ed, edprev, es)
% Computes the stiffness matrix used in energy conserving dynamical
% analysis

% Inputs:
%   ec:     reference coordinates       [3 x 2]
%
%   t:      thickness of the element 
%
%   D:      constitutive matrix         [6 x 6]
%
%   ed:     displacement coordinates    [6 x 1]
%
%   edprev: displacement coordinates    [6 x 1]
%           in the previous step
%
%   es:     mean stress in the time     [3 x 1]
%           in the interval (tn, tn+1)
A02 = det([1 1 1; ec]);

x32 = ec(1, 3) - ec(1, 2);
x13 = ec(1, 1) - ec(1, 3);
x21 = ec(1, 2) - ec(1, 1); 

y23 = ec(2, 2) - ec(2, 3);
y31 = ec(2, 3) - ec(2, 1);
y12 = ec(2, 1) - ec(2, 2);

F = cont2D3ts(ec, ed) - [1 0 0 1]'; 
edbar = (ed + edprev)/2;
Fbar = cont2D3ts(ec, edbar) - [1 0 0 1]';

% A(u)
A = [F(1) 0      F(3) 0                   
     0    F(2)   0    F(4)
     F(2) F(1)   F(4) F(3)];

% A(ubar)
Abar = [Fbar(1) 0       Fbar(3) 0                    
        0       Fbar(2) 0       Fbar(4)
        Fbar(2) Fbar(1) Fbar(4) Fbar(3)];

% Linear part
Bl0 = [y23 0     y31 0     y12 0
       0   x32   0   x13   0   x21
       x32 y23   x13 y31   x21 y12] / A02;

% Nonlinear part
H0 = [y23 0     y31 0     y12 0
      x32 0     x13 0     x21 0
      0   y23   0   y31   0   y12
      0   x32   0   x13   0   x21] / A02;

B0 = Bl0 + A*H0;
B0bar = Bl0 + Abar*H0;

% R0 matrix
R0 = [es(1), es(3);
      es(3), es(2)];
R0 = [R0         zeros(2, 2); 
     zeros(2, 2) R0];

Ke = (B0bar'*D*B0 + H0'*R0*H0)*A02*t/2;
end

