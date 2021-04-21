function  Ke = contAxi3te(ec, t, D, ed, stress)
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
r0=(ec(1,1)+ec(1,2)+ec(1,3))/3;
z0=(ec(2,1)+ec(2,2)+ec(2,3))/3;

r32 = ec(1, 3) - ec(1, 2);
r13 = ec(1, 1) - ec(1, 3);
r21 = ec(1, 2) - ec(1, 1); 

z23 = ec(2, 2) - ec(2, 3);
z31 = ec(2, 3) - ec(2, 1);
z12 = ec(2, 1) - ec(2, 2);

N1bar=(ec(1,2)*ec(2,3) - ec(1,3)*ec(2,2)+(ec(2,2)-ec(2,3))*r0 +(ec(1,3)-ec(1,2))*z0);
N2bar=(ec(1,3)*ec(2,1) - ec(1,1)*ec(2,3)+(ec(2,3)-ec(2,1))*r0 +(ec(1,1)-ec(1,3))*z0);
N3bar=(ec(1,1)*ec(2,2) - ec(1,2)*ec(2,1)+(ec(2,1)-ec(2,2))*r0 +(ec(1,2)-ec(1,1))*z0);

% Deformation gradient
F = contAxi3ts(ec, ed) - [1 0 0 1 1]'; 

% A(u)
A = [F(1)   0      F(3)   0     0;                    
     0      F(2)   0      F(4)  0;
     0      0      0      0     F(5);
     F(2)   F(1)   F(4)   F(3)  0];

% Linear part
Bl0 = [z23        0     z31      0     z12      0;
       0         r32    0       r13  	0       r21;
       N1bar/r0   0   N2bar/r0   0   N3bar/r0   0
       r32       z23   r13      z31   r21       z12] / A02;

% Nonlinear part
H0 = [z23        0     z31      0     z12      0;
      r32        0     r13      0     r21      0;
      0         z23     0      z31     0      z12;
      0         r32     0      r13     0      r21
      N1bar/r0   0   N2bar/r0   0   N3bar/r0   0] / A02;
     
B0 = Bl0 + A*H0;

% R0 matrix
R0 = [stress(1), stress(4);
      stress(4), stress(2)];
R0 = [R0         zeros(2, 2); 
     zeros(2, 2) R0];
R0 =[R0         zeros(4,1)
    zeros(1,4)      stress(3)];


Ke = pi*r0*(B0'*D*B0 + H0'*R0*H0)*A02;

end

