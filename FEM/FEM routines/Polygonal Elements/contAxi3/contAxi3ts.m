function F = contAxi3ts(ec, ed)

% computes the defomration gradient for the 2 dimentional three-node
% element given ec and ed
% 
% ec = [x1 x2 x3
%       y1 y2 y3];
%
% ed = [a1 ... a6];

% Area of the element

A02 = det([1, 1, 1; ec]);

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

H0 = [z23        0     z31      0     z12      0;
      r32        0     r13      0     r21      0;
      0         z23     0      z31     0      z12;
      0         r32     0      r13     0      r21
      N1bar/r0   0   N2bar/r0   0   N3bar/r0   0] / A02;
     
F = H0 * ed + [1 0 0 1 1]';
end

