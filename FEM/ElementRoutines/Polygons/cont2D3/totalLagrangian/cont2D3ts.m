function F = cont2D3ts(ec, ed)
% computes the defomration gradient for the 2 dimentional three-node
% element given ec and ed
% 
% ec = [x1 x2 x3
%       y1 y2 y3];
%
% ed = [a1 ... a6];

% Area of the element
A02 = det([1, 1, 1; ec]);

x32 = ec(1, 3) - ec(1, 2);
x13 = ec(1, 1) - ec(1, 3);
x21 = ec(1, 2) - ec(1, 1); 

y23 = ec(2, 2) - ec(2, 3);
y31 = ec(2, 3) - ec(2, 1);
y12 = ec(2, 1) - ec(2, 2);

H0 = [y23 0     y31 0     y12 0;
      x32 0     x13 0     x21 0;
      0   y23   0   y31   0   y12;
      0   x32   0   x13   0   x21] / A02;
     
F = H0 * ed + [1 0 0 1]';
end