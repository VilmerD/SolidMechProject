function Et = bstiff(ec, ed, k, r)
% Computes the tangent stiffness for the contact model
%
% Inputs:
%   ec:     reference coordinates       [3x2]
%
%   ed:     displacement coordinates    [6x1]
%
%   k:      bar stiffness constant      
%
%   r:      radius of cylinder  
x0 = ec(:, 1) - ec(:, 2);   l0 = norm(x0);
dx = ed(1:3) - ed(4:6);     l = norm(x0 + dx);
lam = l/l0;                 lamc = r / l0;

Et = 0;
if lam < lamc
    Et = k/lam^3;
end
end