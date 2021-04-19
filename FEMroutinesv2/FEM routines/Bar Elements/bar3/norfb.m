function N = norfb(ec, ed, k, r)
% Computes the normal force in the bar element used in the contact model
%
% Inputs:
%   ec:     reference coordinates       [3 x 2]
%
%   ed:     displacement coordinates    [6 x 1]
%
%   k:      bar stiffness constant      
%
%   r:      radius of cylinder
x0 = ec(:, 1) - ec(:, 2);   l0 = norm(x0);
dx = ed(1:3) - ed(4:6);     l = norm(x0 + dx);
lam = l/l0;                 lamc = r / l0;

N = 0;
if lam < lamc
    N = k * (lam/lamc - 1) / lam;
end
end