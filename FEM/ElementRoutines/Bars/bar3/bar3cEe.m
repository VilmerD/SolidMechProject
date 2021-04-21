function S = bar3Ecs(ec, ed, edprev, k, r)
% Computes the normal force in the bar element used in the contact model
%
% Inputs:
%   ec:     reference coordinates       [3 x 2]
%
%   ed:     displacement coordinates    [6 x 1]
%
%   edprev: displacement coordinates    [6 x 1]
%           at previous time step
%
%   k:      bar stiffness constant      
%
%   r:      radius of cylinder
x0 = ec(:, 1) - ec(:, 2);               l0 = norm(x0);
dx = ed(1:3) - ed(4:6);                 l = norm(x0 + dx);
dxprev = edprev(1:3) - edprev(4:6);     lprev = norm(x0 + dxprev);

lam = l/l0;                 strain = (lam^2 - 1)/2; 
lamprev = lprev/l0;         strainprev = (lamprev^2 - 1)/2;
dstrain = (strain - strainprev);

wk = bar3en(ec, ed, k, r)/l0;
wprev = bar3en(ec, edprev, k, r)/l0;  

S = 0;
if norm(dstrain) > 0
    S = (wk - wprev)/dstrain^2;
end
end