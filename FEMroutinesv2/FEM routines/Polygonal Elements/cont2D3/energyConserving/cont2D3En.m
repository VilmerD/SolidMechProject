function [KinE, IntE] = cont2D3En(flag, mpara, ec, t, rho, defgrad, du)
% Computes the Kinetic and Internal energy in the element
% analysis

% Inputs:
%   flag:       material model, 
%               - 1 (st:venant)
%               - 2 (neohookean)
%
%   mpara:      youngs modulous E and nu            [1 x 2]
%
%   ec:         reference coordinates               [3 x 2]
%
%   t:          thickness of the element 
%
%   rho:        density of element
%
%   defgrad:    deformation gradient in the element
%
%   du:         velocity in the element             [6 x 1]         

M = cont2D3m(ec, t, rho);
KinE = du'*M*du/2;

% Strain energy
if flag == 1
    w = wMater2D1(mpara, defgrad);
elseif flag == 2
    w = wMater2D2(mpara, defgrad);
end

A = det([1 1 1; ec])/2;

% w is constant in the element so the integral becomes w * A * t
IntE = w*A*t;
end