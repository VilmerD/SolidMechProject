function w = wMater2D2(mpara, defgrad)
% Computes the strain energy for a neo-hookean material
%
% Inputs:
%   mpara:      material parameters E and nu    [1 x 2]
%
%   defgrad:    deformation gradient in element [4 x 1]

Ey = mpara(1);      nu = mpara(2);
mu = Ey/(1 + nu)/2; kappa = Ey/(1 - 2*nu)/3;

F = reshape(defgrad, 2, 2)';    F = [F, [0; 0]; [0 0] 1];
J = det(F);

C = F' * F;

w = kappa*((J^2 - 1)/2 - log(J))/2 + mu*(J^(-2/3)*trace(C) - 3)/2;
end