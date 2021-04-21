function w = wMater2D1(mpara, defgrad)
% Computes the strain energy for st:venant material
%
% Inputs:
%   mpara:      material parameters E and nu    [1 x 2]
%
%   defgrad:    deformation gradient in element [4 x 1]

F = reshape(defgrad, 2, 2)';    F = [F, [0; 0]; [0 0] 1];
E = (F' * F - eye(3))/2;
E = [E(1, 1) E(2, 2) 2*E(1, 2)]';

D = dMater2D1(1, mpara, defgrad);

w = E'*D*E/2;
end