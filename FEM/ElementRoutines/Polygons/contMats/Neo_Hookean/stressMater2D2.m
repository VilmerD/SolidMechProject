function s = stressMater2D2(flag, mpara, defgrad)
% Computes the stress for a plain strain problem. The flag indicates the
% output
% flag:     1: Second Piola-Kirchoff
%           2: Kirchhoff
%           3: Cauchy

% Using the energy function we get
% S = 2*mu*(E + nu/(1 - 2*nu)*eye_E);
Ey = mpara(1);      nu = mpara(2);
mu = Ey/(1 + nu)/2; kappa = Ey/(1 - 2*nu)/3;

if flag == 1
    % Second Piola Kirchoff or otherwise
    % Reshaping deformation gradient and computing the cauchy-green strain
    % tensor
    J = defgrad(1)*defgrad(4) - defgrad(2)*defgrad(3);
    
    efsq = defgrad.^2;
    C1 = efsq(1) + efsq(3);
    C2 = defgrad(1)*defgrad(2) + defgrad(3)*defgrad(4);
    C4 = efsq(2) + efsq(4);
    detCs = C1*C4 - C2^2;
    
    % Computing the second Piola Kirchoff
    Ci = [C4 C1 -C2]'/detCs;
    s = kappa/2*(J ^ 2 - 1)*Ci + mu*J^(-2/3)*([1 1 0]' - (sum(efsq)+1)/3*Ci);
else
    Fs = reshape(defgrad, 2, 2)';
    F = [Fs, [0; 0]; [0 0] 1];
    % Computing
    s = F*S*F';
    
    % Kirchoff scales by J
    if flag == 3
        s = s/det(F);
    end
    s = [s(1, 1); s(2, 2); s(1, 2)];
end
end

