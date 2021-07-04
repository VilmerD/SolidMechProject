function s = Stress(flag, mpara, defgrad)
% Computes the stress for a plain strain problem. The flag indicates the
% output
% flag:     1: Second Piola-Kirchoff
%           2: Kirchhoff
%           3: Cauchy
J = defgrad(1)*defgrad(4) - defgrad(2)*defgrad(3);
C1 = defgrad(1)^2 + defgrad(3)^2;
C2 = defgrad(1)*defgrad(2) + defgrad(3)*defgrad(4);
C4 = defgrad(2)^2 + defgrad(4)^2;

% Computing the second Piola Kirchoff
k2 = mpara(1)/(1 + mpara(2))/(2*J^(2/3));
k1 = (mpara(1)/(1 - 2*mpara(2))/3/2*(J ^ 2 - 1) - k2*(C1 + C4 + 1)/3)/(C1*C4 - C2^2);
s = [C4 * k1 + k2; C1 * k1 + k2; -C2*k1];
if flag == 1
    % Second Piola Kirchoff or otherwise
    % Reshaping deformation gradient and computing the cauchy-green strain
    % tensor
    return
else
    S = [s(1) s(3);
         s(3) s(2)];
    Fs = reshape(defgrad, 2, 2)';
    F = [Fs, [0; 0]; [0 0] 1];
    % Computing
    s = F*S*F';
    
    % Kirchoff scales by J
    if flag == 3
        s = s/J;
    end
    s = [s(1, 1); s(2, 2); s(1, 2)];
end
end

