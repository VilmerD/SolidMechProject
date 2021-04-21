function s = stressMaterAxi2(flag,mpara,defgrad)
% Computes the stress for a plain strain problem. The flag indicates the
% output
% flag:     1: Second Piola-Kirchoff
%           2: Kirchhoff
%           3: Cauchy
    
    % Using the energy function we get 
    % S = 2*mu*(E + nu/(1 - 2*nu)*eye_E);
    Ey = mpara(1);      nu = mpara(2);
    mu = Ey/(1 + nu)/2; kappa = Ey/(1 - 2*nu)/3;
    
    % Reshaping deformation gradient and computing the cauchy-green strain
    % tensor
    F=defgrad;
    F = [F(1)  F(2) 0
         F(3)  F(4) 0
         0      0   F(5)];
    J = det(F);
    
    C = F' * F;
    cinv = inv(C);

    
    % Computing the second Piola Kirchoff
    S = kappa/2*(J ^ 2 - 1)*cinv + ...
        mu * J^(-2/3) * (eye(3) - trace(C)/3 * cinv);
    
    if flag == 1
        % Second Piola Kirchoff or otherwise
        s = S;
    else
        % Computing
        s = F*S*F';
        
        % Kirchoff scales by J
        if flag == 2
            s = s/det(F);
        end
    end 
    s = [s(1, 1); s(2, 2); s(3,3);s(1, 2)]; 
end

