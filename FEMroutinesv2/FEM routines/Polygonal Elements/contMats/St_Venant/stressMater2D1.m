 function s = stressMater2D1(flag, mpara, defgrad)
% flag:     1: Second Piola-Kirchoff
%           2: Kirchhoff
%           3: Cauchy
    
    % Using the energy function we get 
    % S = 2*mu*(E + nu/(1 - 2*nu)*eye_E);
    Ey = mpara(1);  nu = mpara(2);
    mu = Ey/(1 + nu)/2;
    
    % Reshaping deformation gradient and computing the cauchy-green strain
    % tensor
    F = reshape(defgrad, 2, 2)';    F = [F, [0; 0]; [0 0] 1];
    E = (F' * F - eye(3))/2;
                            
    S = 2*mu*(E + nu/(1 - 2*nu)*trace(E)*eye(3));
    
    if flag ~= 1
        S = F*S*F';
        
        % Kirchoff scales by J
        if flag == 3
            S = S/det(F);
        end
    end
    % Reformating
    s = [S(1, 1); S(2, 2); S(1, 2)];
end