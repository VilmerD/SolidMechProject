function D = dMater2D2(flag, mpara, defgrad)
    Ey = mpara(1);      nu = mpara(2);
    mu = Ey/(1 + nu)/2; kappa = Ey/(1 - 2*nu)/3;
    
    F = reshape(defgrad, 2, 2)';    F = [F, [0; 0]; [0 0] 1];
    J = det(F);
    
    C = F' * F;
    cinv = inv(C);
    
    a1 = kappa*J^2 + 2*mu/9*J^(-2/3)*trace(C);
    a2 = 2*mu/3*J^(-2/3);
    a3 = mu/3*J^(-2/3)*trace(C) - kappa/2 * (J^2 - 1);
    
    Dcomp = @(I, J, K, L) ...
        a1*cinv(I, J)*cinv(K, L) - ...
        a2*((I==J)*cinv(K, L)+(K==L)*cinv(I, J)) + ...
        a3*(cinv(I, K)*cinv(J, L)+cinv(I, L)*cinv(J, K));
    
    % Total Lagrangian
    if flag == 1
        Ds = [Dcomp(1, 1, 1, 1)
              Dcomp(1, 1, 2, 2)
              Dcomp(1, 1, 2, 1)
              Dcomp(2, 2, 2, 2)
              Dcomp(2, 2, 2, 1)
              Dcomp(2, 1, 2, 1)];
        D = [Ds(1) Ds(2) Ds(3)
             Ds(2) Ds(4) Ds(5)
             Ds(3) Ds(5) Ds(6)];
         
    % Updated Lagrangian
    else
        Ds = [];
        index = [1 1 1 1
                 1 1 2 2
                 1 1 2 1
                 1 1 1 2
                 
                 2 2 1 1
                 2 2 2 2
                 2 2 2 1
                 2 2 1 2
                 
                 2 1 1 1
                 1 2 1 1
                 2 1 2 2
                 1 2 2 2
                 2 1 2 1
                 2 1 1 2
                 1 2 1 2
                 1 2 2 1];
        for IN = 1:16
            R = index(IN, :);
            Dtemp = 0;
            for in = 1:16
                r = index(in, :);
                comp = F(R(1), r(1)) * F(R(2), r(2)) * ...
                    Dcomp(r(1), r(2), r(3), r(4)) ...
                    * F(R(3), r(3)) * F(R(4), r(4));
                Dtemp = Dtemp + comp;
            end
            Ds = [Ds Dtemp];
        end
        D = [Ds(1) Ds(2) Ds(3)
             Ds(5) Ds(6) Ds(7)
             Ds(9) Ds(11) Ds(13)]/J^1;
    end  
end

