function D = dMater2D1(flag, mpara, defgrad)
    Ey = mpara(1);  nu = mpara(2);
    mu = Ey/(1 + nu)/2;
    
    % Material tangent
    d = 2*mu/(1 - 2*nu)*...
        [1 - nu,    nu,     0;
        nu,         1 - nu, 0;
        0,          0,      (1/2 - nu)];
    
    % Total lagrangian
    if flag == 1
        D = d;
    
    % Updated lagrangian
    else
        F = reshape(defgrad, 2, 2)';    F = [F, [0; 0]; [0 0] 1];
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
        p = [1 1
             1 2
             1 3
             1 3
             
             2 1
             2 2
             2 3
             2 3
             
             3 1
             3 1
             3 2
             3 2
             3 3
             3 3
             3 3
             3 3];
        for IN = 1:16
            R = index(IN, :);
            Dtemp = 0;
            for in = 1:16
                r = index(in, :);
                comp = F(R(1), r(1)) * F(R(2), r(2)) * ...
                    d(p(in, 1), p(in, 2)) ...
                    * F(R(3), r(3)) * F(R(4), r(4));
                Dtemp = Dtemp + comp;
            end
            Ds = [Ds Dtemp];
        end
        D = [Ds(1) Ds(2) Ds(3)
             Ds(5) Ds(6) Ds(7)
             Ds(9) Ds(11) Ds(13)]/det(F);
    end
    
end

