function D = dMater2D2(flag, mpara, defgrad)
Ey = mpara(1);      nu = mpara(2);
mu = Ey/(1 + nu)/2; kappa = Ey/(1 - 2*nu)/3;

J = defgrad(1)*defgrad(4) - defgrad(2)*defgrad(3);

efsq = defgrad.^2;
C1 = efsq(1) + efsq(3);
C2 = defgrad(1)*defgrad(2) + defgrad(3)*defgrad(4);
C4 = efsq(2) + efsq(4);
detCs = C1*C4 - C2^2;
tc = (sum(efsq) + 1);

a1 = kappa*J^2 + 2*mu/9*J^(-2/3)*tc;
a2 = 2*mu/3*J^(-2/3);
a3 = mu/3*J^(-2/3)*tc - kappa/2 * (J^2 - 1);


% Total Lagrangian
if flag == 1
    a1 = a1/detCs^2;
    a2 = a2/detCs;
    a3 = a3/detCs^2;
    D1 = a1*C4^2    - 2*a2*C4       + 2*a3*C4^2;
    D2 = a1*C4*C1   - a2*(C4 + C1)  + 2*a3*C2^2;
    D3 = -a1*C4*C2  + a2*C2         - 2*a3*C2*C4;
    D4 = a1*C1^2    - 2*a2*C1       + 2*a3*C1^2;
    D5 = -a1*C1*C2  + a2*C2         - 2*a3*C1*C2;
    D6 = a1*C2^2    - 0             + a3*(C1*C4 + C2^2);
    D = [D1 D2 D3
         D2 D4 D5
         D3 D5 D6];

    % Updated Lagrangian
else
    cinv = [C4  -C2 0
            -C2 C1  0
            0   0   detCs]/detCs;
    Dcomp = @(I, J, K, L) ...
    a1*cinv(I, J)*cinv(K, L) - ...
    a2*((I==J)*cinv(K, L)+(K==L)*cinv(I, J)) + ...
    a3*(cinv(I, K)*cinv(J, L)+cinv(I, L)*cinv(J, K));

    F = reshape(defgrad, 2, 2)';
    F = [F, [0; 0]; [0 0] 1];
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

