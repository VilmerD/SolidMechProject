function [P, u] = JFNR(r, u0, p, nmax)
rtol = 1e-6;
utol = 1e-6;
dP = p/nmax;
ep = 1e-4;

P = zeros(length(dP), nmax);
u = [u0, zeros(length(u0), nmax-1)];
for n = 2:nmax
    Pn = P(:, n-1) + dP;
    P(:, n) = Pn;

    un = u(:, n - 1);
    rc = r(un, Pn);

    first_it = 1;
    while first_it || (norm(rc) > rtol && norm(du) > utol)
        K1 = -(rc - r(un + [ep; 0; 0], Pn))/ep;
        K2 = -(rc - r(un + [0; ep; 0], Pn))/ep;
        K3 = -(rc - r(un + [0; 0; ep], Pn))/ep;
        
        [du, ~] = solveq([K1, K2, K3], -rc);
        un = un + du;

        rc = r(un, Pn);

        first_it = 0;
    end
    u(:, n) = un;
end
