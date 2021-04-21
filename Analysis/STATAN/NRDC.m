function [P, u, uk0] = NRDC(solver, K, r, bc, n, u0)
R = [];
atol = 1e-8; 

resk = r(u0, 0);
bc(:, 2) = bc(:, 2)/n;
ndof = length(u0);

zbc = bc;
zbc(:, 2) = 0;

dd = 1:ndof;
dd(bc(:, 1)) = [];

uk = u0;
for k = 1:n
    fprintf('\n(NR) Taking step')
    % Taking initial step
    dupk = bc(:, 2)*k - uk(bc(:, 1));
    [~, duk] = solver.solveq(K(uk), -resk, [bc(:, 1) dupk]); 
    uk = uk + duk;
    
    % Computing residual forces
    resk = r(uk, 0);
    
    % Computing norm of residual forces in equilibrium
    rk = norm(resk(dd));
    R = [R, rk];
    fprintf(' r: %1.2e', rk)
    
    % Iterating untill convergance
    ninner = 0;
    while rk > atol
       if ninner == 0
           fprintf('\n(NR) Correcting...');
       end
       % Computing new estimate, with zero displacement in prescribed nodes
       [~, duk] = solver.solveq(K(uk), -resk, zbc);
       uk = uk + duk;
       resk = r(uk, 0);
       rk = norm(resk(dd));
       
       R = [R rk];
       ninner = ninner + 1;
       fprintf(' r: %1.2e', rk)
    end
    
    % Saving the first solution for the next step
    if k == 1
       uk0 = uk; 
    end
end
solver.statistics.residuals = [solver.statistics.residuals R];
P = r(uk, 0);
u = uk;