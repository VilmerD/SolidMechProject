fprintf('---- Starting MMA ----');
[g0, dg0dx, g1, dg1dx] = objective(xk);

kktnorm = 2*kkttol;
outit = 0;
fprintf('\n --- Starting iteration --- \n');
while kktnorm > kkttol && outit < maxoutit
    outit = outit + 1;
    fprintf('(%i) %s', outit);
    
    % Solving the mma subproblem at the current design
    [xkp1, ymma, zmma, lam, xsi, eta, mu, zet, s, low, upp] = ...
        mmasub(m, n, outit, xk, xmin, xmax, xkm1, xkm2, ...
        g0, dg0dx, g1, dg1dx, low, upp, a0, a, c, d);
    
    % Update iterates
    xkm2 = xkm1;
    xkm1 = xk;
    xk = xkp1;
    
    % Computing new values
    [g0, dg0dx, g1, dg1dx] = objective(xk);
    
    % Checking kkt-condition
    [~, kktnorm, ~] = kktcheck(m, n, xk, ymma, zmma, lam, xsi, eta, ...
        mu, zet, s, xmin, xmax, dg0dx, g1, dg1dx, a0, a, c, d);
    fprintf('\n')
end

clear kktnorm ymma zmma lam xsi eta mu zet s xmin xmax g0 dg0dx g1 ...
    dg1dx a0 a c d epsimin kkttol low n one outit upp xkm1 xkm2 xkp1 ...
    maxoutit m