headingform = '%-4s | %-6s %-10s\n';
printform = '%04i | %06.3f %010.9f\n';
if isempty(who('logfile_optID')); logfile_optID = 1; end

fprintf(logfile_optID, headingform, 'iter', 'g0', 'alph');

X = xk;
outit = 0;
extras = cell(minoutit, 1);
while outit < minoutit || (alph > alphtol && outit < maxoutit)
    
    % Linearizing objective function and constraints
    [g0, dg0dx, g1, dg1dx, extra] = objective(xk, outit);
    
    % Solving the mma subproblem at the current design
    [xkp1, ~, ~, ~, ~, ~, ~, ~, ~, low, upp] = ...
        mmasub(m, n, outit, xk, xmin, xmax, xkm1, xkm2, ...
        g0, dg0dx, g1, dg1dx, low, upp, a0, a, c, d);
    
    % Update iterates
    xkm2 = xkm1;
    xkm1 = xk;
    xk = xkp1;
    X(:, outit+1) = xkp1;
    
    % Checking change in design
    alph = extra.alph;
    
    extras{outit + 1} = extra;
    outit = outit + 1;
    fprintf(logfile_optID, printform, outit-1, g0, alph);
end

fprintf(logfile_optID, 'Done \n\n');