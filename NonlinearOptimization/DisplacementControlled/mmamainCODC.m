logfilename = 'logoptID';
if isempty(who(logfilename)); logoptID = 1; end

headingform = '%-4s | %5s %5s %5s %5s | %8s %8s %8s \n';
printformall = '%4i | %5i %5i %5i %5.0e | %+8.3f %+8.0e %8.0e \n';
printform_init = '%4i | ';
printform_final = ' | %+8.3f %+8.0e %8.0e \n';
fprintf(logoptID, headingform, 'iter', 'nout', 'nin', ...
    'rstrt', 'res', 'g0', 'g1', 'alph');

X = xk;
outit = 0;
extras = cell(minoutit, 1);
while outit < minoutit || (alph > alphtol && outit < maxoutit)
    fprintf(logoptID, printform_init, outit);
    
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
    fprintf(printform_final, g0, g1, alph);
end

fprintf(logfile_optID, 'Done \n\n');