function [obj, solver] = setupNLCompOpt(model, vq, x0, c0)
% SIMP
Emin = 210e4;    Emax = 210e9;
q = 3;
E = @(z) Emin + (Emax - Emin)*z.^q;
dEdz = @(z) q*(Emax - Emin)*z.^(q - 1);

% Filtering
Mf = model.dfilter();

% Some basic model stuff
ndof = model.ndof;
volumes = model.volumes();
Vmax = sum(volumes)*vq;

% Negative sign to force factorization at first step, as the angle between
% zold and znew is below the threshold
zold = -Mf*x0;     

% Model functions and stuff
K = @(z, uk)                model.Kk(E(z), uk);
r = @(z, uk, fext)          model.fintk(E(z), uk) - fext;
x_times_dfdz = @(z, uk, x)  model.x_times_dfdz(dEdz(z), uk, x);

% Setup Displacement control
nn = 4;
dmax = -0.01;            
bc = model.bc;         
xp = bc(:, 2)*dmax;

% Linear Solver
maxits = 4;
nbasis = 8;
cothmax = 1 - 5e-3;
solver = LinearSolver(bc(:, 1), ndof, maxits, nbasis);

% Compute scaling factor
if nargin < 4
    fprintf('Computing scaling factor')
    Kt = @(u) K(x0, u);
    rt = @(u, f) r(x0, u, f);
    [P, u, u0, R] = NRDispCont(solver, Kt, rt, xp/nn, nn);
    solver.statistics.residuals = [solver.statistics.residuals R];
    fprintf('\n')
    c0 = P'*u;
else
    u0 = zeros(ndof, 1);
    u0(bc(:, 1)) = xp/nn;
end
s = 1/(10*c0);

k = 0;

    function [P, u] = solve(z)
        % Checking if the changes in design are small enough to reuse the
        % factorization
        coth = z'*zold/(norm(zold)*norm(z));
        if coth < cothmax
            solver.force_factorization = 1;
        end
        fprintf('. coth: %0.4f', coth)
        
        % Setup stiffness matrix and residual
        Kt = @(u) K(z, u);
        rt = @(u, f) r(z, u, f);
        
        % Solving problem
        [P, u, u0, R] = NRDispCont(solver, Kt, rt, xp/nn, nn, u0);
        solver.statistics.residuals = [solver.statistics.residuals R];
        zold = z;
    end

    function [g0, g0p, g1, g1p] = cmin(z)
        % Filtering design
        zf = Mf*z;
        fprintf('\n\tComputing displacements');
        
        % Computing compliance
        [P, u] = solve(zf);
        c = P'*u;
        g0 = c*s;
        
        % Computing sensitivities
        fprintf('\n\tComputing sensitivities')
        
        % Forcing factorization since this problem needs to be solved
        % exactly
        solver.force_factorization = 1;
        lambda = solver.solveq(K(zf, u), -P, bc(:, 2));
        
        % The residuals here are 0 (i think)
        solver.statistics.residuals = [solver.statistics.residuals, 0];
        
        % The sensitivities are given by lambda*dfdz
        dcdr = x_times_dfdz(zf, u, lambda);
        
        g0p = s*Mf*dcdr;
        
        % Computing the volume constraint and it's sensitivities
        g1p = volumes'*Mf/Vmax;
        g1 = g1p*z - 1;
        
        solver.recordUpdate(zf);
        k = k + 1;
    end
obj = @cmin;
end