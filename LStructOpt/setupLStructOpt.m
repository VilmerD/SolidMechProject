function objective = setupLStructOpt(model, vq, x0)
% SIMP
q = 3;
rho = @(z) z.^q;
drhodz = @(z) q*z.^(q - 1);

% Model stuff
ndof = model.ndof;
bc = model.bc(:, 1);
xp = model.bc(:, 2);
Mf = model.dfilter();
K = @(z) model.K0(rho(z));
dKdz = @(z) model.K0(drhodz(z));
F = model.F;

% Linear Solver
maxits = 5;
solver = LinearSolver(bc, ndof, maxits);
solver.approx_sens = 1;

% Scaling
k = 1;
    function [g0, dg0, g1, dg1] = obj(z)
        zf = Mf*z;
        Ki = model.K0(z);
        [u, ~, l] = solver.solveq(Ki, F, xp);
        c = u'*F;
        g0 = c*k;
        
        
        
    end

objective = @obj;
end
