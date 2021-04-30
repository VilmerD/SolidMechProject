function obj = Copy_of_SetupLC(model, solver, F, vq, x0, c0)


% Setup Load controlled scheme
nmax = 3;
solver.nsteps = nmax;
ndof = model.ndof;
u0 = zeros(ndof, 1);
bc = model.bc;
np = bc(:, 1);
nf = 1:ndof;
nf(np) = [];


% Compute scaling factor
if nargin < 4
    fprintf('Computing scaling factor \n')
    [Pi, ui] = solve(x0);
    c0 = Pi'*ui;
end
s = 100/c0;
k = 0;

    % Computes the function value, the constraints and their derivatives
    % for the current design z
    function [g0, g0p, g1, g1p] = cmin(z)
        % Checking angle between designs
        fprintf('\nComputing displacements');
        solver.checkAngle(z);
        
        % Solving equilibrium equations
        [~, u] = solve(z);
        
        % Computing objective function
        g0 = objective(F, u);
        
        % Computing senisitivies
        fprintf('\nComputing sensitivities')
        [~, l] = solver.solveq(K(z, u), -F, bc, nmax);
        
        % The sensitivities are given by lambda*dfdzf
        g0p = s*dcdz(z, u, l);

        % Computing the volume constraint and it's sensitivities
        g1p = volumes'*Mf/Vmax;
        g1 = g1p*z - 1;
        
        solver.recordUpdate(z);
        k = k + 1;
    end
obj = @cmin;
end

function g0 = objective(F, u)
    Ff = F(nf); uf = u(nf);
    c = Ff'*uf;
    g0 = s*c;
end

% Solves the equilibrium equations for the current design z
function [P, u] = solve(z, u0)
% Setup stiffness matrix and residual
zf = Mf'z;
Kt = @(u) K(z, u);
rt = @(u, f) r(z, u, f);

% Solving problem
if k == 0
    nstart = 1;
else
    nstart = nmax;
end

[Pfull, ufull] = NRLC(solver, Kt, rt, F, nmax, bc, u0, nstart);
P = Pfull(:, end);
u = ufull(:, end);
        
% Detemine which load state to start at next iteration
end

function [r, K] = getNewtonStuff(z)
% SIMP
Emin = 210e4;    Emax = 210e9;
q = 3;
E =         @(rho) Emin + (Emax - Emin)*rho.^q;
dEdrho =    @(rho) q*(Emax - Emin)*rho.^(q - 1);

% Filtering
Mf = model.dfilter();

% Some basic model stuff
volumes = model.volumes()'*Mf;
Vmax = sum(volumes)*vq;  

% Model functions and stuff
K = @(z, uk)              model.Kk(E(Mf*z), uk);
r = @(z, uk, fext)        model.fintk(E(Mf*z), uk) - fext;
dcdz = @(z, uk, x)        Mf'*model.dcdzf(dEdrho(Mf*z), uk, x);
end