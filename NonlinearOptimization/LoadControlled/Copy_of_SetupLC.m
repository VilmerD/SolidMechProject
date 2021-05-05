function obj = Copy_of_SetupLC(model, solver, F, vq, x0, c0)
% Setup Load controlled scheme
ndof = model.ndof;
bc = model.bc;
np = bc(:, 1);
nf = 1:ndof;
nf(np) = [];

% Mapping
mapping_options = struct('filter', model.dfilter());

% Some basic model stuff
volumes = model.volumes()';
Vmax = sum(volumes)*vq;
K = @model.Kk;
r = @model.rk;

% Compute scaling factor
if nargin < 4
    fprintf('Computing scaling factor \n')
    [Pi, ui] = solveFE(x0);
    c0 = Pi'*ui;
end
s = 100/c0;
k = 0;

    % Computes the function value, the constraints and their derivatives
    % for the current design z
    function [g0, dg0dz, g1, dg1dz] = cmin(z)
        % Map the design variables to the densities
        solver.checkAngle(z);
        [E, dEdz] = z_to_E(z, k, mapping_options);       
        
        % Solving equilibrium equations
        [~, u] = solveFE(E, solver, K, r, F, bc, k);
        
        % Computing objective function and senisitivies
        g0 = objective(F, u, s, nf);
        dg0dz = sensitivities(E, dEdz, u);

        % Computing the volume constraint and it's sensitivities
        [g1, dg1dz] = constraints(z, volumes, Vmax, k, mapping_opts);
        
        k = k + 1;
    end
obj = @cmin;
end

% Computes the objective function 
function g0 = objective(F, u, s, nf)
    Ff = F(nf); uf = u(nf);
    c = Ff'*uf;
    g0 = s*c;
end

% Computes the sensitivities
function g0p = sensitivities(E, dEdz, u)
    % Computing adjoint
    [~, l] = solver.solveq(K(E, u), -F, bc, nmax);
        
    % The sensitivities are given by lambda*dfdE
    g0p = s*dcdE(dEdz, u, l);
end

% Solves the equilibrium equations for the current design z
function [P, u] = solveFE(E, solver, K, r, F, bc, k)
% Linear solver options
nmax = 3;

nstart = 1;
if k > 0
    nstart = nmax;
end

% NewtonRaphson options
NR_options = struct('solver', @solver.solveq, ...
                    'solver_options', solver_options, ...
                    'maxiter', 20, ...
                    'nstart', nstart);

K = @(u) K(E, u);
r = @(u, f) r(E, u, f);
[Pfull, ufull] = Copy_of_NRLC(K, r, F, bc, nmax, NR_options);
P = Pfull(:, end);
u = ufull(:, end);
end

% Mapping from the design variables z to the elasticity modulus E
function [E, dEdz] = z_to_E(z, k, mapping_opts)
    % Filtering
    if isfield('filter', mapping_opts)
        zf = mapping_opts.filter*z;
        dzfdz = mapping_opts.filter';
    else
        zf = z;
        dzfdz = 1;
    end
    
    % Projection
    if isfield('projection', mapping_opts)
        rho = mapping_opts.projection(zf, k);
        drhodzf = mapping_opts.projection_sensitivity(zf, k);
    else
        rho = zf;
        drhodzf = 1;
    end
    
    % SIMP
    if isfield('SIMP', mapping_opts)
        q = mapping_opts.SIMP.q;
        Emax = mapping_opts.SIMP.Emax;
        Emin = mapping_opts.SIMP.Emin;
    else
        q = 3;
        Emax = 210e9;
        Emin = 210e5;
    end
    E = Emin + (Emax - Emin)*rho.^q;
    dEdrho = q*(Emax - Emin)*rho.^(q - 1);
    
    dEdz = dzfdz*drhodzf*dEdrho;
end

function [g1, dg1dz] = constraints(z, volumes, Vmax, k, mapping_opts)
    % Filtering
    if isfield('filter', mapping_opts)
        zf = mapping_opts.filter*z;
        dzfdz = mapping_opts.filter';
    else
        zf = z;
        dzfdz = 1;
    end
    
    % Projection
    if isfield('projection', mapping_opts)
        rho = mapping_opts.projection(zf, k);
        drhodzf = mapping_opts.projection_sensitivity(zf, k);
    else
        rho = zf;
        drhodzf = 1;
    end
    
    g1 = volumes*rho/Vmax - 1;
    dg1dz = (dzfdz*drhodzf*volumes/Vmax)';
end