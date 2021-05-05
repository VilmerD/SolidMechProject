function [P, u, N] = Copy_of_NRLC(Kt, r, p, bc, n, options)
% Computes the load controlled response of the system
%
% Inputs
%   Kt:     Stiffness matrix        (function with 1 input)
%
%   f:      Internal force vector   (function with 1 input)
%
%   p:      Final force             [m x 1]
%
%   nmax:   amount of steps
%   
%   bc:     boundary conditions     [k x 2]
%
%   u0:     initial displacement    [m x 1]
%           optional, default is 0

% Setup some standard data
if isfield('rtol', options)
    rtol = options.rtol;
else
    rtol = 1e-9;
end
m = length(p);

if isfield('nstart', options)
    nsteps = (n - (nstart - 1));
else
    nsteps = nmax;
end
dP = p/nsteps;

% Setup boundary condition and free noes
nf = (1:m);
np = bc(:, 1);
if ~isempty(bc)
    nf(np) = [];
end

% Initialize solution vectors
P = zeros(m, n + 1);
u = zeros(m, n + 1);
u(:, nstart) = u0;
NUMBER_OF_ITERATIONS = 0;
for n = nstart:n
    % Take a load step
    Pn = P(:, n) + dP;
    P(:, n + 1) = Pn;
    
    % Initialize displacement vector and residual
    un = u(:, n);
    rn = r(un, Pn);
    NORM_rf = norm(rn(nf));
    
    nit_inner = 0;
    
    % Correction step
    while NORM_rf > rtol
        K = Kt(un);

        % Solve for correction and update
        if isfield('solver', options)
            [~, du] = options.solver(K, -rn, bc, options.solver_options);
        else
            [~, du] = solveq(K, -rn, bc);
        end
        
        % Update quantities
        un = un + du;
        rn = r(un, Pn);
        NORM_rf = norm(rn(nf));

        nit_inner = nit_inner + 1;
        % Check if the maximum number of iterations has been reached, and
        % if so stop the scheme.
        if isfield('maxiter', options)
            if nit_inner == options.maxiter
                error("Newton didn't converge.")
            end
        end
    end
    
    % Accept quantities
    u(:, n + 1) = un;
    NUMBER_OF_ITERATIONS = NUMBER_OF_ITERATIONS + nit_inner;
    
    % Check if the user wants to be updated
    if isfield('Verbose', options)
        v = options.Verbose;
        if v == 1
            fprintf('\n(NR) r: %1.2e', NORM_rf);
        end
    end
end
N = NUMBER_OF_ITERATIONS;
end