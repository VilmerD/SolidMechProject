function [P, u, N] = NRLC(solver, Kt, r, p, nmax, bc, u0, nstart)
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
rtol = 1e-6;
m = length(p);
nsteps = (nmax - (nstart - 1));

if nargin < 7
    nstart = 1;
    if nargin < 6
        u0 = zeros(m, 1);
    end
end
dP = p/nsteps;

% Setup boundary condition and free noes
nf = (1:m);
np = bc(:, 1);
if ~isempty(bc)
    nf(np) = [];
end

% Initialize solution vectors
P = zeros(m, nmax + 1);
u = zeros(m, nmax + 1);
u(:, nstart) = u0;
NUMBER_OF_ITERATIONS = 0;
for n = nstart:nmax
    % Take a load step
    Pn = P(:, n) + dP;
    P(:, n + 1) = Pn;
    
    % Initialize displacement vector and residual
    un = u(:, n);
    rn = r(un, Pn);
    norm_rn_free = norm(rn(nf));
    
    nit_inner = 0;
    
    % Correction step
    while norm_rn_free > rtol
        K = Kt(un);

        % Solve for correction and update
        [~, du] = solver.solveq(K, -rn, bc, n);
        un = un + du;

        rn = r(un, Pn);
        norm_rn_free = norm(rn(nf));

        nit_inner = nit_inner + 1;
        fprintf('\n(NR) r: %1.2e', norm_rn_free);
        if nit_inner > 40
            error("Too large residual")
        end
    end
    
    % Accept quantities
    u(:, n + 1) = un;
    NUMBER_OF_ITERATIONS = NUMBER_OF_ITERATIONS + nit_inner + 1;
    
end
N = NUMBER_OF_ITERATIONS;
end