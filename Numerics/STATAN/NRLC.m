function [P, u, u0, N] = NRLC(solver, Kt, f, p, nmax, bc, u0)
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
rtol = 1e-6;
dP = p/nmax;

if nargin < 6
    u0 = zeros(length(dP), 1);
end
m = length(dP);

nf = (1:m);
np = bc(:, 1);
if ~isempty(bc)
    nf(np) = [];
end

P = zeros(length(dP), nmax + 1);
u = [u0, zeros(length(u0), nmax)];
NUMBER_OF_ITERATIONS = 0;
for n = 2:nmax + 1
    Pn = P(:, n-1) + dP;
    P(:, n) = Pn;

    un = u(:, n - 1);
    rc = f(un, Pn);
    
    nit_inner = 0;
    while nit_inner == 0 || norm(rc(nf)) > rtol
        K = Kt(un);

        [~, du] = solver.solveq(K, -rc, bc);
        un = un + du;

        rc = f(un, Pn);

        nit_inner = nit_inner + 1;
        if nit_inner > 40
            error("Too large residual")
        end
    end
    u(:, n) = un;
    NUMBER_OF_ITERATIONS = NUMBER_OF_ITERATIONS + nit_inner;
    
    if n == 2
        u0 = un;
    end
end
N = NUMBER_OF_ITERATIONS;
u = u(:, end);
P = P(:, end);
end