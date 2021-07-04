function [P, u] = NRDCFAST(K, r, xp, nmax, ndof, options)
% Computes the displacement controlled response of the system
%
% Inputs
%   Kt:     Stiffness matrix        (function with 1 input)
%
%   r:      Residual force vector   (function with 2 inputs)
%
%   xp:     Prescribed displacement [nbc x 1]
%
%   nmax:   amount of steps
%   
%   ndof:   number of degrees of freedom
%
%   options:
%
rMax = 1e12;
if isfield(options, 'rtol')
    rtol = options.rtol;
else

rtol = 1e-9;
end

if isfield(options, 'u0')
    u0 = options.u0;
    n0 = options.n0;
else
    u0 = zeros(ndof, 1);
    n0 = 1;
end

if isfield(options, 'Verbose')
    verbosity = options.Verbose;
else
    verbosity = 1;
end

if isfield(options, 'N_INNER_MAX')
    N_INNER_MAX = options.N_INNER_MAX;
else
    N_INNER_MAX = 10;
end

% Reset warning so illConditionedMatrix warning can be caught
lastwarn('');
s = warning('off', 'MATLAB:illConditionedMatrix');

% The boundary condition used to correct the solution
correct = xp;
correct(:, 2) = 0;
resn = r(u0, 0);

% Free/Prescribed nodes
np = xp(:, 1);
nf = 1:ndof;
nf(np) = [];

% How much to load each step
dup_k = xp(:, 2)/nmax;            

% Initiating quantities
un = u0;
u = zeros(ndof, nmax);
P = zeros(ndof, nmax);
if verbosity
    printHeading();
end
for n = n0:nmax
    if verbosity
        fprintf('\n%12s%9i%9s', 'Taking Step', n, '');
    end
    % Computing how much should be displaced
    dup_n = n*dup_k - un(np);
    load = [xp(:, 1) dup_n];
    
    % Displacement
    [~, s] = options.solver(K(), -resn, load, n);
    
    un = un + s;
    resn = r(un, 0);                % Residual foces
    r_free = norm(resn(nf));        % Norm of residual in free nodes
    r_tot = norm(resn);
    
   % Check for warnings when assembling r
    checkResidualWarnings();
    
    fprintf('% 1.2e', r_free/r_tot)
    
    % Iterating untill convergance
    N_INNER = 0;
    if verbosity
        fprintf('\n%12s%9i%9i%9s', 'Correcting', n, '', '');
    end
    while r_free/r_tot > rtol
        % Computing new estimate, with zero displacement in prescribed nodes
        [~, s] = options.solver(K(), -resn, correct, n);
        un = un + s;
        resn = r(un, 0);
        r_free = norm(resn(nf));
        r_tot = norm(resn);
        
        checkResidualWarnings();
        
        N_INNER = N_INNER + 1;
        % Update user
        if verbosity
            printAction('', n, N_INNER, r_free/r_tot);
        end
        if N_INNER > N_INNER_MAX || r_free > rMax
            errorStruct.message = sprintf(...
                'Newton failed to converge within %i steps.', N_INNER_MAX);
            errorStruct.identifier = 'NR:ConverganceError';
            error(errorStruct)
        end
    end
    u(:, n) = un;
    P(:, n) = resn;
end
end

function printHeading()
fprintf('\n   Action    n_outer  n_inner     r    ');
end

function printAction(action, nouter, ninner, r)
fprintf('\n%12s%9i%9i% 1.2e', action, nouter, ninner, r)
end

function checkResidualWarnings()
    [warnmsg, ~] = lastwarn;
    if ~isempty(warnmsg)
        errorStruct.message = 'Problems with deformation gradient';
        errorStruct.identifier = 'NR:FE_Error';
        error(errorStruct);
    end
end