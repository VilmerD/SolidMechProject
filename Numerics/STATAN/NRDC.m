function [P, u, ef, es] = NRDC(K, r, sfun, xp, nmax, u0)
% Wraps NRDC
solution_found = false;
u0 = [zeros(size(u0, 1), 1) u0];
RESTARTS_MAX = size(u0, 2);
nrestarts = 0;

while ~solution_found && nrestarts < RESTARTS_MAX
    uold = u0(:, RESTARTS_MAX - nrestarts);
    [P, u, ef, es, flag] = NRDC_inner(K, r, sfun, xp, nmax, uold); 
    if flag == 0
        solution_found = true;
    elseif flag == 1 && nrestarts < RESTARTS_MAX
        nrestarts = nrestarts + 1;
    else
    end
end
if ~solution_found
    errorStruct.message = 'Couldnt solve problem';
    error(errorStruct);
end

end

function [P, u, efn, esn, flag] = NRDC_inner(K, r, sfun, xp, nmax, u0)
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
rtol = 1e-9;
N_INNER_MAX = 10;
ndof = size(u0, 1);

% The boundary condition used to correct the solution
load_correction = xp;
load_correction(:, 2) = 0;
[efn, esn] = sfun(u0);
resn = r(efn, esn, 0);

% Free/Prescribed nodes
np = xp(:, 1);
nf = 1:ndof;
nf(np) = [];

% How much to load each step
% Total length to load from u0
dutot = (xp(:, 2) - u0(np));
ltot = norm(dutot);

% Total length to load from 0 (corresponding to nmax steps)
du0 = norm(xp(:, 2));
lstep0 = du0/nmax;

nstep = max(1, ceil(ltot/lstep0)); % Take at least one step
dustep = dutot/nstep;

% Initiating quantities
un = u0;
u = zeros(ndof, nmax);
P = zeros(ndof, nmax);
printHeader();
for n = (nmax - nstep + 1):nmax
    % Computing how much should be displaced
    load = [xp(:, 1) dustep];
    
    % Displacement
    Kn = K(efn, esn);
    dun = msolveq(Kn, -resn, load);
    
    un = un + dun;
    [efn, esn] = sfun(un);
    resn = r(efn, esn, 0);           % Residual foces
    r_free = norm(resn(nf));        % Norm of residual in free nodes
    r_tot = norm(resn);
    
    % Iterating untill convergance
    N_INNER = 0;
    while r_free/r_tot > rtol
        % Computing new estimate, with zero displacement in prescribed nodes
        Kn = K(efn, esn);
        dun = msolveq(Kn, -resn, load_correction);
        un = un + dun;
        [efn, esn] = sfun(un);
        resn = r(efn, esn, 0);
        r_free = norm(resn(nf));
        r_tot = norm(resn);
        
        N_INNER = N_INNER + 1;
        if N_INNER > N_INNER_MAX || r_free > rMax
            flag = 1;
            return
        end
    end
    u(:, n) = un;
    P(:, n) = resn;
    
    % Update user
    printAction(n, N_INNER, r_free/r_tot);
end
flag = 0;
end

function printHeader()
fprintf('%-3s %-2s %-8s\n', 'out', 'in', 'res');
end

function printAction(nouter, ninner, r)
fprintf('%-3i %-2i %8.2e\n', nouter, ninner, r)
end