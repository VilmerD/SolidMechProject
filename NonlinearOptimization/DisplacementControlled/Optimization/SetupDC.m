function [obj, Listener] = SetupDC(model, solver, vq, xp, x0, c0)
% Filtering
Mf = model.dfilter();

% Some basic model stuff
volumes = model.volumes();
Vmax = sum(volumes)*vq;
ndof = model.ndof;
nelm = model.nelm;

% Model functions and stuff
% SIMP
nuMin = 210e4;      nuMax = 210e9;
q = 3;
nu =                @(rho) nuMin + (nuMax - nuMin)*rho.^q;
dnu_drho =          @(rho) q*(nuMax - nuMin)*rho.^(q - 1);
K = @(z, uk)        model.Kk(nu(Mf*z), uk);
r = @(z, uk, fext)  model.fintk(nu(Mf*z), uk) - fext;
dr_dz = @(z, uk)    model.drdz(dnu_drho(Mf*z), uk)*Mf;

NR_OPTIONS.solver = @solver.solveq;
zold = x0;
uold = zeros(ndof, 1);
dzGuess = 1;
dzTol = 1e0;
drdzk = zeros(ndof, nelm);

% Setup Disp. controlled scheme
nmax = 10;
solver.nsteps = nmax;
u0 = zeros(ndof, 1);
bc = model.bc;
cothmax = 1 - 5e-4;

% Free and prescribed nodes
np = bc(:, 1);
nf = (1:ndof)';
nf(np) = [];
f_zero = zeros(ndof, 1);

% Compute scaling factor
k = 0;
if nargin < 6
    fprintf('Computing scaling factor \n')
    [Pi, ui] = solve(x0);
    c0 = Pi(np)'*ui(np);
end
s = -1/c0;
Listener = TOListener();
    
    % Makes an initial guess using information from dz and dfintdz
    function du0 = guess(dz)
        if dzGuess && sqrt(dz'*dz) < dzTol
            b = -drdzk*dz;
            bf = b(nf);
            du0 = zeros(ndof, 1);
            Rold = solver.Rold{nmax};
            du0(nf) = Rold\(Rold'\bf);
        else
            du0 = 0;
        end
    end

    % Wrapper function for solving the equilibrium equations given the design z
    function [Pout, uout] = solve(z)
        dz = z - zold;
        if k > 0
            factorize = solver.forceFactorization;
            solver.forceFactorization = 0;
            nstart = nmax;
            ustart = uold(:, nstart);
            du0 = guess(dz);
            solver.forceFactorization = factorize;
        else
            nstart = 1;
            ustart = u0;
            du0 = 0;
        end
        % Call the inner function
        restarts = 0;
        [P, u] = solveInner(zold, dz, du0);
        
        % If the inner solver fails, solve in two steps with respect to z
        if numel(P) == 1
            restarts = 0;
            nstart = 1;
            ustart = u0;
            [~, u] = solveInner(zold, dz/2, 0);
            
            restarts = 0;
            nstart = 1;
            ustart = u(:, nstart);
            [P, u] = solveInner(zold + dz/2, dz/2, 0);
        end
        
        % Inner function that solves the equilibrium equations, and resets
        % if the iteration doesn't converge
        function [P, u] = solveInner(zold, dz, du0)
            znew = zold + dz;
            try
                NR_OPTIONS.n0 = nstart;
                NR_OPTIONS.u0 = ustart + du0;
                [P, u] = NRDC(@(u) K(znew, u), @(u, f) r(znew, u, f), ...
                    xp, nmax, ndof, NR_OPTIONS);
            catch ME
                % For some reason there was an error and we need to adjust
                % parameters
                restarts = restarts + 1;
                
                % If the number of restarts is more than 2 stop using
                % shortcuts and start from beginning
                if restarts > 2
                    if restarts == 3
                        nstart = 1;
                        ustart = u0;
                        [P, u] = solveInner(zold, dz, 0);
                        return
                    else
                        P = 0;
                        u = 0;
                        return;
                    end
                end
                switch ME.identifier
                    % If NR didnt converge, try forcing factorization
                    case 'NR:ConverganceError'
                        if solver.forceFactorization ~= 1
                            solver.forceFactorization = 1;
                            [P, u] = solveInner(zold, dz, du0);
                            return;
                        elseif nstart > 1
                            nstart = nstart-1;
                            ustart = uold(:, nstart);
                            [P, u] = solveInner(zold, dz, 0);
                            return;
                        end
                        % If there is an error in the assembly of the FE-
                        % components, we need to start closer to the last step
                    case 'NR:FE_Error'
                        if nstart > 1
                            nstart = nstart-1;
                            ustart = uold(:, nstart);
                            [P, u] = solveInner(zold, dz, 0);
                            return;
                        end
                    otherwise
                        rethrow(ME)
                end
                P = 0;
                u = 0;
            end
        end
        % Update the old solution z and u
        uold = u;
        zold = z;
        
        Pout = P(:, end);
        uout = u(:, end);
    end

    function g0 = objective(F, u)
        Fp = F(np); up = u(np);
        c = Fp'*up;
        g0 = s*c;
    end

    % Checks if the a factorization should be forced
    function checkAngle(znew)
        if k > 0
            coth = zold'*znew/(norm(zold)*norm(znew));
            if coth < cothmax
                solver.forceFactorization = 1;
            else
                solver.forceFactorization = 0;
            end
        else
            solver.forceFactorization = 1;
        end
    end

    % Computes the function value, the constraints and their derivatives
    % for the current design z
    function [g0, g0p, g1, g1p] = cmin(z)
        % Checking angle between designs
        checkAngle(z);
        
        % Solving equilibrium equations
        [F, u] = solve(z);
        
        % Computing objective function
        g0 = objective(F, u);
        
        % Computing senisitivies
        [~, l] = solver.solveq(K(z, u), f_zero, xp, nmax);
        a([nf; np]) = [l(nf); u(np)];   a = a(:);
        drdzk = dr_dz(z, u);
        g0p = s*(a'*drdzk)';
        
        % Computing the volume constraint and it's sensitivities
        g1 = volumes'*Mf*z/Vmax - 1;
        g1p = volumes'*Mf/Vmax;
        
        k = k + 1;
        stats = solver.statistics;
        stats.g0 = g0;
        stats.g1 = g1;
        stats.design = z;
        Listener.registerUpdate(stats);
        stats.ncalls = 0;
        stats.factorizations = 0;
        solver.statistics = stats;
    end
obj = @cmin;
end