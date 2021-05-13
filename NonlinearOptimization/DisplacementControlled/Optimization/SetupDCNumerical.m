function obj = SetupDC(model, solver, vq, xp, x0, c0)
% Filtering
Mf = model.dfilter();

% Some basic model stuff
volumes = model.volumes();
Vmax = sum(volumes)*vq;
ndof = model.ndof;
nelm = model.nelm;

% Model functions and stuff
% SIMP
Emin = 210e4;    Emax = 210e9;
q = 3;
E =         @(rho) Emin + (Emax - Emin)*rho.^q;
dEdrho =    @(rho) q*(Emax - Emin)*rho.^(q - 1);
K = @(z, uk)            model.Kk(E(Mf*z), uk);
r = @(z, uk, fext)      model.fintk(E(Mf*z), uk) - fext;
drdz = @(z, uk)         model.drdz(dEdrho(Mf*z), uk)*Mf;

NR_OPTIONS.solver = @solver.solveq;
NR_OPTIONS.N_INNER_MAX = 10;
zold = x0;
uold = zeros(ndof, 1);
drdzk = zeros(ndof, nelm);
Kold = zeros(ndof);

% Setup Disp. controlled scheme
nmax = 10;
solver.nsteps = nmax;
u0 = zeros(ndof, 1);
bc = model.bc;

xzero = bc;
xzero(:, 2) = 0;

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

% Numerical sensitivities
ind = [900 1307];
solver.statistics.del_dc = zeros(2*length(ind), 1);

% Numerical sensitivities at some elements
    function dc = numerical_dcdzf(c, z)
        solver.forceFactorization = 1;
        h = 1e-7;
        dc = [];
        for i = ind
            % f(x + h)
            zpe = z;
            zpe(i) = zpe(i) + h;
            zfpe = Mf*zpe;
            
            [Ppe, upe] = solve(zfpe);
            cpe = (Ppe'*upe);
            
            dce = (cpe - c)/h;
            dc = [dc; dce];
        end
        solver.forceFactorization = 0;
    end

% Makes an initial guess using information from dz and dfintdz
% The computational effort used here should be small?
    function du0 = guess(dz)
        [~, du0] = solver.solveq(Kold, -drdzk*dz, xzero, nmax);
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
                switch ME.identifier
                    % If NR didnt converge, force factorization
                    case 'NR:ConverganceError'
                        if solver.forceFactorization ~= 1
                            solver.forceFactorization = 1;
                            [P, u] = solveInner(zold, dz, du0);
                        else
                            nstart = 1;
                            ustart = uold(:, nstart);
                            [P, u] = solveInner(zold, dz, 0);
                        end
                    % If there is an error in the assembly of the FE-
                    % components, we need to start closer to the last step
                    case 'NR:FE_Error'
                        if restarts < 1
                            n = 2;
                            dzk = dz/n;
                            duk = du0/n;
                            for i = 1:n
                                dui = duk*i;
                                dzi = dzk*i;
                                zoldi = zold + dzk*(i - 1);
                                [P, u] = solveInner(zoldi, dzi, dui);
                                ustart = u(:, nstart);
                            end
                        else
                            nstart = 1;
                            ustart = uold(:, nstart);
                            [P, u] = solveInner(zold, dz, 0);
                        end
                    otherwise
                        rethrow(ME)
                end
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

% Computes the function value, the constraints and their derivatives
% for the current design z
    function [g0, g0p, g1, g1p] = cmin(z)
        % Checking angle between designs
        solver.checkAngle(z);
        
        % Solving equilibrium equations
        [F, u] = solve(z);
        
        % Computing objective function
        g0 = objective(F, u);
        
        % Computing senisitivies
        Kold = K(z, u);
        [~, l] = solver.solveq(Kold, f_zero, xp, nmax);
        a([nf; np]) = [l(nf); u(np)];   a = a(:);
        drdzk = drdz(z, u);
        g0p = s*(a'*drdzk)';
        
        % Computing the volume constraint and it's sensitivities
        g1 = volumes'*Mf*z/Vmax - 1;
        g1p = volumes'*Mf/Vmax;
        
        solver.recordUpdate(z);
        k = k + 1;
    end
obj = @cmin;
end