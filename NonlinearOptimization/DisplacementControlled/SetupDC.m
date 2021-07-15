function [obj, Listener, x0] = SetupDC(model, solver, xp)
% Filtering
DF = DensityFilter(model.ec, model.volumes()/model.t, 15e-3);
kupdate = @(k, b) b + 0.5*~mod(k, 10);

SF = SIMPFilter(1, 1e-5, 3);
HPd = HeavisideProjection(kupdate, 0.8);
DesignF = SF * HPd * DF;

HPv = HeavisideProjection(kupdate, 0.5);
VolumeF = HPv * DF;
    
% Some basic model stuff
vq = 0.3;
volumes_normalized = model.volumes()/(sum(model.volumes())*vq);
g1func = @(z) volumes_normalized' * z - 1;
g1pfunc = @(z) volumes_normalized' * VolumeF.backward(z);

ndof = model.ndof;
nelm = model.nelm;
x0 = ones(nelm, 1)*vq;

% Model functions and stuff
K = @(z)            model.K(z);
r = @(z, uk, fext)  model.fint(uk, z) - fext;
dr_dz = @(z)        model.drdE()*DesignF.backward(z);

NR_OPTIONS.solver = @solver.solveq;
NR_OPTIONS.rtol = 1e-6;
zold = x0;
dzTol = 1e0;
drdzk = 0;
MAXITS_MULTIPLICATION_FREQUENCY = 50;

% Setup Disp. controlled scheme
% If newton doesn't converge at the first optimization step increase nmax!!
nmax = 30;
solver.nsteps = nmax;
u0 = zeros(ndof, 1);
uold = u0;
cothmax = 1 - 1e-3;

% Free and prescribed nodes
np = xp(:, 1);
nf = (1:ndof)';
nf(np) = [];
f_zero = zeros(ndof, 1);

k = 0;
s = 1;
Listener = TOListener();

    % Makes an initial guess using information from dz and dfintdz
    function du0 = guess(dz)
        if sqrt(dz'*dz) < dzTol
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
        
        % Inner function that solves the equilibrium equations, and resets
        % if the iteration doesn't converge
        function [P, u] = solveInner(zold, dz, du0)
            znew = zold + dz;
            try
                NR_OPTIONS.n0 = nstart;
                NR_OPTIONS.u0 = ustart + du0;
                [P, u] = NRDCFAST(@() K(znew), @(u, f) r(znew, u, f), ...
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
        
        if numel(P) == 1
            errStruct.message = 'Newton failed to converge';
            errStruct.id = 'NR:ConverganceError';
            error(errStruct);
        end
        
        % Update the old solution z and u
        uold = u;
        zold = z;
        
        Pout = P(:, end);
        uout = u(:, end);
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
        zDFiltered = DesignF.forward(z);
        zVFiltered = VolumeF.forward(z);
        checkAngle(zDFiltered);
        
        % Solving equilibrium equations
        [F, u] = solve(zDFiltered);
        
        % Computing objective function
        Fp = F(np); up = u(np);
        c = Fp'*up;
        if k == 0
            s = -1/c;
        end
        g0 = s*c;
        
        % Computing senisitivies
        [~, l] = solver.solveq(K(zDFiltered), f_zero, xp, nmax);
        a([nf; np]) = [l(nf); up];
        a = a(:);
        drdzk = dr_dz(z);
        dcdz = (a'*drdzk)';
        g0p = s*dcdz;
                
        % Computing the volume constraint and it's sensitivities
        g1 = g1func(zVFiltered);
        g1p = g1pfunc(z);
        
        k = k + 1;
        if mod(k, MAXITS_MULTIPLICATION_FREQUENCY) == 0
            solver.maxits = min(solver.maxits*2, 64);
        end
        stats = solver.getStats();
        stats.g0 = g0;
        stats.g1 = g1;
        stats.design = zDFiltered;
        Listener.registerUpdate(stats);
        Listener.registerCustom('sensitivities', g0p)
    end
obj = @cmin;
end