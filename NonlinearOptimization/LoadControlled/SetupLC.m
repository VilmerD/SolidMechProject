function obj = SetupLC(model, solver, F, vq, x0, c0)
% SIMP
Emin = 210e4;    Emax = 210e9;
q = 3;
E =         @(rho) Emin + (Emax - Emin)*rho.^q;
dEdrho =    @(rho) q*(Emax - Emin)*rho.^(q - 1);

% Filtering
Mf = model.dfilter();

% Some basic model stuff
volumes = model.volumes();
Vmax = sum(volumes)*vq;

% Model functions and stuff
K = @(z, uk)            model.Kk(E(Mf*z), uk);
r = @(z, uk, fext)      model.fintk(E(Mf*z), uk) - fext;
drdz = @(z, uk)         model.drdz(dEdrho(Mf*z), uk)*Mf;
NR_OPTIONS.solver = @solver.solveq;
zold = x0;
uold = 0;

% Guess correction
BETA_SHRINK = 0.5;
NUMBER_OF_SHRINKAGES_MAX = 2;

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

% Numerical sensitivities
ind = [900 1307];
solver.statistics.del_dc = zeros(2*length(ind), 1);

% Numerical sensitivities at some elements
    function dc = numerical_dcdzf(c, z)
        solver.force_factorization = 1;
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
        solver.force_factorization = 0;
    end

% Makes an initial guess using information from dz and dfintdz
% The computational effort used here should be small?
    function du0 = guess(zold, dz, uold)
        [~, du0] = solver.solveq(K(zold, uold), ...
            -drdz(zold, uold)*dz, bc, nmax);
    end

% Wrapper function for solving the equilibrium equations given the design z
    function [Pout, uout] = solve(z)
        beta = 1;
        if k > 0
            nstart = nmax;
            ustart = uold(:, nstart);
            du0 = guess(zold, dz, ustart);
        else
            nstart = 1;
            ustart = u0;
            du0 = 0;
        end
        % Call the inner function
        [P, u] = solveInner(ustart, nstart);
        
        
        % Inner function that solves the equilibrium equations, and resets
        % if the iteration doesn't converge
        function [P, u] = solveInner(ustart, nstart)
            NR_OPTIONS.u0 = ustart + beta*du0;
            NR_OPTIONS.n0 = nstart;
            
            % Attempt solving problem
            try
                [P, u] = NRLC(@(u) K(z, u), @(u, f) r(z, u, f), ...
                    F, nmax, bc, NR_OPTIONS);
            catch ME
                % Attempt restarting with softer parameters
                switch ME.identifier
                    
                    % If the correction step doesn't converge it is
                    % probably due to CA, and so factorization is forced
                    case 'NRLC:ConverganceError'
                        if solver.force_factorization ~= 1
                            solver.force_factorization = 1;
                            [P, u] = solveInner(ustart, nstart);
                        else
                            disp(['Newton failed to converge', ...
                                'desipe full factorization was used.']);
                            rethrow(ME)
                        end
                    
                    % Otherwise the error probably comes from starting too
                    % far away from equilibrium, and we need to use a new
                    % initial guess
                    case 'NRLC:FE_Error'
                        NUMBER_OF_RESTARTS = log(beta)/log(BETA_SHRINK);
                        if NUMBER_OF_RESTARTS <= NUMBER_OF_SHRINKAGES_MAX
                            % Shrink beta
                            beta = beta*BETA_SHRINK;
                            [P, u] = solveInner(ustart, nstart);
                        elseif beta > 0
                            % Set beta to zero
                            [P, u] = solveInner(ustart, nstart);
                        elseif nstart > 0
                            % Restart the iteration at a previous step
                            beta = 1;
                            nstart = nstart - 1;
                            ustart = uold(:, nstart);
                            du0 = guess(zold, z, ustart);
                            
                            [P, u] = solveInner(ustart, nstart);
                        else
                            % If all attempts to restart fails the
                            % algorithm gives up.
                            disp('Wtf dude.')
                            rethrow(ME)
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
        Ff = F(nf); uf = u(nf);
        c = Ff'*uf;
        g0 = s*c;
    end

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
        g0p = s*(drdz(z, u)'*l);
        
        % Computing the volume constraint and it's sensitivities
        g1 = volumes'*Mf*z/Vmax - 1;
        g1p = volumes'*Mf/Vmax;
        
        solver.recordUpdate(z);
        k = k + 1;
    end
obj = @cmin;
end