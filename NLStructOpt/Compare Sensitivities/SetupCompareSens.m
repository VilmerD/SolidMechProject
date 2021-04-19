function obj = SetupCompareSens(model, solver, vq, x0, c0)
% SIMP
Emin = 210e4;    Emax = 210e9;
q = 3;
E =         @(rho) Emin + (Emax - Emin)*rho.^q;
dEdrho =    @(rho) q*(Emax - Emin)*rho.^(q - 1);

% Filtering
Mf = model.dfilter();
Mf = 1;

% Some basic model stuff
volumes = model.volumes();
Vmax = sum(volumes)*vq;

% Negative sign to force factorization at first step, as the angle between
% zold and znew is below the threshold
zold = -x0/norm(x0);     

% Model functions and stuff
K = @(rho, uk)              model.Kk(E(rho), uk);
r = @(rho, uk, fext)        model.fintk(E(rho), uk) - fext;
dcdzf = @(rho, uk, x)       model.dcdzf(dEdrho(rho), uk, x);

% Setup Displacement control
nn = 1;
dmax = -0.01;            
bc = model.bc;         
bc(:, 2) = bc(:, 2)*dmax;
u0 = zeros(model.ndof, 1);
u0(bc(:, 1)) = 0 + bc(:, 1)/nn*(nn > 1);

% Linear Solver
cothmax = 1 - 5e-3;

% Compute scaling factor
if nargin < 4
    fprintf('Computing scaling factor \n')
    [Pi, ui] = solve(x0);
    c0 = Pi'*ui;
end
s = -100/c0;
k = 0;

% Numerical sensitivities
ind = [900 901 1306 1307];
solver.statistics.del_dc = zeros(2*length(ind), 1);

    % Numerical sensitivities at some elements
    function dc = numerical_dcdzf(c, z)
        solver.force_factorization = 1;
        h = 1e-6;
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
    
    function [P, u] = solve(z)
        % Setup stiffness matrix and residual
        Kt = @(u) K(z, u);
        rt = @(u, f) r(z, u, f);
        
        % Solving problem
        [P, u, u0] = NRDispCont(solver, Kt, rt, bc, nn, u0);
    end

    function [g0, g0p, g1, g1p] = cmin(z)
        % Filtering design
        zf = Mf*z;
        
        fprintf('\nComputing displacements');
        % Checking angle between designs
        nz = norm(z);
        coth = (z/nz)'*zold;
        fprintf('. coth: %0.4f', coth)
        if coth < cothmax
            solver.force_factorization = 1;
        end
        
        % Computing compliance
        [P, u] = solve(zf);
        c = P'*u;
        g0 = c*s;
        
        fprintf('\nComputing sensitivities')
        [~, lambda] = solver.solveq(K(zf, u), -P);
        
        % The sensitivities are given by lambda*dfdzf
        dc = Mf*dcdzf(zf, u, lambda);
        g0p = s*dc;
        
        num_dc = numerical_dcdzf(c, z);
        solver.statistics.del_dc(:, k + 1) = [num_dc; dc(ind)];

        % Computing the volume constraint and it's sensitivities
        g1p = volumes'*Mf/Vmax;
        g1 = g1p*z - 1;
        
        solver.recordUpdate(zf);
        zold = z/nz;
        k = k + 1;
    end
obj = @cmin;
end