function objective = SetupDC(model, vq, alph, NSTEPS)
% Volume constraint
volumes = model.volumes;
Vmax = sum(volumes)*vq;

% Model functions and stuff
S = @(ed) model.defrom(ed);
K = @(ef, es, z) model.K(ef, es, z);
r = @(ef, es, fext, z) model.fint(ef, es, z) - fext;
dr = @(ef, es, z) model.drdE(ef, es, z);

drdzf = 0;
zold = zeros(model.nelm, 1);
zoldn = zold;
RESTARTS_MAX = 4;
uold = zeros(model.ndof, NSTEPS);

% Setup Disp. controlled scheme
% Free and prescribed nodes
bc = model.bc;
np = bc(:, 1);
bc = [np bc(:, 2)*alph];
nf = (1:model.ndof)';
nf(np) = [];

fz = zeros(model.ndof, 1);
s = 1;

% Wrapper function for solving the equilibrium equations given the design z
    function [u, F, es, ef] = solve(z)
        dz = z - zold;
        [u, F, es, ef] = NRDCoptstep(K, r, S, dz, zold, bc, ...
            uold, RESTARTS_MAX, NSTEPS);
        
        % Update uold
        NSOLS = find(normAlong(u, 2, 1) ~= 0, 1, 'first');
        uold(:, NSOLS:end) = u(:, NSOLS:end);
    end

h = 1e-5;
indx = [10 100 666];

    function dg0dz = numSensObj(z)
        dg0dz = zeros(numel(indx), 1);
        zf = filter.forward(z);
        zfp = SF.forward(zf);
        [F, u, ~, ~] = solve(zfp);
        Ffinal = F(:, end);
        ufinal = u(:, end);
        g00 = Ffinal(np)'*ufinal(np);
        
        for i = 1:numel(indx)
            ii = indx(i);
            zph = z;
            zph(ii) = zph(ii) + h;
            zphf = filter.forward(zph);
            zphfp = SF.forward(zphf);
            
            [Fph, uph, ~, ~] = solve(zphfp);
            Fphfinal = Fph(:, end);
            uphfinal = uph(:, end);
            g0ph = Fphfinal(np)'*uphfinal(np);
            
            dg0dz(i) = (g0ph - g00)/h*s;
        end
    end

    function drnum = numSensRes(ef, es, z)
        drnum = zeros(model.ndof, numel(indx));
        zf = filter.forward(z);
        zfp = SF.forward(zf);
        r0 = r(ef, es, zfp, 0);
        
        for i = 1:numel(indx)
            ii = indx(i);
            zph = z;
            zph(ii) = zph(ii) + h;
            zphf = filter.forward(zph);
            zphfp = SF.forward(zphf);
            
            rph = r(ef, es, zphfp, 0);
            drnum(:, i) = (rph - r0)/h;
        end
    end

% Computes the function value, the constraints and their derivatives
% for the current design z
    function [g0, dg0dz, g1, dg1dz, extra] = cmin(z, iter)
        % Checking angle between designs
        model.thresholding.update(iter);
        zf = model.density_filter.forward(z);
        dzfdz = model.density_filter.backward(z);
        
        % Solving equilibrium equations
        [u, F, ef, es] = solve(zf);
        Ffinal = F(:, end);
        ufinal = u(:, end);
        
        % Computing objective function
        Fp = Ffinal(np); up = ufinal(np);
        c = Fp'*up;
        if iter == 0; s = -1/c; end
        g0 = s*c;
        
        % Computing senisitivies
        Kn = K(ef, es, zf);
        l = msolveq(Kn, fz, bc);
        a([nf; np], 1) = [l(nf); up];
        drdzf = dr(ef, es, zf);
        drdz = drdzf*dzfdz;
        dcdz = (a'*drdz)';
        dg0dz = s*dcdz;
        
        % Computing the volume constraint and it's sensitivities
        g1 = volumes'*zf/Vmax - 1;
        dg1dz = volumes'*dzfdz/Vmax;
        
        % Computing change in design
        zn = z/norm(z);
        th = min(1, dot(zn, zoldn));
        alph = sqrt(1 - th^2);
        extra.alph = alph;
        
        zold = z;
        zoldn = zn;
    end
objective = @cmin;
end