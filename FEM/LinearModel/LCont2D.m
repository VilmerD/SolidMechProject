classdef LCont2D < handle
    % CONTMODEL contains routines for the FEM model
    
    properties
        % Geometry
        edof;
        ex;
        ey;
        nelm;
        ndof;
        
        % Matrices
        K0_mats = {0};
        M0_mats = {0};
        
        % Material
        mpara;
        ep;
        rho;
        V0;
        ir = 3;             % Number of integration points, second order
        
        % Problem
        bc;
        F
        
        % Density filter
        w = @(x, r) max(0, 1 - x/r);
        fr = 0;
    end
    
    methods
        function obj = LCont2D(edof, ex, ey, mpara, t, rho, bc, F)
            % Geometrical data
            obj.edof = edof;
            obj.ex = ex;
            obj.ey = ey;
            obj.mpara = mpara;
            [nelm, ~] = size(obj.edof);
            obj.nelm = nelm;
             
            ndof = max(obj.edof(:));
            obj.ndof = ndof;
            
            w = max(ex(:));
            h = max(ey(:));
            obj.V0 = w*h*t;
            
            % Plane stress
            ptype = 1;          
            obj.ep = [ptype, t];
            obj.rho = rho;
            
            % Problem data
            obj.bc = bc;
            obj.F = F;
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------ Utility ------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Computes the compliance for the given design x
        function C = comp(obj, x)
            K = obj.stiffness(x);
            [a, Q] = solve(K, obj.F, obj.bc);
            C = a'*(Q + obj.F);
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------- FE-analysis ------------------------------ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

        % Computes the stiffness matrix for the given x
        function K = stiffness(obj, x)
            if obj.K0_mats{1} == 0
                obj.K0_mats = obj.K0();
            end
            
            K = zeros(obj.ndof, obj.ndof);
            for elm = 1:obj.nelm
                ke = obj.K0_mats{elm}*x(elm);
                
                dofs = obj.edof(elm, 2:end);
                K(dofs, dofs) = K(dofs, dofs) + ke;
            end
            K = sparse(K);
        end
        
        % Computes the mass matrix for the given x 
        function M = massMat(obj, x)
            if obj.M0_mats{1} == 0
                obj.M0_mats = obj.M0();
            end
            
            M = zeros(obj.ndof, obj.ndof);
            for elm = 1:obj.nelm
                Me = obj.M0_mats{elm};
                xe = x(elm);
                
                dofs = obj.edof(elm, 2:end);
                M(dofs, dofs) = M(dofs, dofs) + Me;
            end
            M = sparse(M);
        end
        
        % Computes the "derivative" of K
        function K0 = K0(obj)
            if obj.K0_mats{1} == 0
                K0 = cell(obj.nelm, 1);
                D0 = obj.dmat();

                for elm = 1:obj.nelm
                    exk = obj.ex(elm, :);
                    eyk = obj.ey(elm, :);
                    ke = planqe(exk, eyk, obj.ep, D0);

                    K0{elm} = ke;
                end
                obj.K0_mats = K0;
            else
                K0 = obj.K0_mats;
            end
        end
        
        % Computes the "derivative" of M
        function M0 = M0(obj)
            if obj.M0_mats{1} == 0
                M0 = cell(obj.nelm, 1);
                rho_t = obj.rho * obj.ep(2);

                for elm = 1:obj.nelm
                    exk = obj.ex(elm, :);
                    eyk = obj.ey(elm, :);
                    Me = planei4_m(exk, eyk, rho_t, obj.ir);

                    M0{elm} = Me;
                end
                obj.M0_mats = M0;
            else
                M0 = obj.M0_mats;
            end
        end
            
        % Computes the D-matrix representation of the constitutive tensor
        function D = dmat(obj)
            E = obj.mpara(1);
            nu = obj.mpara(2);
            D = [1      nu  0
                 nu     1   0
                 0      0   (1 - nu)/2];
            D = D*E/(1 - nu^2);
        end
        
        % Computes the areas
        function areas = areas(obj)
            areas = zeros(obj.nelm, 1);
            for k = 1:obj.nelm
                exk = obj.ex(k, :);
                eyk = obj.ey(k, :);
                dx = exk(2) - exk(1);
                dy = eyk(3) - eyk(2);
                areas(k) = dx*dy;
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------ Filter -------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Computes the matrix representation of the density filter
        function M = dfilter(obj)
            if obj.fr > 1e-6
                wr = @(x) obj.w(x, obj.fr);
                a = obj.areas();
                v = a*obj.ep(2);
                xe = [obj.ex(:, 3) + obj.ex(:, 1), ...
                    obj.ey(:, 3) + obj.ey(:, 1)]/2;

                M = zeros(obj.nelm, obj.nelm);
                for i = 1:obj.nelm
                    dxi = xe(i, :) - xe;
                    ndxi = sqrt(dxi(:, 1).^2 + dxi(:, 2).^2);
                    wi = wr(ndxi).*v;
                    wi = wi/sum(wi);
                    M(i, :) = wi;
                end
                M = sparse(M);
            else
                M = 1;
            end
        end
    end
end