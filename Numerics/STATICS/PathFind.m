classdef PathFind < handle
    properties  
        % Stiffness matrix (function) and residual (also function)
        Kt;
        r;
        
        % Problem parameters
        bc;
        p;
        nmax;
        alpha;
        
        % Solution vectors
        u;
        l = [0];
        P;
        
        % Other variables used in the confined method
        c;
        psi;
        d0;
        rtol;
        
        % Which method to use, regular NR or constrained
        method;
        
        % Either fully linearized or semilinearized method to compute
        % D_lambda
        linearization;
    end
    
    methods
        % Initialtes the object
        function obj = PathFind(Kt, r, u0, bc, load_increment, method)
            obj.Kt = Kt;
            obj.r = r;
            obj.bc = bc;
            obj.p = load_increment;
            
            obj.u = [u0];
            obj.psi = 1;
            
            obj.compute_alpha(bc(:, 1), length(u0));
            
            if nargin < 6 
                obj.method = @obj.NR;
            elseif method == "bordering"
                obj.linearization = @obj.bordering;
                obj.method = @obj.ConNR;
            elseif method == "constrain"
                obj.linearization = @obj.constrain;
                obj.method = @obj.ConNR;                
            end
        end
        
        % Sets the load increment
        function set_p(obj, pnew)
            obj.p = pnew;
        end
        
        % Finds the complementay nodoes to the boundary conditions,
        % ie where the residual should be checked
        function compute_alpha(obj, col, l)
           m = [];
           for k = 1:l
               exists = col == k;
               if ~sum(exists)
                   m = [m k];
               end
           end
           obj.alpha = m;
        end
        
        % A constrained NR method
        function ConNR(obj)           
            ctol = 1e-3;
            
            n = 0;
            while n < obj.nmax
                n = n + 1;

                % Displacement u
                un = obj.u(:, end);
                du_tot = zeros(size(un));

                % Load pattern lambda (l)
                ln = obj.l(end);
                dl_tot = 0;

                r0 = obj.r(un, ln * obj.p);        % To prevent taking too many steps
                
                % To ensure at least one step is taken, we check
                % if dl_tot is zero. 
                while dl_tot == 0 || (norm(r0(obj.alpha)) > obj.rtol || norm(obj.c) > ctol)
                    % Computing new stiffness matrix
                    Ktn = obj.Kt(un);

                    % Computing pseudo displacement increments
                    [du_r, ~] = solve(Ktn, -r0, obj.bc);
                    [du_p, ~] = solve(Ktn, obj.p, obj.bc);

                    % Computing load parameter
                    dl = obj.linearization(du_r, du_p, du_tot, dl_tot, n);
                    
                    % Updating increments
                    ln = ln + dl;
                    dl_tot = ln - obj.l(end);

                    un = un + du_r + dl*du_p;
                    du_tot = un - obj.u(:, end);

                    r0 = obj.r(un, ln * obj.p);
                    obj.c = norm(du_tot) ^ 2 + obj.psi * dl_tot ^ 2 *...
                        norm(obj.p) ^ 2 - obj.d0^2 ;
                end  
                obj.u = [obj.u un];
                obj.l = [obj.l ln];
                
                disp(n)
            end
        end
        
        % Sets a new d0 if the step size is too small
        function setd0(obj, dnew)
            obj.d0 = dnew;
        end
        
        % Computes dl if we're using a fully linear model, ie both c and r
        % are linearized
        function dl = bordering(obj, du_r, du_p, du_tot, dl_tot, n)
            if dl_tot == 0
                % If it's the absolute first step s = 1
                if n == 1
                    s = 1;
                else
                    du_previous = obj.u(:, end) - obj.u(:, end-1);
                    dl_previous = obj.l(end) - obj.l(end - 1);
                    s = sign(du_previous' * du_p + dl_previous);
                end
                dl = s *  obj.d0 / (sqrt(norm(du_p) ^ 2 + obj.psi * norm(obj.p) ^ 2));
            else
                dl = - (obj.c + 2*du_tot'*du_r) / ...
                    (2*du_tot'*du_p + 2 * obj.psi * norm(obj.p) ^ 2 * dl_tot);
            end            
        end
        
        % Computes dl if we're using a semi-linear model, ie c is solved
        % exactly
        function dl = constrain(obj, du_r, du_p, du_tot, dl_tot, n)
            % Special case when we're at the first iteration in this step
            if dl_tot == 0
                % If it's the absolute first step s = 1
                if n == 1
                    s = 1;
                else
                    du_previous = obj.u(:, end) - obj.u(:, end-1);
                    dl_previous = obj.l(end) - obj.l(end - 1);
                    s = sign(du_previous' * du_p + dl_previous);
                end
                    dl = s * obj.d0 / (sqrt(norm(du_p) ^ 2 + obj.psi * norm(obj.p) ^ 2));
            % Otherwise compute the a-factors and continue
            else
                a1 = norm(du_p) ^ 2 + obj.psi * norm(obj.p) ^ 2;
                a2 = 2*du_p' * (du_tot + du_r) + 2 * obj.psi * dl_tot * norm(obj.p) ^ 2;
                a3 = norm(du_tot + du_r) ^ 2 + obj.psi * dl_tot ^ 2 * norm(obj.p) ^ 2 - obj.d0 ^ 2;
                dlv = roots([a1 a2 a3]);
                
                % If we have complex solutions we need to decrease the
                % length
                if ~isreal(dlv)
                    d = obj.d0;
                    obj.setd0(obj.d0*2);
                    dl = obj.constrain(du_r, du_p, du_tot, dl_tot, n);
                    obj.setd0(d);
                
                % If we have two distinct solutions, the one which
                % maximizes the cosine of the angle is chosen
                elseif dlv(1) ~= dlv(2)
                    a4 = du_tot' * (du_tot + du_r);
                    a5 = du_tot' * du_p;

                    co = (a4 + a5 * dlv)/obj.d0^2;
                    dl = dlv(co == (max(co)));
                else
                    dl = dlv(1);
                end
            end
                
        end
        
        % Load controlled solver, not used in the project
        function NR(obj)
        rtol = 1e-6;
        utol = 1e-6;
        dP = obj.p/obj.nmax;

        obj.l = zeros(length(dP), obj.nmax);

        for n = 2:obj.nmax
            Pn = obj.l(:, n-1) + dP;
            obj.l(:, n) = Pn;

            un = obj.u(:, n - 1);
            rc = obj.r(un, Pn);

            first_it = 1;
            while first_it || (norm(rc(obj.alpha)) > rtol ...
                    && norm(du(obj.alpha)) > utol)
                K = obj.Kt(un);

                [du, ~] = solveq(K, -rc, obj.bc);
                un = un + du;

                rc = obj.r(un, Pn);

                first_it = 0;
            end
            obj.u(:, n) = un;
        end
        end
        
        % Solves the problem
        function [P, u] = compute(obj, rtol, n, d0, psi)
            obj.u = [obj.u(:, 1)];
            obj.l = [0];
            obj.psi = psi;
            obj.nmax = n;
            obj.rtol = rtol;
            obj.d0 = d0;
            obj.method();
            P = obj.l;
            u = obj.u;
        end
    end
end