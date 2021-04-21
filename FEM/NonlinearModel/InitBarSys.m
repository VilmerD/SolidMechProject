classdef InitBarSys < handle
    properties
        % Setting up geometry        
        Nelm;
        Edof;
        Ndof;      
        ec;                 % In column format, (Ndof x 1)
        
        % Material variables
        ep;
        stress;             % Computes stress given E and stretch
        dmat;               % Computes dmat given E and stretch
        
        % Spring variables, not used in the project
        Nspring;
        Springdof;
        k;
        
        % Solver details
        bc;
        Pm;                 % Load increment
        solver;             % The path finding object
                
    end

    methods
        function obj = InitBarSys(ec, Edof, ndof, bc, Pm, ep, materialmodel) 
            obj.Edof = Edof;
            s = size(Edof);
            obj.Nelm = s(1);
            obj.Ndof = ndof;
            
            obj.ep = ep;
            obj.bc = bc;
            obj.Pm = Pm;
            obj.ec = ec;
            
            obj.setMaterialModel(materialmodel);
        end        
        
        
        % Sets the material model to either linear or nonlinear
        function setMaterialModel(obj, model)            
            if model == "linear"
                obj.stress = @(E, stretch) E/2 * (stretch ^ 2 - 1);
                obj.dmat = @(E, stretch) E;
            else
                obj.stress = @stress1D;
                obj.dmat = @dmat1D;
            end
        end
        
        % Sets the reference configuration, used when a perturbation is
        % introduced
        function setEc(obj, ec)
            obj.ec = ec;
        end
        
        % Sets the load increment, used when a perturbation is
        % introduced 
        function setP(obj, Pnew)
            obj.Pm = Pnew;
        end
        
        % Adds springs to the model, only used in the exercise
        function addSprings(obj, Springdof, k)
           obj.Springdof = Springdof;
           obj.k = k;
           obj.Nspring = length(k);
        end
        
        
        % Computes the path given a method, see the PathFind class
        function [P, u] = computePath(obj, method)
            if nargin > 1
                obj.solver = PathFind(@obj.Stiffness, @obj.response, ...
                    zeros(obj.Ndof, 1), obj.bc, obj.Pm, method);
            else
                obj.solver = PathFind(@obj.Stiffness, @obj.response, ...
                    zeros(obj.Ndof, 1), obj.bc, obj.Pm);
            end
            [P, u] = obj.solver.compute();
        end
        
        function resp = response(obj, ed, P)
        % response(a, P) computes the response to the simple system given
        % the displacements a and the force P
        fext = P;

        fint = zeros(obj.Ndof, 1);

        for elm = 1:obj.Nelm
            % Finds the displacements and reference configuration for
            % this element
            disps = obj.Edof(elm, 2:end)';            
            ed_curr = ed(disps);
            ec_curr = reshape(obj.ec(disps), 3, 2);
            
            % Computes the stretch and stress es
            dx0 = ec_curr(:, 1) - ec_curr(:, 2);
            du = ed_curr(1:3) - ed_curr(4:6);
            es = obj.stress(obj.ep(1), stretch1D(du, dx0));

            % Using the above data the element stiffness matrix can be
            % computed and inserted into the correct position
            felm = bar3gf(ec_curr, ed_curr, es);
            fint(disps) = fint(disps) + felm;
        end
        
        % Computes spring forces, not used in this project, but used in the
        % exercise.
        for spring = 1:obj.Nspring
            ks = obj.k(spring);
            disps = obj.Springdof(spring, 2:end)';
            
            ed_curr = ed(disps);

            fspring = bar1f(ks, ed_curr);

            fint(disps) = fint(disps) + fspring;
        end

        resp = fint - fext;
        end

        function Kt = Stiffness(obj, ed)
        % Stiffness(a, P) computes the stiffness matrix of the simple system
        % given the displacements a and force P
        Kt = zeros(obj.Ndof);

        % Iterates over the elements
        for elm = 1:obj.Nelm
            % Finds the displacements and reference configuration for
            % this element
            disps = obj.Edof(elm, 2:end)';            
            ed_curr = ed(disps);          
            ec_curr = reshape(obj.ec(disps), 3, 2);
            
            % Computes the stretch, material tangent Et and stress es
            dx0 = ec_curr(:, 1) - ec_curr(:, 2);
            du = ed_curr(1:3) - ed_curr(4:6);
            stretch = stretch1D(du, dx0);
            Et = obj.dmat(obj.ep(1), stretch);
            es = obj.stress(obj.ep(1), stretch);    
            
            % Using the above data the element stiffness matrix can be
            % computed
            Kelm = bar3ge(ec_curr, ed_curr, [Et, obj.ep(2)], es);
            
            % Inserting the element stiffness matrix into the correct pos
            disp_column = disps(:);
            Kt(disp_column, disp_column) = Kt(disp_column, disp_column)...
                + Kelm;
        end
        
        % Computes spring stiffness matrices, not used in this project
        % but used in the exercise.
        for spring = 1:obj.Nspring
            disps = obj.Springdof(spring, 2:end)';
            
            ks = obj.k(spring);
            Kspring = bar1e(ks);

            Kt(disps, disps) = Kt(disps, disps)...
                + Kspring;
        end        
        end            
    end
end