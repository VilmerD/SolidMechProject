classdef NLCont2D < handle
    % init2D represents a FE Model, and contains routies for
    % - fint, stiffness matrix for total and updated lagrangian
    % - mass matrix etc
    properties
        % Setting up geometry
        nelm;
        edof;
        ndof;
        ec;
        
        % Preallocation
        I;
        J;
        esm;
        efm;
        matrices;
        
        % Material variables
        material;
        t = 1;
        
        % Element
        element
    end
    
    methods
        function obj = NLCont2D(ec, edof, ndof, t, element, material)
            obj.element = element;
            obj.edof = edof;
            s = size(edof);
            obj.nelm = s(1);
            obj.ndof = ndof;
            
            obj.t = t;
            obj.ec = ec;
            
            obj.matrices = obj.element.precompute(ec);
            obj.preAssemble();
            obj.material = material;
        end
        
        % -----------------------------------------------------------------%
        %%%%%%%%%%%%%%%%%%%%%%%%%% Utility %%%%%%%%%%%%%%%%%%%%%%%%%%
        % -----------------------------------------------------------------%
        
        % Computes the volumes of the elements
        function v = volumes(obj)
            dJ = obj.matrices(:, 1);
            v = zeros(obj.nelm, 1);
            for i = 1:obj.nelm
                v(i) = obj.element.area(dJ{i});
            end
            v = v*obj.t;
        end
        
        % Comptues the row and column vectors used to assemble the
        % stiffness matrix
        function preAssemble(obj)
            nne = (2*obj.element.npoints)^2;
            obj.I = zeros(obj.nelm*nne, 1);
            obj.J = obj.I;
            for elm = 1:obj.nelm
                k0 = ((elm - 1)*nne + 1); ke = (elm*nne);
                ij = meshgrid(obj.edof(elm, 2:end));
                obj.I(k0:ke) = reshape(ij', 1, nne);
                obj.J(k0:ke) = reshape(ij, 1, nne);
            end
        end
        
        % ----------------------------------------------------------------%
        %%%%%%%%%%%%%%%%%%%%%%%% FEM routines %%%%%%%%%%%%%%%%%%%%%%%
        % ----------------------------------------------------------------%
        
        % internal force
        function f = fint(obj, ed, k)
            % fint(ed) computes the internal forces given the displacements
            % ed
            if nargin < 3
                k = ones(obj.nelm, 1);
            end
            
            f = zeros(obj.ndof, 1);
            obj.esm = cell(obj.nelm, 1);
            obj.efm = cell(obj.nelm, 1);
            for elm = 1:obj.nelm
                dofs = obj.edof(elm, 2:end)';
                edk = ed(dofs);
                
                % Computes the defgrad, if updated lagrangian use old ef
                ef = obj.element.defgrad(obj.matrices{elm, 2}, edk);
                
                % Computes stresses
                es = obj.material.es(ef);
                obj.efm{elm} = ef;
                obj.esm{elm} = es;
                
                % Using the above data the element stiffness matrix can be
                % computed and inserted into the correct position
                felm = k(elm)*obj.element.force(obj.matrices{elm, 2}, ...
                    obj.matrices{elm, 3}, obj.matrices{elm, 1}, obj.t, ef, es);
                f(dofs) = f(dofs) + felm;
            end
        end
        
        % stiffness
        function K = K(obj, k)
            % Stiffness() computes the stiffness matrix of the simple system
            % given the displacements a and force P
            nne = (2*obj.element.npoints)^2;
            X = zeros(obj.nelm*nne, 1);
            
            if nargin < 2
                k = ones(obj.nelm, 1);
            end
            
            % Iterates over the elements
            for elm = 1:obj.nelm
                % Computes stresses
                ef = obj.efm{elm};
                es = obj.esm{elm};
                
                % dmat
                D = obj.material.dm(ef);
                
                
                % Using the above data the element stiffness matrix can be
                % computed
                Kelm = k(elm)*obj.element.stiffness(obj.matrices{elm, 2}, ...
                    obj.matrices{elm, 3}, obj.matrices{elm, 1}, obj.t, D, ef, es);
                
                % Inserting the element stiffness matrix into the correct pos
                k0 = ((elm - 1)*nne + 1); ke = (elm*nne);
                X(k0:ke, 1) = Kelm(:);
            end
            K = sparse(obj.I, obj.J, X);
        end
        
        % computes the von-stress for the plane strain problem
        function S = stresses(obj, ed)
            S = zeros(obj.nelm, 1);
            
            for elm = 1:obj.nelm
                dofs = obj.edof(elm, 2:end);
                edk = ed(dofs);
                
                ef = obj.element.defgrad(obj.matrices{elm, 2}, edk);
                
                % Computes stresses
                es = obj.material.es(ef);
                seff = vonMises(es, obj.material.mpara);
                S(elm) = seff;
            end
            
        end
        
        % Computes the sensitivity of the residual
        function y = drdE(obj)
            nne = (obj.element.npoints*2);
            It = zeros(nne*obj.nelm, 1);
            Jt = It;
            X = It;
            for elm = 1:obj.nelm
                dofs = obj.edof(elm, 2:end)';
                
                % Computes the defgrad
                ef = obj.efm{elm};
                es = obj.esm{elm};
                
                % Using the above data the element forces can be computed
                felm = obj.element.force(obj.matrices{elm, 2}, ...
                    obj.matrices{elm, 3}, obj.matrices{elm, 1}, obj.t, ef, es);
                
                k0 = (elm - 1)*nne + 1; ke = elm*nne;
                It(k0:ke, 1) = dofs;
                Jt(k0:ke, 1) = elm;
                X(k0:ke, 1) = felm;
            end
            y = sparse(It, Jt, X);
        end
    end
end