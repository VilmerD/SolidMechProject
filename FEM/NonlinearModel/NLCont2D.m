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
        ex;
        ey;
        
        % Problem parameters
        bc;
        f;
        
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
        function obj = NLCont2D(ex, ey, edof, ndof, t, element, material, ...
                bc, f)
            obj.element = element;
            obj.edof = edof;
            s = size(edof);
            obj.nelm = s(1);
            obj.ndof = ndof;
            
            obj.t = t;
            obj.ec = [ex, ey];
            obj.ex = ex;
            obj.ey = ey;
            
            obj.matrices = obj.element.precompute(obj.ec);
            obj.preAssemble();
            obj.material = material;
            
            obj.bc = bc;
            obj.f = f;
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
            obj.I = reshape(kron(obj.edof(:, 2:end), ...
                ones(2*obj.element.npoints, 1))', [], 1);
            obj.J = reshape(kron(obj.edof(:, 2:end), ...
                ones(1, 2*obj.element.npoints))', [], 1);
        end
        
        % ----------------------------------------------------------------%
        %%%%%%%%%%%%%%%%%%%%%%%% FEM routines %%%%%%%%%%%%%%%%%%%%%%%
        % ----------------------------------------------------------------%
        
        % stresses and strains
        function [ef, es] = defrom(obj, ed)
            ef = cell(obj.nelm, 1);
            es = cell(obj.nelm, 1);
            edelms = ed(obj.edof(:, 2:end))';
            for elm = 1:obj.nelm
                edelm = edelms(:, elm);
                ef{elm} = obj.element.defgrad(obj.matrices{elm, 2}, edelm);
                es{elm} = obj.material.es(ef{elm});
            end
        end
        
        % internal force
        function f = fint(obj, ef, es, x)
            
            f = zeros(obj.ndof, 1);
            for elm = 1:obj.nelm
                % Extract stresses
                efelm = ef{elm};
                eselm = es{elm};
                
                % Using the above data the element stiffness matrix can be
                % computed and inserted into the correct position
                felm = x(elm)*obj.element.force(obj.matrices{elm, 2}, ...
                    obj.matrices{elm, 3}, obj.matrices{elm, 1}, obj.t, ...
                    efelm, eselm);
                % TODO: Fix this horrible indexing to sparse format
                dofs = obj.edof(elm, 2:end)';
                f(dofs) = f(dofs) + felm;
            end
        end
        
        % stiffness
        function K = K(obj, ef, es, x)
            nne = (2*obj.element.npoints)^2;
            X = zeros(obj.nelm*nne, 1);
            
            % Iterates over the elements
            for elm = 1:obj.nelm
                % Extract stresses
                efelm = ef{elm};
                eselm = es{elm};
                
                % dmat
                Delm = obj.material.dm(efelm);
                
                % Using the above data the element stiffness matrix can be
                % computed
                Kelm = x(elm)*obj.element.stiffness(obj.matrices{elm, 2}, ...
                    obj.matrices{elm, 3}, obj.matrices{elm, 1}, obj.t, ...
                    Delm, efelm, eselm);
                
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
        function y = drdE(obj, ef, es, x)
            nne = (obj.element.npoints*2);
            It = reshape(obj.edof(:, 2:end)', [], 1);
            Jt = kron((1:obj.nelm)', ones(2*obj.element.npoints, 1));
            X = zeros(size(It));
            for elm = 1:obj.nelm
                % Computes the defgrad
                efelm = ef{elm};
                eselm = es{elm};
                
                % Using the above data the element forces can be computed
                felm = obj.element.force(obj.matrices{elm, 2}, ...
                    obj.matrices{elm, 3}, obj.matrices{elm, 1}, obj.t, ...
                    efelm, eselm);
                
                k0 = (elm - 1)*nne + 1; ke = k0 + nne - 1;
                X(k0:ke, 1) = x(elm)*felm;
            end
            y = sparse(It, Jt, X);
        end
    end
end