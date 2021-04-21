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
        bc;
        
        % Material variables
        stress;             % Computes stress given E and stretch
        dmat;               % Computes dmat given E and stretch
        mpara;              % E and nu
        t = 1;              % Thickness, used in the plain strain case
        
        % Functions to compute
        defgrad;            % deformation gradient
        force;              % internal force
        stiffness;          % element stiffness
        mass;
        
        fold;
        
        % Number of points per element
        npoints;
        
        % Numerical integration or exact solution?
        numint = 0;         % Nr of numerical integration points if any
        model;
        lagflag;
        stressflag;         % total lagrangian = 1, updated = 3;
        typeflag;           % plane = 1, axi = 2;
        
        % Contact variables
        edofc;              % edof for the contacts
        ecc;                % ec for the contacts
        nelmc;              % number of contact elements
        bforce;
        bstiffness;
        benergy;
        k; 
        r;
        
        % Structural optimization
        w = @(x, r) max(0, 1 - x/r);
        fr;
    end
    
    methods
        function obj = NLCont2D(ec, edof, ndof, mpara, t, eltype, bc, matmod)
            obj.edof = edof;
            s = size(edof);
            obj.nelm = s(1);
            obj.ndof = ndof;
            
            obj.mpara = mpara;
            obj.t = t;
            obj.ec = ec;
            obj.bc = bc;
            
            obj.setElementData(eltype);
            
            if nargin < 8
                matmod = 1;
            end
            obj.setMaterialModel(matmod)
        end
        
        % -----------------------------------------------------------------%
        %%%%%%%%%%%%%%%%%%%%%%%%%% Utility %%%%%%%%%%%%%%%%%%%%%%%%%%
        % -----------------------------------------------------------------%
        % Sets the material model to either linear or nonlinear
        function setMaterialModel(obj, model)
            % which type of stress should be computed?
            stresst = obj.stressflag;
            lag = obj.lagflag;
            m = obj.mpara;
            if model == 1
                if obj.typeflag == 1
                    obj.stress = @(f) stressMater2D1(stresst, m, f);
                    obj.dmat = @(f) dMater2D1(lag, m, f);
                else
                    % Placeholder for implementing material 1
                    obj.setMaterialModel(2)
                    warning("Axisymmetry with S:t Venant stress " + ...
                        "not implemented. Neo-Hookean model set")
                end
            elseif model == 2
                if obj.typeflag == 1
                    obj.stress = @(f) stressMater2D2(stresst, m, f);
                    obj.dmat = @(f) dMater2D2(lag, m, f);
                else
                    obj.stress = @(f) stressMaterAxi2(stresst, m, f);
                    obj.dmat = @(f) dMaterAxi2(lag, m, f);
                end
            end
            obj.model = model;
        end
        
        % Sets the reference configuration, used when a perturbation is
        % introduced
        function setEc(obj, ec)
            obj.ec = ec;
        end
        
        % Computes the volumes of the elements
        function v = volumes(obj)
            areas = zeros(obj.nelm, 1);
            for i = 1:obj.nelm
                exk = obj.ec(i, 1:obj.npoints);
                eyk = obj.ec(i, (obj.npoints + 1):end);
                dx = exk(2) - exk(1);
                dy = eyk(3) - eyk(2);
                areas(i) = dx*dy;
            end
            v = areas*obj.t;
        end
        
        % Fetches the correct function calls for you :)
        function [K, f] = getFuncs(obj, method)
            if strcmpi("total", method)
                K = @obj.K;
                f = @obj.fint;
            elseif strcmpi("updated", method)
                K = @obj.Ku;
                f = @obj.ru;
            elseif strcmpi("energy", method)
                K = @obj.KE;
                f = @obj.fintstar;
            elseif strcmpi("contact", method)
                K = @obj.Kc;
                f = @obj.fintc;
            elseif strcmpi("contact_energy", method)
                K = @obj.KEc;
                f = @obj.fintstarc;
            end
        end
        
        % Gets the element data
        function setElementData(obj, eltype)
            % 3 node triangle total lagrangian
            if eltype == "2D3t"
                obj.defgrad = @cont2D3ts;
                obj.force = @(ec, ed, es) cont2D3tf(ec, obj.t, ed, es);
                obj.stiffness = @(ec, D, ed, es) ...
                    cont2D3te(ec, obj.t, D, ed, es);
                obj.mass = @cont2D3m;
                
                obj.npoints = 3;
                obj.numint = 0;
                
                obj.stressflag = 1;
                obj.typeflag = 1;
                obj.lagflag = 1;
                
                % 3 node triangle updated lagrangian
            elseif eltype == "2D3u"
                obj.defgrad = @cont2D3us;
                obj.force = @(ec, ed, es) cont2D3uf(ec, obj.t, ed, es);
                obj.stiffness = @(ec, D, ed, es) ...
                    cont2D3ue(ec, obj.t, D, ed, es);
                
                obj.npoints = 3;
                obj.numint = 0;
                
                obj.stressflag = 3;
                obj.typeflag = 1;
                obj.lagflag = 2;
                
                obj.fold = cell(obj.nelm, 1);
                for j = 1:obj.nelm
                    obj.fold{j, 1} = [1 0 0 1]';
                end
                
                % Axisymmetry
            elseif eltype == "Axi3"
                obj.defgrad = @contAxi3ts;
                obj.force = @(ec, ed, es) contAxi3tf(ec, obj.t, ed, es);
                obj.stiffness = @(ec, D, ed, es) ...
                    contAxi3te(ec, obj.t, D, ed, es);
                
                obj.npoints = 3;
                obj.numint = 0;
                
                obj.stressflag = 1;
                obj.typeflag = 2;
                obj.lagflag = 1;
                
                % 4 node quadrilateral total lagrangian
            elseif eltype == "2D4t"
                obj.defgrad = @cont2D4ts;
                obj.force = @(ec, ed, es) cont2D4tf(ec, obj.t, ed, es);
                obj.stiffness = @(ec, D, ed, es) ...
                    cont2D4te(ec, obj.t, D, ed, es);
                
                obj.npoints = 4;
                obj.numint = 4;
                
                obj.stressflag = 1;
                obj.typeflag = 1;
                obj.lagflag = 1;
                
                % 6 node triangle updated lagrangian
            elseif eltype == "2D6u"
                obj.defgrad = @cont2D6us;
                obj.force = @(ec, ed, es) cont2D6uf(ec, obj.t, ed, es);
                obj.stiffness = @(ec, D, ed, es) ...
                    cont2D6ue(ec, obj.t, D, ed, es);
                
                obj.npoints = 6;
                obj.numint = 3;
                
                obj.stressflag = 3;
                obj.typeflag = 1;
                obj.lagflag = 2;
                
                % Setting up fold
                obj.fold = cell(obj.nelm, 1);
                for i = 1:obj.nelm
                    foldk = cell(3, 1);
                    for j = 1:3
                        foldk{j, 1} = [1 0 0 1]';
                    end
                    obj.fold{i} = foldk;
                end
                
            else
                error('foo:bar', "Unrecognized element type, " + ...
                    "possible options: \n" + ...
                    "2D3t: Total lagrangian 3-node triangle element. \n" + ...
                    "2D3u: Updated lagrangian 3-node triangle element. \n" + ...
                    "Axi3: Total lagrangian axisymmetric 3-node " + ...
                    "triangle element. \n" + ...
                    "2D4t: Total lagrangian quadrilateral element. \n" + ...
                    "2D6u: Updated lagrangian 6-node triangle element. \n")
                
            end
        end
        
        % Extracts the coordinates, deformations and dofs for the element
        function [eck, edk, dofs] = extract(obj, type, elm, ed)
            % Finds the displacements and reference configuration for
            % this element
            if type == "element"
                Ec = obj.ec;
                Edof = obj.edof;
            elseif type == "contact"
                Ec = obj.ecc;
                Edof = obj.edofc;
            end
            dofs = Edof(elm, 2:end)';
            edk = ed(dofs);
            eck = reshape(Ec(elm, :), obj.npoints, 2)';
        end
        % ----------------------------------------------------------------%
        %%%%%%%%%%%%%%%%%% Total Lagrangian setting %%%%%%%%%%%%%%%%%%
        % ----------------------------------------------------------------%
        % internal force
        function fint = fint(obj, ed)
            % fint(ed) computes the internal forces given the displacements
            % ed
            fint = zeros(obj.ndof, 1);
            for elm = 1:obj.nelm
                [eck, edk, dofs] = obj.extract("element", elm, ed);
                
                % Computes the defgrad, if updated lagrangian use old ef
                ef = obj.defgrad(eck, edk);
                
                % Computes stresses
                if obj.numint ~= 0
                    es = cell(obj.numint, 1);
                    for i = 1:obj.numint
                        es{i} = obj.stress(ef{i});
                    end
                else
                    es = obj.stress(ef);
                end
                
                % Using the above data the element stiffness matrix can be
                % computed and inserted into the correct position
                felm = obj.force(eck, edk, es);
                fint(dofs) = fint(dofs) + felm;
            end
        end
        
        % stiffness
        function Kt = K(obj, ed)
            % Stiffness(a, ) computes the stiffness matrix of the simple system
            % given the displacements a and force P
            Kt = zeros(obj.ndof);
            
            % Iterates over the elements
            for elm = 1:obj.nelm
                [eck, edk, dofs] = obj.extract("element", elm, ed);
                
                % Computes the defgrad, if updated lagrangian use old ef
                ef = obj.defgrad(eck, edk);
                
                % Computes stresses and dmat
                if obj.numint ~= 0
                    es = cell(obj.numint, 1);
                    D = cell(obj.numint, 1);
                    for i = 1:obj.numint
                        es{i} = obj.stress(ef{i});
                        D{i} = obj.dmat(ef{i});
                    end
                else
                    es = obj.stress(ef);
                    D = obj.dmat(ef);
                end
                
                % Using the above data the element stiffness matrix can be
                % computed
                Kelm = obj.stiffness(eck, D, edk, es);
                
                % Inserting the element stiffness matrix into the correct pos
                disp_column = dofs(:);
                Kt(disp_column, disp_column) = Kt(disp_column, disp_column)...
                    + Kelm;
            end
            Kt = sparse(Kt);
            
        end
        
        % computes the von-stress for the plain strain problem
        function S = stresses(obj, ed)
            S = zeros(obj.nelm, 1);
            
            for elm = 1:obj.nelm
                [eck, edk, ~] = obj.extract("element", elm, ed);
                
                % Computes the defgrad, if updated lagrangian use old ef
                ef = obj.defgrad(eck, edk);
                
                % Computes stresses
                es = obj.stress(ef);
                seff = vonMises(es, obj.mpara);
                S(elm) = seff;
            end
            
        end
        % ----------------------------------------------------------------%
        %%%%%%%%%%%%%%%%%%%%% Dynamical analysis %%%%%%%%%%%%%%%%%%%%%
        % ----------------------------------------------------------------%
        % mass matrix
        function M = M(obj, rho)
            M = zeros(obj.ndof);
            for elm = 1:obj.nelm
                dofs = obj.edof(elm, 2:end)';
                eck = reshape(obj.ec(elm, :), obj.npoints, 2)';
                Me = obj.mass(eck, obj.t, rho);
                
                dof_col = dofs(:);
                M(dof_col, dof_col) = M(dof_col, dof_col) + Me;
            end
            M = sparse(M);
        end
        
        % stiffness matrix in energy conserving analysis. Only supports
        % 3- node triangular elements in total lagrangian setting, with
        % st:venant material model
        function Kt = KE(obj, ed, edprev)
            % given the displacements a and force P
            Kt = zeros(obj.ndof);
            
            % Iterates over the elements
            for elm = 1:obj.nelm
                % Finds the displacements and reference configuration for
                % this element
                dofs = obj.edof(elm, 2:end)';
                edk = ed(dofs);
                edprevk = edprev(dofs);
                eck = reshape(obj.ec(elm, :), obj.npoints, 2)';
                
                % Computes the defgrad
                ef = obj.defgrad(eck, edk);
                efprev = obj.defgrad(eck, edprevk);
                es = obj.stress(ef);
                esprev = obj.stress(efprev);
                
                esbar = (es + esprev)/2;
                D = obj.dmat((ef + efprev)/2);
                
                % Using the above data the element stiffness matrix can be
                % computed
                Kelm = cont2D3Ee(eck, obj.t, D, edk, edprevk, esbar);
                
                % Inserting the element stiffness matrix into the correct pos
                Kt(dofs, dofs) = Kt(dofs, dofs) + Kelm;
            end
            Kt = sparse(Kt);
        end
        
        % internal force* in energy conserving dynamical analysis
        function fint = fintstar(obj, ed, edprev)
            fint = zeros(obj.ndof, 1);
            
            for elm = 1:obj.nelm
                % Finds the displacements and reference configuration for
                % this element
                dofs = obj.edof(elm, 2:end)';
                edk = ed(dofs);
                edprevk = edprev(dofs);
                edbar = (edk + edprevk)/2;
                eck = reshape(obj.ec(elm, :), obj.npoints, 2)';
                
                % Computes the defgrad, if updated lagrangian use old ef
                ef = obj.defgrad(eck, edk);
                efprev = obj.defgrad(eck, edprevk);
                es = obj.stress(ef);
                esprev = obj.stress(efprev);
                
                esbar = (es + esprev)/2;
                
                % Using the above data the element stiffness matrix can be
                % computed and inserted into the correct position
                felm = obj.force(eck, edbar, esbar);
                fint(dofs) = fint(dofs) + felm;
            end
        end
        
        % computes the energy for the configuration
        function [KinE, IntE] = energy(obj, rho, ed, du)
            KinE = 0;
            IntE = 0;
            for elm = 1:obj.nelm
                dofs = obj.edof(elm, 2:end)';
                edk = ed(dofs);
                eck = reshape(obj.ec(elm, :), obj.npoints, 2)';
                duk = du(dofs);
                
                % Computes the defgrad
                ef = obj.defgrad(eck, edk);
                
                [KinEe, IntEe] = cont2D3En(obj.model, obj.mpara, eck, ...
                    obj.t, rho, ef, duk);
                KinE = KinE + KinEe;
                IntE = IntE + IntEe;
            end
        end
        % ----------------------------------------------------------------%
        %%%%%%%%%%%%%%%%%%%%%% Contact analysis %%%%%%%%%%%%%%%%%%%%%%
        % ----------------------------------------------------------------%
        % adds contact bars connecting a point with elements
        function addContact(obj, ecc, edofc, nelmc, k, r)
            obj.edofc = edofc;
            obj.ecc = ecc;
            obj.nelmc = nelmc;
            obj.bforce = @(ec, ed) norfb(ec, ed, k, r);
            obj.bstiffness = @(ec, ed) bstiff(ec, ed, k, r);
            obj.benergy = @(ec, ed) bar3en(ec, ed, k, r);
            obj.k = k;
            obj.r = r;
        end
        
        % Stiffness due to contact
        function Kt = Kc(obj, ed)
            Kt = zeros(obj.ndof);
            
            for elm = 1:obj.nelmc
                dofs = obj.edofc(elm, 2:end)';
                edk = ed(dofs);
                eck = obj.ecc{elm};
                
                es = obj.bforce(eck, edk);
                if norm(es) ~= 0
                    Et = obj.bstiffness(eck, edk);
                    Ke = bar3ge(eck, edk, [Et, 1], es);
                    
                    dofcol = dofs(:);
                    Kt(dofcol, dofcol) = Kt(dofcol, dofcol) + Ke;
                end
            end
            Kt = sparse(Kt);
        end
        
        % internal force due to contact
        function fint = fintc(obj, ed)
            fint = zeros(obj.ndof, 1);
            
            for elm = 1:obj.nelmc
                dofs = obj.edofc(elm, 2:end)';
                edk = ed(dofs);
                eck = obj.ecc{elm};
                
                es = obj.bforce(eck, edk);
                
                if norm(es) ~= 0
                    fe = bar3gf(eck, edk, es);
                    dofcol = dofs(:);
                    fint(dofcol) = fint(dofcol) + fe;
                end
            end
        end
        
        % Energy in the bars
        function IntE = energyc(obj, ed)
            IntE = 0;
            
            for elm = 1:obj.nelmc
                dofs = obj.edofc(elm, 2:end)';
                edk = ed(dofs);
                eck = obj.ecc{elm};
                
                IntE = IntE + obj.benergy(eck, edk);
            end
        end
        % ----------------------------------------------------------------%
        %%%%%%%%%%%%%%%%%% Structural optimization %%%%%%%%%%%%%%%%%%
        % ----------------------------------------------------------------%
        % internal force
        function fint = fintk(obj, k, ed)
            % fint(k, ed) computes the internal forces given the displacements
            % ed and premultiplies each element force vector with the
            % constant(s) in k
            fint = zeros(obj.ndof, 1);
            for elm = 1:obj.nelm
                [eck, edk, dofs] = obj.extract("element", elm, ed);
                
                % Computes the defgrad, if updated lagrangian use old ef
                ef = obj.defgrad(eck, edk);
                
                % Computes stresses
                if obj.numint ~= 0
                    es = cell(obj.numint, 1);
                    for i = 1:obj.numint
                        es{i} = obj.stress(ef{i});
                    end
                else
                    es = obj.stress(ef);
                end
                
                % Using the above data the element stiffness matrix can be
                % computed and inserted into the correct position
                felm = k(elm)*obj.force(eck, edk, es);
                fint(dofs) = fint(dofs) + felm;
            end
        end
        
        function Kt = Kk(obj, k, ed)
            % Kk(k, ed) computes the stiffness matrix given the
            % displacements ed and premultiplies each element stiffness
            % with the constant(s) in k
            Kt = zeros(obj.ndof);
            
            % Iterates over the elements
            for elm = 1:obj.nelm
                [eck, edk, dofs] = obj.extract("element", elm, ed);
                
                % Computes the defgrad, if updated lagrangian use old ef
                ef = obj.defgrad(eck, edk);
                
                % Computes stresses and dmat
                if obj.numint ~= 0
                    es = cell(obj.numint, 1);
                    D = cell(obj.numint, 1);
                    for i = 1:obj.numint
                        es{i} = obj.stress(ef{i});
                        D{i} = obj.dmat(ef{i});
                    end
                else
                    es = obj.stress(ef);
                    D = obj.dmat(ef);
                end
                
                % Using the above data the element stiffness matrix can be
                % computed
                Kelm = k(elm)*obj.stiffness(eck, D, edk, es);
                
                % Inserting the element stiffness matrix into the correct pos
                disp_column = dofs(:);
                Kt(disp_column, disp_column) = Kt(disp_column, disp_column)...
                    + Kelm;
            end
            Kt = sparse(Kt);
        end
        
        function y = dcdzf(obj, k, ed, x)
            y = zeros(obj.nelm, 1);
            for elm = 1:obj.nelm
                [eck, edk, dofs] = obj.extract("element", elm, ed);
                
                % Computes the defgrad
                ef = obj.defgrad(eck, edk);
                
                % Computes stresses
                if obj.numint ~= 0
                    es = cell(obj.numint, 1);
                    for i = 1:obj.numint
                        es{i} = obj.stress(ef{i});
                    end
                else
                    es = obj.stress(ef);
                end
                
                % Using the above data the element forces can be computed
                felm = obj.force(eck, edk, es);
                xe = k(elm)*x(dofs)'*felm;
                y(elm) = xe;
            end
        end
        
        % OBS Only works for rectangular elements
        function Mf = dfilter(obj)
            if obj.fr > 1e-6
                wr = @(x) obj.w(x, obj.fr);
                a = obj.volumes();
                v = a*obj.t;
                xe = [obj.ec(:, 3) + obj.ec(:, 1), ...
                    obj.ec(:, 7) + obj.ec(:, 5)]/2;

                M = zeros(obj.nelm, obj.nelm);
                for i = 1:obj.nelm
                    dxi = xe(i, :) - xe;
                    ndxi = sqrt(dxi(:, 1).^2 + dxi(:, 2).^2);
                    wi = wr(ndxi).*v;
                    wi = wi/sum(wi);
                    M(i, :) = wi;
                end
                Mf = sparse(M);
            else
                Mf = 1;
            end
        end
        % ----------------------------------------------------------------%
        %%%%%%%%%%%%%%%%%% Updated Lagrangian setting %%%%%%%%%%%%%%%%%%
        % ----------------------------------------------------------------%
        % residual updated
        function r = ru(obj, edinc, P)
            % response(a, P) computes the response to the simple system given
            % the displacements a and the force P
            fext = P;
            fint = zeros(obj.ndof, 1);
            
            for elm = 1:obj.nelm
                % Finds the displacements and reference configuration for
                % this element
                dofs = obj.edof(elm, 2:end)';
                edinck = edinc(dofs);
                eck = reshape(obj.ec(elm, :), obj.npoints, 2)';
                
                % Computes the defgrad
                ef = obj.defgrad(eck, edinck, obj.fold{elm});
                
                % Computes stresses and dmat
                if obj.numint ~= 0
                    es = cell(obj.numint, 1);
                    for i = 1:obj.numint
                        es{i} = obj.stress(ef{i});
                    end
                else
                    es = obj.stress(ef);
                end
                
                % Using the above data the element stiffness matrix can be
                % computed and inserted into the correct position
                felm = obj.force(eck, obj.t, es);
                fint(dofs) = fint(dofs) + felm;
            end
            
            r = fint - fext;
        end
        
        % stiffness updated
        function Kt = Ku(obj, edinc)
            % Stiffness(a, P) computes the stiffness matrix of the simple system
            % given the displacements a and force P
            Kt = zeros(obj.ndof);
            
            % Iterates over the elements
            for elm = 1:obj.nelm
                % Finds the displacements and reference configuration for
                % this element
                dofs = obj.edof(elm, 2:end)';
                edinck = edinc(dofs);
                eck = reshape(obj.ec(elm, :), obj.npoints, 2)';
                
                % Computes the defgrad, if updated lagrangian use old ef
                ef = obj.defgrad(eck, edinck, obj.fold{elm});
                
                % Computes stresses and dmat
                if obj.numint ~= 0
                    es = cell(obj.numint, 1);
                    D = cell(obj.numint, 1);
                    for i = 1:obj.numint
                        es{i} = obj.stress(ef{i});
                        D{i} = obj.dmat(ef{i});
                    end
                else
                    es = obj.stress(ef);
                    D = obj.dmat(ef);
                end
                
                % Using the above data the element stiffness matrix can be
                % computed
                Kelm = obj.stiffness(eck, obj.t, D, es);
                
                % Inserting the element stiffness matrix into the correct pos
                disp_column = dofs(:);
                Kt(disp_column, disp_column) = Kt(disp_column, disp_column)...
                    + Kelm;
            end
            Kt = sparse(Kt);
        end
        
        % Updates the reference coordinates and the deformation gradient
        function update(obj, edinc)
            for elm = 1:obj.nelm
                dofs = obj.edof(elm, 2:end)';
                edinck = edinc(dofs);
                eck = reshape(obj.ec(elm, :), obj.npoints, 2)';
                obj.fold{elm} = ...
                    obj.defgrad(eck, edinck, obj.fold{elm});
            end
            elementEd = edinc(obj.edof(:, 2:end));
            elementIndex = 1:obj.npoints*2;
            elementIndex = [elementIndex(1:2:end) elementIndex(2:2:end)];
            elementEd = elementEd(:, elementIndex);
            obj.ec = obj.ec + elementEd;
        end
        
    end
end