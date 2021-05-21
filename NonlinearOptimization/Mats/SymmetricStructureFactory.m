classdef SymmetricStructureFactory < handle
    properties
        baseGeometry;
        height;
    end
    
    methods
        function obj = SymmetricStructureFactory(resolution, h)
            if strcmpi(resolution, 'Coarse')
                obj.baseGeometry = 'beamSymCoarse.mat';
            else
                obj.baseGeometry = 'beamSymFine.mat';
            end
            obj.height = h;
        end
        
        function filename = makeStructure(obj, prescribeDisplacement, S)
            addpath(genpath('SolidMechanics/NonlinearOptimization/Mats'))
            load(obj.baseGeometry, 'coord', 'dof', 'edof', 'enod')
            % Basic data
            if nargin < 3
                W  = 1;
                H = 1;
            else
                W = S(1);
                H = S(2);
            end
            [nelm, ~] = size(edof);
            ndof = max(edof(:));
            
            % Scaling
            ycoord = coord(:, 2)*obj.height;
            xcoord = coord(:, 1);
            
            % Reference element
            referenceElement = 1;
            referenceDOFS = enod(referenceElement, 2:end);
            ycoordsReference = ycoord(referenceDOFS);
            elementHeight = ycoordsReference(3) - ycoordsReference(2);
            elementWidth = elementHeight;
            elementsHigh = floor(obj.height/elementHeight);
            elementsWide = floor(nelm/elementsHigh);
            
            width = elementsWide*elementWidth;
            xcoord(:, 1) = xcoord*width;
            
            % Nodes that are welded to walls
            vertLeft = abs(xcoord - 0) < 1e-6;
            dofLeft = dof(vertLeft, 1);   dofLeft = dofLeft(:);
            vertRight = abs(xcoord - width) < 1e-6;
            dofRight = dof(vertRight, :); dofRight = dofRight(:);
            dofsZeroDisplacement = [dofLeft; dofRight];
            z = zeros(size(dofsZeroDisplacement));
            
            % Displaced Nodes
            horizontal = abs(xcoord - 0) < W*elementWidth-1e-6;
            vertical = abs(ycoord - obj.height) < H*elementHeight-1e-6;
            verts_disp = find(horizontal.*vertical);
            dofs_disp = dof(verts_disp, 2);
            nz = ones(size(dofs_disp));
            
            F = zeros(ndof, 1);
            if prescribeDisplacement
                bc = [dofsZeroDisplacement, z; dofs_disp, nz];
            else
                bc = [dofsZeroDisplacement, z];
                F(verts_disp) = 1;
            end
            
            % Element data
            ex = xcoord(enod(:, 2:end));
            ey = ycoord(enod(:, 2:end));
            coord = [xcoord, ycoord];
            
            % Save file
            filename = ['SolidMechanics/NonlinearOptimization/Mats/', ...
                obj.baseGeometry(1:end-4), 'New.mat'];
            save(filename, 'F', 'bc', 'coord', 'dof', 'edof', 'enod', ...
                'nelm', 'ndof', 'ex', 'ey')
        end
        
    end
end