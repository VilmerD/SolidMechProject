classdef StructureFactory    
    properties
        baseGeometry;
        height;
    end
    
    methods
        function obj = StructureFactory(resolution, h)
            if strcmpi(resolution, 'Coarse')
                obj.baseGeometry = ...
                    'beamFullCoarse.mat';
            elseif strcmpi(resolution, 'Fine')
                obj.baseGeometry = ...
                    'beamFullFine.mat';
            end
            obj.height = h;
        end
        
        function filename = makeStructure(obj, prescribe_force, W, H)
            load(obj.baseGeometry, 'coord', 'dof', 'edof', 'enod')
            % Basic data
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
            vert_left = abs(xcoord - 0) < 1e-6;
            dof_left = dof(vert_left, :);   dof_left = dof_left(:);
            vert_right = abs(xcoord - width) < 1e-6;
            dof_right = dof(vert_right, :); dof_right = dof_right(:);
            dofs_zerodisp = [dof_left; dof_right];
            z = zeros(size(dofs_zerodisp));
            
            % Displaced Nodes
            horizontalCenter = abs(xcoord - width/2) < W*elementWidth-1e-6;
            verticalCenter = abs(ycoord - obj.height) < H*elementHeight-1e-6;
            verts_disp = find(horizontalCenter.*verticalCenter);
            dofs_disp = dof(verts_disp, 2);
            nz = ones(size(dofs_disp));
            
            F = zeros(ndof, 1);
            if prescribe_force
                bc = [dofs_zerodisp, z; dofs_disp, nz];
            else
                bc = [dofs_zerodisp, z];
                F(verts_disp) = 1;
            end
            
            % Element data
            ex = xcoord(enod(:, 2:end));
            ey = ycoord(enod(:, 2:end));
            coord = [xcoord, ycoord];
            
            % Save file
            filename = ['NonlinearOptimization\Mats\', ...
                obj.baseGeometry(1:end-4), 'New.mat'];
            save(filename, 'F', 'bc', 'coord', 'dof', 'edof', 'enod', ...
                'nelm', 'ndof', 'ex', 'ey')
        end
    end
end

