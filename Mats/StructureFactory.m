classdef StructureFactory < handle
    properties
        % Structure props
        height;
        width;
        resolution;
        
        % Boundary conditions
        bcf = {};
        
        % Prescribed forces
        Ff = {};
    end
    
    methods
        function obj = StructureFactory(resolution, dimensions)
            obj.resolution = resolution;
            obj.width = dimensions(1);
            obj.height = dimensions(2);
        end
        
        function filename = make(obj, element)
            addpath(genpath('NonlinearOptimization//Mats'))
            generator = obj.getGenerator(element);
            [ex, ey, coord, dof, edof] = generator(obj.resolution);
            [nelm, np] = size(edof);
            ndof = max(edof(:));
            
            % Basic data
            w = obj.width;
            h = obj.height;
            
            xcoord = coord(:, 1)*w;
            ycoord = coord(:, 2)*h;
            ex = ex*w;
            ey = ey*h;
            
            bc = [];
            [nbc, ~] = size(obj.bcf);
            for i = 1:nbc
                fi = obj.bcf{i, 1};
                vertsi = fi(xcoord, ycoord);
                dofsi = reshape(dof(vertsi, obj.bcf{i, 2}), [], 1);
                bc = [bc; dofsi obj.bcf{i, 3}*ones(numel(dofsi), 1)];
            end
            
            F = [];
            [nF, ~] = size(obj.Ff);
            for i = 1:nF
                fi = obj.Ff{i, 1};
                vertsi = fi(xcoord, ycoord);
                dofsi = reshape(dof(vertsi, obj.Ff{i, 2}), [], 1);
                F = [F; dofsi obj.Ff{i, 3}*ones(numel(dofsi), 1)];
            end
            
            % Save file
            sform = 'beamSym%ix%iQ%i.mat';
            name = sprintf(sform, obj.resolution, (np-1)/2);
            filename = ['Mats//', name];
            save(filename, 'F', 'bc', 'edof', 'nelm', 'ndof', 'ex', ...
                'ey')
        end
        
        function addBoundaryCondition(obj, func, coords, val)
            obj.bcf{end + 1, 1} = func;
            obj.bcf{end, 2} = coords;
            obj.bcf{end, 3} = val;
        end
        
        function addPrescribedForce(obj, func, coords, val)
            obj.Ff{end + 1, 1} = func;
            obj.Ff{end, 2} = coords;
            obj.Ff{end, 3} = val;
        end
        
    end
    methods (Static)
        function generator = getGenerator(element)
            if (element.npoints == 8)
                generator = @makeMatQ8;
            elseif (element.npoints == 4)
                generator = @makeMatQ4;
            else
                disp('Factory does not support that element type');
            end
        end
    end
end