classdef StructureFactory < handle
    properties
        % Structure props
        height;
        width;
        resolution;
        
        % Boundary conditions
        bcf = {};
        
        % Prescribed forces
        ff = {};
        
        % Spring elements
        sprf = {};
    end
    
    methods
        function obj = StructureFactory(resolution, dimensions)
            obj.resolution = resolution;
            obj.width = dimensions(1);
            obj.height = dimensions(2);
        end
        
        function make(obj, npoints, filename)
            generator = obj.getGenerator(npoints);
            [ex, ey, coord, dof, edof, enod] = generator(obj.resolution);
            [nelm, ~] = size(edof);
            ndof = max(edof(:));
            
            % Basic data
            xcoord = coord(:, 1)*obj.width;
            ycoord = coord(:, 2)*obj.height;
            ex = ex*obj.width;
            ey = ey*obj.height;
            
            % Add boundary conditions
            bc = [];
            [nbc, ~] = size(obj.bcf);
            for i = 1:nbc
                fi = obj.bcf{i, 1};
                vertsi = fi(xcoord, ycoord);
                dofsi = reshape(dof(vertsi, obj.bcf{i, 2}), [], 1);
                bc = [bc; dofsi obj.bcf{i, 3}*ones(numel(dofsi), 1)];
            end
            
            % Add prescribed forces
            F = zeros(ndof, 1);
            [nF, ~] = size(obj.ff);
            for i = 1:nF
                fi = obj.ff{i, 1};
                vertsi = fi(xcoord, ycoord);
                dofsi = reshape(dof(vertsi, obj.ff{i, 2}), [], 1);
                F(dofsi) = obj.ff{i, 3}*ones(numel(dofsi), 1);
            end
            
            % Add 1d-springs
            nspring_conditions = size(obj.sprf, 1);
            edofs = [];
            for i = 1:nspring_conditions
                % Find vertices and dofs of new springs
                fi = obj.sprf{i, 1};
                vertsi = fi(xcoord, ycoord);
                dofsi = dof(vertsi, obj.sprf{i, 2});
                
                % Add dofs for enpoints of springs
                nsprings = size(dofsi, 1);
                spring_elements = size(edofs, 1) + (1:nsprings)';
                newdofs = ndof + 1:nsprings;
                
                % Add to spring-edof
                edofsnew = [spring_elements dofsi newdofs];
                edofs = [edofs; edofsnew];
                
                % Coordinates
                coordfunci = obj.sprf{i, 3};
                coordsold = [xcoord(vertsi) ycoord(vertsi)];
                coordsnew = coordfunci(coordsold(:, 1), coordsold(:, 2)); 
                
                xcoord = [xcoord; coordsnew(:, 1)];
                ycoord = [ycoord; coordsnew(:, 2)];
                ndof = newdofs(end);
            end
            
            % Save file
            width = obj.width;
            height = obj.height;
            save(filename, 'F', 'bc', 'edof', 'edofs', 'nelm', 'ndof', 'ex', ...
                'ey', 'enod', 'width', 'height')
        end
        
        function addBoundaryCondition(obj, func, coords, val)
            obj.bcf{end + 1, 1} = func;
            obj.bcf{end, 2} = coords;
            obj.bcf{end, 3} = val;
        end
        
        function addPrescribedForce(obj, func, coords, val)
            obj.ff{end + 1, 1} = func;
            obj.ff{end, 2} = coords;
            obj.ff{end, 3} = val;
        end
        
        function addSpring(obj, func, dof, coords)
            obj.sprf{end + 1, 1} = func;
            obj.sprf{end, 2} = dof;
            obj.sprf{end, 3} = coords;
        end
        
    end
    methods (Static)
        function generator = getGenerator(npoints)
            if (npoints == 8)
                generator = @makeMatQ8;
            elseif (npoints == 4)
                generator = @makeMatQ4;
            else
                disp('Factory does not support that element type');
            end
        end
    end
end