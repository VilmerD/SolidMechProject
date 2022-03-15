classdef DensityFilter < Filter
    % DensityFilter filters with a density based approach
    
    properties
        radius;         % filter radius
        Mf;             % filter matrix
    end
    
    methods
        function obj = DensityFilter(ex, ey, areas, radius, kernel)
            % DensityFilter Constructs an instance of this class
            % 
            % Input: 
            %       ex:         element x-coordinates
            %       ey:         element y-coordinates
            %       areas:      element areas
            %       radius:     filter radius
            %       kernel:     filter kernel
            %
            
            % Initialize filter class
            obj@Filter();
            
            % Default kernel is hat function
            if nargin < 5
                kernel = @(x, r) max(0, 1 - x/r);
            end
            obj.Mf = DensityFilter.filterMatrix(ex, ey, areas, ...
                radius, kernel);
            obj.radius = radius;
            
            obj.forward = @(z) obj.F(z);
            obj.backward = @(z) obj.dFdz(z);
        end
        
        function x = F(obj, z)
            x = obj.Mf*z;
        end
        
        function dxt = dFdz(obj, ~)
            dxt = obj.Mf;
        end
    end
    
    % Density filter matrix
    methods (Static)
        function Mf = filterMatrix(ex, ey, areas, radius, kernel)
            % filterMatrix assembles the filter matrix
            %
            % Input:
            %       ex:         element x-coordinates
            %       ey:         element y-coordinates
            %       areas:      element areas
            %       radius:     filter radius
            %       kernel:     filter kernel
            %
            % Output:
            %       Mf:         filter matrix
            %
            
            nelm = size(ex, 1);
            
            % If radius is too small no filtering is done
            if radius < 1e-9
                Mf = speye(nelm);
                return
            end
            w = @(x) kernel(x, radius);
            
            % Midpoints of elements
            element_centers = [mean(ex, 2)'; mean(ey, 2)']';
            
            % Preallocate the I and J vectors
            Im = zeros(nelm*ceil(pi*radius^2/areas(1)), 1);
            Jm = Im;
            X = Im;
            k = 1;
            
            % Iterate through each element
            for elm = 1:nelm
                % Distance from element midpoint to other elements
                % midpoints
                center_dist = element_centers(elm, :) - element_centers;
                
                % Evaluate kernel on center distances
                ndxi = sqrt(center_dist(:, 1).^2 + center_dist(:, 2).^2);
                wi = w(ndxi);
                
                % Extract values which are nonzero for sparse structure
                mask = wi > 0;
                wi_nonzero = wi(mask);
                wi_nonzero = wi_nonzero/sum(wi_nonzero);
                
                % Extract indices of nonzero elements
                [ii, ~] = find(mask);
                n = length(ii);
                
                % Insert nonzero elements
                Im(k:k+n-1) = elm;
                Jm(k:k+n-1) = ii;
                X(k:k+n-1) = wi_nonzero;
                k = k+n;
            end
            % There are possibly some zeros in Im and Jm, so extract
            % nonzero elements (otherwise sparse gets mad)
            mask = Im > 0;
            Mf = sparse(Im(mask), Jm(mask), X(mask));
        end
    end
end