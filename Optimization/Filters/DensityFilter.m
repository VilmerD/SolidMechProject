classdef DensityFilter < Filter
    
    properties
       filter_radius;
       Mf;
    end
    
    methods
        function obj = DensityFilter(ec, areas, filter_radius, kernel)
            obj@Filter();
            if nargin < 4
                kernel = @(x, r) max(0, 1 - x/r);
            end
            obj.Mf = DensityFilter.dfilter(kernel, ec, areas, filter_radius);
            obj.filter_radius = filter_radius;
            
            obj.forward = @(z) obj.F(z);
            obj.backward = @(z) obj.dF_dz(z);
        end
        
        function y = F(obj, z)
            y = obj.Mf*z;
        end
        
        function dy = dF_dz(obj, ~)
            dy = obj.Mf;
        end
    end
    
    % Density filter matrix
    methods (Static)
        function Mf = dfilter(kernel, ec, areas, filter_radius)
            if filter_radius > 1e-6
                w = @(x) kernel(x, filter_radius);
                np = size(ec, 2)/2;
                nelm = size(ec, 1);
                xe = [mean(ec(:, 1:np), 2)'; mean(ec(:, np+1:2*np), 2)']';
                
                Im = zeros(nelm*ceil(2*pi*filter_radius^2/areas(1)), 1);
                Jm = Im;
                X = Im;
                k = 1;
                for elm = 1:nelm
                    dxi = xe(elm, :) - xe;
                    ndxi = sqrt(dxi(:, 1).^2 + dxi(:, 2).^2);
                    wi = w(ndxi).*areas;
                    wi = wi/sum(wi);
                    
                    mask = wi > 0;
                    [ii, ~] = find(mask);
                    n = length(ii);
                    Im(k:k+n-1) = ii;
                    Jm(k:k+n-1) = elm;
                    X(k:k+n-1) = wi(mask);
                    k = k+n;
                end
                mask = Im > 0;
                Mf = sparse(Im(mask), Jm(mask), X(mask));
            else
                Mf = 1;
            end
        end
    end
end