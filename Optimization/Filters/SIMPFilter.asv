classdef SIMPFilter < Filter
    
    properties
        q;
    end
    
    methods
        function obj = SIMPFilter(E, e_prop, q)
            if nargin < 3
                q = 3;
                if nargin < 3
                    e_prop = 1e-5;
                end
            end
            forward = @(z) SIMPFilter.SIMPForward(z, E, e_prop, q);
            backward = @(z) SIMPFilter.SIMPBackward(z, E, e_prop, q);
            obj@Filter(forward, backward);
            obj.q = q;
        end
    end
    
    methods (Static)
        function rho = SIMPForward(z, Emax, e_prop, q)
            rho = Emax*(e_prop + (1 - e_prop) * z.^q);
        end
        function rhoprime = SIMPBackward(z, Emax, e_prop, q)
            nz = length(z);
            rhoprime = Emax*(q * (1 - e_prop) * z.^(q - 1));
            rhoprime = spdiags(rhoprime, 0, nz, nz);
        end
    end
end