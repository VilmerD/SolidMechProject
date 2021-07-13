classdef SIMPFilter < Filter
    
    properties
        E;
        e_prop;
        q;
    end
    
    methods
        function obj = SIMPFilter(E, e_prop, q)
            obj@Filter();
            if nargin == 0
                E = 1;
                e_prop = 1e-5;
                q = 3;
            end
            obj.E = E;
            obj.e_prop = e_prop;
            obj.q = q;
            
            obj.forward = @(z) obj.F(z);
            obj.backward = @(z) obj.dF_dz(z);
        end
        
        function y = F(obj, z)
            y = SIMPFilter.SIMPForward(z, obj.E, obj.e_prop, obj.q);
        end
        
        function dy = dF_dz(obj, z)
            dy = SIMPFilter.SIMPBackward(z, obj.E, obj.e_prop, obj.q);
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