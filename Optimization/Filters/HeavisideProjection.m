classdef HeavisideProjection < Filter
    
    properties
        k
        b
        bupdate
        eta
    end
    
    methods
        function obj = HeavisideProjection(bupdate, eta)
            % The constructor of the superclass must be called before
            % using any properties
            obj@Filter();
            
            % Then properties can be initiated
            obj.b = 1e-1;
            obj.k = 0;
            obj.eta = eta;
            
            obj.bupdate = @() bupdate(obj.k, obj.b);
            obj.forward = @(z) obj.F(z);
            obj.backward = @(z) obj.dF_dz(z);
        end
        
        % Updates k and b according to given rule
        function update(obj)
            obj.k = obj.k + 1;
            obj.b = obj.bupdate();
        end
        
        function y = F(obj, z)
            obj.update();
            y = HeavisideProjection.H(z, obj.b, obj.eta);
        end
        
        function dy = dF_dz(obj, z)
            nz = length(z);
            d = HeavisideProjection.dH_dx(z, obj.b, obj.eta);
            dy = spdiags(d, 0, nz, nz);
        end
    end
    
    methods (Static)
        function y = H(x, b, eta)
            y = (tanh(b*eta) + tanh(b*(x - eta)))./...
                (tanh(b*eta) + tanh(b*(1 - eta)));
        end
        
        function y = dH_dx(x, b, eta)
            y = b*cosh(b*(x - eta)).^(-2)./...
                (tanh(b*eta) + tanh(b*(1 - eta)));
        end
    end
end