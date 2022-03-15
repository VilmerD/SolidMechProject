classdef ExponentialProjection < Filter
    % ExponentialProjection is a min-max operator
    
    properties
        k
        b
        b0 = 0;
        bupdate
        
        type;
    end
    
    methods
        function obj = ExponentialProjection(params)
            obj.b = obj.b0;
            obj.k = 0;
            
            obj.bupdate = @() ExponentialProjection.bupdateLinear(...
                obj.b, obj.k, params{1:end-1});
            
            
            obj.type = params{end};
            if strcmpi(type, 'dial')
                obj.forward = @(x) obj.F(x);
                obj.backward = @(x) obj.dF(x);
            elseif strcmpi(type, 'ero')
                obj.forward = @(x) 1 - obj.F(1 - x);
                obj.backward = @(x) dF(1 - x);
            end
        end
        
        function xt = F(obj, x)
            xt = ExponentialProjection.MME(x, obj.b);
        end
        
        function dxt = dF(obj, x)
            dxt = ExponentialProjection.dMME(x, obj.b);
        end
            
    end
    
    methods (Static)
        function xt = MME(x, b)
            xt = 1 - exp(-b*x) + x*exp(-b);
        end
        
        function dxt = dMME(x, b)
            dxt = b*exp(-b*x) + exp(-b);
            
            n = length(x);
            dxt = spdiags(dxt, 0, n, n);
        end
        
        function bnew = bupdateLinear(b, k, k0, bmax, bstep, bincrement)
            bupatenow = ~mod(k - k0, bstep);
            bthreshold = (k >= k0);
            bstepk = min(bincrement, bmax - b);
            bnew = b + bstepk*bupatenow*bthreshold;
        end
    end
end

