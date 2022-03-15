classdef HeavisideProjection < Filter
    % HeavisideProjection projects onto 0-1 space 
    
    properties
        k           % uteration number used for ramping the projection
        b           % beta-value is the steepness of the projection
        b0 = 1e-9;  % starting b-value
        bupdate     % update function
        eta         % eta-value is the midpoint of projection
    end
    
    methods
        function obj = HeavisideProjection(params)
            % HeavisideProjection creates an instance of the class
            %
            % Input:
            %       bupdate:        update function                       1
            %       eta:            eta-value                             1
            %
            
            obj.b = obj.b0;         % Starting value of b
            obj.k = 0;              % Iteration number
            obj.eta = params{end};          % The midpoint of the projection
            
            obj.bupdate = @() HeavisideProjection.bupdateLinear(...
                obj.b, obj.k, params{1:end-1});
            obj.forward = @(x) obj.F(x);
            obj.backward = @(x) obj.dF(x);
        end
        
        % Updates k and b according to given rule
        function update(obj, k)
            if k == 0
                obj.b = obj.b0;
                obj.k = 0;
            else
                obj.k = k;
                obj.b = obj.bupdate();
            end
        end
        
        function xt = F(obj, x)
            xt = HeavisideProjection.H(x, obj.b, obj.eta);
        end
        
        function dxt = dF(obj, x)
            dxt = HeavisideProjection.dH(x, obj.b, obj.eta);
        end
    end
    
    methods (Static)
        function xt = H(x, b, eta)
            % H is an approximation of the heaviside step function
            % 
            % Input:
            %       x:      values to project                       (n x 1)
            %       b:      projection steepness                          1
            %       eta:    projection midpoint                           1
            % 
            % Output:
            %       xt:     projected values                        (n x 1)
            
            xt = (tanh(b*eta) + tanh(b*(x - eta)))./...
                (tanh(b*eta) + tanh(b*(1 - eta)));
        end
        
        function dxt = dH(x, b, eta)
            % dH is an approximation of the heaviside step function
            % 
            % Input:
            %       x:      values to project                       (n x 1)
            %       b:      projection steepness                          1
            %       eta:    projection midpoint                           1
            % 
            % Output:
            %       dxt:    jacobian                                (n x n)
            
            dxt = b*cosh(b*(x - eta)).^(-2)./...
                (tanh(b*eta) + tanh(b*(1 - eta)));
            
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