classdef RAMPFilter < Filter
    % RAMPFilter penalizes using the RAMP method
    
    properties
        xmin;   % minimum x-value
        q;      % penalty exponent
    end
    
    methods
        function obj = RAMPFilter(xmin, q)
            % RAMPFilter Constructs an instance of this class
            % 
            % Input:
            %       xmin:       minimum x-value
            %       q:          penalty exponent
            % 
            
            % Default values
            if nargin == 0
                xmin = 1e-4;
                q = 8;
            end
            obj.Emin = xmin;
            obj.q = q;
            
            obj.forward = @(x) obj.F(x);
            obj.backward = @(x) obj.dF(x);
        end
        
        function xt = F(obj, x)
            xt = RAMPFilter.RAMPForward(x, obj.xmin, obj.q);
        end
        
        function dxt = dF(obj, x)
            dxt = RAMPFilter.RAMPBackward(x, obj.xmin, obj.q);
        end
    end
    
    methods (Static)
        function xt = RAMPForward(x, xmin, q)
            % RAMPForward penalizes x using RAMP
            % 
            % Input:
            %       x:      values to penalize                      (n x 1)
            %       xmin:   minimal x-value                               1
            %       q:      penalty exponent                              1
            %
            % Output:
            %       xt:     penalized values  
            
            xt = xmin + x*(1 - xmin)./(1 + q*(1 - x));
        end
        
        function dxt = RAMPBackward(x, xmin, q)
            % RAMPBackward is the derivative of the forward penalization
            % 
            % Input:
            %       x:      values to penalize                      (n x 1)
            %       xmin:   minimal x-value                               1
            %       q:      penalty exponent                              1
            %
            % Output:
            %       dxt:    jacobian                                (n x n)
            
            dxt = (1 - xmin)*(1 + q)./(1 + q*(1 - x)).^2;
            n = length(x);
            dxt = spdiags(dxt, 0, n, n);
        end
    end
end

