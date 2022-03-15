classdef SIMPFilter < Filter
    % SIMPFilter penalizes using the SIMP method
    
    properties
        xmin;       % Minimal x-value
        q;          % Penalty exponent
    end
    
    methods
        function obj = SIMPFilter(xmin, q)
            % SIMPFilter Constructs an instance of this class
            % 
            % Input:
            %       xmin:   Minimum x-value                               1
            %       q:      Penalty exponent                              1
            % 
            
            % Default value
            if nargin == 0
                xmin = 1e-4;
                q = 3;
            end
            obj.xmin = xmin;
            obj.q = q;
            
            obj.forward = @(x) obj.F(x);
            obj.backward = @(x) obj.dF(x);
        end
        
        function xt = F(obj, x)
            xt = SIMPFilter.SIMPForward(x, obj.xmin, obj.q);
        end
        
        function dxt = dF(obj, x)
            dxt = SIMPFilter.SIMPBackward(x, obj.xmin, obj.q);
        end
    end
    
    methods (Static)
        function xt = SIMPForward(x, xmin, q)
            % SIMPForward penalizes x using SIMP 
            % 
            % Input:
            %       x:      Values to penalize                      (n x 1)
            %       xmin:   Minimal x-value                               1
            %       q:      Penalty exponent                              1
            %
            % Output:
            %       xt:     Penalized values                        (n x 1)
            
            xt = xmin + (1 - xmin)*x.^q;
        end
        
        function dxt = SIMPBackward(x, xmin, q)
            % SIMPBackward is the derivative of the forward penalization
            % 
            % Input:
            %       x:      Values to penalize                      (n x 1)
            %       xmin:   Minimal x-value                               1
            %       q:      Penalty exponent                              1
            %
            % Output:
            %       dxt:    Jacobian                                (n x n)
            
            dxt = (1 - xmin)*q*x.^(q - 1);
            
            n = length(x);
            dxt = spdiags(dxt, 0, n, n);
        end
    end
end