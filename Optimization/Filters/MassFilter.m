classdef MassFilter < Filter
    % MassFilter penalizes low density regions
    
    methods
        function obj = MassFilter()
            % MassFilter Construct an instance of this class
            
            xmin = 1e-1;                    % threshold penalization value
            P1 = [1 0];                     % x
            P2 = [-5e6, 6e5, zeros(1, 6)];  % -5e6x^7 + 6e5x^6
            obj.forward = @(x) MassFilter.MassForward(x, P1, P2, xmin);
            
            dP1 = P1(1:end-1).*((length(P1)-1):-1:1);
            dP2 = P2(1:end-1).*((length(P2)-1):-1:1);
            obj.backward = @(x) MassFilter.MassBackward(x, dP1, dP2, xmin);
        end
        
        function xt = F(obj, x)
            xt = MassFilter.MassForward(x, obj.xmin, obj.q);
        end
        
        function dxt = dF(obj, x)
            dxt = MassFilter.MassBackward(x, obj.xmin, obj.q);
        end
    end
    
    methods (Static)
        function xt = MassForward(x, P1, P2, xmin)
            % MassForward penalizes low density regions
            %
            % Input:
            %       x:      values to penalize                      (n x 1)
            %
            % Output:
            %       xt      penalized values                        (n x 1)
            
            %
            xt = polyval(P1, x).*(x >= xmin) + polyval(P2, x).*(x < xmin);
        end
        
        function dxt = MassBackward(x, dP1, dP2, xmin)
            % MassBackward is the derivative of the forward penalization
            % Input:
            %       x:      values to penalize                      (n x 1)
            %
            % Output:
            %       dxt:    jacobian                                (n x n)
            %
            dxt = polyval(dP1, x).*(x >= xmin) + polyval(dP2, x).*(x < xmin);
            
            n = length(x);
            dxt = spdiags(dxt, 0, n, n);
        end
    end
end

