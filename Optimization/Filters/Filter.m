classdef Filter < handle
    %Filter base class. Used to construct specific filters.
    
    properties
        forward
        backward
    end
    
    methods
        function obj = Filter(forward, backward)
            if nargin == 0
                forward = @(x) x;
                backward = @(x) 1;
            end
            obj.forward = forward;
            obj.backward = backward;
        end
        
        % In reality the filters are composed with the order obj2(obj1) ie
        % the composition is read right to left
        function newobj = times(obj1, obj2)
            f_new = @(z) obj1.forward(obj2.forward(z));
            b_new = @(z) obj1.backward(obj2.forward(z)) * obj2.backward(z);
            newobj = Filter(f_new, b_new);
        end
        
        function newobj = mtimes(obj1, obj2)
            newobj = times(obj1, obj2);
        end
        
        function newobj = plus(obj1, obj2)
            newobj = times(obj1, obj2);
        end
    end
    
    methods (Static)
        function obj = Identity()
            obj = Filter.Scalar(1);
        end
        
        function obj = Scalar(a)
           obj = Filter(@(x) a*x, @(x) a);
        end
    end
end