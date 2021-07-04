classdef NHAxi < abstractMaterial
    methods
        function obj = NHAxi(flag, mpara)
            obj.flag = flag;
            obj.mpara = mpara;
        end
    end
    
    methods (Static)
        D = DMatrix(flag, mpara, defgrad);
        S = Stress(flag, mpara, defgrad);
    end
end

