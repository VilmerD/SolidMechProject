classdef NHCont < abstractMaterial
    methods
        function obj = NHCont(flag, mpara)
            obj.flag = flag;
            obj.mpara = mpara;
        end
    end
    
    methods (Static)
        D = DMatrix(flag, mpara, defgrad);
        S = Stress(flag, mpara, defgrad);
        E = Energy(mpara, defgrad);
    end
end

