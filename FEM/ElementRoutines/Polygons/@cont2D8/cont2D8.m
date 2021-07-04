classdef cont2D8 < abstractElement
    
    methods
        function obj = cont2D8()
            obj.npoints = 8;
            obj.numint = 9;
        end
    end
    
    methods (Static)
        A = area(dJ);
        [dJ, H, B] = matrices(ec, dN);
        ef = defgrad(H, ed);
        felm = force(H, B, dJ, t, defgrad, stress);
        kelm = stiffness(H, B, dJ, t, D, defgrad, stress);
        dN = baseFuncs();
    end
end

