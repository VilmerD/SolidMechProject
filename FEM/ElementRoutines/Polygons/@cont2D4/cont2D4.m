classdef cont2D4 < abstractElement
    
    methods
        function obj = cont2D4()
            obj.npoints = 4;
            obj.numint = 4;
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

