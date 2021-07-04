classdef abstractElement < handle
% Element interface that is used as reference in the nonlinear model
    
    properties
        npoints;
        numint;
    end
    
    methods (Abstract, Static)
        baseFuncs(); 
        area(dJ);
        matrices(ec, dN);
        defgrad(H, ed);
        force(H, B, dJ, t, defgrad, stress);
        stiffness(H, B, dJ, t, D, defgrad, stress);
    end
    
    methods
        function M = precompute(obj, ec)
            [nelm, ~] = size(ec);
            M = cell(nelm, 3);
            dN = obj.baseFuncs();
            for elm = 1:nelm
                eck = reshape(ec(elm, :), obj.npoints, 2)';
                [dJ, H, B] = obj.matrices(eck, dN);
                M{elm, 1} = dJ;
                M{elm, 2} = H;
                M{elm, 3} = B;
            end
        end
    end
end

