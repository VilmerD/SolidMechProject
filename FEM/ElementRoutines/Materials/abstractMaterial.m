classdef abstractMaterial < handle
    
    properties
        mpara;
        flag;
    end
    
    methods (Abstract, Static)
        DMatrix(flag, mpara, defgrad)
        Stress(flag, mpara, defgrad)
    end
    
    methods
        function D = dm(obj, defgrad)
            if iscell(defgrad)
                np = numel(defgrad);
                D = cell(np, 1);
                for i = 1:np
                    D{i} = obj.DMatrix(obj.flag, obj.mpara, defgrad{i});
                end
                return
            else
                D = obj.DMatrix(obj.flag, obj.mpara, defgrad{i});
                return
            end
        end
        
        function S = es(obj, defgrad)
            if iscell(defgrad)
                np = numel(defgrad);
                S = cell(np, 1);
                for i = 1:np
                    S{i} = obj.Stress(obj.flag, obj.mpara, defgrad{i});
                end
                return
            else
                S = obj.Stress(obj.flag, obj.mpara, defgrad{i});
                return
            end
        end
    end
end

