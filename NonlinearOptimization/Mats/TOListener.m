classdef TOListener < handle
    
    properties
        statistics
    end
    
    methods
        function obj = TOListener()
            obj.statistics.g0 = [];
            obj.statistics.g1 = [];
            obj.statistics.designs = [];
            obj.statistics.factorizations = [];
            obj.statistics.lineqs = [];
        end
        
        function registerUpdate(obj, statistics)
            obj.statistics.g0(end + 1) = statistics.g0;
            obj.statistics.g1(end + 1) = statistics.g1;
            obj.statistics.designs(:, end+1) = statistics.design;
            obj.statistics.factorizations(end + 1) = statistics.factorizations;
            obj.statistics.lineqs(end + 1) = statistics.ncalls; 
        end
        
        function registerCustom(obj, name, value)
            try
                obj.statistics.(name)(:, end + 1) = value;
            catch e
                obj.statistics.(name) = [];
                obj.statistics.(name)(:, 1) = value;
            end
        end
    end
end

