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
    end
end

